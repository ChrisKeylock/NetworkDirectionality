function [betatau,rhotau,epstau,Ztau]=BetaTau(A,n1,n2)

%Takes a matrix A, obtains epstau for it and then uses this value to find
%betatau. As part of thix process, the full Ztau contour is defined and
%this can also be returned.

%In the context of our work (Phys. Rev. E, 2026), A should be the DGL.

%n1 is the number of points used on the transect
%n2 is the discretisation of the real line

%The defaults for n1 and n2 are for high precision and are therefore slow
%to run: reducing these values speeds up the code significantly

if nargin<3
    n2=201;
end

if nargin<2
    n1=71;
end

sizer=length(A);

%Look for the maximum between the eigenvalues on the real line
eA=eig(A);
ReA=sort(real(eA));
min_e=ReA(1)-(eps*1000);
max_e=ReA(sizer)+(eps*1000);

pos=linspace(min_e,max_e,n1);
for loop1=1:length(pos)
    Transect_tau(loop1,1)=pse_Point(A,pos(loop1));
end
[~,loc]=max(Transect_tau);

%Interpolate the maximum
pp = csapi(real(pos),Transect_tau);
fprime = fnder(pp,1);
fprime2 = fnder(pp,2);
count=1;
for loop1=loc-1:loc
    xint=fprime.breaks(loop1+1)-fprime.breaks(loop1);
    deriv(count,1)=fprime.coefs(loop1,3);
    deriv(count,2)=fprime.coefs(loop1,1)*xint^2+fprime.coefs(loop1,2)*xint+fprime.coefs(loop1,3);
    deriv(count,3)=fprime2.coefs(loop1,2);
    deriv(count,4)=fprime2.coefs(loop1,1)*xint+fprime.coefs(loop1,2);
    count=count+1;
end
%Extract with change in sign on deriv1 and negative on deriv2
temp=find(sign(deriv(:,1))~=sign(deriv(:,2)) & sign(deriv(:,3))<0 & sign(deriv(:,4))<0);
%Now solve
if temp==1
    temp=loc-1;
else
    temp=loc;
end
xloc(1)=(-fprime.coefs(temp,2)+sqrt(fprime.coefs(temp,2)^2-4*fprime.coefs(temp,1)*fprime.coefs(temp,3) ))/(2*fprime.coefs(temp,1));
xloc(2)=(-fprime.coefs(temp,2)-sqrt(fprime.coefs(temp,2)^2-4*fprime.coefs(temp,1)*fprime.coefs(temp,3) ))/(2*fprime.coefs(temp,1));
xint=fprime.breaks(temp+1)-fprime.breaks(temp);

for loop1=1:2
    yloc(loop1)=pp.coefs(temp,1)*xloc(loop1)^3+pp.coefs(temp,2)*xloc(loop1)^2+pp.coefs(temp,3)*xloc(loop1)+pp.coefs(temp,4);
end

epstau=max(yloc)+(eps*10000);
imagtol = 1e-10*max(norm(A),epstau);

theta(:,1)=linspace(0, pi,n2)';
O=zeros(sizer);
mE=-epstau*eye(sizer);

for j=1:n2
   % I*e^(i*theta(j))
   eith = (cos(theta(j)) + 1i*sin(theta(j))) * eye(sizer);
   % I*e^(-i*theta(j))
   emith = (cos(theta(j)) - 1i*sin(theta(j))) * eye(sizer);
   
   % compute the generalized eigenvalues for the pencil H - lambda*K
   H = [mE A; A' mE];
   K = [O eith; emith O];
   eHK = eig(H,K);
         
   if min(abs(imag(eHK))) <= imagtol % check if pencil H-lambda*K has a real generalized eigenvalue
      indx = find(abs(imag(eHK)) <= imagtol);  % extract such eigenvalues
      ZCont = real(max(real(eHK(indx))));
      Ztau(j,1)=ZCont*cos(theta(j)) + 1i*ZCont*sin(theta(j));
   end 
   
end

% flag=0;
% while flag==0
    %loop until there are no outlier distances
    clear diff newang newvec
    %Look for large changes in the coordinates and add extra angles in here
    diff(:,1)=Ztau;
    diff(:,2)=circshift(Ztau,[1 0]);
    diff(:,3)=abs(diff(:,2)-diff(:,1));
    diff(1,3)=0;
    typicaldist=median(diff(2:length(diff),3));
    MAD=1.4826*median(abs(diff(2:length(diff),3)-typicaldist));
    
    thresh=sqrt(2*log(length(diff)))*MAD;

    %find the locations exceeding the threshold and then in-fill an appropriate
    %number of points
    clear temp
    temp=find(diff(:,3)-typicaldist>thresh);
    if isempty(temp)
        flag=1;
    else
        multiples=floor((diff(temp,3)-typicaldist)/thresh);
        newvec=zeros(sum(multiples),1);
        count=1;
        for loop1=1:length(temp)
            vectorR=linspace(real(diff(temp(loop1),1)),real(diff(temp(loop1),2)),multiples(loop1)+2);
            vectorI=linspace(imag(diff(temp(loop1),1)),imag(diff(temp(loop1),2)),multiples(loop1)+2);
            newvec(count:count+length(vectorR)-3,1)=vectorR(2:length(vectorR)-1)+1i*vectorI(2:length(vectorR)-1);
            count=count+length(vectorR)-2;
        end
        newang=angle(newvec);

        for j=1:length(newang)
            % I*e^(i*theta(j))
            eith = (cos(newang(j)) + 1i*sin(newang(j))) * eye(sizer);
            % I*e^(-i*theta(j))
            emith = (cos(newang(j)) - i*sin(newang(j))) * eye(sizer);
   
            % compute the generalized eigenvalues for the pencil H - lambda*K
            H = [mE A; A' mE];
            K = [O eith; emith O];
            eHK = eig(H,K);
         
            if min(abs(imag(eHK))) <= imagtol % check if pencil H-lambda*K has a real generalized eigenvalue
                indx = find(abs(imag(eHK)) <= imagtol);  % extract such eigenvalues
                ZCont = real(max(real(eHK(indx))));
                Ztemp(j,1)=ZCont*cos(newang(j)) + 1i*ZCont*sin(newang(j));
            end 
        end

        theta=[theta;newang]; [theta,pos]=sort(theta);
        Ztau=[Ztau;Ztemp]; Ztau=Ztau(pos);
    end
% end

        
[rhotau(1),rhotau(2)]=max(abs(Ztau));
rhotau(2)=theta(rhotau(2));
[betatau(1),betatau(2)]=max(imag(Ztau(:,1)));
Ztau(:,2)=theta;

end


%%%%%%%%%%%%%%%%%%%%%%
function Z=pse_Point(A,coord);

%Trefethen's 1999 inverse Lanczos code is adopted here

%Make Schur form
[U,T]=schur(A,'complex');
T=triu(T);
eA=diag(T);

%Rotate the Schur form with loop1 the current eigenvalue
%and loop2 all of the others, so we make lots of 2x2 matrices and rotate
%them into a subspace
T2=T;
for loop1 = 1:length(eA)
   for loop2 = length(eA)-1:-1:loop1
      G([2 1],[2 1]) = planerot([T2(loop2,loop2+1) T2(loop2,loop2)-T2(loop2+1,loop2+1)]')';
      T2(:,[loop2:loop2+1]) = T2(:,[loop2:loop2+1])*G; T2([loop2:loop2+1],:) = G'*T2([loop2:loop2+1],:);
   end
end
T2 = triu(T2(1:length(eA),1:length(eA)));

[L,U,P] = lu(T2-coord*eye(length(eA)));
L1 = L'; U1 = U';

sigold = 0; qold = zeros(length(eA),1); beta = 0; H = [];
q = randn(length(eA),1) + sqrt(-1)*randn(length(eA),1); q = q/norm(q);
for loop1=1:99
     v = L1\(U1\(U\(L\q))) - beta*qold;
     %v = T1\(Tt\q)-beta*qold;
     alpha = real(q'*v); 
     v = v - alpha*q;
     beta = norm(v); 
     qold = q; 
     q = v/beta;
     H(loop1+1,loop1) = beta; 
     H(loop1,loop1+1) = beta; 
     H(loop1,loop1) = alpha;
     try
        sig = max(eig(H(1:loop1, 1:loop1)));
     catch
        sig=1e308; break;
     end;
     if (abs(sigold/sig-1)<.001), break; end;
     sigold = sig; 
end 
Z = 1/sqrt(sig);
end

