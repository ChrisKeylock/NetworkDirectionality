function [betatau,rhotau,Ztau]=BetaTauImposeEpsTau(A,epstau,n2)

%Takes a matrix A and a value for epstau and then uses this value to find
%betatau. As part of this process, the full Ztau contour is defined and
%this can also be returned.

%In terms of our directionality work (Phys. Rev. E, 2026), A should be the
%DGL

%n2 is the discretisation for theta. Reducing this value will speed up the
%code significantly

if nargin<3
    n2=201;
end


sizer=length(A);

imagtol = 1e-10*max([norm(A),epstau]);

%Apply epstau to obtain the Zcontour
theta=linspace(0,pi,n2)';
O=zeros(sizer);
mE=-epstau*eye(sizer);

for j=1:n2
   % I*e^(i*theta(j))
   eith = (cos(theta(j)) + 1i*sin(theta(j))) * eye(sizer);
   % I*e^(-i*theta(j))
   emith = (cos(theta(j)) - i*sin(theta(j))) * eye(sizer);
   
   % compute the generalized eigenvalues for the pencil H - lambda*K
   H = [mE A; A' mE];
   K = [O eith; emith O];
   eHK = eig(H,K);
         
   if min(abs(imag(eHK))) <= imagtol % check if pencil H-lambda*K has a real generalized eigenvalue
      indx = find(abs(imag(eHK)) <= imagtol);  % extract such eigenvalues
      ZCont = abs(max(real(eHK(indx))));
      Ztau(j,1)=ZCont*cos(theta(j)) + 1i*ZCont*sin(theta(j));
   end 
   
end

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



[rhotau(1),rhotau(2)]=max(abs(Ztau));
rhotau(2)=theta(rhotau(2));
[betatau(1),betatau(2)]=max(imag(Ztau(:,1)));
Ztau(:,2)=theta;

end