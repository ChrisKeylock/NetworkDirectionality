function output=DGL(A)

sizer=size(A);
if sizer(1)~=sizer(2)
    return
end

outD=sum(A')';
D=eye(size(A)).*outD;
P=inv(D)*A;
[~,eP,LEP]=eig(P); eP=real(diag(eP));
temp=find(abs(eP)>=1-(eps*1000));
StatVec=LEP(:,temp);
if sum(StatVec)<0
    StatVec=-StatVec;
end
StatVec=StatVec./sum(StatVec);
Stat=eye(length(A)).*((StatVec));

Pval=(Stat^0.5*P*Stat^-0.5);

output = eye(sizer(1))-((Stat^0.5)*P*(Stat^(-0.5)));
