function [indices,output]=NetworkLaplacianDirectionality(Adjacency,n1,n2)

% This is the primary function for the analysis of the directed graph 
% Laplacian (DGL) described in Keylock and Carbone (2026, Phys. Rev. E)

% It takes a square adjacency matrix, calculates the DGL following 
% Li and Zhang (2012), displaces the origin to improve the accuracy of 
% calculations by angle given the eigenvalue at the origin for the DGL 
% (see paper for discussion), obtains the beta contour for the DGL and for
% the normal variant (referred to as SchurForm, below), and then obtains
% d_epsilon, the distance between these contours by projecting back on to
% the beta contour for SchurForm using the pseudoeigenvalue surface for 
% the original DGL. These final variables are called Ztau1atprofile and
% then d_eps

%Summary indices are also returned including the Li and Zhang (2012) 
% statistic and the variant of this described in the paper. The naming for
% these indices follows that in the paper.

%Each of the variable names follows closely to the notation in the paper
%and they are all returned as part of output. E.g. output.d_eps is the
%value for d_eps.

%The computational cost increases with large values for the discretisation
%terms, n1 and n2

%After running this code, PlottingToolsUtility.m can be passed the output
% variable from this function to generate some relevant graphical output.

% An example test case to run this code: 
% Adjacency = magic(11).*(1-eye(11));
% [indices,output]=NetworkLaplacianDirectionality(Adjacency,100,70)

if nargin<3
    n2=151;
end
if nargin<2
    n1=81;
end

sizer=size(Adjacency);
if size(1)~=size(2)
   disp('Input is not valid as the adjacency matrix must be square');
   return
end

temp=find(diag(Adjacency)~=0);
if ~isempty(temp)
    disp('Input is not valid as there are self-loops in the adjacency matrix');
   return
end

temp=find(Adjacency(:)<0);
if ~isempty(temp)
    disp('Input is not a valid adjacency matrix');
   return
end


dgl=DGL(Adjacency);

[U,T]=schur(dgl,'complex');
N=triu(T,1);
L=T-N;
SchurForm=U*L*U';
LiZhang=(dgl-dgl')/2;


dgl=dgl-eye(sizer);
SchurForm=SchurForm-eye(sizer);

[betatau{1},rhotau{1},epstau,Ztau{1}]=BetaTau(dgl,n1,n2);
[betatau{2},rhotau{2},Ztau{2}]=BetaTauImposeEpsTau(SchurForm,epstau(1),n2);

%We work radially to avoid problem of multiple values when the Schur 
% form tends back to the origin at gaps between eigenvalues

%Interpolate back onto a regular mesh
angs=linspace(0,pi,2*n2);
xVals=real(Ztau{2}(:,2));
yVals=abs(Ztau{2}(:,1));
ychip=pchip(xVals,yVals,angs);
zcoords=ychip.*(cos(angs)+1i*sin(angs));
for loop1=1:length(ychip)
    zmat=(zcoords(1,loop1).*eye(length(dgl))-dgl);
    [~,Ex,~]=svd(zmat);
    Ztau1atprofile(loop1,1)=min(diag(Ex));
end
d_eps(:,1)=(epstau-Ztau1atprofile)/epstau;
%Force to zero in case numerical/discretisation results in negative values
temp=find(d_eps(:,1)<0);
d_eps(temp,1)=0;
d_eps(:,2)=angs;

output.d_eps=d_eps;
output.epstau=epstau;
output.Ztau=Ztau;
output.dgl=dgl;
output.SchurForm=SchurForm;

indices.epstau=epstau;
indices.Dir1=betatau{1}(1);
indices.RelDir1=(betatau{1}(1)-betatau{2}(1))/betatau{2}(1);
indices.RelDir2=median(d_eps(:,1));

[~,EL,~]=svd(LiZhang); E(:,1)=diag(EL);
[~,EN,~]=svd(N); E(:,2)=diag(EN);
indices.Dir2=max(E(:,1));
indices.Dir2b=max(E(:,2))*sqrt(0.5);
