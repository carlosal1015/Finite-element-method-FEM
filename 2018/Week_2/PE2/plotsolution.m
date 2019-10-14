function plotsolution(name,mesh,x)
%PLOTSOLUTION Summary of this function goes here
%   Detailed explanation goes here
figure('Name',name);       % create a matlab figure

[N,DIM] = size(mesh);      % mesh nodes and dimension
M=N^(1/DIM);               % nodes in every direction

h=1/(M-1);                 % grid size

if (DIM==2)
    MX = reshape(mesh(:,1),M,M);
    MY = reshape(mesh(:,2),M,M);
    Z  = reshape(x,M,M);
    surf(MX,MY,Z);
elseif (DIM==3)
    [MX,MY,MZ]=meshgrid([0:h:1]);
    Z=zeros(N+1,N+1,N+1);
    Z(2:N,2:N,2:N)=reshape(x,N-1,N-1,N-1);
    xslice=[0.3,0.7];  yslice=0.5; zslice=0.3;
    slice(MX,MY,MZ,Z,xslice,yslice,zslice);
end
end
