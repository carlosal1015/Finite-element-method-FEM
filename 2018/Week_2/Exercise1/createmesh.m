function [mesh] = createmesh(DIM,M)
%CREATEMESH creates a regular finite element mesh in DIM dimensions
%   M is the number of elements in every mesh direction
%   therefore, h=1/M is the mesh size
%   
%   The mesh has (M+1) points in every line
%     0 on the boundary, 1,2,...,M-1 innerpoint and M on the boundary
%   The total number of points is N = (M+1)^DIM
%
%   The mesh is a N x DIM array.
%   The coordinate of point i is: mesh(i,0), mesh(i,1)            (in 2d)
%                            and: mesh(i,0), mesh(i,1), mesh(i,2) (in 3d)

h=1/M;                % compute mesh size
N=(M+1)^DIM;          % number of points
mesh=zeros(N,DIM);
if DIM==1
    mesh = h*[0:1:M];
elseif DIM==2
    mesh(:,1) = repmat(h*[0:1:M]',M+1,1);
    mesh(:,2) = repelem(h*[0:1:M]',M+1,1); % h*kron([0:1:M]', ones(M+1,1));
elseif DIM==3
    mesh(:,1) = repmat(h*[0:1:M]',(M+1)*(M+1),1);
    mesh(:,2) = repmat(repelem(h*[0:1:M]',(M+1),1),(M+1),1);
    mesh(:,3) = repelem(h*[0:1:M]',(M+1)*(M+1),1); % h*kron([0:1:M]', ones(M+1,1));
end



