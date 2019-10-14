function [A] = assemblematrix2d(mesh)
%Assemble the Finite Element Matrix in 2D
%The matris given by its stencil, a 3x3 matrix
%
[N,DIM] = size(mesh);
assert(DIM==2);           % check that we are in 2d

M = N^(1/DIM);            % number of points in every direction
h = 1/(M-1);              % element size

S = [ -1/3, -1/3, -1/3 ; -1/3,  8/3, -1/3 ;  -1/3, -1/3, -1/3];


% First, we create an empty sparse matrix
A=sparse(N,N);

% Write matrix stencil into Matrix A, row by row

for my=1:M                  % row of the mesh
    for mx=1:M              % column of the mesh 
        
        % index of the row in the matrix belonging to meshpoint my,mx
        rowA = (my-1)*M + mx;
        
        % enter the diagonal
        A(rowA,rowA) = S(2,2);
        % enter left and right - if not mx=1 or mx=M
        if (mx>1) A(rowA,rowA-1) = S(2,1); end
        if (mx<M) A(rowA,rowA+1) = S(2,3); end
        % enter up and down - if not my=1 or my=M
        if (my>1) A(rowA,rowA-M) = S(1,2); end
        if (my<M) A(rowA,rowA+M) = S(3,2); end
        % enter diagonals - if not ...
        if (mx>1) && (my>1) A(rowA,rowA-1-M) = S(1,1); end
        if (mx<M) && (my>1) A(rowA,rowA+1-M) = S(1,3); end
        if (mx>1) && (my<M) A(rowA,rowA-1+M) = S(3,1); end
        if (mx<M) && (my<M) A(rowA,rowA+1+M) = S(3,3); end
    end
end

%%% set dirichlet values
for m=1:M        %
    % bottom boundary
    A(m,:) = 0;
    A(m,m) = 1;
    % left boundary
    A((m-1)*M+1,:) = 0;
    A((m-1)*M+1,(m-1)*M+1) = 1;
    % top boundary
    A((M-1)*M+m,:) = 0;
    A((M-1)*M+m,(M-1)*M+m) = 1;
    % right boundary
    A(M*m,:) = 0;
    A(M*m,M*m) = 1;
end
    
end