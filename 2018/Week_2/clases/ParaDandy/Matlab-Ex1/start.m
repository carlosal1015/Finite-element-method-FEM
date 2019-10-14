function start()
% M is the number of mesh elements in every direction
% N=(M+1)*(M+1) is the total number of points 
M = 3;

% create the mesh
mesh = createmesh(2,M);

% assemble the matrix
A = assemblematrix2d(mesh);

% set the right hand side
b = assemblerhs2d(mesh);

% solve the linear system
x = A\b;

% show the solution
plotsolution('solution',mesh,x);
end

