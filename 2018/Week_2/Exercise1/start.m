function start()
% M is the number of mesh elements in every direction
% N=(M+1)*(M+1) is the total number of points 
M = 400; %M = 20 Increase the number of mesh

% el numero de nodos es M+1

% create the mesh
mesh = createmesh(2,M);
%mesh
%pause

% assemble the matrix
A = assemblematrix2d(mesh);

% set the right hand side
b = assemblerhs2d(mesh);

% solve the linear system
x = A\b; % u = A\b, u es la solución, b es rhs.

% show the solution
%x = zeros((M+1) * (M+1), 1);
%x(36,1) = 1;
%plotsolution('solution',mesh,x);
plotsolution('rhs',mesh,x); % Visualizar la función
end
% righthandside.m
