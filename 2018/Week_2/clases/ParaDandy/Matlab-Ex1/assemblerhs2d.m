function [rhs] = assemblerhs(mesh)

[N,DIM] = size(mesh);     % nr. of nodes and dimension

assert(DIM==2);           % make sure that we are in 2d

N1d=N^(1/DIM);            % points in every direction
M = N1d-1;                % elements in every direction
h=1/M;                    % grid size

rhs=zeros(N,1);           % vector of zero's to carry right hand side

% Gauss Points & Guass Weights
gp = h*[0.5-0.5/sqrt(3.0),0.5-0.5/sqrt(3.0);
        0.5+0.5/sqrt(3.0),0.5-0.5/sqrt(3.0);
        0.5-0.5/sqrt(3.0),0.5+0.5/sqrt(3.0);
        0.5+0.5/sqrt(3.0),0.5+0.5/sqrt(3.0)];
gw = h*h*[0.25,0.25,0.25,0.25];

% The for test functions
PHI = [(1-gp(:,1)/h).*(1-gp(:,2)/h) ,(gp(:,1)/h).*(1-gp(:,2)/h),(1-gp(:,1)/h).*(  gp(:,2)/h) ,(  gp(:,1)/h).*(  gp(:,2)/h)];


% loop over all mesh elements in y- and x-direction
% numerical quadrature with a 2x2-point Gauss rule
for my=1:M                  % row of the mesh element
    for mx=1:M              % column of the mesh element
        
        row = (my-1)*N1d + mx;  % index in lower/left corner
        x=mesh(row,1);     % x-coordinate of this point
        y=mesh(row,2);     % y-coordinate of this point
        
        % Assemble local right hand side vector
        %   for the 4 test functions in the element
        bloc=zeros(4,1);  
        for j=1:4         % loop over the four test functions of the element
            for k=1:4         % loop over the gauss points
                X=x+gp(k,1);
                Y=y+gp(k,2);
                bloc(j) = bloc(j)+gw(k)*righthandside(X,Y)*PHI(j,k); % <<<--- here's the rhs
            end
        end
        
        % enter local vector to global vector
        rhs(row)       = rhs(row)       + bloc(1);
        rhs(row+1)     = rhs(row+1)     + bloc(2);
        rhs(row+N1d)   = rhs(row+N1d)   + bloc(3);
        rhs(row+N1d+1) = rhs(row+N1d+1) + bloc(4);
    end
end
    
% set Dirichlet values
for m=1:N1d
    % 
    rhs(m) = 0;     % <<<--- here's the Dirichlet values
    rhs((m-1)*N1d+1) = 0;
    rhs(m*N1d) = 0;
    rhs((N1d-1)*N1d+m) = 0;
end

end

