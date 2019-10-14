function [error] = h1error(mesh,u)

[N,DIM] = size(mesh);     % nr. of nodes and dimension

assert(DIM==2);           % make sure that we are in 2d

N1d=N^(1/DIM);            % points in every direction
M = N1d-1;                % elements in every direction
h=1/M;                    % grid size

% Gauss Points & Guass Weights
gp = h*[0.5-0.5/sqrt(3.0),0.5-0.5/sqrt(3.0);
        0.5+0.5/sqrt(3.0),0.5-0.5/sqrt(3.0);
        0.5-0.5/sqrt(3.0),0.5+0.5/sqrt(3.0);
        0.5+0.5/sqrt(3.0),0.5+0.5/sqrt(3.0)];
gw = h*h*[0.25,0.25,0.25,0.25];

% The for test functions
PHI_X = [(-1/h).*(1-gp(:,2)/h) ,(1/h).*(1-gp(:,2)/h),(-1/h).*(  gp(:,2)/h) ,(1/h).*(  gp(:,2)/h)];
PHI_Y = [(1-gp(:,1)/h).*(-1/h) ,(gp(:,1)/h).*(-1/h),(1-gp(:,1)/h).*(1/h) ,(  gp(:,1)/h).*(1/h)];

error = 0;

% loop over all mesh elements in y- and x-direction
% numerical quadrature with a 2x2-point Gauss rule
for my=1:M                  % row of the mesh element
    for mx=1:M              % column of the mesh element
        
        row = (my-1)*N1d + mx;  % index in lower/left corner
        x0=mesh(row,1);     % x-coordinate of this point
        y0=mesh(row,2);     % y-coordinate of this point

        % Finite element solution on this element
        UH = [u(row),u(row+1),u(row+N1d),u(row+N1d+1)];
        
        % Assemble local right hand side vector
        %   for the 4 test functions in the element
        
        for k=1:4         % loop over the gauss points
            x=x0+gp(k,1);
            y=y0+gp(k,2);
            [uex,dx_uex,dy_uex] = exactsolution(x,y);
            uh_x = PHI_X(k,:)*UH';
            uh_y = PHI_Y(k,:)*UH';
            error = error + gw(k)*( (dx_uex-uh_x)^2 + (dy_uex-uh_y)^2 );
        end
    end
end

error = sqrt(error);
end

