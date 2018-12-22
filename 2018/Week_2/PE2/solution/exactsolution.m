function [u,dx_u,dy_u] = exactsolution(x,y)
% returns the exact solution u(x,y) and the derivatives
% d_x u(x,y) and d_y(x,y)
u = 1-((x-0.5)^2+(y-0.5)^2)^0.25;
dx_u =0.25*(1-2*x)*((x-0.5)^2+(y-0.5)^2)^(-0.75);
dy_u =0.25*(1-2*y)*((x-0.5)^2+(y-0.5)^2)^(-0.75);

u = sin(pi*x)*sin(pi*y);
dx_u = pi*cos(pi*x)*sin(pi*y);
dy_u = pi*sin(pi*x)*cos(pi*y);
end

