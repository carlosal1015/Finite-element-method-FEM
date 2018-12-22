function [u,dx_u,dy_u] = exactsolution(x,y)
% returns the exact solution u(x,y) and the derivatives
% d_x u(x,y) and d_y(x,y)
u = sin(4.5*pi*x)*sin(2.5*pi*y);
dx_u = 4.5*pi*cos(4.5*pi*x)*sin(2.5*pi*y);
dy_u = 2.5*pi*sin(4.5*pi*x)*cos(2.5*pi*y);
end

