% Grabe con: run_poisson1D.m
clear all; clc
%%%%%%%%%%%%%%%%%%%%%% Dirichlet Homogeneo
% a = 0; b = 1; % ua = 0; ub = 0;
% n = 10;
% f=@(x) pi^2*sin(pi*x);
% u=@(x) sin(pi*x);
%%%%%%%%%%%%%%% Dirichlet no homogéno
a = 1; b = 2;
ua = -1; ub = 4;
n = 10;
f=@(x) - 2*cos(pi*x)+pi*pi*x.^2*cos(pi*x)+4*pi*x.*sin(pi*x);
u=@(x) x.^2.*cos(pi*x);
[x, poi] = poisson1D(f, a, b, n, ua, ub);
plot(x, poi, 'LineWidth', 2)
hold on
t = linspace(a, b, n+1)';
plot(t, u(t), '-r', 'LineWidth',1)
hold off
shg
disp()