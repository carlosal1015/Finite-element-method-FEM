function [x, poi] = poisson1D(f, a, b, n, ua, ub)
M = zeros(n+1,n+1); c = zeros(n+1,1);
h = (b-a)/n; % Diferencia de los intervalos
x = linspace(a,b,n+1); % Habra n+1 puntos

for k = 1:n
    M_local = (1/h) *[1 -1; -1 1];
    b_local = (h/2) * [f(x(k)); f(x(k+1))];
    M(k:k+1, k:k+1) = M_local + M(k:k+1, k:k+1);
    c(k:k+1,1) = c(k:k+1,1) + b_local;
end
M(1,:) = 0; M(1,1) = 1; M(n+1,:) = 0; M(n+1,n+1) = 1;
c(1) = ua; c(n+1) = ub;
poi = M\C;
end