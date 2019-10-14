function [x, proy] = PL2(f, a, b, n)
M = zeros(n+1, n+1);
c = zeros(n+1, 1);
h = (b-a)/n; % Diferencia de los intervalos
x = linspace(a,b,n+1); % Habra n+1 puntos

for k = 1:n
    M_local = h*[1/3 1/6; 1/6 1/3];
    b_local = (h/2) * [f(x(k)); f(x(k+1))];
    M(k:k+1, k:k+1) = M_local + M(k:k+1, k:k+1);
    c(k:k+1) = c(k:k+1) + b_local;
end

proy = M\c; % Resolviendo el sistema lineal
% display(M)
display(c)
end
% Crear con las otras cuadraturas, trapecia, simpson, revisar whitakker.