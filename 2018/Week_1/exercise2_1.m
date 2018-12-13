% Corrida de la Proyeccion L2 de una funcion
% Grable: run_PL2.m
a = -1;
b = 1;
n = 10; % 100, 250
% f= @(x) 2*x.*sin(2*pi*x)+3;
% f= @(x) 0*x+1;
% f= @(x) x.^3.*(x-1).*(1-2*x);
epsilon = 0.01;
f= @(x) abs(x);
%f=@(x) atan((x-0.5)/epsilon);
[x, proy] = PL2(f,a,b,n); % x va los nodos

plot(x, proy, 'LineWidth', 2)
hold on
t = linspace(a, b, 250);
plot(t,f(t), '-r', 'LineWidth', 1)
shg
hold off

% Un script es un conjunto de instrucciones
% intervalo desde 0 hasta 1
% Es conveniente copiar el encabezado.
% LinwWidth da el espesor de la línea.