a = [2,3,4,7,4,3]
format compact
z = a(1)
v =a(3:end) % a partir del 3 va a ser válido
a = [2, a, [9,1]]
d= max(a)
length(a)
c = [2,3; 4,5;7,1]
[nuf,numc] = size(c)

% Cuando hacemos las réplicas el cálculo se puede hacer lento.
% Matlab es un compilador, no es un traductor, no ejecuta línea por línea,
% C es más veloz.
% Hay caja de herramientas para pasar al lenguaje C
% Primero se trabajo a C.
% Vamos a crear un conjunto de índices
c([1,3], 2)
c([3,1], 2)
A = [2 6 7 9
    3 7 2 5
    4 -2 1 6
    1 9 8 2];
% También un una sola línea

b = [2 -1 3 4]';
% Matlab trata de emplear el mejor método que tiene
% En matlab todas las variables son dinámicas
% Ax=b
% x = A\b

A x= b \implies A^T\cdot Ax= A^T \cdot b^\ast
A*x-b
% el asterisco es el error del cálculo

norm(A*x-b) % debe dar un valor muy cercano a cero.
% en matlab no es eficiente calcular la inversa
% el metodo de choleski se recomienda si la matriz es definida positiva
% el meodo de gauss jordan es el metodo directo.
% metodo de householder es util cuando la matriz
%metodo de gibbens radio y otros medtoodo iterativos como el metodo de jacobi, metodo sor
% Buscar los métodos
cond(A)

%> help("cond")
% cond   Condition number with respect to inversion.
%    cond(X) returns the 2-norm condition number (the ratio of the
%    largest singular value of X to the smallest).  Large condition
%    numbers indicate a nearly singular matrix.
% 
%    cond(X,P) returns the condition number of X in P-norm:
% 
%       NORM(X,P) * NORM(INV(X),P). 
% 
%    where P = 1, 2, inf, or 'fro'. 
% 
%    Class support for input X:
%       float: double, single
% 
%    See also rcond, condest, condeig, norm, normest.
%
%    Reference page for cond
%    Other functions named cond
% La mtriz de hilbert es condicional
%x = linspace()


t = 0:0,3:2
t_1 = linspace(0,0.3,5)
t_2 = 0:00.1:2;
y = t_1.*sin(t_1);
plot(x,y)

format short

t = 0:0.1:6;

plot(t,t.*sin(t),t,t.^2)

t = 0:0.1:2
y1 = t.*sin(t);
plot(t,y1,)
shg % Para qu'e sirve
%help("shg")
% shg    Show graph window.
%    shg brings the current figure window forward.
x = 3:0.1:5;
y = x.*exp(x);

plot(x,y,'r', 'LineWidth',3)
shg
xlabel('Horas')
% Los apostrofes en Matlab signifca una cadena en matlab, lo que sigue es
% el espeso en la línea.
ylabel('Temperatura(°C)')
title('Figura ...')

axis equal% Para desactivar usar axis off
shg
x_2 = 2:0.1:3;
plot(x,x,':ok',x,x.^2,'--xr',x,x.^3,'-b')
legend('y=x','y=x^{2}','y=x^{3}')