b = [8, -7, 6, 3, 4]
% b()
% Parte a

% Parte b

for j=1:5
    j^2
end

i = 1:0.33:4

n = 100;
s = 0; % inicializa
for k = 1:n
    s = s+1/(k^2+1)
end
%calibri, una biblioteca
% Depurador en Matlab
comb(3,4)% k = 80, n = 100
% factorial de 100 no lo calcula.
home
k = 88; n = 100
comb(k,n)

% Para código largo se escribe en el editor
% En la consola se escriba código corto.
% En el function tiene una entrada y una salida
% En el editor no tiene entrada.

% 0 es falso, 1 es verdadero, distinto de 1 y no cero también es verdadero.

4 < 6
3 > 8
4 == 4

a = 1, b = 5, n = 10, s='g', %j =1
% Vamos a graficar phi_1
grafica_hat(1,a,b,n,s)% phi_1
hold on
grafica_hat(2,a,b,n,s)% phi_2

figure
hold on
pause(5)
% Tarea: Ponerle los colores.
for j = 1:n
    grafica_hat(j,a,b,n,s)% phi_1,phi_2, phi_3
end