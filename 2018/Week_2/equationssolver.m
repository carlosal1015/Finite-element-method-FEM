% Crear su gráfica de puntos de n datos vs el E(n), el error después de n
% iteraciones.

% Calcule el máximo error absoluto que se puede cometer al redondear una
% cantidad a n cifras decimales, considerando que si la siguiente cifra
% después de la n-ésima es mayor o igual a 5 el redondeo es hacia arriba.

%    $e^x - x = 2$
% Calcule la cota máxima del error que se comente en la última
% aproximación.
% Cómo hallar fórmulas de recurrencia.
x =1; z = 0;
prompt = 'Ingrese el valor de n:';
n = input(prompt);
for k = 0:n
    z = z + 1/factorial(k);
    for i = 0:10
        x = x - (z - x - 2)/(z -1);
    end
end
z
x
fprintf("El método de aproximaciones sucesivas no es eficiente.\n Hay que usar el método de paso.")