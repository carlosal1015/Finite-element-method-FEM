% Crear un programa en MATLAB que pida ingresar la cantidad de términos de
% la expansión de Taylor de e, imprima el error absoluto y el error
% relativo.
z = 0;
prompt = 'Ingrese el valor de n:';
n = input(prompt);
for k = 0:n
  z = z + 1/factorial(k);
end
fprintf("El valor de z es");
z
fprintf("El error relativo es")
error_absolute = abs(exp(1) - z)
fprintf("El error absoluto es")
error_relative = error_absolute/ exp(1)