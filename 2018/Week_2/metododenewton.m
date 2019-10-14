% Resolver la ecuación dada por 2cos x - exp(x) = 0
x = 1;
for n = 0: 40
    x = x - (2*cos(x) - exp(x) )/(-2*sin(x)- exp(x));
end
x
% Comparar con la función solve de matlab.