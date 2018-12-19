function [g] = dirichlet(x,y)
% Ejemplo inicial.
%g=0;

% Problema 2
%g = sin(4.5* pi * x) * sin(2.7* pi * y);
% No colocar multiplos de 2pi porque se anula.

%    if x < y
%        g = 1
%    else
%        g = -1
%    end

%g = x ^ 2 + y ^2;
g = ( (x - 0.5) ^ 2 + (y - 0.5) ^ 2 ) ^ (0.25);

end