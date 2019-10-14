function y = evalfun(x)
% evalfun calcula el valor $x\mapsto x^2 -1$
% https://www.mathworks.com/matlabcentral/answers/25568-defining-functions
    if ~isnumeric(x)
        error('Debe ingresar un número real.')
    end
y = x^2 - 1;
end