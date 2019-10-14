function [g] = dirichlet(x,y)
%g=0;
    if (x + y) < 1/2
        g = 1;
    else
        g = 0;
    end
end