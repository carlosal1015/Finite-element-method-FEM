function [g] = dirichlet(x,y)
g = 1-((x-0.5)^2+(y-0.5)^2)^0.25;   % dirichlet boundary value in the point (x,y)
g = sin(pi*x)*sin(pi*y);
end

