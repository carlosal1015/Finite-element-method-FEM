function rhs = righthandside(x,y)
% evalates the right hand side f(x,y)
%rhs = 0.25*((x-0.5)^2+(y-0.5)^2)^(-0.75);
rhs = sin(pi*x)*sin(pi*y)*2*pi*pi;
%rhs = 1;

end

