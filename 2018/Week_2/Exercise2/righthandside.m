function rhs = righthandside(x,y)
% evalates the right hand side f(x,y)
%rhs = 1;
%rhs = -4;
% evalates the right hand side f(x,y)
%  rhs=1
%[x, y] = meshgrid(-1:.1:1); 
%k = 5;
%u = sin(k * pi * x) * sin(k * pi * y);
%rhs = 2 * k * k * pi * pi * u;

%if (1/3 < x < 2/3 && 1/3 < y < 2/3)
%    rhs =-();
%else rhs = 0;

rhs = -0.25 * ( (x - 0.5) ^ 2 + (y - 0.5) ^ 2 ) ^ (-0.75);

end