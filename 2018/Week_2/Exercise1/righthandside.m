function rhs = righthandside(x,y)
% evalates the right hand side f(x,y)
%  rhs=1
%[x, y] = meshgrid(-1:.1:1); 
k = 5;
u = sin(k * pi * x) * sin(k * pi * y);
rhs = 2 * k * k * pi * pi * u;

%if (1/3 < x < 2/3 && 1/3 < y < 2/3)
%    rhs =-();
%else rhs = 0;
end


%if [x >= 1/3 & x <= 2/3, y >= 1/3 & y <= 2/3]
%    rhs = 1
%else
%    rhs = 0
%end

% Definir la función f(x,y)=
% -1, 1/3 < x,y < 2/3
% 0, caso contrario