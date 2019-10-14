function grafica_hat(j,a,b,n,s)
% n = numero de segmentos
% 
x = linspace(a,b,n+1);
if j == 1
    xx = [x(1), x(2), x(n+1)];
    yy = [1, 0, 0];
else
    if j == n
        xx = [x(1), x(2), x(n+1)];
        yy = [0, 0, 1];
    else % i<j<n
        xx = [x(1), x(j-1), x(j), x(j+1), x(n+1)];
        yy = [0,0,1,0,0];
    end
end
plot(xx,yy,s)
shg
end