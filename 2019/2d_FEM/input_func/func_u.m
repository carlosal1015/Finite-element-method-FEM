function val=func_u(x,y,flag)

switch flag
    case 0
        val=ones(size(x));
    case 1
        val=x+y;
    case 2
        val=x.^2+y.^2+x.*y+1;
    case 3
        val=x.^3+y.^3;   
    case 4
        val=(1-2*x.^2+x.^4).*(1-2*y.^2+y.^4); % zero grad and zero func value on boundary
    case 5
        val=sin(x.*y);
        
        
    case 7
        val=(exp(x)+(x+1).*exp(y)).*x.^2.*y.^2.*(1-x).^2.*(1-y).^2;
    case 8 % adaptivity bsp
%         u = chi*u0
        
        x0=[1;1]/2; alpha=1/2;        
        
        r=sqrt( (x-x0(1)).^2 + (y-x0(2)).^2 );
        valchi=x.*(1-x).*y.*(1-y); % cut-off function        
        valu0=r.^alpha;
        val = valchi.*valu0;

    case 9 % Lshape example
        [theta,rho]= cart2pol(x,y);
        theta(theta<=0)=theta(theta<=0)+2*pi; 
        
        val=rho.^(2/3).*sin(2/3*theta-pi/3);
    case 10
        val=sin(x).*cos(y);
        
end


end