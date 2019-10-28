function val=func_grad_u(x,y,flag)

switch flag
    case 0
        val(:,1)=zeros(length(x),1);
        val(:,2)=zeros(length(x),1);
    case 1
        val(:,1)=ones(length(x),1);
        val(:,2)=ones(length(x),1);
    case 2
        val(:,1)=2*x+y;
        val(:,2)=2*y+x;
    case 3
        val(:,1)=3*x.^2;
        val(:,2)=3*y.^2;
    case 4
        val(:,1)=4*x.*(x.^2 - 1).*(y.^2 - 1).^2; %4*x.*(x.^2 - 1).*(y.^2 - 1).^2;
        val(:,2)=4*y.*(x.^2 - 1).^2.*(y.^2 - 1); %4*y.*(x.^2 - 1).^2.*(y.^2 - 1);
    case 5
        xy=x.*y;
        val(:,1)=y.*cos(xy);
        val(:,2)=x.*cos(xy);
        
        
    case 7
        exp_x=exp(x); exp_y=exp(y); y2=y.^2; x2=x.^2;
        yexp_y=y.*exp_y; y2exp_y=y2.*exp_y; xexp_y=x.*exp_y;
        xM1yM1=(x - 1).*(y - 1);
        
        val(:,1)=x.*y2.*xM1yM1.*(y - 1).*(x2.*exp_x - 2.*exp_y - 2.*exp_x + 5.*x2.*exp_y + 3.*x.*exp_x + xexp_y);
        val(:,2)=x2.*y.*(x - 1).*xM1yM1.*(y2exp_y - 2.*exp_y - 2.*exp_x - 2.*xexp_y + 4.*y.*exp_x + 3.*yexp_y + 3.*x.*yexp_y + x.*y2exp_y);
    case 8 % adaptivity bsp
        
        % f=-Delta(chi*u0)
        x0=[1;1]/2; alpha=1/2;               
        
        r=sqrt( (x-x0(1)).^2 + (y-x0(2)).^2 );
        valchi=x.*(1-x).*y.*(1-y); % cut-off function
        gradchi1=(1-2*x).*y.*(1-y);
        gradchi2=x.*(1-x).*(1-2*y);

        valu0=r.^alpha;
        tmp=alpha*r.^(alpha-2);
        gradu01=(x-x0(1)).*tmp;
        gradu02=(y-x0(2)).*tmp;
        
        val(:,1)=gradchi1.*valu0 + valchi.*gradu01;
        val(:,2)=gradchi2.*valu0 + valchi.*gradu02;
        
        val(r==0,1)=0;        
        val(r==0,2)=0;        

    case 9 % Lshape example
        [theta,rho]= cart2pol(x,y);
        theta(theta<=0)=theta(theta<=0)+2*pi;
        ur = 2/3*rho.^(-1/3).*sin(2/3*theta-pi/3);
        ut = 2/3*rho.^(2/3).*cos(2/3*theta-pi/3);
        
        val(:,1)=ur.*x./rho-ut.*y./rho.^2;
        val(:,2)=ur.*y./rho+ut.*x./rho.^2;        
    case 10
        val(:,1)=cos(x).*cos(y);
        val(:,2)=-sin(x).*sin(y);
      
end

end