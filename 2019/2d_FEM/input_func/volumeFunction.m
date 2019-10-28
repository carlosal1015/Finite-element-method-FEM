function val=volumeFunction(x,y,flag)
% laplace u = u_xx + u_yy

switch flag
    case 0
        val=zeros(size(x));
    case 1
        val=zeros(size(x));
    case 2
        val=4*ones(size(x));
    case 3
        val=6*x+6*y;
    case 4
        val=4*(x.^2 - 1).*(y.^2 - 1).^2 + 4*(x.^2 - 1).^2.*(y.^2 - 1) + 8*x.^2.*(y.^2 - 1).^2 + 8*y.^2.*(x.^2 - 1).^2  ;% (12*y.^2 - 4).*(x.^4 - 2*x.^2 + 1) + (12*x.^2 - 4).*(y.^4 - 2*y.^2 + 1);
    case 5
        val=-sin(x.*y).*(x.^2 + y.^2);
        
        
    case 7
        exp_x=exp(x); exp_y=exp(y); y2=y.^2; x2=x.^2;
        
        yexp_y=y.*exp_y; y2exp_y=y2.*exp_y;
        xexp_y=x.*exp_y; x2exp_x=x2.*exp_x; yexp_x=y.*exp_x;
        x2exp_y=x2.*exp_y; xexp_x=x.*exp_x;
        xyexp_y=x.*yexp_y; xy2exp_y=x.*y2exp_y;
        xM1yM1=(x - 1).*(y - 1);
        
        val=x2.*(x - 1).*xM1yM1.*(y2exp_y - 2*exp_y - 2*exp_x - 2*xexp_y + 4*yexp_x + 3*yexp_y + 3*xyexp_y + xy2exp_y) + x.*y2.*(y - 1).^2.*(x2exp_x - 2*exp_y - 2*exp_x + 5*x2exp_y + 3*xexp_x + xexp_y) + y2.*xM1yM1.*(y - 1).*(x2exp_x - 2*exp_y - 2*exp_x + 5*x2exp_y + 3*xexp_x + xexp_y) + x2.*y.*(x - 1).^2.*(y2exp_y - 2.*exp_y - 2.*exp_x - 2.*xexp_y + 4.*yexp_x + 3.*yexp_y + 3.*xyexp_y + xy2exp_y) + x.*y2.*xM1yM1.*(y - 1).*(exp_x + exp_y + x2exp_x + 5.*xexp_x + 10.*xexp_y) + x2.*y.*(x - 1).*xM1yM1.*(4.*exp_x + exp_y + y2exp_y + xexp_y + 5.*yexp_y + 5.*xyexp_y + xy2exp_y);
    case 8 % adaptivity bsp, no solution
        
        % f=-Delta(chi*u0)
        x0=[1;1]/2; alpha=1/2;                
        
        r=sqrt( (x-x0(1)).^2 + (y-x0(2)).^2 );
        tmpy=y.*(1-y);
        tmpx=x.*(1-x);
        valchi=tmpx.*tmpy; % cut-off function
        gradchi1=(1-2*x).*tmpy;
        gradchi2=tmpx.*(1-2*y);
        deltachi=-2*(tmpy+tmpx);                
%         valchi=x.*(1-x).*y.*(1-y); % cut-off function
%         gradchi1=(1-2*x).*y.*(1-y);
%         gradchi2=x.*(1-x).*(1-2*y);
%         deltachi=-2*y.*(1-y)-2*x.*(1-x);        
        
        valu0=r.^alpha;
        tmp=alpha*r.^(alpha-2);
        gradu01=(x-x0(1)).*tmp;
        gradu02=(y-x0(2)).*tmp;
        % deltau0=alpha^2*r.^(alpha-2);
        deltau0=alpha*tmp;
        
        val=(deltachi.*valu0+2*(gradchi1.*gradu01+gradchi2.*gradu02)+valchi.*deltau0);
        val(r==0)=0;       
        
    case 9 % Lshape example
        val=zeros(size(x));    
    case 10     
        val=-2*sin(x).*cos(y);
        
end
val = -val;

end