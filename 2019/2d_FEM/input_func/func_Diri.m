function val=func_Diri(x,y,flag)

switch flag
    case 9 % Lshape example
        [theta,rho]= cart2pol(x,y);
        theta(theta<=0)=theta(theta<=0)+2*pi;
        
        if all(x==0) || all(y==0)
            val=zeros(size(x));
        else
            val=rho.^(2/3).*sin(2/3*theta-pi/3);
        end
    otherwise
        val=func_u(x,y,flag);
        
end


end