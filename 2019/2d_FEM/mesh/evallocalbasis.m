function tab=evallocalbasis(ptyp,x,p)
% size(tab)= polydegreee X nodes

switch ptyp
    
    case 98 % triangle, p=1
        n1=size(x,1);
        if p==1
            tab=zeros(3,n1);
            tab(1,:)=1-x(:,1)-x(:,2);
            tab(2,:)=x(:,1);
            tab(3,:)=x(:,2);
        elseif p==2
            tab=zeros(6,n1);
            
            xy=x(:,1).*x(:,2);
            x2=x(:,1).^2;
            y2=x(:,2).^2;
            
            tab(1,:)=1 -3*x(:,1) -3*x(:,2) +2*x2 +4*xy +2*y2;
            tab(2,:)=  -1*x(:,1)           +2*x2            ; 
            tab(3,:)=            -1*x(:,2)             +2*y2; 
            
            tab(4,:)=  +4*x(:,1)           -4*x2 -4*xy      ;            
            tab(5,:)=                            +4*xy      ;             
            tab(6,:)=            +4*x(:,2)       -4*xy -4*y2; 
        end
    otherwise
        error('Unknown ptyp %i in evallocalbasis.',ptyp)
end

end