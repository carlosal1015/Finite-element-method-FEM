function tab=evallocalbasis_prime(ptyp,x,p)
% size(tab)= polydegreee X nodes

switch ptyp
    case 98 %  triangle, p=1
        n1=size(x,1);
        if p==1
            tabx=zeros(3,n1);
            tabx(1,:)=-1;
            tabx(2,:)=1;
            tabx(3,:)=0;
            
            taby=zeros(3,n1);
            taby(1,:)=-1;
            taby(2,:)=0;
            taby(3,:)=1;
        elseif p==2
            tabx=zeros(6,n1);
            taby=zeros(6,n1);
             
            tabx(1,:)= -3  +4*x(:,1) +4*x(:,2) ;
            tabx(2,:)=  -1          +4*x(:,1)            ; 
            tabx(3,:)=  0; 
            tabx(4,:)=  +4           -8*x(:,1) -4*x(:,2)      ;            
            tabx(5,:)=                            +4*x(:,2)      ;             
            tabx(6,:)=                   -4*x(:,2)  ; 
            
            
            taby(1,:)= -3   +4*x(:,1) +4*x(:,2);
            taby(2,:)=            0          ; 
            taby(3,:)=            -1             +4*x(:,2); 
            taby(4,:)=             -4*x(:,1)     ;            
            taby(5,:)=                            +4*x(:,1)      ;             
            taby(6,:)=            +4       -4*x(:,1) -8*x(:,2);             
        end
        tab{1}=tabx;
        tab{2}=taby;
        
    otherwise
        error('Unknown ptyp %i in evallocalbasis_prime.',ptyp)
end

end