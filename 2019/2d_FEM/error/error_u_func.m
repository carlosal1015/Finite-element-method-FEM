function [L2Error,H1Error]=error_u_func(element2connectivity,element2vertex,vertex_coor,u,func_u,func_grad_u,gqnf,p)

[gqk, gqw]=dunavant_triangle_quad(p+gqnf);

basis=evallocalbasis(98,gqk,p);
basisPrime=evallocalbasis_prime(98,gqk,p);
basisPrimeX=basisPrime{1};
basisPrimeY=basisPrime{2};

L2Error=0;
H1Error=0;

for i=1:size(element2vertex,1)
    p1=vertex_coor(element2vertex(i,1),:)';
    p2=vertex_coor(element2vertex(i,2),:)';
    p3=vertex_coor(element2vertex(i,3),:)';
    
    a  =p1;
    a1 =p2-p1;
    a2 =p3-p1;
    det_jac_Ableitung=a1(1)*a2(2)-a1(2)*a2(1);
    det_jac=abs(det_jac_Ableitung);
    %     det_jac=abs(a1(1)*a2(2)-a1(2)*a2(1)); % verzerrungsfaktor
        
    val_u= func_u(a(1)+a1(1)*gqk(:,1)+a2(1)*gqk(:,2),a(2)+a1(2)*gqk(:,1)+a2(2)*gqk(:,2));
    val_grad_u= func_grad_u(a(1)+a1(1)*gqk(:,1)+a2(1)*gqk(:,2),a(2)+a1(2)*gqk(:,1)+a2(2)*gqk(:,2));    
    
    tmp=u(element2connectivity(i,:))';
    val_uh=tmp*basis;
    val_uh_prime_X=tmp*( a2(2)*basisPrimeX-a1(2)*basisPrimeY) ;
    val_uh_prime_Y=tmp*(-a2(1)*basisPrimeX+a1(1)*basisPrimeY);
    val_uh_prime_X=val_uh_prime_X/det_jac_Ableitung;
    val_uh_prime_Y=val_uh_prime_Y/det_jac_Ableitung;
    
    L2Error=L2Error+ ( (val_u'-val_uh).^2 )*gqw*det_jac/2; % /2 wegen der dunavant_triangle_quad
    H1Error=H1Error+ ( (val_grad_u(:,1)'-val_uh_prime_X).^2 + (val_grad_u(:,2)'-val_uh_prime_Y).^2 )*gqw*det_jac/2; % /2 wegen der dunavant_triangle_quad
end

H1Error=sqrt(H1Error + L2Error);
L2Error=sqrt(L2Error);
end