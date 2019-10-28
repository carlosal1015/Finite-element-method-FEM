function fig_handle=drawSolution_gradient(element2connectivity,element2vertex,vertex_coor,u,nrRefinements,scaling,p)

% u= u_all
fig_handle=figure;
hold on;

h1=subplot(1,2,1);
h2=subplot(1,2,2);


points=[0 0; 1 0; 0 1];
tri=1:3;

j=0;
for l=1:nrRefinements % number of refinements
    for k=1:4^(l-1)
        j=j+1;
                
        p_neu1=sum(points(tri(j,1:2),:),1)/2;
        p_neu2=sum(points(tri(j,2:3),:),1)/2;
        p_neu3=sum(points(tri(j,[1 3]),:),1)/2;
        
        n=size(points,1);
        points=[points;p_neu1;p_neu2;p_neu3];        
        
        tri(end+1,:)=[tri(j,1) n+1 n+3];
        tri(end+1,:)=[n+1 tri(j,2) n+2];
        tri(end+1,:)=[n+3 n+2 tri(j,3)];
        tri(end+1,:)=[n+1 n+2 n+3];
        
    end
end
tri=tri(j+1:end,:);

% gqn=size(points,1);
basisPrime=evallocalbasis_prime(98,points,p);
basisPrimeX=basisPrime{1};
basisPrimeY=basisPrime{2};

for i=1:size(element2vertex,1)
    p1=vertex_coor(element2vertex(i,1),:)';
    p2=vertex_coor(element2vertex(i,2),:)';
    p3=vertex_coor(element2vertex(i,3),:)';
        
    a  =p1;
    a1 =p2-p1;
    a2 =p3-p1;
    det_jac_Ableitung=a1(1)*a2(2)-a1(2)*a2(1);
    
    X=a(1)+a1(1)*points(:,1)+a2(1)*points(:,2);
    Y=a(2)+a1(2)*points(:,1)+a2(2)*points(:,2);
    
    % eval grad u_h    
%     val_uh_prime_X=zeros(1,gqn); val_uh_prime_Y=zeros(1,gqn);
    val_uh_prime_X=u(element2connectivity(i,:))'*( a2(2)*basisPrimeX-a1(2)*basisPrimeY) ;
    val_uh_prime_Y=u(element2connectivity(i,:))'*(-a2(1)*basisPrimeX+a1(1)*basisPrimeY);
    val_uh_prime_X=val_uh_prime_X/det_jac_Ableitung;
    val_uh_prime_Y=val_uh_prime_Y/det_jac_Ableitung;

%     for j=1:3*p
%         G=element(i).connectivity(j);
%         val_uh_prime_X=val_uh_prime_X+u(G)*( a2(2)*basisPrimeX(j,:)-a1(2)*basisPrimeY(j,:)) ;
%         val_uh_prime_Y=val_uh_prime_Y+u(G)*(-a2(1)*basisPrimeX(j,:)+a1(1)*basisPrimeY(j,:));
%     end
%     val_uh_prime_X=val_uh_prime_X/det_jac_Ableitung;
%     val_uh_prime_Y=val_uh_prime_Y/det_jac_Ableitung;
    
    
    subplot(h1)
    h3=trisurf(tri,X,Y,val_uh_prime_X);
    set(h3, 'edgecolor','none');
    hold all
    
    subplot(h2)
    h4=trisurf(tri,X,Y,val_uh_prime_Y);
    set(h4, 'edgecolor','none');
    hold all
    
end

subplot(h1)
xlabel('x')
ylabel('y')
view([-16 36])

subplot(h2)
xlabel('x')
ylabel('y')
view([-16 36])
% axis([-1.1*scaling 1.1*scaling -1.1*scaling 1.1*scaling -0.3 2.2])
hold off;


end