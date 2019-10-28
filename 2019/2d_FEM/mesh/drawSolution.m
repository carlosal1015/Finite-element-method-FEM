function fig_handle=drawSolution(element2connectivity,element2vertex,vertex_coor,u,nrRefinements,p)
% u= u_all
fig_handle=figure; 
hold on;

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

basis=evallocalbasis(98,points,p);

for i=1:size(element2vertex,1)
    p1=vertex_coor(element2vertex(i,1),:)';
    p2=vertex_coor(element2vertex(i,2),:)';
    p3=vertex_coor(element2vertex(i,3),:)';
    
    a  =p1;
    a1 =p2-p1;
    a2 =p3-p1;
    
    X=a(1)+a1(1)*points(:,1)+a2(1)*points(:,2);
    Y=a(2)+a1(2)*points(:,1)+a2(2)*points(:,2);
    
    val = u(element2connectivity(i,:))'*basis;
       
    h=trisurf(tri,X,Y,val);
    set(h, 'edgecolor','none');
    
end
view([-16 36])

hold off;
end


