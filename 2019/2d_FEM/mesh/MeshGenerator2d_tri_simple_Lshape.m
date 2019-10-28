function [element2vertex,vertex_coor,vertex_bc]=MeshGenerator2d_tri_simple_Lshape(n,scaling,bc_type)
% bctype=1 is Dirichlet, =0 inner, =2 is Neumann

h=1/n;
vertex_coor=zeros((2*n+1)*(n+1)+(n+1)*n,2); %alle knoten allokiiert
vertex_bc=zeros((2*n+1)*(n+1)+(n+1)*n,1); 

% set vertex.coor
for i=0:n %�ber alle knoten in y-richtung
    y=-1+i*h;
    for j=0:2*n
        x=-1+j*h;
        cnt=i*(2*n+1)+j+1;
        vertex_coor(cnt,1:2)=[x;y]*scaling;
        if (i==n && j>=n)
            vertex_bc(cnt)=1; % dirichlet node
        elseif i==0 || j==0 ||j==2*n
            vertex_bc(cnt)=bc_type; % Neumann node       
        else
            vertex_bc(cnt)=0; % inner node
        end
    end
end

for i=1:n %�ber alle knoten in y-richtung
    y=0+i*h;
    for j=0:n %%�ber alle knoten in x-richtung
        x=-1+j*h;
        cnt=(2*n+1)*(n+1)+(i-1)*(n+1)+j+1;
        vertex_coor(cnt,1:2)=[x;y]*scaling;
        if j==n
            vertex_bc(cnt)=1; % dirichlet node
        elseif i==n || j==0
            vertex_bc(cnt)=bc_type; % Neumann node  
        else
            vertex_bc(cnt)=0; % inner node
        end
    end
end



element2vertex=zeros(2*3*n^2,3); %alle elemente allokiiert

for i=1:n %�ber all elemente in y-richtung
    a=(i-1)*(2*n+1);
    f=(i-1)*(2*n+2*n+1);
    for j=1:2*n % �ber alle elemente in x-richtung
        a=a+1;
        f=f+1;
        cnt=(i-1)*2*n+j;
        element2vertex(2*cnt-1,1:3)=[ a a+1 a+2*n+2];
        element2vertex(2*cnt,1:3)=[ a a+2*n+1 a+2*n+2];
    end
end

for i=1:n %�ber all elemente in y-richtung
    if i==1
        f=2*n*n+(2*n+1)*n;
        a=(2*n+1)*n;
    else
        f=2*n*n+(2*n+1)*n + 2*n+n+1 + (i-2)*(n+n+1);
        a=(2*n+1)*(n+1) + (i-2)*(n+1);
    end
    
    for j=1:n % �ber alle elemente in x-richtung
        a=a+1;
        f=f+1;
        cnt= 2*n*n + (i-1)*n+j;
        
        if i==1
            element2vertex(2*cnt-1,1:3)=[ a a+1  a+2*n+2];
            element2vertex(2*cnt,1:3)=[ a a+2*n+1 a+2*n+2];
        else
            element2vertex(2*cnt-1,1:3)=[ a a+1  a+n+2];
            element2vertex(2*cnt,1:3)=[ a a+n+1 a+n+2];
        end
        
    end
end


end