function [element2vertex,element2face,face2vertex,face_bc,face2neighbor,face_normal,vertex_coor,vertex_bc]=MeshGenerator2d_tri(n,scaling,bc_type)
% bctype=1 is Dirichlet, =0 inner, =2 is Neumann


h=2/n;

% set vertex.coor
vertex_coor=zeros((n+1)^2,2); %alle knoten allokiiert
vertex_bc=zeros((n+1)^2,1);

cnt=1;
for i=0:n %ueber alle knoten in y-richtung
    y=-1+i*h;
    for j=0:n %ueber alle knoten in x-richtung
        x=-1+j*h;
        vertex_coor(cnt,1:2)=[x;y]*scaling;
        if i==0                       % south edge
            vertex_bc(cnt)=1;     % dirichlet node
        elseif i==n || j==0 || j==n   % west, north, east edges are Neumann boundaries
            vertex_bc(cnt)=bc_type;     % Neumann node            
        else
            vertex_bc(cnt)=0;     % inner node
        end
        cnt=cnt+1;
    end
end

% set element
element2vertex = zeros(2*n^2,3);
element2face   = zeros(2*n^2,3);

cnt=1;
face_matrix=zeros(5*n^2,2); cnt_face=0; % face_matrix=[];
for i=1:n %ueber all square elemente in y-richtung
    a=(i-1)*(n+1);
    for j=1:n % ueber alle square elemente in x-richtung
        a=a+1;
        
        face1_vertex=[a a+1];       % unten
        face2_vertex=[a a+n+1];     % links
        face3_vertex=[a+1 a+n+2];   % rechts
        face4_vertex=[a+n+1 a+n+2]; % oben
        face5_vertex=[a a+n+2];     % diagonal
                
        lower_element_neigh=cnt-2*n+1;
        if lower_element_neigh>0
            face_nr1=element2face(lower_element_neigh,2);
        else
            face_nr1=cnt_face+1; %size(face_matrix,1)+1;
            face_matrix(face_nr1,1:2)=face1_vertex;
            cnt_face=cnt_face+1;
        end
        
        
        left_element_neigh=cnt-2;
        if j>1
            face_nr2=element2face(left_element_neigh,2);
        else
            face_nr2=cnt_face+1; cnt_face=cnt_face+1;
            face_matrix(face_nr2,1:2)=face2_vertex;
        end

        
        % immer neu
        face_nr3=cnt_face+1; %size(face_matrix,1)+1;
        face_nr4=face_nr3+1;
        face_nr5=face_nr4+1;
        
        face_matrix(face_nr3:face_nr5,1:2)=[face3_vertex;face4_vertex;face5_vertex];
        cnt_face=cnt_face+3;
 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        element_vertex          = [a a+1 a+n+2];
        element2vertex(cnt,1:3) = element_vertex;
        element2face(cnt,1:3)   = [face_nr1 face_nr3 face_nr5];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        element_vertex            = [a a+n+1 a+n+2];
        element2vertex(cnt+1,1:3) = element_vertex;
        element2face(cnt+1,1:3)   = [face_nr2 face_nr4 face_nr5];
                
        cnt=cnt+2;
    end
end

% set face
% face2vertex = zeros(cnt_face,2);
face2neighbor = zeros(cnt_face,2);
face_bc = zeros(cnt_face,1);
face_normal = zeros(cnt_face,2);
face2vertex=face_matrix(1:cnt_face,1:2);

for i=1:cnt_face
    
    if vertex_bc(face_matrix(i,1))==1 && vertex_bc(face_matrix(i,2))==1 ...
            && ( vertex_coor(face_matrix(i,1),1)==vertex_coor(face_matrix(i,2),1) ||  vertex_coor(face_matrix(i,1),2)==vertex_coor(face_matrix(i,2),2) )
        face_bc(i)=1;
        
        if vertex_coor(face_matrix(i,1),2)==-scaling && vertex_coor(face_matrix(i,2),2)==-scaling
            face_normal(i,1:2)=[0;-1]; % south
        elseif vertex_coor(face_matrix(i,1),2)==+scaling && vertex_coor(face_matrix(i,2),2)==+scaling
            face_normal(i,1:2)=[0;+1]; % north
        elseif vertex_coor(face_matrix(i,1),1)==+scaling && vertex_coor(face_matrix(i,2),1)==+scaling
            face_normal(i,1:2)=[+1;0]; % east
        else
            face_normal(i,1:2)=[-1;0]; % west
        end
        
    elseif vertex_bc(face_matrix(i,1))>=1 && vertex_bc(face_matrix(i,2))>=1 ...
            && ( vertex_coor(face_matrix(i,1),1)==vertex_coor(face_matrix(i,2),1) ||  vertex_coor(face_matrix(i,1),2)==vertex_coor(face_matrix(i,2),2) )
        face_bc(i)=2;
        
        if vertex_coor(face_matrix(i,1),2)==-scaling && vertex_coor(face_matrix(i,2),2)==-scaling
            face_normal(i,1:2)=[0;-1]; % south
        elseif vertex_coor(face_matrix(i,1),2)==+scaling && vertex_coor(face_matrix(i,2),2)==+scaling
            face_normal(i,1:2)=[0;+1]; % north
        elseif vertex_coor(face_matrix(i,1),1)==+scaling && vertex_coor(face_matrix(i,2),1)==+scaling
            face_normal(i,1:2)=[+1;0]; % east
        else
            face_normal(i,1:2)=[-1;0]; % west
        end
        
    else
        face_bc(i)=0;        
        tang=vertex_coor(face_matrix(i,2),:)-vertex_coor(face_matrix(i,1),:);
        face_normal(i,1:2)=[tang(2);-tang(1)]/norm(tang);
    end
end


for i=1:size(element2vertex,1)
    for j=1:3
        tmp=face2neighbor(element2face(i,j),:);
        pos=find(tmp==0,1,'first');
        face2neighbor(element2face(i,j),pos)=i;
        
%         face(element(i).face(j)).neighbor=[face(element(i).face(j)).neighbor i];
    end
end


end

