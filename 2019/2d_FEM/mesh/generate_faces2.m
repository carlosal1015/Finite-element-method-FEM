function [element2face,face2vertex,face_bc,face2neighbor,face_normal]=generate_faces2(element2vertex,vertex_coor,vertex_bc)
% ich habe element.vertex, vertex.coor und vertex.bctype

n=size(element2vertex,1);
element2face=zeros(n,3);

face_vertex_number_matrix=zeros(3*n,3);
face_vertex_number_matrix(:,3)=1:3*n;

node2face=cell(length(vertex_bc),1);

cnt=0;
for i=1:n
    % 1. face
    face_vertex=[element2vertex(i,1) element2vertex(i,2)];
        
    face_nr=node2face{element2vertex(i,1)}(ismember(node2face{element2vertex(i,1)},node2face{element2vertex(i,2)}));
    if isempty(face_nr)
        face_nr=cnt+1;
        cnt=cnt+1;
        node2face{element2vertex(i,1)}(end+1)=face_nr;
        node2face{element2vertex(i,2)}(end+1)=face_nr;        
    end
    
    face_vertex_number_matrix(face_nr,1:2)=face_vertex;   
    element2face(i,1)=face_nr;
    
    
    % 2. face
    face_vertex=[element2vertex(i,2) element2vertex(i,3)];
    
    face_nr=node2face{element2vertex(i,2)}(ismember(node2face{element2vertex(i,2)},node2face{element2vertex(i,3)}));
    if isempty(face_nr)
        face_nr=cnt+1;
        cnt=cnt+1;
        node2face{element2vertex(i,2)}(end+1)=face_nr;
        node2face{element2vertex(i,3)}(end+1)=face_nr;        
    end
    
    face_vertex_number_matrix(face_nr,1:2)=face_vertex;    
    element2face(i,2)=face_nr;
    
    % 3. face
    face_vertex=[element2vertex(i,3) element2vertex(i,1)];
    
    face_nr=node2face{element2vertex(i,3)}(ismember(node2face{element2vertex(i,3)},node2face{element2vertex(i,1)}) );
    if isempty(face_nr)
        face_nr=cnt+1;
        cnt=cnt+1;
        node2face{element2vertex(i,3)}(end+1)=face_nr;
        node2face{element2vertex(i,1)}(end+1)=face_nr;        
    end
    
    face_vertex_number_matrix(face_nr,1:2)=face_vertex;
    element2face(i,3)=face_nr;
end

% generate face
% face2vertex = zeros(cnt,2);
face2vertex=face_vertex_number_matrix(1:cnt,1:2);
face_bc     = zeros(cnt,1);
face2neighbor = zeros(cnt,2);
face_normal = zeros(cnt,2);

for i=1:cnt % ueber alle faces
    
    if vertex_bc(face2vertex(i,1))==1 && vertex_bc(face2vertex(i,2))==1 && min(abs(vertex_coor(face2vertex(i,1),:) - vertex_coor(face2vertex(i,2),:)) )==0 
        face_bc(i)=1; % dirichlet face
    elseif vertex_bc(face2vertex(i,1))>=1 && vertex_bc(face2vertex(i,2))>=1 && min(abs(vertex_coor(face2vertex(i,1),:) - vertex_coor(face2vertex(i,2),:)) )==0 
        face_bc(i)=2; % neumann face
    else
        face_bc(i)=0; % inner face
    end
    
    tang=vertex_coor(face_vertex_number_matrix(i,2),:)-vertex_coor(face_vertex_number_matrix(i,1),:);
    face_normal(i,1:2)=[tang(2);-tang(1)]/norm(tang);    
end

for i=1:n
    for j=1:3
        tmp=face2neighbor(element2face(i,j),:);
        pos=find(tmp==0,1,'first');
        face2neighbor(element2face(i,j),pos)=i;
    end
end


% ensure that normal is outer normal on boundary

for i=1:cnt
   if face2neighbor(i,2)==0 % aussenkante
       face_midpoint=(vertex_coor(face2vertex(i,1),:) + vertex_coor(face2vertex(i,2),:))/2;
       
       el_nr=face2neighbor(i,1);
       mid_point=(vertex_coor(element2vertex(el_nr,1),:) + vertex_coor(element2vertex(el_nr,2),:) + vertex_coor(element2vertex(el_nr,3),:) )/3;
       
%        Frage: auf welcher Seite der Halbebene liegt der mittelpunkt
       
       d=face_midpoint*face_normal(i,:)';
       dist=mid_point*face_normal(i,:)' - d;
       
       if dist>=0
          face_normal(i,:)=-face_normal(i,:);
       end
       
   end
end

end



