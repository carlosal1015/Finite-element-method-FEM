function [element2connectivity,innerDofs,totaldof]=generate_connectivity_data(element2vertex,element2face,face_bc,vertex_bc,p)
% only regular mesh


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                H1 p=1 oder p=2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vertex2dof=zeros(size(vertex_bc));
if p==2
    face2dof=zeros(size(face_bc));
end

% setze dof in nicht Dirichlet Knoten
cnt_reg=0;
for i=1:length(vertex_bc)
    if vertex_bc(i)~=1 %inner knoten
        vertex2dof(i)=cnt_reg+1;
        cnt_reg=cnt_reg+1;
    end
end

if p==2
    for i=1:length(face_bc)
        if face_bc(i)~=1 %inner kante
            face2dof(i)=cnt_reg+1;
            cnt_reg=cnt_reg+1;
        end
    end
end

innerDofs=cnt_reg;

% setze Dirichlet Knoten dofs
for i=1:length(vertex_bc)
    if vertex_bc(i)==1 %Dirichlet knoten
        vertex2dof(i)=cnt_reg+1;
        cnt_reg=cnt_reg+1;
    end
end

if p==2
    for i=1:length(face_bc)
        if face_bc(i)==1 %Dirichlet kante
            face2dof(i)=cnt_reg+1;
            cnt_reg=cnt_reg+1;
        end
    end
end

totaldof=cnt_reg;

element2connectivity=zeros(size(element2vertex,1),3*p);

% setze connectivity matrix Kodierung
if p==1
    for i=1:size(element2vertex,1)
        element2connectivity(i,1:3)=[vertex2dof(element2vertex(i,1)) vertex2dof(element2vertex(i,2)) vertex2dof(element2vertex(i,3)) ];
    end
elseif p==2
    for i=1:size(element2vertex,1)
        element2connectivity(i,1:6)=[vertex2dof(element2vertex(i,1)) vertex2dof(element2vertex(i,2)) vertex2dof(element2vertex(i,3)) ...
            face2dof(element2face(i,1)) face2dof(element2face(i,2)) face2dof(element2face(i,3)) ];
    end
end

end