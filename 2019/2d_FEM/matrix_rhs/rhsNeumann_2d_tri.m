function b=rhsNeumann_2d_tri(element2connectivity,element2vertex,face2vertex,face_bc,face2neighbor,face_normal,vertex_coor,Neuman_func,gqn,innerDofs,p)

d=innerDofs; %innere regulaere dofs
b=zeros(d,1); %lege rhs-vector an

[gqk, gqw]=gaussvalues(p+gqn);
gqw=gqw(:);

for l=1:length(face_bc)
    if face_bc(l)==2 % ist neumann kante
        k=face2neighbor(l,1); % hat nur einen nachbarn
        p1=vertex_coor(face2vertex(l,1),:)'; %linker kanten knoten
        p2=vertex_coor(face2vertex(l,2),:)'; %rechter kanten knoten
        lengthOfFace=norm(p2-p1); %laenge der Kante
              
        espalte = element2connectivity(k,:);
        ezeile  = 1:length(espalte);
        eval    = ones(size(espalte));
          
        belem=compute_local_neumannRHS(element2vertex,vertex_coor,k,Neuman_func,lengthOfFace,p1,p2,gqk,gqw,face_normal(l,:),p);
        for i=1:length(ezeile) % test function
            if espalte(i)<=innerDofs
                b(espalte(i))=b(espalte(i))+eval(i)*belem(ezeile(i));
            end
        end
                
    end
end


end


function belem=compute_local_neumannRHS(element2vertex,vertex_coor,elementNr,Neuman_func,lengthOfFace,p1,p2,gqk,gqw,normal,p)

p1E=vertex_coor(element2vertex(elementNr,1),:)';
p2E=vertex_coor(element2vertex(elementNr,2),:)';
p3E=vertex_coor(element2vertex(elementNr,3),:)';
aE  =p1E;
a1E =p2E-p1E;
a2E =p3E-p1E;

a = (p1+p2)/2;
b = (p2-p1)/2;

tmp=[a(1)+b(1)*gqk
     a(2)+b(2)*gqk];

gqkE=[a1E a2E]\(tmp-aE);
gqkE=gqkE';
basis=evallocalbasis(98,gqkE,p);

tmpN=Neuman_func(tmp(1,:),tmp(2,:),normal).*gqw;
belem=basis*tmpN*lengthOfFace/2;
        
end

