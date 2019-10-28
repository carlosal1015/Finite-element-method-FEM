function b=dirichletdaten_interpolation_u(element2connectivity,element2vertex,face2vertex,face_bc,face2neighbor,vertex_coor,diriFunc,innerDofs, total_dofs, p)

d=total_dofs; %inner + dirichlet dofs
b=zeros(d,1);

for i=1:length(face_bc)
    if face_bc(i)==1 % ist dirichlet kante
        elementNr=face2neighbor(i,1); % hat nur einen nachbarn   
        p1=vertex_coor(face2vertex(i,1),:)';%linker kanten knoten
        p2=vertex_coor(face2vertex(i,2),:)'; %rechter kanten knoten
        if p==1
            gqk=[p1 p2];
        elseif p==2
            gqk=[p1 p2 (p1+p2)/2];
        end
         
        espalte = element2connectivity(elementNr,:);
        ezeile  = 1:length(espalte);
%         eval    = ones(size(espalte));

        [belem,ei]=compute_local_interpolation(element2vertex,vertex_coor,elementNr,diriFunc,p,gqk);
        for j=ei                 
            if espalte(j)>innerDofs
                b(espalte(j))=belem(ezeile(j));
            end
        end
    end
end

end

function [belem,ei]=compute_local_interpolation(element2vertex,vertex_coor,elementNr,diriFunc,p,gqk)

% nur fuer diese h-version
p1=vertex_coor(element2vertex(elementNr,1),:)';
p2=vertex_coor(element2vertex(elementNr,2),:)';
p3=vertex_coor(element2vertex(elementNr,3),:)';

a  =p1;
a1 =p2-p1;
a2 =p3-p1;

gqk=[a1 a2]\(gqk-a);
gqk=gqk';

basis=evallocalbasis(98,gqk,p);
tmp=abs(basis-1)<=10^-12;
[ei,ej,~]=find(tmp);

tmp1= diriFunc(a(1)+a1(1)*gqk(:,1)+a2(1)*gqk(:,2),a(2)+a1(2)*gqk(:,1)+a2(2)*gqk(:,2));
belem=zeros(3*p,1);
belem(ei)=tmp1(ej);
end