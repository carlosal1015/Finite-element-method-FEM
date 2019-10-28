function b=rhsVolume_2d_tri(element2connectivity,element2vertex,vertex_coor,f,gqn,innerDofs,p)

d=innerDofs; %innere regulaere dofs
b=zeros(d,1); %lege rhs-vector an

[gqk, gqw]=dunavant_triangle_quad(p+gqn);
basis=evallocalbasis(98,gqk,p);

for k=1:size(element2connectivity,1) %over all possible elements
    espalte = element2connectivity(k,:);
    ezeile  = 1:length(espalte);
    eval    = ones(size(espalte));

    belem=compLocalVolumeRHS(element2vertex,vertex_coor,f,k,gqk,gqw,basis);
    for i=1:length(ezeile) % test function
        if espalte(i)<=innerDofs
            b(espalte(i))=b(espalte(i))+eval(i)*belem(ezeile(i));
        end
    end
end

end

function belem=compLocalVolumeRHS(element2vertex,vertex_coor,f,elementNr,gqk,gqw,basis)

p1=vertex_coor(element2vertex(elementNr,1),:);
p2=vertex_coor(element2vertex(elementNr,2),:);
p3=vertex_coor(element2vertex(elementNr,3),:);

a  =p1;
a1 =p2-p1;
a2 =p3-p1;
det_jac=abs(a1(1)*a2(2)-a1(2)*a2(1)); % verzerrungsfaktor

tmp1= f(a(1)+a1(1)*gqk(:,1)+a2(1)*gqk(:,2),a(2)+a1(2)*gqk(:,1)+a2(2)*gqk(:,2)).*gqw;
belem=basis*tmp1*det_jac/2;

end

