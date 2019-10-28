function [A,b]=stiffnessMatrix_2d_tri_opt(element2connectivity,element2vertex,vertex_coor,innerDofs, b,G, p)

d=innerDofs;

[gqk, gqw]=dunavant_triangle_quad(2);
basisPrime=evallocalbasis_prime(98,gqk,p);
basisPrimeX=basisPrime{1};
basisPrimeY=basisPrime{2};

cnt=0;
for k=1:size(element2connectivity,1) %over all possible elements
    espalte=element2connectivity(k,:);
    cnt=cnt+sum(espalte<=innerDofs)^2;
end

e1=zeros(1,cnt); e2=zeros(1,cnt); e3=zeros(1,cnt);


cnt=1;
for k=1:size(element2connectivity,1) %over all possible elements
    aelem=compLocalStiffnessMatrix2d(element2vertex,vertex_coor,k,p, gqw,basisPrimeX,basisPrimeY);
    
    espalte = element2connectivity(k,:);
    ezeile  = 1:length(espalte);
    eval    = ones(size(espalte));
    
    for i=1:length(ezeile) % test function
        for j=1:length(ezeile) % trial function
            if espalte(i)<=innerDofs && espalte(j)<=innerDofs
                e1(cnt)=espalte(i);
                e2(cnt)=espalte(j);
                e3(cnt)=eval(i)*eval(j)*aelem(ezeile(i),ezeile(j));
                cnt=cnt+1;
            end
            
            if espalte(i)>innerDofs && espalte(j)<=innerDofs
                b(espalte(j))=b(espalte(j))-G(espalte(i))*eval(i)*eval(j)*aelem(ezeile(i),ezeile(j));
            end
            
        end
    end
end

A=sparse(e1(1:cnt-1),e2(1:cnt-1),e3(1:cnt-1),d,d);
end


function aelem=compLocalStiffnessMatrix2d(element2vertex,vertex_coor,elementNr,p, gqw,basisPrimeX,basisPrimeY)

p1=vertex_coor(element2vertex(elementNr,1),:)';
p2=vertex_coor(element2vertex(elementNr,2),:)';
p3=vertex_coor(element2vertex(elementNr,3),:)';

% a  =p1;
a1 =p2-p1;
a2 =p3-p1;
a2a2=a2'*a2;
a1a1=a1'*a1;
a1a2=a1'*a2;
% det_jac_Ableitung=a1(1)*a2(2)-a1(2)*a2(1);
% det_jac=abs(det_jac_Ableitung);
det_jac=abs(a1(1)*a2(2)-a1(2)*a2(1)); % verzerrungsfaktor

if p==1
    aelem=[a2a2+a1a1-2*a1a2 -a2a2+a1a2 -a1a1+a1a2
        -a2a2+a1a2 a2a2 -a1a2
        -a1a1+a1a2  -a1a2  a1a1 ]/det_jac/2;
elseif p==2
    aelem=zeros(6,6);
    
    for i=1:6 % test
        for j=i:6 % ansatz
%             aelem(i,j)=sum( (a2a2*basisPrimeX(i,:).*basisPrimeX(j,:) + a1a1*basisPrimeY(i,:).*basisPrimeY(j,:) - a1a2*(basisPrimeX(i,:).*basisPrimeY(j,:) +basisPrimeY(i,:).*basisPrimeX(j,:) )   ).*gqw' );
            aelem(i,j)=(a2a2*basisPrimeX(i,:).*basisPrimeX(j,:) + a1a1*basisPrimeY(i,:).*basisPrimeY(j,:) - a1a2*(basisPrimeX(i,:).*basisPrimeY(j,:) +basisPrimeY(i,:).*basisPrimeX(j,:) )  )*gqw;
            aelem(j,i)=aelem(i,j); % symmetry
        end
    end
    aelem=aelem/det_jac/2;  % addition /2 wegen dieser dunavant_triangle_quad
end



end
