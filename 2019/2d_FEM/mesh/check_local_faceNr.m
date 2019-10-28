function element2face=check_local_faceNr(element2vertex,element2face,face2vertex)

for i=1:size(element2face,1)
    tmp=element2face(i,1:3);
    
    tmpa=element2vertex(i,1:2);
    tmpa1=tmpa(end:-1:1);
    
    tmpb=element2vertex(i,2:3);
    tmpb1=tmpb(end:-1:1);
    
    pos=zeros(3,1);
    for j=1:3
        tmp1=face2vertex(tmp(j),1:2);
        if all(tmp1==tmpa) || all(tmp1==tmpa1)
          pos(j)=1;
        elseif all(tmp1==tmpb) || all(tmp1==tmpb1)
           pos(j)=2;
       else
           pos(j)=3;
       end
    end
    element2face(i,1:3)=tmp(pos);
end

end