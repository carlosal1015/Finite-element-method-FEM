function fig_handle=drawMesh(face2vertex,vertex_coor)

fig_handle=figure;
axes('FontSize',10); 
axis equal; 
xlabel(''); ylabel('');
hold on


for i=1:size(face2vertex)
    x(1)=vertex_coor(face2vertex(i,1),1); y(1)=vertex_coor(face2vertex(i,1),2);
    x(2)=vertex_coor(face2vertex(i,2),1); y(2)=vertex_coor(face2vertex(i,2),2);    
    plot(x,y,'-','Color','k');   
  
% %      text( (x(1)+x(2))/2,(y(1)+y(2))/2,num2str(i))
%     
%     if normal_flag==1
%        normal=face(i).normal/10*norm([x(1);y(1)]-[x(2);y(2)]);
%        m=[x(1)+x(2); y(1)+y(2)]/2;       
% %        plot([m(1);m(1)+normal(1)],[m(2);m(2)+normal(2)],'-','Color','r');  
%        quiver(m(1),m(2),normal(1),normal(2),'color','r')
%     end
end


% elNr=2;
% text(vertex(element(elNr).vertex(1)).coor(1),vertex(element(elNr).vertex(1)).coor(2),'1')
% text(vertex(element(elNr).vertex(2)).coor(1),vertex(element(elNr).vertex(2)).coor(2),'2')
% text(vertex(element(elNr).vertex(3)).coor(1),vertex(element(elNr).vertex(3)).coor(2),'3')


hold off


end