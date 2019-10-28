function fig_handle=drawMesh_lshape_zoom(element,face,vertex)

rad=1/8;
xmin=-rad; xmax=rad; ymin=-rad; ymax=rad;

fig_handle=figure('PaperPosition',[0 0 22 22],'PaperSize',[22 22],'Position',[200 400 750 500]);
axes('FontSize',10); axis equal; axis([xmin xmax ymin ymax]) 
xlabel(''); ylabel('');
hold on


for i=1:length(face)
    x(1)=vertex(face(i).vertex(1)).coor(1); y(1)=vertex(face(i).vertex(1)).coor(2);
    x(2)=vertex(face(i).vertex(2)).coor(1); y(2)=vertex(face(i).vertex(2)).coor(2); 
    if norm([x(1) y(1)],'inf')<=rad+10^-10 || norm([x(2) y(2)],'inf')<=rad+10^-10
        plot(x,y,'-','Color','k');   
    end
    
end


hold off


end