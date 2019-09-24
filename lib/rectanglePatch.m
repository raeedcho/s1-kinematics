function rectanglePatch(pos,facecolor)
% This function plots a rectangle with face color using patch (so we can use the colormap if we want)
    x_verts = [pos(1) pos(1)+pos(3) pos(1)+pos(3) pos(1)];
    y_verts = [pos(2) pos(2) pos(2)+pos(4) pos(2)+pos(4)];

    patch(x_verts,y_verts,facecolor)
end
