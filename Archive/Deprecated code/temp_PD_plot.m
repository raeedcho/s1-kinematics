%% plot PDs
for i=1:11
    subplot(4,5,i)
    polar(pi,1,'.')
    hold on
    polar([PD_simple(i) PD_simple(i)]*pi/180,[0 moddepth(i)/max(moddepth)],'b')
    
    [x1,y1]=pol2cart(PD_simple(i)*pi/180,moddepth(i)/max(moddepth)); % needed to fill up the space between the two CI
    [x2,y2]=pol2cart(CI_boot(1,i)*pi/180,moddepth(i)/max(moddepth));
    [x3,y3]=pol2cart(CI_boot(2,i)*pi/180,moddepth(i)/max(moddepth));

    %     jbfill(x1,y1,y2,'b','b',1,0.5);
    x_fill = [0 x2 x1 x3];
    y_fill = [0 y2 y1 y3];

    % fill(x_fill,y_fill,'r');
    patch(x_fill,y_fill,'b','facealpha',0.3);
    
    title(['Chan ' num2str(ul(i,1)) ', Unit ' num2str(ul(i,2))])
end