function ax=ant_plot_4(x,y,ax,A,ca,colorset,colabel,seri,axble,colorbarback,colorbaruni)

edgx=0.05;
edgy=0.05;
 if ~exist('colorbarback','var')
     % third parameter does not exist, so default it to something
      colorbarback = nan;
 end
 
  if ~exist('colorbaruni','var')
     % third parameter does not exist, so default it to something
      colorbaruni = nan;
  end
 
  
if axble==1
    f=figure()
    u = f.Units;
    f.Units = 'centimeters';
    f.Position=[0 0 16 15];
    
    
    nseri='a)';
    ax = subplot(2,2,1)
    ax.Position = [edgx edgy+0.5 0.5-edgx 0.5-edgx]
elseif axble==2
    ax = subplot(2,2,2)
    ax.Position = [0.5 edgy+0.5 0.5-edgx 0.5-edgx]
    nseri='b)';

elseif axble==3
    ax = subplot(2,2,3)
    ax.Position = [edgx edgy/1.5 0.5-edgx 0.5-edgx]
    nseri='c)';
elseif axble==4
    ax = subplot(2,2,4)
    ax.Position = [0.5 edgy/1.5 0.5-edgx 0.5-edgx]
    nseri='d)';
   
end
    
la=-85:5:-50;
lon=0:1:360;
la2=-90:5:-50;
lon2=0:30:360;
[la,lon]=meshgrid(la,lon);
[lon2,la2]=meshgrid(lon2,la2);

[xt,yt]=ll2ps(la,lon);
[xt2,yt2]=ll2ps(la2,lon2);
grayColor = [.7 .7 .7];



h=imagesc(x,y,(A))
caxis([ca(1),ca(2)])

set(gca,'YDir','normal')
set(h, 'AlphaData', 1-isnan((A)))
hold on
plot(xt,yt,'Color',  grayColor);
hold on
plot(xt2,yt2,'Color', grayColor);

draw_coastline
axis equal

% 
xlim([-3330000 3330000]) 
ylim([-3330000 3330000]) 

if axble==1
    xticks([-2000000  0  2000000])
    xticklabels([])
    yticks([-2000000  0  2000000])
    yticklabels({'120^{\circ}W','90^{\circ}W','60^{\circ}W'})
    ytickangle(90)
    ax.XRuler.TickLabelGapOffset = -2;  % default = +2
    set(gca,'TickLength',[0 0])
    
elseif axble==2
    
    xticks([-2000000  0  2000000])
    xticklabels([])
    yticks([-2000000  0  2000000])
    yticklabels({'120^{\circ}E','90^{\circ}E','60^{\circ}E'})
    set(gca,'YAxisLocation','right')
    ytickangle(270)
    ax.XRuler.TickLabelGapOffset = -2;  % default = +2
    set(gca,'TickLength',[0 0])
    
elseif axble==3
    xticks([-2000000  0  2000000])
    xticklabels({'150^{\circ}W','180^{\circ}','150^{\circ}E'})
    yticks([-2000000  0  2000000])
    yticklabels({'120^{\circ}W','90^{\circ}W','60^{\circ}W'})
    ytickangle(90)
    ax.XRuler.TickLabelGapOffset = -2;  % default = +2
    set(gca,'TickLength',[0 0])
    
elseif axble==4
    xticks([-2000000  0  2000000])
    xticklabels({'150^{\circ}W','180^{\circ}','150^{\circ}E'})
    yticks([-2000000  0  2000000])
    yticklabels({'120^{\circ}E','90^{\circ}E','60^{\circ}E'})
    set(gca,'YAxisLocation','right')
    ytickangle(270)
    ax.XRuler.TickLabelGapOffset = -2;  % default = +2
    set(gca,'TickLength',[0 0])
end

ss=50000;
hold on
plot([2500000 3000000],[3000000.9-ss 3000000.9-ss],'k-','linewidth',1.8); 
hold on
plot([2000000 2500000],[2950000.9-ss 2950000.9-ss],'k-','linewidth',1.8); 
hold on
plot([2000000 3000000],[2920000.9-ss 2920000.9-ss],'k-','linewidth',0.2); 
hold on
plot([2000000 3000000],[3025000.9-ss 3025000.9-ss],'k-','linewidth',0.2); 
hold on
plot([2000000 2000000],[2920000.9-ss 3025000.9-ss],'k-','linewidth',0.2); 
hold on
plot([3000000 3000000],[2920000.9-ss 3025000.9-ss],'k-','linewidth',0.2); 


text(2500000,3350000,'1000 km','horiz','center','vert','top'); 

pos2 = plotboxpos(gca);
pos3= [pos2(1),pos2(2)+pos2(4)-0.07/2,0.07/2,0.07/2-0.001];


a=colorbar('west','Position',...
    [pos2(1)+0.03 pos2(2)+0.02 0.01 0.383989501312336/2],...
    'AxisLocation','in');

if isnan(colorset)
else
    [map, descriptorname, description] = colorcet(colorset);
    colormap(ax,map);
end

if isnan(colorbarback)
    a.Label.String = colabel;
    a.Label.Position(1) = -3;
    a.TickLength = .05;
else
    fh = gcf();
    cbh = findall(fh, 'Type','ColorBar');
    
    cbh=cbh(1);
    
    ax2 = axes('position',[pos2(1)+0.03 pos2(2)+0.02 0.01 0.383989501312336/2],'XTick',[],'YTick',[],'Box','on');
    ax2.Position(1) = pos2(1);
    ax2.Position(2) = pos2(2);
    ax2.Position(3) =ax2.Position(3)+.08;
    ax2.Position(4) = ax2.Position(4)+.03;
    
    % a=colorbar('west','Position',...
    %     [pos2(1)+0.03 pos2(2)+0.02 0.01 0.383989501312336/2],...
    %     'AxisLocation','in');
    
    cb2 = colorbar(ax2);
    cb2.Position = cbh.Position;
    cb2.FontSize = cbh.FontSize;
    cb2.Limits = cbh.Limits;
    cb2.LimitsMode = cbh.LimitsMode;
    cb2.TickDirection = cbh.TickDirection;
    cb2.YAxisLocation = cbh.YAxisLocation;
    caxis(ax2,cbh.Limits)
    
    if isnan(colorset)
    else
        [map, descriptorname, description] = colorcet(colorset);
        colormap(ax2,map);
    end
    ax2.Color = [1 1 1 colorbarback];
    
    cb2.Label.String = colabel;
    cb2.Label.Position(1) = -3;
    cb2.TickLength = .05;
end   
 
if isnan(colorbaruni)
    
elseif axble<4
    colorbar off
elseif axble==4
    colorbar off
    a=colorbar('Northoutside','Position',...
        [pos2(1)-0.23 pos2(2)+0.49 0.41 0.02],...
        'AxisLocation','in');
    a.Label.String = colabel;
    a.Label.Position(1) = a.Limits(2)+(a.Limits(2)-a.Limits(1))/15;
    a.Label.Position(2) = 1.5;
    a.TickLength = .04;
    if isnan(colorset)
    else
        [map, descriptorname, description] = colorcet(colorset);
        colormap(a,map);
    end
end

    
%%



annotation('rectangle',pos3,'FaceColor','white')

if isnan(seri)
    annotation('textbox',pos3,'String',seri,'FitBoxToText','off');
else
    annotation('textbox',pos3,'String',nseri,'FitBoxToText','off');
end

