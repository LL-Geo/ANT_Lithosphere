
ANT_913=load('Ant_20_913_att_adiabte_F_priorden.txt');
ANT_913=ANT_913(415390:end,:);


Vs=reshape(ANT_913(:,4),167,167,71);

x_in=reshape(ANT_913(:,1),167,167,71);
y_in=reshape(ANT_913(:,2),167,167,71);
z_in=reshape(ANT_913(:,3),167,167,71);

den_a=readtable("g_den_a.csv");

den_all=readtable("g_den_all_0905.csv");

Den_all=table2array(den_all);
Den_a=table2array(den_a);

% figure()
% for i=1:12
%     subplot(3,4,i)
%     imagesc(flipud(reshape(Den_all(:,i*3),661,661)))
% end
% 
% figure()
% for i=1:40
%     subplot(4,10,i)
%     imagesc(flipud(reshape(Den_a(:,i),661,661)))
%     colorcet('D01A')
%     caxis([-40,40])
% end


xtar=-3.3e6:1e4:3.3e6;
ytar=-3.3e6:1e4:3.3e6;
ztar=-4e5:10e3:-10e3;

[xtar,ytar,ztar]=meshgrid(xtar,ytar,ztar);

VS_mat=griddata(x_in,y_in,z_in,Vs,xtar,ytar,ztar);

Den_all_in=Den_all(:,1:40);
Vs_t=[xtar(:),ytar(:),ztar(:),VS_mat(:),Den_all_in(:)];

save ANT20_grid_PostInversion.txt Vs_t -ascii

Den_all_in=Den_all(:,1:40)-Den_a(:,1:40);
Vs_t=[xtar(:),ytar(:),ztar(:),VS_mat(:),Den_all_in(:)];

save ANT20_grid_PriorInversion.txt Vs_t -ascii


