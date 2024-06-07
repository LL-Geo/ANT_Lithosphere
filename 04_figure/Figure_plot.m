Data_p=load('C:\Users\00103168\LocalData\Todo\GHF_ANT\Lithosphere\code_clean\05_postprocess\submit_version\ANT20_DTE_priorden.txt');

Data=load('C:\Users\00103168\LocalData\Todo\GHF_ANT\Lithosphere\code_clean\05_postprocess\submit_version\ANT20_DTE_postden.txt');

Data_PUM=load('C:\Users\00103168\LocalData\Todo\GHF_ANT\Lithosphere\code_clean\05_postprocess\submit_version\Ant_20_PUM_att_adiabte_F_priorden.txt');
Data_PUM=Data_PUM(415390:end,:);

T=reshape(Data_PUM(:,5),167,167,71);
Visco=reshape(Data_PUM(:,7),167,167,71);
Den=reshape(Data_PUM(:,6),167,167,71);
x_in=reshape(Data_PUM(:,1),167,167,71);
y_in=reshape(Data_PUM(:,2),167,167,71);
z_in=reshape(Data_PUM(:,3),167,167,71);

xtar=-3.3e6:1e4:3.3e6;
ytar=-3.3e6:1e4:3.3e6;
ztar=-4e5:10e3:-10e3;

[xtar,ytar,ztar]=meshgrid(xtar,ytar,ztar);

Visco_mat=griddata(x_in,y_in,z_in,Visco,xtar,ytar,ztar);
T_mat=griddata(x_in,y_in,z_in,T,xtar,ytar,ztar);
Den_mat=griddata(x_in,y_in,z_in,Den,xtar,ytar,ztar);

datablock=reshape(Data,661,661,[],9);
datablock_p=reshape(Data_p,661,661,[],9);

depth=unique(Data(:,3));

load("lab_final.mat")


LAB_old=readtable('C:\Users\00103168\LocalData\Todo\GHF_ANT\Lithosphere\code_clean\00_input\LAB_F_903_priorden.csv');
LAB_old=table2array(LAB_old);
LAB_oldx=reshape(LAB_old(:,1),661,661);
LAB_oldy=reshape(LAB_old(:,2),661,661);

LAB_old=reshape(LAB_old(:,3),661,661);

LAB_new=reshape(LAB_new,661,661);



datablock=reshape(Data,661,661,[],9);

depth=unique(Data(:,3));

datablock(:,:,31,6);

x=-3300000:10000:3300000;
y=-3300000:10000:3300000;
[X,Y]=meshgrid(x,y);

lab_d=[-100,-150,-200,-250];
for i=1:4
ant_plot_4(x,y,1,datablock(:,:,31-(i-1)*5,6),[3250 3450],'L17',' kg m^{-3}','b)',i,nan,1)

hold on

 contour(LAB_oldx,LAB_oldy,-LAB_new/1e3,[lab_d(i),lab_d(i)],'LineColor','b');
end

print(gcf,"Fig4.png",'-dpng','-r300')

for i=1:4
ant_plot_4(x,y,1,datablock(:,:,31-(i-1)*5,5),[400 1600],'L16','^{o}C','b)',i,nan,1)
hold on

 contour(LAB_oldx,LAB_oldy,-LAB_new/1e3,[lab_d(i),lab_d(i)],'LineColor','b');

end


print(gcf,"Fig5.png",'-dpng','-r300')



x=ncread('C:\Users\00103168\LocalData\Antarctica\data_file\Bedmachine\BedMachineAntarctica-v3.nc','x');
y=ncread('C:\Users\00103168\LocalData\Antarctica\data_file\Bedmachine\BedMachineAntarctica-v3.nc','y');

bed=ncread('C:\Users\00103168\LocalData\Antarctica\data_file\Bedmachine\BedMachineAntarctica-v3.nc','bed');

x=x(1:20:end);
y=y(1:20:end);
bed=bed(1:20:end,1:20:end);
[X,Y]=meshgrid(x,y);

x=-3330000:10000:3330000;
y=-3330000:10000:3330000;


for i=1:4
ant_plot_4(x,y,1,datablock(:,:,31-(i-1)*5,9),[89.3 94],'L17','Mg#','b)',i,nan,1)
hold on
 contour(X,Y,bed',[0,1000],'LineColor','w');
draw_coastline
end
print(gcf,"Fig6.png",'-dpng','-r300')


for i=1:4
ant_plot_4(x,y,1,log10(datablock(:,:,31-(i-1)*5,7)),[18 23],'L16','         log_{10} Pa s','b)',i,nan,1)
end

print(gcf,"Fig8.png",'-dpng','-r300')


%% Supporting figure
lab_d=[-100,-150,-200,-250];
for i=1:4
    ant_plot_4(x,y,1,datablock_p(:,:,31-(i-1)*5,6),[3250 3450],'L17',' kg m^{-3}','b)',i,nan,1)
    hold on
    contour(LAB_oldx,LAB_oldy,LAB_old/1e3,[lab_d(i),lab_d(i)],'LineColor','b');
end
print(gcf,"Fig_S3_Initial_Den.png",'-dpng','-r300')

for i=1:4
ant_plot_4(x,y,1,datablock_p(:,:,31-(i-1)*5,5),[400 1600],'L16','^{o}C','b)',i,nan,1)
hold on
contour(LAB_oldx,LAB_oldy,LAB_old/1e3,[lab_d(i),lab_d(i)],'LineColor','b');
end
print(gcf,"Fig_S4_Initial_T.png",'-dpng','-r300')

for i=1:4
ant_plot_4(x,y,1,datablock(:,:,31-(i-1)*5,5)-T_mat(:,:,31-(i-1)*5),[-200 200],'D01A','^{o}C','b)',i,nan,1)
end
print(gcf,"Fig_S5_DT.png",'-dpng','-r300')


for i=1:4
ant_plot_4(x,y,1,log10(datablock(:,:,31-(i-1)*5,7))-real(log10(Visco_mat(:,:,31-(i-1)*5))),[-10 10],'D01A','         log_{10} Pa s','b)',i,nan,1)
end
print(gcf,"Fig_S8_Vis_After.png",'-dpng','-r300')

for i=1:4
ant_plot_4(x,y,1,real(log10(Visco_mat(:,:,31-(i-1)*5))),[18 30],'L16','         log_{10} Pa s','b)',i,nan,1)
end
print(gcf,"Fig_S7_Vis_Initial.png",'-dpng','-r300')

for i=1:4
ant_plot_4(x,y,1,log10(datablock(:,:,31-(i-1)*5,7)),[18 30],'L16','         log_{10} Pa s','b)',i,nan,1)
end
print(gcf,"Fig_S6_Vis.png",'-dpng','-r300')
