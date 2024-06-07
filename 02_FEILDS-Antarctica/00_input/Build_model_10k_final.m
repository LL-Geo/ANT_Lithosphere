clc
clear 

Bedmachine_file='C:\Users\00103168\LocalData\Antarctica\data_file\Bedmachine\BedMachineAntarctica-v3.nc';


% Load Data
x=ncread(Bedmachine_file,'x');
y=ncread(Bedmachine_file,'y');
mask=ncread(Bedmachine_file,'mask');
thickness=ncread(Bedmachine_file,'thickness');
bed=ncread(Bedmachine_file,'bed');
surface=ncread(Bedmachine_file,'surface');
geoid=ncread(Bedmachine_file,'geoid');

kernel = ones(20, 20) / (20*20);
Bed = imfilter(bed, kernel, 'same', 'replicate');
Bed=double(Bed(67:20:end-66,67:20:end-66));

Ice = imfilter(thickness, kernel, 'same', 'replicate');
Ice(thickness==0)=0;
Ice=double(Ice(67:20:end-66,67:20:end-66));

water=surface-thickness-bed;
Water = imfilter(water, kernel, 'same', 'replicate');
Water(water==0)=0;
Water=double(Water(67:20:end-66,67:20:end-66));


Sx=double(x(7:40:end));
Sy=double(y(7:40:end));
[Sx, Sy] = meshgrid(Sx, Sy);


x=-3300000:10000:3300000;
y=-3300000:10000:3300000;

[X, Y] = meshgrid(x,y);


Ice(Ice<=0)=0;
Bed=flipud(Bed');
Ice=flipud(Ice');
Water=flipud(Water');
Density=Ice*0.917+Water.*1.03;


load('C:\Users\00103168\LocalData\Todo\GHF_ANT\Lithosphere\code_clean\00_input\crust\ST2_output_f.mat')


SSB_10k = interp2(Sx, Sy, flipud(MeanSSB_th), X, Y);
SSB_Den_10k = interp2(Sx, Sy, flipud(MeanSSB_den), X, Y);

Moho_10k = interp2(Sx, Sy, flipud(MeanMoho), X, Y);
Crust_den_10k = interp2(Sx, Sy, flipud(MeanCrust_den), X, Y);


Upper_Crust_den_10k = interp2(Sx, Sy, flipud(Mean_Hete_den(:,:,6)), X, Y);
Lower_Crust_den_10k = interp2(Sx, Sy, flipud(Mean_Hete_den(:,:,7)), X, Y);


Base_SSB = Bed- SSB_10k; 

Suf_den=SSB_Den_10k;
Suf_den(isnan(SSB_Den_10k))=Upper_Crust_den_10k(isnan(SSB_Den_10k));



FX=reshape(X,661*661,1);
FY=reshape(Y,661*661,1);
FMoho=reshape(Moho_10k,661*661,1);
FBase_SSB=reshape(Base_SSB,661*661,1);
SDen=reshape(Suf_den,661*661,1);

uDen=reshape(Upper_Crust_den_10k,661*661,1);

lDen=reshape(Lower_Crust_den_10k,661*661,1);
Topo_den=reshape(Density,661*661,1);
FBed=reshape(Bed,661*661,1);


Layer_Surface=[FX FY FBed SDen Topo_den];

Layer_Moho=[FX FY FMoho lDen lDen-lDen+3.3];

% Build layer for interpolation
% X Y Ztop Dtop ZBottom DBottom
Layer_DenModel=[FX FY FBed SDen FBase_SSB SDen FBase_SSB uDen FBase_SSB-(FBase_SSB-FMoho)/2 uDen FBase_SSB-(FBase_SSB-FMoho)/2 lDen FMoho lDen];


cHeader = {'X' 'Y' 'Topography' 'Roh_c' 'M_surf'}; %dummy header
textHeader = strjoin(cHeader, ',');
%write header to file
fid = fopen('Surface_10k.csv','w'); 
fprintf(fid,'%s\n',textHeader)
fclose(fid)
%write data to end of file
dlmwrite('Surface_10k.csv',Layer_Surface,'-append');


cHeader = {'X' 'Y' 'Moho' 'Roh_c' 'Roh_M'}; %dummy header
textHeader = strjoin(cHeader, ',');
%write header to file
fid = fopen('Moho_10k.csv','w'); 
fprintf(fid,'%s\n',textHeader)
fclose(fid)
%write data to end of file
dlmwrite('Moho_10k.csv',Layer_Moho,'-append');


cHeader = {'X' 'Y' 'z1_t' 'm1_t' 'z1_b' 'm1_b' 'z2_t' 'm2_t' 'z2_b' 'm2_b' 'z3_t' 'm3_t' 'z3_b' 'm3_b'}; %dummy header
textHeader = strjoin(cHeader, ',');
%write header to file
fid = fopen('Crust_density.csv','w'); 
fprintf(fid,'%s\n',textHeader)
fclose(fid)
%write data to end of file
dlmwrite('Crust_density.csv',Layer_DenModel,'-append');




[AP,APP]= geotiffread('FAA_R.tif');
Faa=double(reshape(flipud((AP)),661*661,1));

Gravity=[FX FY FX FY Faa Faa-Faa+10000];

cHeader = {'Longitude' 'Latitude' 'X' 'Y' 'FAA' 'Height'}; %dummy header
textHeader = strjoin(cHeader, ',');
%write header to file
fid = fopen('Gravity.csv','w'); 
fprintf(fid,'%s\n',textHeader)
fclose(fid)
%write data to end of file
dlmwrite('Gravity.csv',Gravity,'-append');
