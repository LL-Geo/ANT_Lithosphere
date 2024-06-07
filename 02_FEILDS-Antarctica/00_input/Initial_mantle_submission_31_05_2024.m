clc
clear
ANT_PUM=load('C:\Users\00103168\LocalData\Todo\GHF_ANT\Lithosphere\code_clean\00_input\mantle\Ant_20_PUM_att_adiabte_F_priorden.txt');
ANT_PUM=ANT_PUM(192278:end,:);

ANT_913=load('C:\Users\00103168\LocalData\Todo\GHF_ANT\Lithosphere\code_clean\00_input\mantle\Ant_20_913_att_adiabte_F_priorden.txt');
ANT_913=ANT_913(192278:end,:);


Depth=-unique(-ANT_PUM(:,3));


Depth_box=Depth';
Depth_box=repmat(Depth_box,167*167,1);

Result_vs=reshape(ANT_913(:,8),167*167,79);



Result_T=reshape(ANT_PUM(:,5),167,167,79);

Sx=reshape(ANT_PUM(:,1),167,167,79);
Sy=reshape(ANT_PUM(:,2),167,167,79);

load('LAB_submission_31_05_2024.mat')


LAB_box=reshape(-LAB_new,167*167,1);
LAB_box=repmat(LAB_box,1,79);

mask=zeros(167*167,79);

mask(Depth_box>=LAB_box)=1;


Den_pum=reshape(ANT_PUM(:,6),167*167,79);

Den_913=reshape(ANT_913(:,6),167*167,79);

Den_in=Den_913.*mask+Den_pum.*(1-mask);

max_den_in=max(Den_in);
max_den_in=repmat(max_den_in,27889,1);

Den_in(Result_vs<0)=max_den_in(Result_vs<0);

Den_out=[reshape(Sx(:,:,1),167*167,1),reshape(Sy(:,:,1),167*167,1)];

Den_out(:,3:2:79*2+1)=Depth_box;
Den_out(:,4:2:79*2+2)=Den_in;


cHeader= cell(1,160);

cHeader(1:2) = {'X' 'Y' }; %dummy header

for i=1:79
    cHeader(1+i*2)={append('z',num2str(i))};
end

for i=1:79
    cHeader(2+i*2)={append('d',num2str(i))};
end

textHeader = strjoin(cHeader, ',');
%write header to file
fid = fopen('Mantle_submission.csv','w'); 
fprintf(fid,'%s\n',textHeader)
fclose(fid)
%write data to end of file
dlmwrite('Mantle_submission.csv',Den_out,'-append');


Den_mean=mean(Den_in,2);
figure()
scatter(Den_out(:,1),Den_out(:,2),[],Den_mean,'filled')

for i=1:1:16
    subplot(4,4,i)
    scatter(Den_out(:,1),Den_out(:,2),[],-Den_in(:,i*5-1))
    caxis([-3550,-3250])
    colormap('hot')
end

Result_vs=reshape(ANT_913(:,8),167*167,79);

for i=1:1:16
    subplot(4,4,i)
    scatter(Den_out(:,1),Den_out(:,2),[],Result_vs(:,i*5-1))
    colormap('hot')
end

figure()
scatter(Den_out(:,1),Den_out(:,2),[],Den_in(:,49),'filled')
 hold on
 contour(reshape(Den_out(:,1),167,167),reshape(Den_out(:,2),167,167),reshape(LAB_new,167,167)/1e3,[250,250],'LineColor','k')
 caxis([3250,3450])
