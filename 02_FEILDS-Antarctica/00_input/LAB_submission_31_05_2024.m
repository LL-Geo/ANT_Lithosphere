ANT_PUM=load('C:\Users\00103168\LocalData\Todo\GHF_ANT\Lithosphere\code_clean\00_input\mantle\Ant_20_PUM_att_adiabte_F_priorden.txt');

ANT_PUM=ANT_PUM(192278:end,:);

Depth=unique(abs(ANT_PUM(:,3)));

Result_T=reshape(ANT_PUM(:,5),167*167,79);


Gradient_T=gradient(Result_T);

Result_T(Result_T<0)=nan;


load('C:\Users\00103168\LocalData\Todo\GHF_ANT\Crust_model\Crust_Model\Ant_Crust.mat')

Moho=MeanMoho;

xt=-3330000:20000:3330000;
yt=-3330000:20000:3330000;
Moho_interp=interp2(xt,yt,Moho,ANT_PUM(1:167*167,1),ANT_PUM(1:167*167,2));

LAB=zeros(167*167,1);
for i=1:167*167
    Target=0:500:400000;
    Depths=Depth;
    Depths(Depth<abs(Moho_interp(i)))=nan;
    Depths(Depth<40000)=nan;
    if abs(Moho_interp(i))>30000
        Depths(Depth<50000)=nan;
    end
    s=[Depths,Result_T(i,:)'];
    s=s(sum(isnan(s),2)==0,:);
    if isempty(s)
        continue
    else
        s=[0,0;s];
        s([diff(s(:,2))<0;logical(0)],2)=nan;
        s=s(sum(isnan(s),2)==0,:);
        ST=interp1(s(:,1)',s(:,2)',Target,'linear');
%         ST(gradient(ST)<0)=nan;
        if abs(Moho_interp(i))<20000
            ST(Target>300000)=nan;
        end
        Target(isnan(ST))=[];
        ST(isnan(ST))=[];
        ssk = sign(ST-1200);
        positiveGoingIndexes = strfind(ssk, [-1, 1]);
        if isempty(positiveGoingIndexes)
            continue
        else
        LAB(i)=Target(positiveGoingIndexes(1));
        end
    end
end

LAB_new=LAB;
LAB_new(LAB_new==0)=35000;




x=-3300000:10000:3300000;
y=-3300000:10000:3300000;

[X, Y] = meshgrid(x,y);

LAB_fie = griddata(ANT_PUM(1:167*167,1), ANT_PUM(1:167*167,2), -(LAB_new), X, Y);

LAB_file=[X(:),Y(:),LAB_fie(:)];

cHeader = {'X' 'Y' 'LAB'}; %dummy header
textHeader = strjoin(cHeader, ',');
%write header to file
fid = fopen('LAB_Submission_Run.csv','w'); 
fprintf(fid,'%s\n',textHeader)
fclose(fid)
%write data to end of file
dlmwrite('LAB_Submission_Run.csv',LAB_file,'-append');



