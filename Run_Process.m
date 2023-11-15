clear all
%addpath(genpath('D:\Box\Box\Grants\Cardiac DCE for 7 mins viability scan\Data\code\image-registration-master'));
%setup_image_registration_repository()
%% load Dicom
% subjectfolder='D:\Box\Box\Grants\Cardiac DCE for 7 mins viability scan\Data\code\Data\SofiaD7';
subjectfolder='D:\Box\Box\Grants\Cardiac DCE for 7 mins viability scan\Data\code\RawData\Lisbon_D7';
PreconFoldername=[subjectfolder,'\PreconT1'];
PostconFoldername=[subjectfolder,'\PostconT1'];

listPost=dir(PostconFoldername);
for n=1:length(listPost)-2
    DataPost(n)=loaddicom([PostconFoldername,filesep,listPost(n+2).name]);
    timeinstr{n}=DataPost(n).info.AcquisitionTime;
    timein(n)=str2num(timeinstr{n}(1:2))*60*60+str2num(timeinstr{n}(3:4))*60+str2num(timeinstr{n}(5:end));%s
end
listPre=dir(PreconFoldername);
for n=1:length(listPre)-2
    DataPre(n)=loaddicom([PreconFoldername,filesep,listPre(n+2).name]);
end

%% Assign scan time
timein=timein-min(timein)+60;

clear PostT1Reg
%% Image registration
[optimizer,metric] = imregconfig('multimodal');
fixed=double(DataPre.img);
for n=1:length(DataPost)
    moving=double(DataPost(n).img);
    PostT1(:,:,n) = moving;
    PostT1Reg(:,:,n) = double(imregister(moving,fixed,'affine',optimizer,metric));
end
implay(PostT1Reg/1000)
%% Sort time and image
[timesort,To]=sort(timein);
PostT1Reg=PostT1Reg(:,:,To);
timein=timesort;

%% concentration
dR1=1./PostT1Reg-1./fixed;
gdrelaxivity=0.4;
Gdcon=dR1./gdrelaxivity;
%% Draw Reference ROI

figure;imagesc(fixed);title(['blood pool'])
Blood=roipoly;
%%
title(['Heart mask'])
HeartMask=roipoly;
%%
title(['Myo mask_endo'])
LVendo=roipoly;
title(['Myo mask_epi'])
LVepi=roipoly;
Myomask=~LVendo.*LVepi;

%% Pixel wise fitting
clear tempdR1 Cp Ctoi clear time
%fitting paramaters
fittimelimit=200;
fitpoints=zeros(size(timein));
fitpoints(1:4)=1;
fitpointind=find(fitpoints);
for n=1:length(fitpointind);
    tempdR1=dR1(:,:,fitpointind(n));
    Cp(:,n+1)=mean(mean(tempdR1(Blood)));
    Ctoi(:,n+1)=tempdR1(HeartMask);
    time(n+1)=timein(n);

end

fitresults=repmat(zeros(size(fixed)),[1 1 5]);


%% DCE fit tofts
%Parameters for starting 1min
%{
x0=[0.1 0.01 0];
lb=zeros(3,1);
ub=[2 2 10];
%%%%
clear fitresultsvectorDEC;
for c=1:sum(HeartMask(:))
   sigin=Ctoi(c,:);
   tempfitDCE= fitdcemri(sigin',Cp',time',x0,lb,ub,'Tofts');
   
   fitresultsvectorDEC(c,:)=tempfitDCE;
   %results=a*exp(-bt)+d*exp(-et) and gof
end
for n=1:3;
    temp=zeros(size(fixed));
    temp(HeartMask)=fitresultsvectorDEC(:,n);
    fitresultsDCE(:,:,n)=temp;
end

figure;imagesc(fitresultsDCE(:,:,1));axis equal
figure;imagesc(fitresultsDCE(:,:,2));axis equal
figure;imagesc(fitresultsDCE(:,:,3));axis equal
 temp2=fitresultsDCE(:,:,2);
temp1=fitresultsDCE(:,:,1);
figure;plot(temp1(:),temp2(:),'o')
%}
%% ETK
x0=[0.1 0.08 0.2]+abs(rand(1,3))*0.1;
lb=[-1 -1 -1];
ub=[3 3 1];
%%%%
clear fitresultsvectorDEC;
for c=1:sum(HeartMask(:))
   sigin=Ctoi(c,:);
   tempfitDCE= fitdcemri_etk(sigin',Cp',time',x0,lb,ub,'etk');
   
   fitresultsvectorDECetk(c,:)=tempfitDCE;
   %results=a*exp(-bt)+d*exp(-et) and gof
end
for n=1:4;
    temp=zeros(size(fixed));
    temp(HeartMask)=fitresultsvectorDECetk(:,n);
    fitresultsDCEetk(:,:,n)=temp;
end

figure;imagesc(fitresultsDCEetk(:,:,1).*HeartMask);axis equal;title(['Ktrans Con#=',num2str(length(fitpointind))])
figure;imagesc(fitresultsDCEetk(:,:,2));axis equal;title(['Kep Con#=',num2str(length(fitpointind))])
figure;imagesc(fitresultsDCEetk(:,:,3).*HeartMask);axis equal;title(['Vp Con#=',num2str(length(fitpointind))])
figure;imagesc(fitresultsDCEetk(:,:,1)./fitresultsDCEetk(:,:,2).*HeartMask);axis equal;title(['Ve Con#=',num2str(length(fitpointind))])
figure;imagesc(fitresultsDCEetk(:,:,4));axis equal;title(['Rsq Con#=',num2str(length(fitpointind))])
temp3=fitresultsDCEetk(:,:,3);
temp3=temp3(Myomask==1);
temp2=fitresultsDCEetk(:,:,2);
temp2=temp2(Myomask==1);
temp1=fitresultsDCEetk(:,:,1);
temp1=temp1(Myomask==1);
figure;scatter3(temp1(:),temp2(:),temp3(:),'o')

%% cxm

x0=[0.1 0.01 0.5 0.5]+abs(rand(1,4))*0.1;
lb=zeros(4,1);
ub=[2 2 1 1];
%%%%
clear fitresultsvectorDEC;
for c=1:sum(HeartMask(:))
   sigin=Ctoi(c,:);
   tempfitDCE= fitdcemri_etk(sigin',Cp',time',x0,lb,ub,'cxm');
   
   fitresultsvectorDECcxm(c,:)=tempfitDCE;
   %results=a*exp(-bt)+d*exp(-et) and gof
end
for n=1:5;
    temp=zeros(size(fixed));
    temp(HeartMask)=fitresultsvectorDECcxm(:,n);
    fitresultsDCEcxm(:,:,n)=temp;
end
m1=fitresultsDCEcxm(:,:,1);
m2=fitresultsDCEcxm(:,:,2);
F=fitresultsDCEcxm(:,:,3);
B=fitresultsDCEcxm(:,:,4);
c=m2-B.*(m2-m1);
b=m1.*m2./-c;
vp=c./F;
a=m1+m2-b;
PS=a.*vp-F;
ve=PS./b;
figure;imagesc(F);axis equal;title(['F Con#=',num2str(length(fitpointind))])
figure;imagesc(vp);axis equal;title(['vp Con#=',num2str(length(fitpointind))])
figure;imagesc(PS);axis equal;title(['PS Con#=',num2str(length(fitpointind))])
%figure;imagesc(fitresultsDCE(:,:,1)./fitresultsDCE(:,:,2).*HeartMask);axis equal;title(['Ve Con#=',num2str(length(fitpointind))])
figure;imagesc(ve);axis equal;title(['ve Con#=',num2str(length(fitpointind))])
figure;imagesc(fitresultsDCEcxm(:,:,5));axis equal;title(['Rsq Con#=',num2str(length(fitpointind))])
temp4=fitresultsDCEcxm(:,:,4);
temp3=fitresultsDCEcxm(:,:,3);
temp2=fitresultsDCEcxm(:,:,2);
temp1=fitresultsDCEcxm(:,:,1);
figure;scatter3(temp1(:),temp2(:),temp3(:),'o')
figure;scatter3(temp1(:),temp3(:),temp4(:),'o')

%% ECV
Posttemp=PostT1Reg(:,:,find(abs(time-15*60)==min(abs(time-15*60))));
ECV=(100-40)*(1./Posttemp-1./fixed)/(mean(1./Posttemp(Blood))-mean(1./fixed(Blood)));
figure;imagesc(ECV.*HeartMask);title(['ECV',num2str(length(fitpointind))])
% %%
% clear sigin;
% %bi exponential fit
% for c=1:sum(HeartMask(:))
%    sigin=Ctoi(c,:);
%    [tempfitexp2,tempfitresultsvector5]=fit(time(fitpoints)',sigin(fitpoints)','exp2');
%    
%    fitresultsvector(c,:)=[tempfitexp2.a,tempfitexp2.b,tempfitexp2.c,tempfitexp2.d,tempfitresultsvector5.rsquare]';
%    %results=a*exp(-bt)+d*exp(-et) and gof
% end
% for n=1:5;
%     temp=zeros(size(fixed));
%     temp(HeartMask)=fitresultsvector(:,n);
%     fitresults(:,:,n)=temp;
% end
% 
% %% simulated LGE
% 
% LGEt=7000;
% EEt=10;
% LGER1=HeartMask.*fitresults(:,:,1).*exp(-fitresults(:,:,2)*LGEt)+fitresults(:,:,3).*exp(-fitresults(:,:,4)*LGEt);
% EER1=HeartMask.*fitresults(:,:,1).*exp(-fitresults(:,:,2)*EEt)+fitresults(:,:,3).*exp(-fitresults(:,:,4)*EEt);
% figure;imagesc(LGER1);imcontrast;
% figure;imagesc(EER1);imcontrast;
% %clim([0.005 0.055])
%%
SimLate_R1
%% Clustering
clear ktemp
for n=1:4;
    temp=fitresultsDCEcxm(:,:,n);
    ktemp(n,:)=temp(HeartMask);
end

[idx,C]=kmeans(ktemp(1:3,:)',6);
idxim=zeros(size(temp));
idxim(HeartMask)=idx;
figure;imagesc(idxim)

%% save
mkdir([subjectfolder,'/DCE/Images']);
save([subjectfolder,'/DCE/workspace']) 

mkdir([subjectfolder,'/DCE']);
save([subjectfolder,'/DCE/workspace']) 
%%
function data = loaddicom(path)
list=dir(path);
if ~(strcmp('IMA',list(3).name(end-2:end)))

  fprintf('no Dicom in the directory')  
else
data.img=dicomread([path,filesep,list(3).name]);
data.info=dicominfo([path,filesep,list(3).name]);

end
end
%%

