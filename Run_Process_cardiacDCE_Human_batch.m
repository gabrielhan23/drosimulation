clear all
close all
%addpath(genpath('D:\Box\Box\Grants\Cardiac DCE for 7 mins viability scan\Data\code\image-registration-master'));
%setup_image_registration_repository()
%% load Dicom
%subjectfolder='D:\Box\Box\Grants\Cardiac DCE for 7 mins viability scan\Data\code\RawData\Human\NCKU_May6\8 min T1 map\DCE\Slice3';
subjectfolder_in='D:\Box\Box\Yang Lab\Projects\Cardiac Delayed Phase DCE\Data\Structured data\Human\NCKU\10';
subjectfolder=subjectfolder_in;
Load_Dicom_Human
subjectfolder=subjectfolder_in;
Nslice=size(timeinos{1},2)
%%
for ns=1:Nslice
%% Assign scan time
clear Blood
for c=1:size(timeinos,2)
timeino(c)=timeinos{c}(ns);
end
timein=timeino-min(timeino);


%% Image registration
clear PostT1Reg
[optimizer,metric] = imregconfig('multimodal');
fixed=double(DataPre.img{ns});
%fixed=double(DataPost(length(DataPost)).img);
for n=1:length(DataPost)
    moving=double(DataPost(n).img{ns});
    PostT1(:,:,n) = moving;
    PostT1Rego(:,:,n) = double(imregister(moving,fixed,'affine',optimizer,metric));
end
%% Nonrigid image registration



%% Sort time and image
[timesort,To]=sort(timein);
PostT1Reg=PostT1Rego(:,:,To);
PreT1=double(DataPre.img{ns});
timein=timesort;
%implay(PostT1Reg/1000)
%% concentration
dR1=1./PostT1Reg-1./PreT1;
gdrelaxivity=0.4;
Gdcon=dR1./gdrelaxivity;
%% Draw Reference ROI

if ~exist('Blood','var')
figure;imagesc(fixed);title(['blood pool'])
Blood=roipoly;
%
title(['Heart mask'])
HeartMask=roipoly;
%
title(['Myo mask_endo'])
LVendo=roipoly;
title(['Myo mask_epi'])
LVepi=roipoly;
Myomask=~LVendo.*LVepi;
%title(['MI'])
%MIMask=roipoly;
end

%% mask avoid singularity from registration error
HeartMask= HeartMask.*(min(dR1,[],3)>=0).*(max(dR1,[],3)<Inf);
HeartMask=HeartMask>0;

%% Pixel wise fitting
clear tempdR1 Cp Ctoi clear time
%fitting paramaters
firstscantime=0;%mins
%fittimelimit=8;

fitpointsset={[2:length(timein)]};
for fc=1:length(fitpointsset)
    fitpoints=zeros(size(timein));
    fitpoints(fitpointsset{fc})=1;
    fitpointind=find(fitpoints);
    clear Cp Ctoi time
for n=1:length(fitpointind);
    tempdR1=dR1(:,:,fitpointind(n));
    Cp(:,n+1)=mean(mean(tempdR1(Blood)));
    Ctoi(:,n+1)=tempdR1(HeartMask);
    time(n+1)=timein(fitpointind(n))+firstscantime*60;

end
fitresults=repmat(zeros(size(fixed)),[1 1 5]);

%% dR1 curve
for n=1:size(dR1,3)
    temp=dR1(:,:,n);
    sigcurve_blood(n)=mean(temp(Blood));
    sigcurve_myo(n)=mean(temp(Myomask==1));
%    sigcurve_MI(n)=mean(temp(MIMask==1));
end
figure;plot(timein/60,sigcurve_blood,'.-');hold on;
plot(timein/60,sigcurve_myo,'.-')
%plot(timein/60,sigcurve_MI,'.-')
%% ETK
NityFivpercent=0;
l=0;
x0=[0.06 0.005 0.15];
x0c=repmat(x0,[sum(HeartMask(:)) 1]);
Rsqc=-1*ones([sum(HeartMask(:)) 1]);
lb=[0 0 0];
ub=[1 1 1];
%figure(10);plot(Rsqc);hold on;
iter=0;
clear fitresultsvectorDECetk fitresultsDCEetk tempfitDCE;
while(NityFivpercent<0.9 & iter<10) % initial guess optimization

for c=1:sum(HeartMask(:))
    if Rsqc(c)<0.9
   sigin=Ctoi(c,:);
   tempfitDCE= fitdcemri_etk(sigin',Cp',time',x0,lb,ub,'etk');
   %tempfitDCE= fitdcemri_etk(sigin',Cp',time',x0c(c,:),lb,ub,'etk');
        if tempfitDCE(end)> Rsqc(c)%initial guess potimization
            fitresultsvectorDECetk(c,:)=tempfitDCE;
            %results=a*exp(-bt)+d*exp(-et) and gof
            %x0c(c,:)=x0;
            %x0c(c,:)=tempfitDCE(1:end-1);
            Rsqc(c)=tempfitDCE(end);
            x0c(c,:)=x0;

        end
%       x0=rand(1,3).*ub;
          x0= max(lb, x0c(c,:).*((1-(rand(1,3)-0.5)*1.5)));
          x0=    min(ub, x0);
        
        end
end

%
%{
for c=2230%=1:sum(HeartMask(:))
    if Rsqc(c)<0.9
        x0c(c,:)
        %x0i=[x0c(c,:); max(lb, x0c(c,:).*((1-(rand(20,3)-0.5)*1.5)))];
        ub=[1 1 1 ]
        x0i=[x0c(c,:); max(lb, (rand(100,length(x0)).*repmat((ub-lb),100,1)+lb))];
        %x0i=    min(ub, x0i);
           sigin=Ctoi(c,:)*100;
        for k=1:size(x0i,1)
               tempfitDCE(k,:)= fitdcemri_etk(sigin',Cp',time',x0i(k,:),lb,ub,'etk')
        end
        MaxInd=find(tempfitDCE(:,end)==max(tempfitDCE(:,end)));
   %tempfitDCE= fitdcemri_etk(sigin',Cp',time',x0c(c,:),lb,ub,'etk');
        if tempfitDCE(end)> Rsqc(c)%initial guess potimization
            fitresultsvectorDECetk(c,:)=tempfitDCE(MinInd,:);
            %results=a*exp(-bt)+d*exp(-et) and gof
            %x0c(c,:)=x0;
            %x0c(c,:)=tempfitDCE(1:end-1);
            Rsqc(c)=fitresultsvectorDECetk(c,end);
            x0c(c,:)=x0i(MinInd,:);

        end

        
        end
end
%}
%
%figure(10);plot(Rsqc);

for n=1:length(tempfitDCE);
    temp=zeros(size(fixed));
    temp(HeartMask)=fitresultsvectorDECetk(:,n);
    fitresultsDCEetk(:,:,n)=temp;
end
fitQ= fitresultsDCEetk(:,:,end);
Myomaskadj=Myomask.*(fitQ>0);
f=sort(fitQ(Myomaskadj==1));
sum(Rsqc<0.9)/length(Rsqc)
NityFivpercent=f(floor(sum(Myomaskadj(:))*0.05))% get the 95% R2

%figure;imagesc(fitresultsDCEetk(:,:,4));axis equal;title(['95%=',num2str(NityFivpercent),'x0=',num2str(x0)])

iter=iter+1
end
%%
tek_Ktrans=fitresultsDCEetk(:,:,1);
tek_Kep=fitresultsDCEetk(:,:,2);
tek_Vp=fitresultsDCEetk(:,:,3);
tek_Ve=fitresultsDCEetk(:,:,1)./fitresultsDCEetk(:,:,2).*HeartMask;
tek_Rsq=fitresultsDCEetk(:,:,end);
%%
figure;imagesc(tek_Ktrans);axis equal;title(['Ktrans Con#=',num2str(length(fitpointind))]);caxis([0 1])
figure;imagesc(tek_Kep);axis equal;title(['Kep Con#=',num2str(length(fitpointind))]);caxis([0 0.05])
figure;imagesc(tek_Vp);axis equal;title(['Vp Con#=',num2str(length(fitpointind))]);caxis([0 1])
figure;imagesc(tek_Ve);axis equal;title(['Ve Con#=',num2str(length(fitpointind))]);caxis([0 100])
figure;imagesc(tek_Rsq);axis equal;title(['Rsq Con#=',num2str(length(fitpointind))])
temp3=fitresultsDCEetk(:,:,3);
temp3=temp3(Myomaskadj==1);
temp2=fitresultsDCEetk(:,:,2);
temp2=temp2(Myomaskadj==1);
temp1=fitresultsDCEetk(:,:,1);
temp1=temp1(Myomaskadj==1);
figure;scatter3(temp1(:),temp2(:),temp3(:),'o')

%% cxm
NityFivpercent=0;
l=0;
x0=[0.06 0.01 0.4 0.5];
x0c=repmat(x0,[sum(HeartMask(:)) 1]);
Rsqc=-1*ones([sum(HeartMask(:)) 1]);
lb=[-3 -3 -3 -3];
ub=[5 5 5 5];
iter=0;
clear fitresultsvectorDECcxm fitresultsDCEcxm;
while(NityFivpercent<0.85 &iter<10) % initial guess optimization


for c=1:sum(HeartMask(:))
   %results=a*exp(-bt)+d*exp(-et) and gof
   if Rsqc(c)<0.9   
   % x0=x0c(c,:).*((1-(rand(1,4)-0.5)*1.5));
   sigin=Ctoi(c,:);
   tempfitDCE= fitdcemri_etk(sigin',Cp',time',x0,lb,ub,'cxm');
   %tempfitDCE= fitdcemri_etk(sigin',Cp',time',x0c(c,:),lb,ub,'cxm');
   fitresultsvectorDECcxm(c,:)=tempfitDCE;
   Rsqc(c)=tempfitDCE(end);
       if tempfitDCE(end)> Rsqc(c)%initial guess potimization
        %results=a*exp(-bt)+d*exp(-et) and gof
        %x0c(c,:)=tempfitDCE(1:end-2);
        c
        x0
        x0c(c,:)=x0;
       end
        x0=x0c(c,:).*((1-(rand(1,length(x0))-0.5)*1.5));
    end
end
for n=1:length(tempfitDCE);
    temp=zeros(size(fixed));
    temp(HeartMask)=fitresultsvectorDECcxm(:,n);
    fitresultsDCEcxm(:,:,n)=temp;
end
fitQ= fitresultsDCEetk(:,:,end);
Myomaskadj=Myomask.*(fitQ>0);
f=sort(fitQ(Myomaskadj==1));
sum(Rsqc<0.9)/length(Rsqc)
NityFivpercent=f(floor(sum(Myomaskadj(:))*0.05))% get the 95% R2
%figure;imagesc(fitQ);axis equal;title(['95%=',num2str(NityFivpercent),'x0=',num2str(x0)])
iter=iter+1
end
%%
m1=fitresultsDCEcxm(:,:,1);
m2=fitresultsDCEcxm(:,:,2);
F=fitresultsDCEcxm(:,:,3);
B=fitresultsDCEcxm(:,:,4);
c=m2-B.*(m2-m1);
b=m1.*m2./c;
vp=c./F;
a=m1+m2-b;
PS=a.*vp-F;
ve=PS./b;
%%
figure;imagesc(F);axis equal;title(['F Con#=',num2str(length(fitpointind))]);caxis([0 1])
figure;imagesc(vp);axis equal;title(['vp Con#=',num2str(length(fitpointind))]);caxis([0 1])
figure;imagesc(PS);axis equal;title(['PS Con#=',num2str(length(fitpointind))]);caxis([-1 1])
%figure;imagesc(fitresultsDCE(:,:,1)./fitresultsDCE(:,:,2).*HeartMask);axis equal;title(['Ve Con#=',num2str(length(fitpointind))])
figure;imagesc(ve);axis equal;title(['ve Con#=',num2str(length(fitpointind))]);caxis([-20 20])
figure;imagesc(fitresultsDCEcxm(:,:,end));axis equal;title(['Rsq Con#=',num2str(length(fitpointind))])
temp4=fitresultsDCEcxm(:,:,4);
temp3=fitresultsDCEcxm(:,:,3);
temp2=fitresultsDCEcxm(:,:,2);
temp1=fitresultsDCEcxm(:,:,1);
figure;scatter3(temp1(:),temp2(:),temp3(:),'o')
figure;scatter3(temp1(:),temp3(:),temp4(:),'o')
%% DCE fit tofts
%Parameters for starting 1min

x0=[0.1 0.01 0];
lb=zeros(3,1);
ub=[2 2 10];
%%%%
clear fitresultsvectorDEC;
for c=1:sum(HeartMask(:))
   sigin=Ctoi(c,:);
   tempfitDCE= fitdcemri_etk(sigin',Cp',time',x0,lb,ub,'Tofts');
   
   fitresultsvectorDEC(c,:)=tempfitDCE;
   %results=a*exp(-bt)+d*exp(-et) and gof
end
for n=1:3;
    temp=zeros(size(fixed));
    temp(HeartMask)=fitresultsvectorDEC(:,n);
    fitresultsDCE(:,:,n)=temp;
end
%{
figure;imagesc(fitresultsDCE(:,:,1));axis equal;title('Ktrans');caxis([0 1])
figure;imagesc(fitresultsDCE(:,:,2));axis equal;title('Kep');caxis([0 0.01])
figure;imagesc(fitresultsDCE(:,:,3));axis equal;title('Toft Rsq')
 temp2=fitresultsDCE(:,:,2);
temp1=fitresultsDCE(:,:,1);
figure;plot(temp1(:),temp2(:),'o')
%}
%% ECV
Posttemp=PostT1Reg(:,:,find(abs(time-15*60)==min(abs(time-15*60))));
ECV=(100-40)*(1./Posttemp-1./fixed)/(mean(1./Posttemp(Blood))-mean(1./fixed(Blood)));
figure;imagesc(ECV.*HeartMask);title(['ECV',num2str(length(fitpointind))]);axis equal;caxis([20 100])

%%
SimLate_R1
%% Clustering
%{
clear ktemp
for n=1:4;
    temp=fitresultsDCEcxm(:,:,n);
    ktemp(n,:)=temp(HeartMask);%/max(temp(HeartMask));
end
[idx,C]=(fitresultsDCEcxm,4);
%[idx,C]=kmeans(ktemp',4);
idximcxm=zeros(size(temp));
idximcxm(HeartMask)=idx;
figure;imagesc(idximcxm);axis equal
%%
clear ktemp
for n=1:3;
    temp=fitresultsDCEetk(:,:,n);
    ktemp(n,:)=temp(HeartMask)/max(temp(HeartMask));
end
[idx,C]=kmeans(ktemp',4);
idximetk=zeros(size(temp));
idximetk(HeartMask)=idx;
figure;imagesc(idximetk);axis equal
%}
%% save
Save_Files_DCE_Human
end

%%
end
