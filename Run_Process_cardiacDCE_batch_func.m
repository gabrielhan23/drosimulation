%clear all
%close hidden all
%addpath(genpath('D:\Box\Box\Grants\Cardiac DCE for 7 mins viability scan\Data\code\image-registration-master'));
%setup_image_registration_repository()
%% load Dicom
%subjectfolder='D:\Box\Box\Grants\Cardiac DCE for 7 mins viability scan\Data\code\RawData\Human\NCKU_May6\8 min T1 map\DCE\Slice3';
%subjectfolder_in='C:\Users\yanghj\Box\Yang Lab\Projects\Cardiac Delayed Phase DCE\Data\Structured data\Animals\Dave\Dave_D6';
function []=Run_Process_cardiacDCE_batch_func(subjectfolder_in)  
subjectfolder=subjectfolder_in;
mkdir([subjectfolder,'\DCE']);
Load_Dicom
subjectfolder=subjectfolder_in;
SeriesNum=0;
%% Assign scan time
timein=timeino-min(timeino);


%% Image registration
clear PostT1Reg
[optimizer,metric] = imregconfig('multimodal');
fixed=double(DataPre.img);
%fixed=double(DataPost(length(DataPost)).img);
for n=1:length(DataPost)
    moving=double(DataPost(n).img);
    PostT1(:,:,n) = moving;
    PostT1Rego(:,:,n) = double(imregister(moving,fixed,'affine',optimizer,metric));
end
%% Nonrigid image registration



%% Sort time and image
[timesort,To]=sort(timein);
PostT1Reg=PostT1Rego(:,:,To);
PreT1=double(DataPre.img);
timein=timesort;
%implay(PostT1Reg/1000)
%% concentration
dR1=(1./PostT1Reg-1./PreT1)*1e3; %s^-1 
gdrelaxivity=0.4; %s^-1*mM^-1
Gdcon=dR1./gdrelaxivity;
%% Draw Reference ROI

if ~exist('Blood','var')
figure;imagesc(PostT1Rego(:,:,end));title(['blood pool'])
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
title(['MI'])
MIMask=roipoly;
title(['MVO'])
MVOMask=roipoly;
title(['Remote'])
RemoteMask=roipoly;

end

%% mask avoid singularity from registration error
%{
HeartMask=HeartMask.*min(dR1,[],3)>=0;
HeartMask=HeartMask.*max(dR1,[],3)<Inf;
%}
HeartMask= HeartMask.*(min(dR1,[],3)>=0).*(max(dR1,[],3)<Inf);
HeartMask=HeartMask.*~isnan(HeartMask);
HeartMask=HeartMask>0;
Blood=Blood.*(min(dR1,[],3)>=0).*(max(dR1,[],3)<Inf);
Blood=Blood.*~isnan(Blood);

Blood=(Blood>0);
%% Pixel wise fitting
clear tempdR1 Cp Ctoi time
%fitting paramaters
firstscantime=1;%mins
%fittimelimit=8;

fitpointsset={[1:3], [2:4], [2:11], [1:length(timein)],[2 4 6]};
%fitpointsset={[1:3],[2:4], [2:11]};
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
timepoints=uint16(ceil(time(2:end)/60))

%% dR1 curve
clear sigcurve_blood sigcurve_myo sigcurve_MI sigcurve_MVO sigcurve_Remote
for n=1:size(dR1,3)
    temp=dR1(:,:,n);
    sigcurve_blood(n)=mean(temp(Blood));
    sigcurve_myo(n)=mean(temp(Myomask==1));
    sigcurve_MI(n)=mean(temp(MIMask==1));
    sigcurve_MVO(n)=mean(temp(MVOMask==1));
    sigcurve_Remote(n)=mean(temp(RemoteMask==1));
end
%{
T_interp=min(timein(:)):20/60:max(timein(:));
B_interp=interp1(timein,sigcurve_blood,T_interp,'spline');
T_interp=T_interp+firstscantime;
figure;hold on;
%plot(timein(fitpointind)/60,sigcurve_blood(fitpointind),'.b','MarkerSize',20);hold on;
%plot(T_interp/60,B_interp,'-b','LineWidth',3)
%plot(timein/60,sigcurve_myo,'.-')
MI_interp=interp1(timein,sigcurve_MI,T_interp,'spline');
plot(timein(fitpointind)/60,sigcurve_MI(fitpointind),'.r','MarkerSize',20)
plot(T_interp/60,MI_interp,'-r','LineWidth',3)
MVO_interp=interp1(timein,sigcurve_MVO,T_interp,'spline');
plot(timein(fitpointind)/60,sigcurve_MVO(fitpointind),'.g','MarkerSize',20)
plot(T_interp/60,MVO_interp,'-g','LineWidth',3)
Remote_interp=interp1(timein,sigcurve_Remote,T_interp,'spline');
plot(timein(fitpointind)/60,sigcurve_Remote(fitpointind),'.black','MarkerSize',20)
plot(T_interp/60,Remote_interp,'-black','LineWidth',3)
%}
%% ETK
NityFivpercent=0;
l=0;

x0=[0.2 0.05 0.15];
%x0=[0.06 0.005 0.15];
x0c=repmat(x0,[sum(HeartMask(:)) 1]);
x0ctemp=x0c;
Rsqc=-1*ones([sum(HeartMask(:)) 1]);
lb=[0 0 0];
ub=[1 0.5 1];
%figure(10);plot(Rsqc);hold on;
iter=0;
clear fitresultsvectorDECetk fitresultsDCEetk tempfitDCE;
while(NityFivpercent<0.85 & iter<1) % initial guess optimization

for c=1:sum(HeartMask(:))
    if Rsqc(c)<0.9
   sigin=Ctoi(c,:);
   
   tempfitDCE= fitdcemri_etk(sigin',Cp',time',x0ctemp(c,:),lb,ub,'etk');
   %tempfitDCE= fitdcemri_etk(sigin',Cp',time',x0,lb,ub,'etk');
        
         if tempfitDCE(end)> Rsqc(c)%initial guess optimization
       
            fitresultsvectorDECetk(c,:)=tempfitDCE;
            %results=a*exp(-bt)+d*exp(-et) and gof
            %x0c(c,:)=x0;
            %x0c(c,:)=tempfitDCE(1:end-1);
            Rsqc(c)=tempfitDCE(end);
            x0c(c,:)=x0ctemp(c,:);

        end
%       x0=rand(1,3).*ub;
          x0ctemp(c,:)= max(lb, x0c(c,:).*((1-(rand(1,3)-0.5)*1.5)));
          x0ctemp(c,:)=    min(ub, x0);
       
        end
end



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
%
tek_Ktrans=fitresultsDCEetk(:,:,1);%0-1
tek_Kep=fitresultsDCEetk(:,:,2); %0-0.5 
tek_Vp=fitresultsDCEetk(:,:,3)*100;% 0-100 percentage
tek_Ve=fitresultsDCEetk(:,:,1)./fitresultsDCEetk(:,:,2).*HeartMask; % 0-100 percentage
tek_Rsq=fitresultsDCEetk(:,:,end);
%%
%2CXM bound
NityFivpercent=0;
l=0;
x0=[0.5 5 10 20];
x0c=repmat(x0,[sum(HeartMask(:)) 1]);
Rsqc=-1*ones([sum(HeartMask(:)) 1]);
%lb=[-3 -3 -3 -3];
lb=[0 -1 0 0];
ub=[5 1 100 100];
iter=0;
clear fitresultsvectorDECcxm fitresultsDCEcxm;
while(NityFivpercent<0.85 &iter<10) % initial guess optimization


for c=1:sum(HeartMask(:))
   %results=a*exp(-bt)+d*exp(-et) and gof
   if Rsqc(c)<0.9   
   % x0=x0c(c,:).*((1-(rand(1,4)-0.5)*1.5));
   sigin=Ctoi(c,:);
   tempfitDCE= fitdcemri_etk(sigin',Cp',time',x0,lb,ub,'cxm_bound');
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
F_cxm=fitresultsDCEcxm(:,:,1);%  0-5
PS_cxm=fitresultsDCEcxm(:,:,2);% -1 - 1
vp_cxm=fitresultsDCEcxm(:,:,3);% 0-100
ve_cxm=fitresultsDCEcxm(:,:,4);% 0-100
Rsq_cxm=fitresultsDCEcxm(:,:,end);
%% Plot figures
figure;title(['T',num2str(timepoints),'min'])
subplot(5,2,1)
imagesc(tek_Ktrans);axis equal;title(['Ktrans Con#=T',num2str(timepoints),'min']);caxis([0 0.5]);xlim([1 size(fitresultsDCEcxm,1)]);ylim([1 size(fitresultsDCEcxm,2)]);
subplot(5,2,2);imagesc(tek_Kep);axis equal;title(['Kep Con#=',num2str(length(fitpointind))]);caxis([0 0.015]);xlim([1 size(fitresultsDCEcxm,1)]);ylim([1 size(fitresultsDCEcxm,2)]);
subplot(5,2,3);imagesc(tek_Vp);axis equal;title(['Vp Con#=',num2str(length(fitpointind))]);caxis([0 0.6]);xlim([1 size(fitresultsDCEcxm,1)]);ylim([1 size(fitresultsDCEcxm,2)]);
subplot(5,2,4);imagesc(tek_Ve);axis equal;title(['Ve Con#=',num2str(length(fitpointind))]);caxis([0 100]);xlim([1 size(fitresultsDCEcxm,1)]);ylim([1 size(fitresultsDCEcxm,2)]);
subplot(5,2,5);imagesc(tek_Rsq);axis equal;title(['Rsq Con#=',num2str(length(fitpointind))]);xlim([1 size(fitresultsDCEcxm,1)]);ylim([1 size(fitresultsDCEcxm,2)]);

subplot(5,2,6);imagesc(F_cxm);axis equal;title(['F Con#=',num2str(length(fitpointind))]);caxis([0 1]);xlim([1 size(fitresultsDCEcxm,1)]);ylim([1 size(fitresultsDCEcxm,2)]);
subplot(5,2,7);imagesc(vp_cxm);axis equal;title(['vp Con#=',num2str(length(fitpointind))]);caxis([0 50]);xlim([1 size(fitresultsDCEcxm,1)]);ylim([1 size(fitresultsDCEcxm,2)]);
subplot(5,2,8);imagesc(PS_cxm);axis equal;title(['PS Con#=',num2str(length(fitpointind))]);caxis([-1 1]);xlim([1 size(fitresultsDCEcxm,1)]);ylim([1 size(fitresultsDCEcxm,2)]);
%figure;imagesc(fitresultsDCE(:,:,1)./fitresultsDCE(:,:,2).*HeartMask);axis equal;title(['Ve Con#=',num2str(length(fitpointind))])
subplot(5,2,9);imagesc(ve_cxm);axis equal;title(['ve Con#=',num2str(length(fitpointind))]);caxis([0 100]);xlim([1 size(fitresultsDCEcxm,1)]);ylim([1 size(fitresultsDCEcxm,2)]);
subplot(5,2,10);imagesc(Rsq_cxm);axis equal;title(['Rsq Con#=',num2str(length(fitpointind))]);xlim([1 size(fitresultsDCEcxm,1)]);ylim([1 size(fitresultsDCEcxm,2)]);
saveas(gcf,[subjectfolder,'\DCE\DCEmaps_T',num2str(timepoints),'min.png'])
saveas(gcf,[subjectfolder,'\DCE\DCEmaps_T',num2str(timepoints),'min.fig'])
% temp4=fitresultsDCEcxm(:,:,4);
% temp3=fitresultsDCEcxm(:,:,3);
% temp2=fitresultsDCEcxm(:,:,2);
% temp1=fitresultsDCEcxm(:,:,1);
% figure;scatter3(temp1(:),temp2(:),temp3(:),'o')
% figure;scatter3(temp1(:),temp3(:),temp4(:),'o')
%% Scatterplot
%{
Scatterinput(:,:,1)=tek_Ktrans;
Scatterinput(:,:,2)=tek_Kep;
Scatterinput(:,:,3)=tek_Vp;
DCE_Scatterplot(Scatterinput,RemoteMask,MIMask,MVOMask)
title(['Scatter Con#=',num2str(length(fitpointind))])
xlabel('etk_Ktrans') 
ylabel('etk_Kep') 
zlabel('etk_Vp')
legend({'Remote','MI','MVO'})
saveas(gcf,[subjectfolder,'\DCE\DCEmaps_T',num2str(timepoints),'min_scatterETK.fig'])

Scatterinput(:,:,1)=F_cxm;
Scatterinput(:,:,2)=vp_cxm;
Scatterinput(:,:,3)=PS_cxm;
DCE_Scatterplot(Scatterinput,RemoteMask,MIMask,MVOMask)
title(['Scatter Con#=',num2str(length(fitpointind))])
xlabel('cxm_F') 
ylabel('cxm_Vp') 
zlabel('cxm_PS')
legend({'Remote','MI','MVO'})
saveas(gcf,[subjectfolder,'\DCE\DCEmaps_T',num2str(timepoints),'min_scatterCXM1.fig'])

Scatterinput(:,:,1)=F_cxm;
Scatterinput(:,:,2)=vp_cxm;
Scatterinput(:,:,3)=ve_cxm;
DCE_Scatterplot(Scatterinput,RemoteMask,MIMask,MVOMask)
title(['Scatter Con#=',num2str(length(fitpointind))])
xlabel('cxm_F') 
ylabel('cxm_Vp') 
zlabel('cxm_Ve')
legend({'Remote','MI','MVO'})
saveas(gcf,[subjectfolder,'\DCE\DCEmaps_T',num2str(timepoints),'min_scatterCXM2.fig'])
%}

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
   tempfitDCE= fitdcemri_etk(sigin',Cp',time',x0,lb,ub,'Tofts');
   
   fitresultsvectorDEC(c,:)=tempfitDCE;
   %results=a*exp(-bt)+d*exp(-et) and gof
end
for n=1:3;
    temp=zeros(size(fixed));
    temp(HeartMask)=fitresultsvectorDEC(:,n);
    fitresultsDCE(:,:,n)=temp;
end

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
%figure;imagesc(ECV.*HeartMask);title(['ECV',num2str(length(fitpointind))]);axis equal;caxis([20 100])

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
DerivedDCEMaps.F_cxm=F_cxm;
DerivedDCEMaps.vp_cxm=vp_cxm;
DerivedDCEMaps.PS_cxm=PS_cxm;
DerivedDCEMaps.ve_cxm=ve_cxm;
DerivedDCEMaps.Rsq_cxm=Rsq_cxm;
DerivedDCEMaps.simLateR1cxm=simLateR1cxm;
DerivedDCEMaps.tek_Ktrans=tek_Ktrans;
DerivedDCEMaps.tek_Kep=tek_Kep;
DerivedDCEMaps.tek_Vp=tek_Vp;
DerivedDCEMaps.tek_Ve=tek_Ve;
DerivedDCEMaps.tek_Rsq=tek_Rsq;
DerivedDCEMaps.simLateR1etk=simLateR1etk;
%
Save_Files_DCE_batch
%%
end
