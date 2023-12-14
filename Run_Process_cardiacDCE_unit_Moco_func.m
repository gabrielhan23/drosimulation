function []=Run_Process_cardiacDCE_unit_Moco_func(subjectfolder_in, raw, label, subject)  
subjectfolder=subjectfolder_in;
Load_Dicom_Moco
if raw == 1
    folder = ['DCE_Raw' filesep label];
    mkdir(fullfile(subjectfolder, folder));
else
    folder = ['DCE' filesep label];
    mkdir(fullfile(subjectfolder, folder));
    
end

SeriesNum=0;
%% Assign scan time
timein=timeino_moco-min(timeino_moco);


%% Image registration
% since we have already performed the registration, no need of this
clear PostT1Reg
% [optimizer,metric] = imregconfig('multimodal');
fixed = double(DataPre.img);
moving = double(DataPost_moco(1).img);
% corrected_pre = double(imregister(moving,fixed,'affine',optimizer,metric));
% DataPre.img = corrected_pre;
for n = 1:length(DataPost_moco)
    PostT1Rego(:,:,n) = double(DataPost_moco(n).img);
end

%% Sort time and image
[timesort,To]=sort(timein);
%To=circshift(To,-1);%debug use;
PostT1Reg=PostT1Rego(:,:,To);
PreT1=double(DataPre.img);
timein=timesort;
%implay(PostT1Reg/1000)
%% concentration
dR1=(1./PostT1Reg-1./PreT1)*1e3; %s^-1 
gdrelaxivity=4.97; %s^-1*mM^-1
Gdcon=dR1./gdrelaxivity; %mM
%% Draw Reference ROI

if ~exist('Blood','var')
figure;imagesc(PostT1Rego(:,:,1), [100 , 2000]);
title(['blood pool'])
Blood=roipoly;

title(['Heart mask'])
HeartMask=roipoly;

title(['Myo mask_endo'])
LVendo=roipoly;

title(['Myo mask_epi'])
LVepi=roipoly;
Myomask=~LVendo.*LVepi;

title(['MI'])
MIMask=roipoly;

% lesion

% title(['MVO'])
% MVOMask=roipoly;
% save(fullfile(subjectfolder, 'ROI_MVOMask.mat'))
title(['Remote'])
RemoteMask=roipoly;

mkdir(fullfile(subjectfolder, 'DCE', label, 'ROI'))
save(fullfile(subjectfolder, 'DCE', label, 'ROI', 'ROI.mat'),'Blood', 'HeartMask', 'LVendo', 'Myomask', 'MIMask', 'RemoteMask')
end

%% mask avoid singularity from registration error
%{
HeartMask=HeartMask.*min(dR1,[],3)>=0;
HeartMask=HeartMask.*max(dR1,[],3)<Inf;
%}
HeartMask= HeartMask.*(min(dR1,[],3)>=0).*(max(dR1,[],3)<Inf);
HeartMask=HeartMask.*~isnan(HeartMask);
HeartMask=HeartMask>0;

% Mona: dilate the blood pool region
% se = strel('square',10);
% Blood = imdilate(Blood, se);
% figure, imagesc(Blood)
Blood=Blood.*(min(dR1,[],3)>=0).*(max(dR1,[],3)<Inf);
Blood=Blood.*~isnan(Blood);

Blood=(Blood>0);
%% Pixel wise fitting
clear tempdR1 Cp Ctoi time
%fitting paramaters

firstscantime=2;%1mins
% firstscantime=1;%1mins

% fitpointsset={[1 2], [1 2 3], [1 2 3 4] };
fitpointsset={[1 2 3], [1 2 3 4] };

% 
%fitpointsset={[1:3], [2:4], [1 3 4], [2:length(timein)], [1:length(timein)],[1 3 5],[1 4 6]};
for fc=1:length(fitpointsset)
    fitpoints=zeros(size(timein));
    fitpoints(fitpointsset{fc})=1;
    fitpointind=find(fitpoints);
    clear Cp Ctoi time
for n=1:length(fitpointind);
    tempdGdcon=Gdcon(:,:,fitpointind(n));
    Cp(:,n+1)=mean(mean(tempdGdcon(Blood)));
    Ctoi(:,n+1)=tempdGdcon(HeartMask);
    time(n+1)=timein(fitpointind(n))+firstscantime*60;

end
fitresults=repmat(zeros(size(fixed)),[1 1 5]);
timepoints=uint16(ceil(time(2:end)/60))

%% dR1 curve
clear sigcurve_blood sigcurve_myo sigcurve_MI sigcurve_MVO sigcurve_Remote
for n=1:size(dR1,3)
    temp=Gdcon(:,:,n);
    sigcurve_blood(n)=mean(temp(Blood));
    sigcurve_myo(n)=mean(temp(Myomask==1));
    sigcurve_MI(n)=mean(temp(MIMask==1));
    % sigcurve_MVO(n)=mean(temp(MVOMask==1));
    sigcurve_Remote(n)=mean(temp(RemoteMask==1));
end
%% plot dR1 curve
T_interp=min(timein(:)):20/60:max(timein(:));
B_interp=interp1(timein,sigcurve_blood,T_interp,'spline');
T_interp=T_interp+firstscantime;
figure;hold on;
plot(timein(fitpointind)/60,sigcurve_blood(fitpointind),'.b','MarkerSize',20);hold on;
plot(T_interp/60,B_interp,'-b','LineWidth',3)
plot(timein/60,sigcurve_myo,'.-')
MI_interp=interp1(timein,sigcurve_MI,T_interp,'spline');
plot(timein(fitpointind)/60,sigcurve_MI(fitpointind),'.r','MarkerSize',20)
plot(T_interp/60,MI_interp,'-r','LineWidth',3)
% if sum(MVOMask(:))>1
% MVO_interp=interp1(timein,sigcurve_MVO,T_interp,'spline');
% plot(timein(fitpointind)/60,sigcurve_MVO(fitpointind),'.g','MarkerSize',20)
% plot(T_interp/60,MVO_interp,'-g','LineWidth',3)
% end
Remote_interp=interp1(timein,sigcurve_Remote,T_interp,'spline');
plot(timein(fitpointind)/60,sigcurve_Remote(fitpointind),'.black','MarkerSize',20)
plot(T_interp/60,Remote_interp,'-black','LineWidth',3)
%%
T_interp=min(timein(:)):20/60:max(timein(:));
B_interp=interp1(timein,sigcurve_blood,T_interp,'spline');
T_interp=T_interp+firstscantime;
figure;hold on;
plot(timein(fitpointind)/60,sigcurve_blood(fitpointind)/sigcurve_blood(1),'.b','MarkerSize',20);hold on;
plot(T_interp/60,B_interp/sigcurve_blood(1),'-b','LineWidth',3)
plot(timein/60,sigcurve_myo,'.-')
MI_interp=interp1(timein,sigcurve_MI,T_interp,'spline');
plot(timein(fitpointind)/60,sigcurve_MI(fitpointind)/sigcurve_MI(1),'.r','MarkerSize',20)
plot(T_interp/60,MI_interp/sigcurve_MI(1),'-r','LineWidth',3)
% if sum(MVOMask(:))>1
% MVO_interp=interp1(timein,sigcurve_MVO,T_interp/sigcurve_MVO(1),'spline');
% plot(timein(fitpointind)/60,sigcurve_MVO(fitpointind)/sigcurve_MVO(1),'.g','MarkerSize',20)
% plot(T_interp/60,MVO_interp,'-g','LineWidth',3)
% end
Remote_interp=interp1(timein,sigcurve_Remote,T_interp,'spline');
plot(timein(fitpointind)/60,sigcurve_Remote(fitpointind)/Remote_interp(1),'.black','MarkerSize',20)
plot(T_interp/60,Remote_interp/Remote_interp(1),'-black','LineWidth',3)
%% ETK
NityFivpercent=0;
l=0;

x0=[0.2 0.05 0.15];
%x0=[0.06 0.005 0.15];
x0c=repmat(x0,[sum(HeartMask(:)) 1]);
x0ctemp=x0c;
Rsqc=-1*ones([sum(HeartMask(:)) 1]);
lb=[0 0 0];
ub=[1 1 1];
%figure(10);plot(Rsqc);hold on;
iter=0;
clear fitresultsvectorDECetk fitresultsDCEetk tempfitDCE;
while(NityFivpercent<0.9 & iter<5) % initial guess optimization

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
disp('ETK')
end
%
tek_Ktrans=fitresultsDCEetk(:,:,1);%0-1
tek_Kep=fitresultsDCEetk(:,:,2); %0-0.5 
%tek_Ve=fitresultsDCEetk(:,:,2); %0-0.5 
tek_Vp=fitresultsDCEetk(:,:,3)*100;% 0-100 percentage
tek_Ve=fitresultsDCEetk(:,:,1)./fitresultsDCEetk(:,:,2).*HeartMask; % 0-100 percentage
%tek_Kep=fitresultsDCEetk(:,:,1)./fitresultsDCEetk(:,:,2).*HeartMask; % 0-100 percentage
tek_Rsq=fitresultsDCEetk(:,:,end);
%%
figure;
subplot(3,2,1);imagesc(tek_Ktrans);axis equal;title(['tek_Ktrans Con#=T',num2str(timepoints),'min']);caxis([0 1]);xlim([1 size(fitresultsDCEetk,1)]);ylim([1 size(fitresultsDCEetk,2)]);
subplot(3,2,2);imagesc(tek_Kep);axis equal;title(['tek_Kep Con#=T',num2str(timepoints),'min']);caxis([0 0.01]);xlim([1 size(fitresultsDCEetk,1)]);ylim([1 size(fitresultsDCEetk,2)]);
subplot(3,2,3);imagesc(tek_Ve);axis equal;title(['tek_Ve Con#=T',num2str(timepoints),'min']);caxis([0 100]);xlim([1 size(fitresultsDCEetk,1)]);ylim([1 size(fitresultsDCEetk,2)]);
subplot(3,2,4);imagesc(tek_Vp);axis equal;title(['tek_Vp Con#=T',num2str(timepoints),'min']);caxis([0 100]);xlim([1 size(fitresultsDCEetk,1)]);ylim([1 size(fitresultsDCEetk,2)]);
subplot(3,2,5);imagesc(tek_Rsq);axis equal;title(['tek_Rsq Con#=T',num2str(timepoints),'min']);caxis([0 1]);xlim([1 size(fitresultsDCEetk,1)]);ylim([1 size(fitresultsDCEetk,2)]);

%% Delayed CXM


NityFivpercent=0;
l=0;
x0=[0.1 0.5 0.2];
x0c=repmat(x0,[sum(HeartMask(:)) 1]);
Rsqc=-1*ones([sum(HeartMask(:)) 1]);
%lb=[-3 -3 -3 -3];
lb=[0 0 0];
ub=[10 1 1];
iter=0;
clear fitresultsvectorDECcxmD fitresultsDCEcxmD;
while(NityFivpercent<0.85 &iter<1) % initial guess optimization

for c=1:sum(HeartMask(:))
   %results=a*exp(-bt)+d*exp(-et) and gof
   if Rsqc(c)<0.9   
   % x0=x0c(c,:).*((1-(rand(1,4)-0.5)*1.5));
   sigin=Ctoi(c,:);
   tempfitDCE= fitdcemri_etk(sigin',Cp',time',x0,lb,ub,'cxmD');
   %tempfitDCE= fitdcemri_etk(sigin',Cp',time',x0c(c,:),lb,ub,'cxmD');
   fitresultsvectorDECcxmD(c,:)=tempfitDCE;
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
    temp(HeartMask)=fitresultsvectorDECcxmD(:,n);
    fitresultsDCEcxmD(:,:,n)=temp;
end
fitQ= fitresultsDCEcxmD(:,:,end);
Myomaskadj=Myomask.*(fitQ>0);
f=sort(fitQ(Myomaskadj==1));
sum(Rsqc<0.9)/length(Rsqc)
NityFivpercent=f(floor(sum(Myomaskadj(:))*0.05))% get the 95% R2
%figure;imagesc(fitQ);axis equal;title(['95%=',num2str(NityFivpercent),'x0=',num2str(x0)])
iter=iter+1
disp('CXMD')
end

PS_cxmD=fitresultsDCEcxmD(:,:,1);% -1 - 1
vp_cxmD=fitresultsDCEcxmD(:,:,2);% 0-100
ve_cxmD=fitresultsDCEcxmD(:,:,3);% 0-100
Rsq_cxmD=fitresultsDCEcxmD(:,:,end);


%%
figure;
subplot(2,2,1);imagesc(PS_cxmD);axis equal;title(['PS_cxmD Con#=T',num2str(timepoints),'min']);caxis([0 0.1]);xlim([1 size(PS_cxmD,1)]);ylim([1 size(PS_cxmD,2)]);
subplot(2,2,2);imagesc(vp_cxmD);axis equal;title(['vp_cxmD Con#=T',num2str(timepoints),'min']);caxis([0 1]);xlim([1 size(PS_cxmD,1)]);ylim([1 size(PS_cxmD,2)]);
subplot(2,2,3);imagesc(ve_cxmD);axis equal;title(['ve_cxmD Con#=T',num2str(timepoints),'min']);caxis([0 1]);xlim([1 size(PS_cxmD,1)]);ylim([1 size(PS_cxmD,2)]);
subplot(2,2,4);imagesc(Rsq_cxmD);axis equal;title(['Rsq_cxmD Con#=T',num2str(timepoints),'min']);caxis([0 1]);xlim([1 size(PS_cxmD,1)]);ylim([1 size(PS_cxmD,2)]);


%% Delayed CXMDiff

if length(Cp)>4
NityFivpercent=0;
l=0;
x0=[0.001 0.05 0.2];
x0c=repmat(x0,[sum(HeartMask(:)) 1]);
Rsqc=-1*ones([sum(HeartMask(:)) 1]);
%lb=[-3 -3 -3 -3];
lb=[0 0 0];
ub=[2 1 100];
iter=0;
clear fitresultsvectorDECcxmDiff fitresultsDCEcxmDiff;
while(NityFivpercent<0.85 &iter<5) % initial guess optimization

for c=1:sum(HeartMask(:))
   %results=a*exp(-bt)+d*exp(-et) and gof
   if Rsqc(c)<0.9   
   % x0=x0c(c,:).*((1-(rand(1,4)-0.5)*1.5));
   sigin=Ctoi(c,:);
   tempfitDCE= fitdcemri_etk(sigin',Cp',time',x0,lb,ub,'cxmDiff');
   %tempfitDCE= fitdcemri_etk(sigin',Cp',time',x0c(c,:),lb,ub,'cxmD');
   fitresultsvectorDECcxmDiff(c,:)=tempfitDCE;
   Rsqc(c)=tempfitDCE(end);
       if tempfitDCE(end)> Rsqc(c)%initial guess potimization
        %results=a*exp(-bt)+d*exp(-et) and gof
        %x0c(c,:)=tempfitDCE(1:end-2);
        x0c(c,:)=x0;
       end
        x0=x0c(c,:).*((1-(rand(1,length(x0))-0.5)*1.5));
    end
end
for n=1:length(tempfitDCE);
    temp=zeros(size(fixed));
    temp(HeartMask)=fitresultsvectorDECcxmDiff(:,n);
    fitresultsDCEcxmDiff(:,:,n)=temp;
end
fitQ= fitresultsDCEcxmDiff(:,:,end);
Myomaskadj=Myomask.*(fitQ>0);
f=sort(fitQ(Myomaskadj==1));
sum(Rsqc<0.9)/length(Rsqc)
NityFivpercent=f(floor(sum(Myomaskadj(:))*0.05))% get the 95% R2
%figure;imagesc(fitQ);axis equal;title(['95%=',num2str(NityFivpercent),'x0=',num2str(x0)])
iter=iter+1
disp('CXMDiff')
end

PS_cxmDiff=fitresultsDCEcxmDiff(:,:,1);% -1 - 1
vp_cxmDiff=fitresultsDCEcxmDiff(:,:,2);% 0-100
ve_cxmDiff=fitresultsDCEcxmDiff(:,:,3);% 0-100
Rsq_cxmDiff=fitresultsDCEcxmDiff(:,:,end);


%%
figure;
subplot(2,2,1);imagesc(PS_cxmDiff);axis equal;title(['PS_cxmDiff Con#=T',num2str(timepoints),'min']);caxis([0 2]);xlim([1 size(PS_cxmD,1)]);ylim([1 size(PS_cxmD,2)]);
subplot(2,2,2);imagesc(vp_cxmDiff);axis equal;title(['vp_cxmDiff Con#=T',num2str(timepoints),'min']);caxis([0 1]);xlim([1 size(PS_cxmD,1)]);ylim([1 size(PS_cxmD,2)]);
subplot(2,2,3);imagesc(ve_cxmDiff);axis equal;title(['ve_cxmDiff Con#=T',num2str(timepoints),'min']);caxis([0 1]);xlim([1 size(PS_cxmD,1)]);ylim([1 size(PS_cxmD,2)]);
subplot(2,2,4);imagesc(Rsq_cxmDiff);axis equal;title(['Rsq_cxmDiff Con#=T',num2str(timepoints),'min']);caxis([0 1]);xlim([1 size(PS_cxmD,1)]);ylim([1 size(PS_cxmD,2)]);

end


%%
%2CXM bound
NityFivpercent=0;
l=0;
x0=[0.5 0.1 10 20];
x0c=repmat(x0,[sum(HeartMask(:)) 1]);
Rsqc=-1*ones([sum(HeartMask(:)) 1]);
%lb=[-3 -3 -3 -3];
lb=[0 -1 0 0];
ub=[5 1 100 100];
iter=0;
clear fitresultsvectorDECcxm fitresultsDCEcxm;
while(NityFivpercent<0.85 &iter<5) % initial guess optimization

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
fitQ= fitresultsDCEcxm(:,:,end);
Myomaskadj=Myomask.*(fitQ>0);
f=sort(fitQ(Myomaskadj==1));
sum(Rsqc<0.9)/length(Rsqc)
NityFivpercent=f(floor(sum(Myomaskadj(:))*0.05))% get the 95% R2
%figure;imagesc(fitQ);axis equal;title(['95%=',num2str(NityFivpercent),'x0=',num2str(x0)])
iter=iter+1
disp('CXM')
end
F_cxm=fitresultsDCEcxm(:,:,1);%  0-5
PS_cxm=fitresultsDCEcxm(:,:,2);% 0 - 100
vp_cxm=fitresultsDCEcxm(:,:,3);% 0-100
ve_cxm=fitresultsDCEcxm(:,:,4);% 0-100
Rsq_cxm=fitresultsDCEcxm(:,:,end);
%% Plot figures
figure;title(['T',num2str(timepoints),'min'])
subplot(3,2,1);imagesc(F_cxm);axis equal;title(['F_cxm Con#=T',num2str(timepoints),'min']);caxis([0 1]);xlim([1 size(fitresultsDCEcxm,1)]);ylim([1 size(fitresultsDCEcxm,2)]);
subplot(3,2,2);imagesc(PS_cxm);axis equal;title(['PS_cxm Con#=T',num2str(timepoints),'min']);caxis([0 1]);xlim([1 size(fitresultsDCEcxm,1)]);ylim([1 size(fitresultsDCEcxm,2)]);
subplot(3,2,3);imagesc(ve_cxm);axis equal;title(['ve_cxm Con#=T',num2str(timepoints),'min']);caxis([0 50]);xlim([1 size(fitresultsDCEcxm,1)]);ylim([1 size(fitresultsDCEcxm,2)]);
subplot(3,2,4);imagesc(vp_cxm);axis equal;title(['vp_cxm Con#=T',num2str(timepoints),'min']);caxis([0 10]);xlim([1 size(fitresultsDCEcxm,1)]);ylim([1 size(fitresultsDCEcxm,2)]);

subplot(3,2,5);imagesc(Rsq_cxm);axis equal;title(['Rsq_cxm Con#=T',num2str(timepoints),'min']);caxis([0 1]);xlim([1 size(fitresultsDCEcxm,1)]);ylim([1 size(fitresultsDCEcxm,2)]);
%% mixed fitting 2CXM

% identify area BP (VP~=1)
% in BP, PS=inf, F=inf, Vp=1 , Ve=0
BP_F_thre=otsuthresh(F_cxm(HeartMask));
BPmask=F_cxm>BP_F_thre;


% identify area MVO (VP~=0)
%in MVO F~=0, PS=inf 
MVO_Vp_thre=0.01;
MVOmask=(tek_Vp<MVO_Vp_thre)&HeartMask;

% Mix model: set number for illpost values
% ETK model for Vp small
% ETK model for F inf
% ETK model for PS inf
% 2CMX model for the rest

tekmix_Ktrans=tek_Ktrans;%0-1
tekmix_Kep=tek_Kep; %0-0.5 
tekmix_Vp=tek_Vp;% 0-100 percentage
tekmix_Ve=tek_Ve;
tekmix_Rsq=tek_Rsq;
F_cxmmix=F_cxm;%  0-5
PS_cxmmix=PS_cxm;% -1 - 100
vp_cxmmix=vp_cxm;% 0-100
ve_cxmmix=ve_cxm;% 0-100


F_cxmmix(MVOmask)=0;%  0-5
PS_cxmmix(MVOmask)=0;% -1 - 100
vp_cxmmix(MVOmask)=tek_Vp(MVOmask);% 0-100
ve_cxmmix(MVOmask)=tek_Ve(MVOmask);% 0-100
Rsq_cxmmix(MVOmask)=tek_Rsq(MVOmask);% 0-100

F_cxmmix(BPmask)=5;%  0-5
PS_cxmmix(BPmask)=1;% -1 - 100
vp_cxmmix(BPmask)=100;% 0-100
ve_cxmmix(BPmask)=0;% 0-100
Rsq_cxmmix=Rsq_cxm;
%% Plot figures
figure;title(['T',num2str(timepoints),'min'])
subplot(3,2,1);imagesc(F_cxmmix);axis equal;title(['F_cxmmix Con#=T',num2str(timepoints),'min']);caxis([0 1]);xlim([1 size(fitresultsDCEcxm,1)]);ylim([1 size(fitresultsDCEcxm,2)]);
subplot(3,2,2);imagesc(PS_cxmmix);axis equal;title(['PS_cxmmix Con#=T',num2str(timepoints),'min']);caxis([0 0.1]);xlim([1 size(fitresultsDCEcxm,1)]);ylim([1 size(fitresultsDCEcxm,2)]);
subplot(3,2,3);imagesc(ve_cxmmix);axis equal;title(['ve_cxmmix Con#=T',num2str(timepoints),'min']);caxis([0 50]);xlim([1 size(fitresultsDCEcxm,1)]);ylim([1 size(fitresultsDCEcxm,2)]);
subplot(3,2,4);imagesc(vp_cxmmix);axis equal;title(['vp_cxmmix Con#=T',num2str(timepoints),'min']);caxis([0 10]);xlim([1 size(fitresultsDCEcxm,1)]);ylim([1 size(fitresultsDCEcxm,2)]);

subplot(3,2,5);imagesc(Rsq_cxmmix);axis equal;title(['Rsq_cxm Con#=T',num2str(timepoints),'min']);caxis([0 1]);xlim([1 size(fitresultsDCEcxm,1)]);ylim([1 size(fitresultsDCEcxm,2)]);


%% ECV
Posttemp=PostT1Reg(:,:,find(abs(time-15*60)==min(abs(time-15*60))));
ECV=(100-40)*(1./Posttemp-1./fixed)/(mean(1./Posttemp(Blood))-mean(1./fixed(Blood)));
%figure;imagesc(ECV.*HeartMask);title(['ECV',num2str(length(fitpointind))]);axis equal;caxis([20 100])

%%
% try
    SimLate_R1
    simLateR1cxmmix=simLateR1cxm;
    MVOmask_sim=repmat(MVOmask,[1 1 30]);
    simLateR1cxmmix(MVOmask_sim)=simLateR1etk(MVOmask_sim);
    BPmask_sim=repmat(BPmask,[1 1 30]);
    
    simLateR1cxmmix(BPmask_sim)==simLateR1etk(BPmask_sim);
% catch
%     disp("The SimLate fitting failed")
%     simLateR1cxmmix = 0;
%     simLateR1etk = 0;
%     simLateR1cxm = 0;
%     simLateR1cxm_Masked = 0;
% end
%mixed



%% Blood pool mask from CXM f maps
% use remote F value to perform F map thresholding mean+5SD.
Remote_Fmean=mean(F_cxm(RemoteMask));
Remote_Fsd=std(F_cxm(RemoteMask));
Bloodthreshold=Remote_Fmean+Remote_Fsd*10;
%Bloodthreshold=mean(F_cxm(Blood))-5*std(F_cxm(Blood));
CXM_BloodMask=F_cxm>Bloodthreshold;

% Use the top threshold to create blood pool mask 
simLateR1cxm_Masked=(simLateR1cxm).*(-1*repmat(CXM_BloodMask,...
    [1 1 size(simLateR1cxm,3)])+0.5).*(simLateR1cxm>0).*repmat(HeartMask,[1 1 size(simLateR1cxm,3)]);


%% save
DerivedDCEMaps.F_cxm=F_cxm;
DerivedDCEMaps.vp_cxm=vp_cxm;
DerivedDCEMaps.PS_cxm=PS_cxm;
DerivedDCEMaps.ve_cxm=ve_cxm;
DerivedDCEMaps.Rsq_cxm=Rsq_cxm;
DerivedDCEMaps.simLateR1cxm=simLateR1cxm;
DerivedDCEMaps.simLateR1cxm_Masked=simLateR1cxm_Masked;
DerivedDCEMaps.tek_Ktrans=tek_Ktrans;
DerivedDCEMaps.tek_Kep=tek_Kep;
DerivedDCEMaps.tek_Vp=tek_Vp;
DerivedDCEMaps.tek_Ve=tek_Ve;
DerivedDCEMaps.tek_Rsq=tek_Rsq;
DerivedDCEMaps.simLateR1etk=simLateR1etk;

% Mona, July 3, 2023
DerivedDCEMaps.vp_cxmD=vp_cxmD;
DerivedDCEMaps.PS_cxmD=PS_cxmD;
DerivedDCEMaps.ve_cxmD=ve_cxmD;
DerivedDCEMaps.Rsq_cxmD=Rsq_cxmD;
DerivedDCEMaps.simLateR1cxmD=0;
%{
DerivedDCEMaps.vp_cxmD=vp_cxmD;
DerivedDCEMaps.PS_cxmD=PS_cxmD;
DerivedDCEMaps.ve_cxmD=ve_cxmD;
DerivedDCEMaps.Rsq_cxmD=Rsq_cxmD;
DerivedDCEMaps.simLateR1cxmD=simLateR1cxmD;
%}
DerivedDCEMaps.BloodMask_cxm=CXM_BloodMask;
%{
DerivedDCEMaps.tekb_Ktrans=tekb_Ktrans;
DerivedDCEMaps.tekb_Kep=tekb_Kep;
DerivedDCEMaps.tekb_Vp=tekb_Vp;
DerivedDCEMaps.tekb_Ve=tekb_Ve;
DerivedDCEMaps.tekb_Rsq=tekb_Rsq;
DerivedDCEMaps.simLateR1etk=simLateR1etk;
%}
DerivedDCEMaps.F_cxmmix=F_cxmmix;
DerivedDCEMaps.vp_cxmmix=vp_cxmmix;
DerivedDCEMaps.PS_cxmmix=PS_cxmmix;
DerivedDCEMaps.ve_cxmmix=ve_cxmmix;
DerivedDCEMaps.Rsq_cxmmix=Rsq_cxmmix;
DerivedDCEMaps.simLateR1cxmmix=simLateR1cxmmix;

%
Save_Files_DCE_batch_Moco
%%
end
