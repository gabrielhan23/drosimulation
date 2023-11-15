%% Simulate LGE
clear simLateR1vectorDEC simLateR1cxm simLateR1etk simLateR1cxmD
tempsimLateR1vectorDEC=zeros(size(fixed));
%fit AIF
%[tempCpexp2,tempfitresultsvector]=fit(time(2:end)',Cp(2:end)','exp2');%Fit Cp
[tempCpexp2,tempfitresultsvector]=fit([time(2:end)'; 6*60*60],[Cp(2:end)';0],'exp2')%Fit Cp adding 0 as the last
timesim=0:60:1800;
Cpsim=tempCpexp2(timesim);%Cp at 1000sec
Cpsim(1)=0;
%%
for c=1:sum(HeartMask(:))
   sigin=Ctoi(c,:);
   tempsimLateR1= etk(fitresultsvectorDECetk(c,:)',[timesim',Cpsim]);
   
   simLateR1vectorDEC(c,:)=tempsimLateR1;
   %results=a*exp(-bt)+d*exp(-et) and gof
end
%
for n=2:size(simLateR1vectorDEC,2)
    tempsimLateR1vectorDEC(HeartMask)=simLateR1vectorDEC(:,n);
simLateR1etk(:,:,n-1)=tempsimLateR1vectorDEC;
end
%handle=implay(simLateR1etk*1000/7);
%set(handle.Parent, 'Name', 'etk')
%{
figure;
for s=1:15
subplot(5,3,s)
temp=simLateR1etk(:,:,s*2);
Mn=max(mean(temp(RemoteMask>0&~isnan(temp))),0);
SD=std(temp((RemoteMask>0)&~isnan(temp)));
imagesc(temp);colormap gray;caxis([0 Mn+5*SD])
end

title(['ETKsim_T',num2str(timepoints),'min'])
saveas(gcf,[subjectfolder,'\DCE\ETKsim_T',num2str(timepoints),'min.fig'])
%}

%% plot ETK curve
%{

for n=1:size(simLateR1etk,3)
    temp=simLateR1etk(:,:,n);
    sigcurve_bloodetk(n)=mean(temp(Blood));
    sigcurve_myoetk(n)=mean(temp(Myomask==1));
    sigcurve_MIetk(n)=mean(temp(MIMask==1));
    sigcurve_MVOetk(n)=mean(temp(MVOMask==1));
    sigcurve_Remoteetk(n)=mean(temp(RemoteMask==1));
end
T_interp=1:5:1800;
B_interp=interp1(timesim(2:end),Cpsim(2:end),T_interp,'spline');
figure;hold on
plot(timesim(2:end)/60,Cpsim(2:end),'.b','MarkerSize',20);hold on;
plot(T_interp/60,B_interp,'-b','LineWidth',3)
MI_interp=interp1(timesim(2:end),sigcurve_MIetk,T_interp,'spline');
plot(timesim(2:end)/60,sigcurve_MIetk,'.r','MarkerSize',20)
plot(T_interp/60,MI_interp,'-r','LineWidth',3)
MVO_interp=interp1(timesim(2:end),sigcurve_MVOetk,T_interp,'spline');
plot(timesim(2:end)/60,sigcurve_MVOetk,'.g','MarkerSize',20)
plot(T_interp/60,MVO_interp,'-g','LineWidth',3)
Remote_interp=interp1(timesim(2:end),sigcurve_Remoteetk,T_interp,'spline');
plot(timesim(2:end)/60,sigcurve_Remoteetk,'.black','MarkerSize',20)
plot(T_interp/60,Remote_interp,'-black','LineWidth',3)
xlim([1 30])
title(['ETK curve'])
%}
%%
%%
for c=1:sum(HeartMask(:))
   sigin=Ctoi(c,:);
   %tempsimLateR1= cxm(fitresultsvectorDECcxm(c,:)',[timesim',Cpsim]);
   tempsimLateR1= cxm_b(fitresultsvectorDECcxm(c,:)',[timesim',Cpsim]);
   
   simLateR1vectorDEC(c,:)=tempsimLateR1;
   %results=a*exp(-bt)+d*exp(-et) and gof
end
%
for n=2:size(simLateR1vectorDEC,2)
    tempsimLateR1vectorDEC(HeartMask)=simLateR1vectorDEC(:,n);
simLateR1cxm(:,:,n-1)=tempsimLateR1vectorDEC;
end
%figure;imagesc(simLateR1(:,:,end));axis equal
%handle=implay(simLateR1cxm*500/7);
%set(handle.Parent, 'Name', 'cxm')
%{
figure;
for s=1:15
subplot(5,3,s)
temp=simLateR1cxm(:,:,s*2);
Mn=max(mean(temp(RemoteMask>0&~isnan(temp))),0);
SD=std(temp((RemoteMask>0)&~isnan(temp)));
imagesc(temp);colormap gray;caxis([0 Mn+5*SD])
end
title(['CXMsim_T',num2str(timepoints),'min'])
saveas(gcf,[subjectfolder,'\DCE\CXMsim_T',num2str(timepoints),'min.fig'])
%}

%% cxmD
%{
for c=1:sum(HeartMask(:))
   sigin=Ctoi(c,:);
   %tempsimLateR1= cxm(fitresultsvectorDECcxm(c,:)',[timesim',Cpsim]);
   tempsimLateR1D= cxmD(fitresultsvectorDECcxmD(c,:)',[timesim',Cpsim]);
   
   simLateR1vectorDECD(c,2:length(tempsimLateR1D)+1)=tempsimLateR1D;
   %results=a*exp(-bt)+d*exp(-et) and gof
end
%
for n=2:size(simLateR1vectorDECD,2)
    tempsimLateR1vectorDEC(HeartMask)=simLateR1vectorDECD(:,n)+Ctoi(:,2);
simLateR1cxmD(:,:,n-1)=tempsimLateR1vectorDEC;
end
%}
%figure;imagesc(simLateR1(:,:,end));axis equal
%handle=implay(simLateR1cxm*500/7);
%set(handle.Parent, 'Name', 'cxm')
%{
figure;
for s=1:15
subplot(5,3,s)
temp=simLateR1cxmD(:,:,s*2);
Mn=max(mean(temp(RemoteMask>0&~isnan(temp))),0);
SD=std(temp((RemoteMask>0)&~isnan(temp)));
imagesc(temp);colormap gray;caxis([Mn-2*SD Mn+9*SD])
end
title(['CXMsimD_T',num2str(timepoints),'min'])
saveas(gcf,[subjectfolder,'\DCE\CXMDsimD_T',num2str(timepoints),'min.fig'])
%}

%% plot cxm curve
%{
for n=1:size(simLateR1etk,3)
    temp=simLateR1cxm(:,:,n);
    sigcurve_bloodcxm(n)=mean(temp(Blood));
    sigcurve_myocxm(n)=mean(temp(Myomask==1));
    sigcurve_MIcxm(n)=mean(temp(MIMask==1));
    sigcurve_MVOcxm(n)=mean(temp(MVOMask==1));
    sigcurve_Remotecxm(n)=mean(temp(RemoteMask==1));
end
T_interp=1:5:1800;
B_interp=interp1(timesim(2:end),Cpsim(2:end),T_interp,'spline');
figure;plot(timesim(2:end)/60,Cpsim(2:end),'.b','MarkerSize',20);hold on;
plot(T_interp/60,B_interp,'-b','LineWidth',3)
MI_interp=interp1(timesim(2:end),sigcurve_MIcxm,T_interp,'spline');
plot(timesim(2:end)/60,sigcurve_MIcxm,'.r','MarkerSize',20)
plot(T_interp/60,MI_interp,'-r','LineWidth',3)
MVO_interp=interp1(timesim(2:end),sigcurve_MVOcxm,T_interp,'spline');
plot(timesim(2:end)/60,sigcurve_MVOcxm,'.g','MarkerSize',20)
plot(T_interp/60,MVO_interp,'-g','LineWidth',3)
Remote_interp=interp1(timesim(2:end),sigcurve_Remotecxm,T_interp,'spline');
plot(timesim(2:end)/60,sigcurve_Remotecxm,'.black','MarkerSize',20)
plot(T_interp/60,Remote_interp,'-black','LineWidth',3)
xlim([1 30])
%ylim([0 0.016])
title(['CXM curve'])
%}
%%
function C_toi = cxm(beta,X)
    
    time=X(:,1);
    Cp = X(:,2);
    
    m1 = beta(1);
    m2 = beta(2);
    F=beta(3);
    B=beta(4);
    C_toi = F*convolution(Cp,B*exp(-m1*time)+(1-B)*exp(-m2*time));

end

function C_toi = cxm_b(beta,X)
    
    time=X(:,1);
    Cp = X(:,2);
    
    F = beta(1);
    PS = beta(2);
    vp=beta(3);
    ve=beta(4);
    a=(F+PS)/vp;
    b=PS/ve;
    c=F/vp;
    M1=0.5*(a+b+sqrt((a+b).^2-4*b*c));
    M2=0.5*(a+b-sqrt((a+b).^2-4*b*c));
    B=(M2-c)/(M2-M1);
    C_toi = F*convolution(Cp,B*exp(-M1*time)+(1-B)*exp(-M2*time));

end

function C_toi = etk(beta,X)
    
    time=X(:,1);
    Cp = X(:,2);
    
    Ktrans = beta(1);
    kepTOI = beta(2);
    Vp=beta(3);
    
    C_toi = Vp*Cp+Ktrans*convolution(Cp,exp(-kepTOI.*time));

end
function C_toi = cxmD(beta,X)
    
    time=X(:,1);
    Cp = X(:,2);
%     [tempCpexp2,tempfitresultsvector]=fit([time(2:end)'; 6*60*60],[Cp(2:end)';0.1],'exp2')%Fit Cp adding 0.1 as noise level
%     a1=tempCpexp2.a;
%     m1=-tempCpexp2.b;
%     a2=tempCpexp2.c;
%     m2=-tempCpexp2.d;
%    C_toi=PS*(a1/(ms-PS)*(exp(-PS*time)-exp(-m2*time)))
 %   C_toi = F*convolution(H,Cp);
 %   C_toid = vp*(Cpd)+ve*(PS*(convolution(Cp,exp(-PS/ve*time))-convolution(Cp,exp(-PS/ve*time)));
 %   Ced=PS*(convolution(Cpd,exp(-PS/ve*timed)))/ve;

     PS = beta(1); %1/t=sec^-1
    vp=beta(2); % 1/100
    ve=beta(3); %1/100

    Ce=(convolution(Cp,exp(-PS/ve*time)));
   
    C_toi = vp*Cp+ve*Ce;
 
 
 % original cxmD
 %{
 
            timed=time(2:end)-time(2);
            Cpd=Cp(2:end)-Cp(2);
    
    PS = beta(1); %1/t=sec^-1
    vp=beta(2); % 1/100
    ve=beta(3); %1/100

    Ce=(convolution(Cp,exp(-PS*time)));
   
    Ced=Ce(2:end)-Ce(2);
    C_toi = vp*Cpd+ve*Ced;
 
%}
end

function c = convolution(a, b, shape)
%CONVOLUTION MODIFIED BY JULIO CARDENAS, MAY OF 2011.
%   SAME THAN CONV BUT RETURN A VECTOR FROM a(1) to a(end), not the central
%   section as described for the usual convolution function.
%  

if ~isvector(a) || ~isvector(b)
  error(message('MATLAB:conv:AorBNotVector'));
end

if nargin < 3
    shape = 'full';
end

if ~ischar(shape)
  error(message('MATLAB:conv:unknownShapeParameter'));
end

% compute as if both inputs are column vectors
[rows,~]=size(a);
c = conv2(a(:),b(:),shape);
c=c(1:rows);

% restore orientation
if shape(1) == 'f'
    if length(a) > length(b)
        if size(a,1) == 1 %row vector
            c = c.';
        end
    else
        if size(b,1) == 1 %row vector
            c = c.';
        end
    end
else
    if size(a,1) == 1 %row vector
        c = c.';
    end
end

end