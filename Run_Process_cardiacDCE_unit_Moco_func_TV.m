% For DCE with cxmmix model, using the total variation
% Mona, March 25, 2024

function [] = Run_Process_cardiacDCE_unit_Moco_func_TV(subjectfolder, basefolder, raw, label, subject, tips)
 % load the data
folder = ['DCE_' raw filesep label '_' tips];

outputfolder = fullfile(subjectfolder, folder);
mkdir(outputfolder);

Load_Dicom_Moco

fixed=double(DataPre.img);
for n = 1:length(DataPost_moco)
    PostT1Rego(:,:,n) = double(DataPost_moco(n).img);
end
% Assign scan time
SeriesNum=0;
timein=timeino_moco-min(timeino_moco);
% Sort time and image
[timesort,To]=sort(timein);
PostT1Reg=PostT1Rego(:,:,To);
PreT1=double(DataPre.img);
timein=timesort;

%% concentration
dR1=(1./PostT1Reg-1./PreT1)*1e3; %s^-1 
gdrelaxivity=4.97; %s^-1*mM^-1
Gdcon=dR1./gdrelaxivity; %mM

%% Draw Reference ROI

if ~exist('Blood','var')
figure;imagesc(PostT1Rego(:,:,2), [100 , 1000]);
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

mkdir(fullfile(subjectfolder, 'MASK'))
save(fullfile(subjectfolder, 'MASK', [label, '_ROI.mat']),'Blood', 'HeartMask', 'LVendo', 'Myomask', 'MIMask', 'RemoteMask')
end

%% mask avoid singularity from registration error
HeartMask= HeartMask.*(min(dR1,[],3)>=0).*(max(dR1,[],3)<Inf);
HeartMask=HeartMask.*~isnan(HeartMask);
HeartMask=HeartMask>0;

Blood=Blood.*(min(dR1,[],3)>=0).*(max(dR1,[],3)<Inf);
Blood=Blood.*~isnan(Blood);

Blood=(Blood>0);

%% Pixel wise fitting
clear tempdR1 Cp Ctoi time
firstscantime=2;
fitpointsset={[1 2 3], [1 2 3 4]};
for fc=1:length(fitpointsset)
    fitpoints=zeros(size(timein));
    fitpoints(fitpointsset{fc})=1;
    fitpointind=find(fitpoints);
    clear Cp Ctoi time
    for n=1:length(fitpointind)
        tempdGdcon=Gdcon(:,:,fitpointind(n));
        Cp(:,n+1)=mean(mean(tempdGdcon(Blood)));
        Ctoi(:,n+1)=tempdGdcon(HeartMask);
        time(n+1)=timein(fitpointind(n))+firstscantime*60;
    
    end
    fitresults=repmat(zeros(size(fixed)),[1 1 5]);
    timepoints=uint16(ceil(time(2:end)/60))
    
    %% plot dR1 curve
    clear sigcurve_blood sigcurve_myo sigcurve_MI sigcurve_MVO sigcurve_Remote
    for n=1:size(dR1,3)
        temp=Gdcon(:,:,n);
        sigcurve_blood(n)=mean(temp(Blood));
        sigcurve_myo(n)=mean(temp(Myomask==1));
        sigcurve_MI(n)=mean(temp(MIMask==1));
        % sigcurve_MVO(n)=mean(temp(MVOMask==1));
        sigcurve_Remote(n)=mean(temp(RemoteMask==1));
    end
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
    Remote_interp=interp1(timein,sigcurve_Remote,T_interp,'spline');
    plot(timein(fitpointind)/60,sigcurve_Remote(fitpointind),'.black','MarkerSize',20)
    plot(T_interp/60,Remote_interp,'-black','LineWidth',3)
    %
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
    Remote_interp=interp1(timein,sigcurve_Remote,T_interp,'spline');
    plot(timein(fitpointind)/60,sigcurve_Remote(fitpointind)/Remote_interp(1),'.black','MarkerSize',20)
    plot(T_interp/60,Remote_interp/Remote_interp(1),'-black','LineWidth',3)
    
    %% 2CXM bound model
    iters = 10;
    [Nx, Ny] = size(HeartMask);
    x0_initial = [0.5 0.1 10 20];
    x0c = repmat(x0_initial, [sum(HeartMask(:)) 1, iters]);
    %lb=[-3 -3 -3 -3];
    lb = [0 -1 0 0];
    ub = [5 1 100 100];
    fitting_results = zeros(Nx, Ny, 6);
    for tvit = 1:iters
        for i = 1:sum(HeartMask(:))
            sigin = Ctoi(i,:);
            x0 = x0c(i, :, tvit);
            try
                x0 = x0c(i, :, tvit);
                [F, Ps, Vp, Ve, pred, rmse, Rsq] = fit_cxm_bound(sigin',Cp',time',x0,lb,ub);
            catch
                x0 = x0c(i, :, 1);
                [F, Ps, Vp, Ve, pred, rmse, Rsq] = fit_cxm_bound(sigin',Cp',time',x0,lb,ub);
            end
            fitresultsvectorDECcxm(i,:) = [F, Ps, Vp, Ve, rmse, Rsq];
        end
        for n = 1:size(fitting_results, 3)
            temp = zeros(size(fixed));
            temp(HeartMask) = fitresultsvectorDECcxm(:,n);
            fitting_results(:,:,n)=temp;
        end
    
        %% Total Variation
        lambda = 0.001;
        F_cxm = reduceTV(fitting_results(:, :, 1), lambda, 100, 0.01);
        Ps_cxm = reduceTV(fitting_results(:, :, 2), lambda, 100, 0.01);
        Vp_cxm = reduceTV(fitting_results(:, :, 3), lambda, 100, 0.01);
        Ve_cxm = reduceTV(fitting_results(:, :, 4), lambda, 100, 0.01);
        Rsq_cxm = fitting_results(:, :, 5);
        x0c(:, :, tvit+1) = [F_cxm(HeartMask) Ps_cxm(HeartMask) Vp_cxm(HeartMask) Ve_cxm(HeartMask)];

        figure;title(['T',num2str(timepoints),'min Iter', num2str(tvit)])
        subplot(3,2,1);imagesc(F_cxm); axis off equal; title(['F_cxm Con#=T',num2str(timepoints),'min']); clim([0 1]);
        subplot(3,2,2);imagesc(Ps_cxm); axis off equal; title(['PS_cxm Con#=T',num2str(timepoints),'min']);clim([0 1]);
        subplot(3,2,3);imagesc(Vp_cxm); axis off equal; title(['ve_cxm Con#=T',num2str(timepoints),'min']);clim([0 50]);
        subplot(3,2,4);imagesc(Ve_cxm);axis off equal; title(['vp_cxm Con#=T',num2str(timepoints),'min']);clim([0 10]);
        
        subplot(3,2,5);imagesc(Rsq_cxm);axis off equal;title(['Rsq_cxm Con#=T',num2str(timepoints),'min']);clim([0 1]);
        savefig(gcf, fullfile(outputfolder, ['T',num2str(timepoints),'min Iter', num2str(tvit), '.fig']))
        saveas(gcf, fullfile(outputfolder, ['T',num2str(timepoints),'min Iter', num2str(tvit), '.png']))
    end

    % simulate the R1
    timesim=0:60:1800;
    try
    [tempCpexp2,tempfitresultsvector]=fit([time(2:end)'; 6*60*60],[Cp(2:end)';0],'exp2')
    catch
        disp('The fitting failed, try add starting points')
        [tempCpexp2,tempfitresultsvector]=fit([time(2:end)'; 6*60*60],[Cp(2:end)';0],'exp2', 'StartPoint',[0,0,0,0])%Fit Cp adding 0 as the last
    end
    Cpsim=tempCpexp2(timesim);%Cp at 1000sec
    Cpsim(1)=0;
    tempsimLateR1vectorDEC = zeros(size(fixed));
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

    %% Blood pool mask from CXM f maps
    Remote_Fmean=mean(F_cxm(RemoteMask));
    Remote_Fsd=std(F_cxm(RemoteMask));
    Bloodthreshold=Remote_Fmean+Remote_Fsd*10;
    %Bloodthreshold=mean(F_cxm(Blood))-5*std(F_cxm(Blood));
    CXM_BloodMask=F_cxm>Bloodthreshold;
    
    % Use the top threshold to create blood pool mask 
    simLateR1cxm_Masked=(simLateR1cxm).*(-1*repmat(CXM_BloodMask,...
        [1 1 size(simLateR1cxm,3)])+0.5).*(simLateR1cxm>0).*repmat(HeartMask,[1 1 size(simLateR1cxm,3)]);
    
    % save the results
    DerivedDCEMaps.F_cxm=F_cxm.*HeartMask;
    DerivedDCEMaps.vp_cxm=Vp_cxm.*HeartMask;
    DerivedDCEMaps.PS_cxm=Ps_cxm.*HeartMask;
    DerivedDCEMaps.ve_cxm=Ve_cxm.*HeartMask;
    DerivedDCEMaps.Rsq_cxm=Rsq_cxm.*HeartMask;
    DerivedDCEMaps.simLateR1cxm=simLateR1cxm;
    DerivedDCEMaps.simLateR1cxm_Masked=simLateR1cxm_Masked;

    timepoints=uint16(ceil(time(2:end)/60))
    imagedir=[outputfolder filesep 'T', num2str(timepoints),'min' filesep 'Images'];
    
    mkdir(imagedir);
    dicominfo = DataPost_moco(1).info
    %
    if SeriesNum<10
    % SeriesNum=savedicom(ECV,imagedir,dicominfo,SeriesNum,['ECV'])
    temp=dR1(:,:,1);
    SeriesNum=savedicom(dR1/median(temp(Blood))*100,imagedir,dicominfo,SeriesNum,['dR1'])
    end
    vp_cxm_w=DerivedDCEMaps.vp_cxm;
    F_cxm_w=DerivedDCEMaps.F_cxm*100;
    PS_cxm_w=DerivedDCEMaps.PS_cxm*100;
    ve_cxm_w=DerivedDCEMaps.ve_cxm;
    simLateR1cxm_w=DerivedDCEMaps.simLateR1cxm*10;
    
    SeriesNum=savedicom(F_cxm_w,imagedir,dicominfo,SeriesNum,['cxm_F_T',num2str(timepoints),'min'])
    SeriesNum=savedicom(vp_cxm_w,imagedir,dicominfo,SeriesNum,['cxm_Vp_T',num2str(timepoints),'min'])
    SeriesNum=savedicom(PS_cxm_w,imagedir,dicominfo,SeriesNum,['cxm_PS_T',num2str(timepoints),'min'])
    SeriesNum=savedicom(ve_cxm_w,imagedir,dicominfo,SeriesNum,['cxm_Ve_T',num2str(timepoints),'min'])
    SeriesNum=savedicom(Rsq_cxm,imagedir,dicominfo,SeriesNum,['cxm_Rsq_T',num2str(timepoints),'min'])
    SeriesNum=savedicom(simLateR1cxm_w,imagedir,dicominfo,SeriesNum,['SimR1cxm_T',num2str(timepoints),'min'])%in sec
    %savedicom(idximcxm,imagedir,dicominfo,'Kmeancluster_cxm')

    save([imagedir filesep 'workspace T',num2str(timepoints),'min']) 
end
end


function [t1_map_reduced] = reduceTV(t1_map, lambda, nIterations, dt)
% reduceTV performs Total Variation reduction on a 2D T1 map
%
% Inputs:
%   t1_map       - The original 2D T1 map (Ny by Nx matrix)
%   lambda       - Regularization parameter controlling the strength of TV reduction
%   nIterations  - Number of iterations for the gradient descent process
%   dt           - Time step for the gradient descent update

% Initialize the reduced T1 map with the original map
t1_map_reduced = t1_map;

[Ny, Nx] = size(t1_map); % Get the dimensions of the T1 map

for iter = 1:nIterations
    % Calculate gradients
    gradX = [diff(t1_map_reduced, 1, 2), zeros(Ny, 1)]; % Gradient along x-axis
    gradY = [diff(t1_map_reduced, 1, 1); zeros(1, Nx)]; % Gradient along y-axis

    % Calculate divergence from the gradients
    divX = [zeros(Ny, 1), -diff(gradX, 1, 2)];
    divY = [zeros(1, Nx); -diff(gradY, 1, 1)];
    div = divX + divY;

    % Calculate the gradient of the TV term
    g = -div;

    % Update the T1 map using gradient descent
    t1_map_reduced = t1_map_reduced - dt * lambda * g;
end
end

function [F, Ps, Vp, Ve, pred, rmse, Rsq] = fit_cxm_bound(toi, Cp, time, x0, lb, ub)

    % Make sure all vectors are column, not row
    if size(toi,1) == 1
        toi = toi';
    end
    if size(Cp,1) == 1
        Cp = Cp';
    end
    if size(time,1) == 1
        time = time';
    end

    % Make sure data vectors are same size
    if size(toi,1) ~= size(time,1) || size(toi,1) ~= size(Cp,1) || size(Cp,1) ~= size(time,1)
       error('dce_mri_fit:toi,Cp,and time vectors must all be same size'); 
    end
    
    % Set options for non-linear LSQ fitting
    S=optimset; S.Algorithm='trust-region-reflective'; S.Display='off';
    %S.TolFun=1e-3; S.TolX=1e-3;% RY change fitting threshold
    S.TolFun=1e-10; S.TolX=1e-10;% RY change fitting threshold
    S.MaxIter=1000;

    % Perform the fitting (the function "rrm" is in the nested functions below)
    [B,~,~] = lsqcurvefit(@cxm_b,x0,[time,Cp],toi,lb,ub,S);

    % Store each parameter (see reference)
    F = B(1);    %F
    Ps = B(2);    %PS
    Vp = B(3);    %Vp
    Ve = B(4);    %Ve
    pred = cxm_b(B,[time,Cp]);
    %pars(5,1)=rsquare(toi,pred);%ry
    [Rsq, rmse] = rsquare(toi,pred);%ry
end

function C_toi = cxm_b(beta,X)
    
    time=X(:,1);
    Cp = X(:,2);
    
    F = beta(1); %liter/(sec*ml)
    PS = beta(2); %1/t=sec^-1
    vp=beta(3); % 1/100
    ve=beta(4); %1/100
    a=(F+PS)/vp;
    b=PS/ve;
    c=F/vp;
    M1=0.5*(a+b+sqrt((a+b).^2-4*b*c));
    M2=0.5*(a+b-sqrt((a+b).^2-4*b*c));
    B=(M2-c)/(M2-M1);
    H=B*exp(-M1*time)+(1-B)*exp(-M2*time);
    %C_toi = F*convolution(H,Cp);
    C_toi = F*convolution(H,Cp);
end

function [r2, rmse] = rsquare(y,f,varargin)
% Compute coefficient of determination of data fit model and RMSE
%
% [r2 rmse] = rsquare(y,f)
% [r2 rmse] = rsquare(y,f,c)
%
% RSQUARE computes the coefficient of determination (R-square) value from
% actual data Y and model data F. The code uses a general version of 
% R-square, based on comparing the variability of the estimation errors 
% with the variability of the original values. RSQUARE also outputs the
% root mean squared error (RMSE) for the user's convenience.
%
% Note: RSQUARE ignores comparisons involving NaN values.
% 
% INPUTS
%   Y       : Actual data
%   F       : Model fit
%
% OPTION
%   C       : Constant term in model
%             R-square may be a questionable measure of fit when no
%             constant term is included in the model.
%   [DEFAULT] TRUE : Use traditional R-square computation
%            FALSE : Uses alternate R-square computation for model
%                    without constant term [R2 = 1 - NORM(Y-F)/NORM(Y)]
%
% OUTPUT 
%   R2      : Coefficient of determination
%   RMSE    : Root mean squared error
%
% EXAMPLE
%   x = 0:0.1:10;
%   y = 2.*x + 1 + randn(size(x));
%   p = polyfit(x,y,1);
%   f = polyval(p,x);
%   [r2 rmse] = rsquare(y,f);
%   figure; plot(x,y,'b-');
%   hold on; plot(x,f,'r-');
%   title(strcat(['R2 = ' num2str(r2) '; RMSE = ' num2str(rmse)]))
%   
% Jered R Wells
% 11/17/11
% jered [dot] wells [at] duke [dot] edu
%
% v1.2 (02/14/2012)
%
% Thanks to John D'Errico for useful comments and insight which has helped
% to improve this code. His code POLYFITN was consulted in the inclusion of
% the C-option (REF. File ID: #34765).

if isempty(varargin); c = true; 
elseif length(varargin)>1; error 'Too many input arguments';
elseif ~islogical(varargin{1}); error 'C must be logical (TRUE||FALSE)'
else c = varargin{1}; 
end

% Compare inputs
if ~all(size(y)==size(f)); error 'Y and F must be the same size'; end

% Check for NaN
tmp = ~or(isnan(y),isnan(f));
y = y(tmp);
f = f(tmp);

if c; r2 = max(0,1 - sum((y(:)-f(:)).^2)/sum((y(:)-mean(y(:))).^2));
else r2 = 1 - sum((y(:)-f(:)).^2)/sum((y(:)).^2);
    if r2<0
    % http://web.maths.unsw.edu.au/~adelle/Garvan/Assays/GoodnessOfFit.html
        warning('Consider adding a constant term to your model') %#ok<WNTAG>
        r2 = 0;
    end
end

rmse = sqrt(mean((y(:) - f(:)).^2));

end