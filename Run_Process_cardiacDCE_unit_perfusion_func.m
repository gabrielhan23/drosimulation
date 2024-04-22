% For DCE with cxmmix model, using the total variation
% Mona, March 25, 2024

% function [] = Run_Process_cardiacDCE_unit_perfusion_func(subjectfolder, output, subject, label)
 % load the data
folder = [subject filesep 'DCE' filesep label];

outputfolder = fullfile(output, folder);
mkdir(outputfolder);

Load_Dicom_Perfusion

%% concentration
PostT1Reg(PostT1Reg == 0) = PostT1Reg(PostT1Reg == 0) + 1e-3; % zero exist
PostT1Reg_aif(PostT1Reg_aif == 0) = PostT1Reg_aif(PostT1Reg_aif == 0) + 1e-3;
PreT1(PreT1 == 0) = PreT1(PreT1 == 0) + 1e-3; % zero exist
PreT1_aif(PreT1_aif == 0) = PreT1_aif(PreT1_aif == 0) + 1e-3; % zero exist
dR1 = (PostT1Reg-PreT1)./PreT1 .*1e3; %s^-1 
gdrelaxivity = 4.97; %s^-1*mM^-1
Gdcon = dR1./gdrelaxivity; %mM

scale = 1/10;
dR1_aif = (PostT1Reg_aif-PreT1_aif)./PreT1_aif .*1e3;
Gdcon_aif = scale .* dR1_aif ./ gdrelaxivity; %mM
%% Draw Reference ROI

if ~exist('Blood','var')
figure;imagesc(PreT1, [0 , 200]);
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

figure;imagesc(PostT1Reg_aif(:,:,1), [0 , 500]);
title(['blood pool'])
AIF=roipoly;

mkdir(fullfile(outputfolder, 'MASK'))
save(fullfile(outputfolder, 'MASK', [label, '_ROI.mat']), 'AIF', 'Blood', 'HeartMask', 'LVendo', 'Myomask', 'MIMask', 'RemoteMask')
end
%% mask avoid singularity from registration error
% HeartMask= HeartMask.*(min(dR1,[],3)>=0).*(max(dR1,[],3)<Inf);
% HeartMask=HeartMask.*~isnan(HeartMask);
% HeartMask=HeartMask>0;
% 
% Blood=Blood.*(min(dR1,[],3)>=0).*(max(dR1,[],3)<Inf);
% Blood=Blood.*~isnan(Blood);
% 
% Blood=(Blood>0);

%% Plot the timepoints
figure, plotdR1vsTime_aif(Gdcon_aif, AIF, timein);
figure, plotdR1vsTime(Gdcon, Blood, Myomask, timein);

%% Create pairs
temporal_ratio = timein(2) - timein(1);
index_interval_1 = round(10 / temporal_ratio);
index_interval_2 = round(20 / temporal_ratio);

t_start_indexs = 1:index_interval_1:round(3*60/temporal_ratio);
t_end_indexs = round(3*60/temporal_ratio):index_interval_2:round(10*60/temporal_ratio);


for i = 1:length(t_start_indexs)
    for j = 1:length(t_end_indexs)
        t_start_index = t_start_indexs(i);
        t_end_index = t_end_indexs(j);
        indexs = t_start_index:index_interval_1:t_end_index;
        output_params = fullfile(outputfolder, sprintf('t_%d_%d', t_start_index, t_end_index));
        mkdir(output_params)
        if exist(fullfile(output_params, 'PS.nii'), 'file')
            disp("The file existed.")
        else
            
            aif = mean(Gdcon_aif(:, :, indexs), 1:2);
            myo_curves = Gdcon(:, :, indexs);
            x0=[0.5 0.1 10 20];
            lb=[0 -1 0 0];
            ub=[5 1 100 100];
            timepoints = timein(indexs);
            if length(timepoints) < 4
                disp("Time range too small")
                continue
            end
            
            [fitresultsDCEcxm, simLateR1] = Double_cxm(squeeze(aif), myo_curves, HeartMask, Myomask, timepoints, x0, lb, ub, 0.85, 1);
            
            saveParameters(fitresultsDCEcxm, output_params)
    
            mat2Nifti(simLateR1, fullfile(output, 'simLateR1_cxm_bound.nii'), [1, 1, 1])
        end

    end
end

function [] = saveParameters(params, output)
    PS_cxm_bound=params(:,:,1);% -1 - 1
    vp_cxm_bound=params(:,:,2);% 0-100
    ve_cxm_bound=params(:,:,3);% 0-100
    Rsq_cxm_bound=params(:,:,end);

    mat2Nifti(PS_cxm_bound, fullfile(output, 'PS.nii'), [1, 1, 1]);
    mat2Nifti(vp_cxm_bound, fullfile(output, 'vp.nii'), [1, 1, 1]);
    mat2Nifti(ve_cxm_bound, fullfile(output, 've.nii'), [1, 1, 1]);
    mat2Nifti(Rsq_cxm_bound, fullfile(output, 'Rssq.nii'), [1, 1, 1]);
end

function [] = mat2Nifti(volume, savepath, voxelSize)
% save Nifti
% reference https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image
temp_nii = make_nii(volume);
temp_nii.hdr.dime.pixdim(2:4) = voxelSize;

save_nii(temp_nii, savepath);
disp(["Suceess to save the nifti files" + savepath])
end

function [fitresultsDCEcxm, simLateR1] = Double_cxm(aif, myo_curves, mask, Myomask, time, x0, lb, ub, threshold, iters)
NityFivpercent=0;

[rows, cols, len] = size(myo_curves);

fitresultsDCEcxm = zeros(rows, cols, 6);
Rsqc=-1*ones([rows*cols 1]);
x0c=repmat(x0,[rows*cols 1]);
% cp: blood pool concentration mean
% sigin: ctoi, concentration within the heartmask
iter = 0;
while(NityFivpercent<threshold & iter<iters) % initial guess optimization
    for i = 1:rows
        for j = 1:cols
            if mask(i, j) == 1
                c = (i-1) * rows + j;
                if Rsqc(c)<0.9   
                    sigin = squeeze(myo_curves(i, j, :));
                    cp = aif;
                    tempfitDCE = fitdcemri_etk(sigin, cp,time,x0,lb,ub,'cxm_bound');
                    fitresultsDCEcxm(i, j, :) = tempfitDCE;
                    Rsqc(c)=tempfitDCE(end);
                    if tempfitDCE(end) > Rsqc(c)%initial guess otimization
                        x0c(c,:)=x0;
                        Rsqc(c)=tempfitDCE(end);
                    end
                    x0 = x0c(c,:).*((1-(rand(1,length(x0))-0.5)*1.5));
                end
            else
                fitresultsDCEcxm(i, j, :) = 0;
            end

        end
    end
    fitQ = fitresultsDCEcxm(:,:,end);
    maskadj= Myomask.*(fitQ>0);
    f = sort(fitQ(maskadj==1));
    NityFivpercent = f(floor(sum(maskadj(:))*0.05));% get the 95% R2
    fprintf("Iter %d ---- 95 percentage R2 in iter: %.3f\n", iter, NityFivpercent)
    iter=iter+1;
end

simLateR1 = simulatedR1(time, mask, aif, myo_curves, fitresultsDCEcxm);
end

function [simLateR1] = simulatedR1(time, mask, Cp, Ctoi, fitResults)
% simulate the R1
timesim=0:60:1800;
try
[tempCpexp2,tempfitresultsvector]=fit([time(1:end)'; 6*60*60],[Cp(1:end);0],'exp2');
catch
    disp('The fitting failed, try add starting points')
    [tempCpexp2,tempfitresultsvector]=fit([time(1:end)'; 6*60*60],[Cp(1:end);0],'exp2', 'StartPoint',[0,0,0,0]);%Fit Cp adding 0 as the last
end
Cpsim = tempCpexp2(timesim); %Cp at 1000sec

[x, y, t] = size(fitResults);
simLateR1 = zeros(x, y, size(Cpsim, 1));
for i = 1:x
    for j = 1:y
        if mask(i, j) == 1
           sigin=Ctoi(i ,j,:);
           tempsimLateR1= cxm_b(squeeze(fitResults(i, j,:)), [timesim',Cpsim]);
           
           simLateR1(i, j,:) = tempsimLateR1;
           %results=a*exp(-bt)+d*exp(-et) and gof
        end
    end
end
end

function [t] = plotdR1vsTime(data, Blood, Myo, timepoints)
Blood = repmat(Blood, [1, 1, size(data, 3)]);
avg_Blood = squeeze(mean(data .* Blood, 1:2));

Myo = repmat(Myo, [1, 1, size(data, 3)]);
avg_Myo = squeeze(mean(data .* Myo, 1:2));


% figure, errorbar(timepoints, avg, stdev)
t = tiledlayout(2,1);
nexttile
h1 = plot(timepoints, avg_Blood); hold on
leg_all = legend([h1], 'Blood','fontsize',12, 'FontName', 'Arial', ...
    'Location', 'northeast', 'FontWeight', 'bold');
nexttile
h2 = plot(timepoints, avg_Myo);
leg_all = legend([h2], 'Myo','fontsize',12, 'FontName', 'Arial', ...
    'Location', 'northeast', 'FontWeight', 'bold');
t.Padding = 'none';

end

function [t] = plotdR1vsTime_aif(data, AIF, timepoints)
AIF = repmat(AIF, [1, 1, size(data, 3)]);
avg_AIF = squeeze(mean(data .* AIF, 1:2));



% figure, errorbar(timepoints, avg, stdev)
t = tiledlayout(1,1);
nexttile
h1 = plot(timepoints, avg_AIF); hold on
leg_all = legend([h1], 'AIF','fontsize',12, 'FontName', 'Arial', ...
    'Location', 'northeast', 'FontWeight', 'bold');

t.Padding = 'none';
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