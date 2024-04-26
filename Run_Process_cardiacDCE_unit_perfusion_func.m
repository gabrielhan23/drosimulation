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


dR1_aif = (PostT1Reg_aif-PreT1_aif)./PreT1_aif .*1e3;
Gdcon_aif = dR1_aif ./ gdrelaxivity; %mM
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


avg_Blood = squeeze(mean(Gdcon .* repmat(Blood, [1, 1, size(Gdcon, 3)]), 1:2));

avg_AIF = squeeze(mean(Gdcon_aif .* repmat(AIF, [1, 1, size(Gdcon_aif, 3)]), 1:2));
scale = avg_Blood(100)/avg_AIF(100);
Gdcon_aif = Gdcon_aif .* scale;
%% Plot the timepoints

avg_AIF_rescale = squeeze(mean(Gdcon_aif .* repmat(AIF, [1, 1, size(Gdcon_aif, 3)]), 1:2));
figure,
h1 = plot(timein, avg_Blood, '-bo'); hold on
h2 = plot(timein, avg_AIF_rescale, '-go');
h3 = plot(timein, avg_AIF, '-ro');
legend([h1, h2, h3], 'HR Blood', 'LR AIF Rescale', 'LR AIF Original');
%% Create pairs
temporal_ratio = timein(2) - timein(1);
index_interval_1 = round(10 / temporal_ratio);
index_interval_2 = round(20 / temporal_ratio);

t_start_indexs = [1, 36, 54, 73];
t_end_indexs = [116, 149, 187, 255, 262, 301, 340, 378, 384];

CXMB = CXM_bound();
CXMB.debug = true;
CXMB.x0 = [0.5 5 500 500];
CXMB.lb = [0 -10 0 0];
CXMB.ub = [50 10 1000 1000];
for i = 1:length(t_start_indexs)
    for j = 1:length(t_end_indexs)
        t_start_index = t_start_indexs(i);
        t_end_index = t_end_indexs(j);
        indexs = t_start_index:index_interval_1:t_end_index;
        output_params = fullfile(outputfolder, sprintf('t_%d_%dsec', round(timein(t_start_index)), round(timein(t_end_index))));
        mkdir(output_params)
        if exist(fullfile(output_params, 'PS.nii'), 'file')
            disp("The file existed.")
        else
            AIF_mask = repmat(AIF, [1, 1, size(Gdcon_aif(:,:,indexs), 3)]);
            aif = squeeze(mean(Gdcon_aif(:,:,indexs) .* AIF_mask, 1:2));
            myo_curves = Gdcon(:, :, indexs) .* repmat(HeartMask, [1, 1, size(Gdcon(:,:,indexs), 3)]);
            
            timepoints = timein(indexs)';

            if length(timepoints) < 4
                disp("Time range too small")
                continue
            end
            [fitresultsDCEcxm, simulatedCXM] = CXMB.SOLVER_CXM_BOUND_ITERATIONS( ...
                squeeze(aif), permute(myo_curves, [3, 1, 2]), HeartMask, Myomask, timepoints);
            saveParameters(fitresultsDCEcxm, output_params)
    
            mat2Nifti(simulatedCXM, fullfile(output_params, 'simLateR1_cxm_bound.nii'), [1, 1, 1])

            save(fullfile(output_params, 'workspace'))
        end

    end
end

function [] = saveParameters(params, output)
    Flow_cxm_bound=params(:,:,1);
    PS_cxm_bound=params(:,:,2);
    Vp_cxm_bound=params(:,:,3);
    Ve_cxm_bound=params(:,:,4);
    rmse_cxm_bound=params(:,:,5);
    Rsq_cxm_bound=params(:,:,6);
    
    mat2Nifti(Flow_cxm_bound, fullfile(output, 'Flow.nii'), [1, 1, 1]);
    mat2Nifti(PS_cxm_bound, fullfile(output, 'PS.nii'), [1, 1, 1]);
    mat2Nifti(Vp_cxm_bound, fullfile(output, 'vp.nii'), [1, 1, 1]);
    mat2Nifti(Ve_cxm_bound, fullfile(output, 've.nii'), [1, 1, 1]);
    mat2Nifti(rmse_cxm_bound, fullfile(output, 'Rmse.nii'), [1, 1, 1]);
    mat2Nifti(Rsq_cxm_bound, fullfile(output, 'Rssq.nii'), [1, 1, 1]);
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
