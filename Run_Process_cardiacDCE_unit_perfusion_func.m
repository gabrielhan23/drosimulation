% For DCE with cxmmix model, using the total variation
% Mona, March 25, 2024

% function [] = Run_Process_cardiacDCE_unit_perfusion_func(subjectfolder, output, subject, label)
% load the data
folder = [subject, filesep, 'DCE', filesep, label];

outputfolder = fullfile(output, folder);
mkdir(outputfolder);

Load_Dicom_Perfusion

%% concentration
PostT1Reg(PostT1Reg == 0) = PostT1Reg(PostT1Reg == 0) + 1e-3; % zero exist
PostT1Reg_aif(PostT1Reg_aif == 0) = PostT1Reg_aif(PostT1Reg_aif == 0) + 1e-3;
PreT1(PreT1 == 0) = PreT1(PreT1 == 0) + 1e-3; % zero exist
PreT1_aif(PreT1_aif == 0) = PreT1_aif(PreT1_aif == 0) + 1e-3; % zero exist
dR1 = (PostT1Reg - PreT1) ./ PreT1 .* 1e3; %s^-1
gdrelaxivity = 4.97; %s^-1*mM^-1
Gdcon = dR1 ./ gdrelaxivity; %mM


dR1_aif = (PostT1Reg_aif - PreT1_aif) ./ PreT1_aif .* 1e3;
Gdcon_aif = dR1_aif ./ gdrelaxivity; %mM

%% Draw Reference ROI

if ~exist('Blood', 'var')
    figure;
    imagesc(PreT1, [0, 200]);
    title(['blood pool'])
    Blood = roipoly;

    title(['Heart mask'])
    HeartMask = roipoly;

    title(['Myo mask_endo'])
    LVendo = roipoly;

    title(['Myo mask_epi'])
    LVepi = roipoly;
    Myomask = ~LVendo .* LVepi;

    title(['MI'])
    MIMask = roipoly;

    % lesion

    % title(['MVO'])
    % MVOMask=roipoly;
    % save(fullfile(subjectfolder, 'ROI_MVOMask.mat'))
    title(['Remote'])
    RemoteMask = roipoly;

    figure;
    imagesc(PostT1Reg_aif(:, :, 1), [0, 500]);
    title(['blood pool'])
    AIF = roipoly;

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

Blood_Gdcon = Gdcon;
Blood_Gdcon(repmat(Blood, [1, 1, size(Gdcon, 3)]) == 0) = NaN;
avg_Blood = squeeze(mean(Blood_Gdcon, 1:2, "omitnan"));
Gdcon_aif(repmat(AIF, [1, 1, size(Gdcon_aif, 3)]) == 0) = NaN;
avg_AIF = squeeze(mean(Gdcon_aif, 1:2, "omitnan"));


% add 100 and 600 
range = 10;
ind = interp1(timein,1:length(timein),100,'nearest');
avg_Blood_100 = mean(avg_Blood(ind-range:ind+range));
avg_Blood_600 = mean(avg_Blood(end-2*range:end));

ratio = (timein(end)-100) / (avg_Blood_600 - avg_Blood_100);
Gdcon_aif = Gdcon_aif .* abs(ratio);

% avg_AIF = squeeze(mean(Gdcon_aif, 1:2, "omitnan"));
% ratio2 = mean(avg_AIF(187-range:300+range)) ./ mean(avg_Blood(187-range:300+range));
% Gdcon_aif = Gdcon_aif ./ ratio2;

% scale = avg_Blood(100) / avg_AIF(100);
% Gdcon_aif = Gdcon_aif .* scale;

%% Plot the timepoints

Gdcon_aif(repmat(AIF, [1, 1, size(Gdcon_aif, 3)]) == 0) = NaN;
avg_AIF_rescale = squeeze(mean(Gdcon_aif, 1:2, "omitnan"));

Myo_Gdcon = Gdcon;
Myo_Gdcon(repmat(Myomask, [1, 1, size(Gdcon, 3)]) == 0) = NaN;
avg_Myo = squeeze(mean(Myo_Gdcon, 1:2, "omitnan"));

figure,
h1 = plot(timein, avg_Blood, '-bo'); hold on
h2 = plot(timein, avg_AIF_rescale, '-go'); hold on
h3 = plot(timein, avg_AIF, '-ro'); hold on
h4 = plot(timein, avg_Myo, '-yo');
legend([h1,h2, h3,h4], 'HR Blood', 'LR AIF Rescale', 'LR AIF Original', 'Myo');
saveas(gcf, fullfile(outputfolder, 'all_aif_myo.png'))
%% Create pairs
temporal_ratio = timein(2) - timein(1);
index_interval_1 = round(10 / temporal_ratio);
index_interval_2 = round(20 / temporal_ratio);

t_start_indexs = 1;
t_end_indexs = ind:30:length(avg_Blood);

CXMB = CXM_bound();
CXMB.debug = false;
CXMB.x0 = [0.13, 0.01, 0.017, 0.21];
CXMB.lb = [0, 0, 0, 0];
CXMB.ub = [0.5, 0.5, 0.5, 2];
for i = 1:length(t_start_indexs)
    for j = 1:length(t_end_indexs)
        t_start_index = t_start_indexs(i);
        t_end_index = t_end_indexs(j);
        indexs = t_start_index:index_interval_1:t_end_index;
        output_params = fullfile(outputfolder, sprintf('t_%d_%dsec', round(timein(t_start_index)), round(timein(t_end_index))));
        mkdir(output_params)
        if exist(fullfile(output_params, 'PS.nii'), 'file')  && resume
            disp("The file existed.")
        else

            aif = squeeze(mean(Gdcon_aif(:, :, indexs), 1:2, 'omitmissing'));
            myo_curves = Gdcon(:, :, indexs) .* repmat(HeartMask, [1, 1, size(Gdcon(:, :, indexs), 3)]);

            timepoints = timein(indexs)' ./ 60; % convert to mins

            if length(timepoints) < 4
                disp("Time range too small")
                continue
            end

            scale_x = 500;

            [fitresultsDCEcxm, simulatedCXM] = CXMB.SOLVER_CXM_BOUND_ITERATIONS( ...
                squeeze(aif)./scale_x, permute(myo_curves, [3, 1, 2])./scale_x, HeartMask, Myomask, timepoints);
            
            simulatedCXM = simulatedCXM .* scale_x;
            
            tmp_simulatedCXM = permute(simulatedCXM, [2, 3, 1]);
            tmp_simulatedCXM(repmat(Myomask, [1, 1, size(tmp_simulatedCXM, 3)]) == 0) = NaN;
            avg_tmp_simulatedCXM = squeeze(mean(tmp_simulatedCXM, 1:2, "omitnan"));
            
            myo_curves(repmat(Myomask, [1, 1, size(myo_curves, 3)]) == 0) = NaN;
            avg_Myo = squeeze(mean(myo_curves, 1:2, "omitnan"));

            figure,
            h1 = plot(timein(indexs), aif', 'bo'); hold on
            h2 = plot(timein(indexs), avg_tmp_simulatedCXM, '-g'); hold on
            h4 = plot(timein(indexs), avg_Myo, 'o');
            legend([h1, h2, h4], 'AIF', 'Simulated Curve', 'Myo');
            % savefig(gcf, fullfile(output_params, 'curve_fit.png'))
            % demonstrate the results
            CMRmap=[0 0 0;.15 .15 .5;.3 .15 .75;.6 .2 .50;1 .25 .15;.9 .5 0;.9 .75 .1;.9 .9 .5;1 1 1];
    
            figure, t = tiledlayout(2,2);
            width = 400;    % Replace with your desired figure width
            height = 400;   % Replace with your desired figure height
            fig = gcf;      % Get the current figure handle
            set(fig, 'Position', [100, 100, width, height]);
            labels = {'Flow', 'PS', 'Vp', 'Ve'};
            vmax = [0.1, 0.1, 1, 1];
            for idx = 1:4
                pred = fitresultsDCEcxm(:,:,idx)';
                nexttile
                h1 = imshow(pred, [0, vmax(idx)]); colormap(CMRmap);
                colorbar
                title(labels{idx})
                t.Padding = 'none';
                t.TileSpacing = 'tight';
            end

            CXMB.mat2Nifti(fitresultsDCEcxm(:, :, 1), fullfile(output_params, 'Flow.nii'), [1, 1, 1]);
            CXMB.mat2Nifti(fitresultsDCEcxm(:, :, 2), fullfile(output_params, 'PS.nii'), [1, 1, 1]);
            CXMB.mat2Nifti(fitresultsDCEcxm(:, :, 3), fullfile(output_params, 'vp.nii'), [1, 1, 1]);
            CXMB.mat2Nifti(fitresultsDCEcxm(:, :, 4), fullfile(output_params, 've.nii'), [1, 1, 1]);
            CXMB.mat2Nifti(fitresultsDCEcxm(:, :, 5), fullfile(output_params, 'Rmse.nii'), [1, 1, 1]);
            CXMB.mat2Nifti(fitresultsDCEcxm(:, :, 6), fullfile(output_params, 'Rssq.nii'), [1, 1, 1]);

            CXMB.mat2Nifti(simulatedCXM, fullfile(output_params, 'simLateR1_cxm_bound.nii'), [1, 1, 1])
            
            save(fullfile(output_params, 'workspace'))
        end

    end
end