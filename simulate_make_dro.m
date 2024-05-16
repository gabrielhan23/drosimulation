% Simulate the four parameters, vp, ve, ps, flow
% depend on Classdef SimulatedDRO and CXM_bound
% Mona, April 2024
clear all
addpath(genpath('.'))
output = '/Users/mona/Library/CloudStorage/Dropbox/0.MAC-SYNC/0.PROJECT/DCE/Data/SimulatedDRO';
dims = [40, 120, 3];

% simulation
t_start = (0:60:600);
DRO = SimulatedDRO();
DRO.snr = 17.5;
DRO.t_start = 0;
DRO.t_end = 1200./60;
DRO.t_intval = 0.02; % mins
DRO.ps_slice = 2;
[timepoints_all, aif_all, myo_curves_all, myo_pars] = DRO.SIMULATE_2CXM(dims);
for s_idx = 1:length(t_start)
    t_end = (t_start(s_idx)+90:30:1200);
    for e_idx = 1:length(t_end)
        start_index = t_start(s_idx)/60/DRO.t_intval;
        end_index = t_end(e_idx)/60/DRO.t_intval;
        timepoints = timepoints_all(start_index+1:end_index);
        aif = aif_all(start_index+1:end_index);
        myo_curves = myo_curves_all(start_index+1:end_index,:,:);

        name = sprintf('start%.1f_end%.1f_interval%.2f_snr%.2f', t_start(s_idx)/60, t_end(e_idx)/60, DRO.t_intval, DRO.snr);
        fprintf("Processing simulation %s...\n", name)
        if exist(fullfile(output, [name, '_gt_pred.png']), "file")
            fprintf("Already calculated\n")
            continue
        end
        % plot the aif and myo_curves
        figure,
        h1 = plot(timepoints, aif, '-o'); hold on
        h2 = plot(timepoints, myo_curves(:, end/2, end/2), '-o'); hold on
        leg_all = legend([h1 h2], 'AIF', 'Myo curves', 'fontsize',12, 'FontName', 'Arial', ...
            'Location', 'northeast', 'FontWeight', 'bold');
        
        saveas(gcf, fullfile(output, [name, '_aif_myo_curve.png']))
        %
        CXMB = CXM_bound();
        mask = ones(dims(1), dims(2));
        [fitresultsDCEcxm, simulatedCXM] = CXMB.SOLVER_CXM_BOUND_ITERATIONS(aif, myo_curves, mask, mask, timepoints);
        % plot the results
        F_cxm = fitresultsDCEcxm(:,:,1);
        PS_cxm = fitresultsDCEcxm(:,:,2);
        Vp_cxm = fitresultsDCEcxm(:,:,3);
        Ve_cxm = fitresultsDCEcxm(:,:,4);
        Rmse_cxm = fitresultsDCEcxm(:,:,5);
        Rsq_cxm = fitresultsDCEcxm(:,:,6);


        myo_pars = squeeze(myo_pars);

        % demonstrate the results
        CMRmap=[0 0 0;.15 .15 .5;.3 .15 .75;.6 .2 .50;1 .25 .15;.9 .5 0;.9 .75 .1;.9 .9 .5;1 1 1];

        t = tiledlayout(4,2);
        width = 300;    % Replace with your desired figure width
        height = 800;   % Replace with your desired figure height
        fig = gcf;      % Get the current figure handle
        set(fig, 'Position', [100, 100, width, height]);
        title(t, sprintf('start%.1f end%.1f interval%.2f snr%.2f', DRO.t_start, DRO.t_end, DRO.t_intval, DRO.snr))
        labels = {'Flow', 'PS', 'Vp', 'Ve'};
        vmax = [3.0, 4.0, 0.24, 0.60];
        for i = 1:4
            pred = fitresultsDCEcxm(:,:,i)';
            gt = myo_pars(:, :, i)';
            MSE = immse(gt, pred);
            SSIM = ssim(pred, gt);
            nexttile
            h1 = imshow(gt, [0, vmax(i)]); colormap(CMRmap);
            colorbar
            title(labels{i})
            nexttile
            h2 = imshow(pred, [0, vmax(i)]); colormap(CMRmap);
            title(sprintf("mse: %.3f\nssim: %.3f", MSE, SSIM))
            colorbar
            t.Padding = 'none';
            t.TileSpacing = 'tight';
        end
        saveas(gcf, fullfile(output, [name, '_gt_pred.png']))

        % Save the simulated DRO and predicted DRO
        mkdir(fullfile(output, name))
        CXMB.mat2Nifti(permute(myo_curves, [2, 3, 1]), fullfile(output, name, 'SimulationDRO.nii'), [1, 1, 1]);
        CXMB.mat2Nifti(permute(simulatedCXM, [2, 3, 1]), fullfile(output, name, 'PredictedDRO.nii'), [1, 1, 1]);
        close all
    end
end






