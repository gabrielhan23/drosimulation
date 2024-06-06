% Test different parameter
clear all
addpath(genpath('.'))
% output = '/Users/mona/Library/CloudStorage/Dropbox/0.MAC-SYNC/0.PROJECT/DCE/Data/SimulatedDRO/Tryparameters';
% mkdir(output)
resume = false;

dims = [40, 120, 3];

% simulation
t_start = (0:60:600);
DRO = SimulatedDRO();
DRO.snr = 0;
DRO.t_start = 0;
DRO.t_end = 1200./60;
DRO.t_intval = 0.05; % mins
DRO.ps_slice = 1;

% large vp, low ve
flow = 0.03;
PS = 0.001;
vp = 0.02;
ve = 1e-6;

for i = 1:length(ve)
DRO.flow(1) = flow;
DRO.ps(1) = PS;
DRO.vp(1) = vp;
DRO.ve(1) = ve(i);
[timepoints_all, aif_all, myo_curves_all, myo_pars] = DRO.SIMULATE_2CXM(dims);
figure,
h1 = plot(timepoints_all, aif_all, '-o'); hold on
h2 = plot(timepoints_all, myo_curves_all(:, 1, 1), '-o'); hold on
leg_all = legend([h1 h2], 'AIF', 'Myo curves', 'fontsize',12, 'FontName', 'Arial', ...
    'Location', 'northeast', 'FontWeight', 'bold');

title(sprintf("Flow:%.3f, PS:%.3f, Vp:%.3f, Ve:%.3f", myo_pars(1, 1, 1), myo_pars(1, 1, 2), myo_pars(1, 1, 3), myo_pars(1, 1, 4)));
end
saveas(gcf, fullfile(output, sprintf("Flow%.3fPS%.3fVp%.4fVe%.4f_aif_myo.png", myo_pars(1, 1, 1), myo_pars(1, 1, 2), myo_pars(1, 1, 3), myo_pars(1, 1, 4))))
%%
dims = [1, 1, 1];

t_start = (0:60:600);
DRO.t_start = 0;
DRO.t_end = 1200./60;
DRO.t_intval = 0.05; % mins

for s_idx = 2:length(t_start)
    t_end = (t_start(s_idx)+90:60:1200);
    for e_idx = length(t_end)
        start_index = t_start(s_idx)/60/DRO.t_intval;
        end_index = t_end(e_idx)/60/DRO.t_intval;
        timepoints = timepoints_all(start_index+1:end_index);
        aif = aif_all(start_index+1:end_index);
        myo_curves = myo_curves_all(start_index+1:end_index,1,1);
        % myo_curves = reshape(myo_curves, [length(myo_curves), 1, 1]);

        name = sprintf('start%.1f_end%.1f_interval%.2f_snr%.2f', t_start(s_idx)/60, t_end(e_idx)/60, DRO.t_intval, DRO.snr);
        fprintf("Processing simulation %s...\n", name)
        if exist(fullfile(output, [name, '_gt_pred.png']), "file") && resume
            fprintf("Already calculated\n")
            continue
        end
        % plot the aif and myo_curves
        figure,
        h1 = plot(timepoints, aif, '-o'); hold on
        h2 = plot(timepoints, myo_curves(:), '-o'); hold on
        leg_all = legend([h1 h2], 'AIF', 'Myo curves', 'fontsize',12, 'FontName', 'Arial', ...
            'Location', 'northeast', 'FontWeight', 'bold');
        title(sprintf("Flow:%.1f, PS:%.1f, Vp:%.2f, Ve:%.2f", myo_pars(1, 1, 1), myo_pars(1, 1, 2), myo_pars(1, 1, 3), myo_pars(1, 1, 4)));
        
        % saveas(gcf, fullfile(output, [name, '_aif_myo_curve.png']))
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
        
        MSE = immse(myo_curves, simulatedCXM);
        SSIM = ssim(myo_curves, simulatedCXM);
        figure,
        width = 800;    % Replace with your desired figure width
        height = 300;   % Replace with your desired figure height
        fig = gcf;      % Get the current figure handle
        set(fig, 'Position', [100, 100, width, height]);
        h1 = plot(timepoints, aif, 'o'); hold on
        h2 = plot(timepoints, myo_curves(:), 'o'); hold on
        h3 = plot(timepoints, simulatedCXM(:), '-'); hold on
        leg_all = legend([h1 h2 h3], 'AIF', 'Myo curves', 'Simulated Myo curves', 'fontsize',12, 'FontName', 'Arial', ...
            'Location', 'northeast', 'FontWeight', 'bold');
        title(sprintf("Flow:%.4f, PS:%.4f, Vp:%.4f, Ve:%.4f, Rsq:%.4f, MSE:%4f, SSIM:%4f\nFlow:%.4f, PS:%.4f, Vp:%.4f, Ve:%.4f", ...
            F_cxm, PS_cxm, Vp_cxm, Ve_cxm, Rsq_cxm, MSE, SSIM, myo_pars(1, 1, 1), myo_pars(1, 1, 2), myo_pars(1, 1, 3), myo_pars(1, 1, 4)));
        saveas(gcf, fullfile(output, [name, '_aif_myo_simulated_curve.png']))


        myo_pars = squeeze(myo_pars);

        % demonstrate the results
        % CMRmap=[0 0 0;.15 .15 .5;.3 .15 .75;.6 .2 .50;1 .25 .15;.9 .5 0;.9 .75 .1;.9 .9 .5;1 1 1];
        % 
        % t = tiledlayout(4,2);
        % width = 300;    % Replace with your desired figure width
        % height = 800;   % Replace with your desired figure height
        % fig = gcf;      % Get the current figure handle
        % set(fig, 'Position', [100, 100, width, height]);
        % title(t, sprintf('start%.1f end%.1f interval%.2f snr%.2f', DRO.t_start, DRO.t_end, DRO.t_intval, DRO.snr))
        % labels = {'Flow', 'PS', 'Vp', 'Ve'};
        % vmax = [0.01, 0.05, 1, 1];
        % for i = 1:4
        %     pred = fitresultsDCEcxm(:,:,i)';
        %     gt = myo_pars(1, 1, i)';
        %     MSE = immse(gt, pred);
        %     SSIM = ssim(pred, gt);
        %     nexttile
        %     h1 = imshow(gt, [0, vmax(i)]); colormap(CMRmap);
        %     colorbar
        %     title(labels{i})
        %     nexttile
        %     h2 = imshow(pred, [0, vmax(i)]); colormap(CMRmap);
        %     title(sprintf("mse: %.3f\nssim: %.3f", MSE, SSIM))
        %     colorbar
        %     t.Padding = 'none';
        %     t.TileSpacing = 'tight';
        % end
        % saveas(gcf, fullfile(output, [name, '_gt_pred.png']))

        % Save the simulated DRO and predicted DRO
        mkdir(fullfile(output, name))
        CXMB.mat2Nifti(permute(myo_curves, [2, 3, 1]), fullfile(output, name, 'SimulationDRO.nii'), [1, 1, 1]);
        CXMB.mat2Nifti(permute(simulatedCXM, [2, 3, 1]), fullfile(output, name, 'PredictedDRO.nii'), [1, 1, 1]);
        save(fullfile(output, name, 'workspace.mat'))
        close all
    end
end                     