% Simulate the four parameters, vp, ve, ps, flow
% depend on Classdef SimulatedDRO and CXM_bound
% Mona, April 2024
clear all

output = '/Users/mona/Library/CloudStorage/Dropbox/0.MAC-SYNC/0.PROJECT/DCE/Data/SimulatedDRO';
dims = [40, 120, 3];


DRO = SimulatedDRO();
DRO.snr = 17.5;
DRO.t_start = 0;
DRO.t_end = 2;
DRO.t_intval = 0.02;
DRO.ps_slice = 2;
[timepoints, aif, myo_curves, myo_pars] = DRO.SIMULATE_2CXM(dims);

% plot the aif and myo_curves
figure, 
h1 = plot(timepoints, aif, '-o'); hold on
h2 = plot(timepoints, myo_curves(:, end/2, end/2), '-o'); hold on 
leg_all = legend([h1 h2], 'AIF', 'Myo curves', 'fontsize',12, 'FontName', 'Arial', ...
    'Location', 'northeast', 'FontWeight', 'bold');

%%
CXMB = CXM_bound();
mask = ones(dims(1), dims(2));
[fitresultsDCEcxm, simulatedCXM] = CXMB.SOLVER_CXM_BOUND_ITERATIONS(aif, myo_curves, mask, mask, timepoints);
%% plot the results
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
saveas(gcf, fullfile(output, 'gt_pred.png'))

% Save the simulated DRO and predicted DRO

CXMB.mat2Nifti(permute(myo_curves, [2, 3, 1]), fullfile(output, 'SimulationDRO.nii'), [1, 1, 1]);
CXMB.mat2Nifti(permute(simulatedCXM, [2, 3, 1]), fullfile(output, 'PredictedDRO.nii'), [1, 1, 1]);



