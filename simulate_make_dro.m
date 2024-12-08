% Simulate the four parameters, vp, ve, ps, flow
% depend on Classdef SimulatedDRO and CXM_bound
% Mona, April 2024
function results = simulation(DRO, output, num)
    results = zeros(num, 4, 2);
    resume = true;
    reimage = true;
    t_start = (0:15:120); % previously was 
    t_start = [0];

    DRO.snr = 0;
    DRO.t_start = 0;
    DRO.t_end = 1200/60;
    DRO.t_intval = 0.05; % mins
    DRO.ps_slice = 1;
    [dims, timepoints_all, aif_all, myo_curves_all, myo_pars] = DRO.SIMULATE_2CXM();
    figure('visible','off');
    h1 = plot(timepoints_all, aif_all, '-o'); hold on
    h2 = plot(timepoints_all, myo_curves_all(:, 1, 1), '-o'); hold on
    leg_all = legend([h1 h2], 'AIF', 'Myo curves', 'fontsize', 12, 'FontName', 'Arial', ...
        'Location', 'northeast', 'FontWeight', 'bold');
    saveas(gcf, fullfile(output, 'alltimepoints_aif_myo_curve.png'))
    for s_idx = 1:length(t_start)
        t_end = (t_start(s_idx)+90:60:1200);
        t_end = t_end(1:num);
        for e_idx = 1:length(t_end)
            start_index = t_start(s_idx)/60/DRO.t_intval;
            end_index = t_end(e_idx)/60/DRO.t_intval;
            timepoints = timepoints_all(start_index+1:end_index);
            aif = aif_all(start_index+1:end_index);
            myo_curves = myo_curves_all(start_index+1:end_index,:,:);
    
            name = sprintf('start%.1f end%.1f interval%.2f snr%.2f', t_start(s_idx)/60, t_end(e_idx)/60, DRO.t_intval, DRO.snr);
            filename = sprintf('start%.1f_end%.1f_interval%.2f_snr%.2f', t_start(s_idx)/60, t_end(e_idx)/60, DRO.t_intval, DRO.snr);

            CXMB = CXM_bound();
            mask = ones(dims(1), dims(2));
            if exist(fullfile(output, [filename, '_gt_pred.png']), "file") && resume
                fprintf("Already calculated %s in %s\n", name, output)
                S = load(fullfile(output, [filename, '_data.mat']));
                results = S.results;
                if reimage == false
                    return
                end
                fitresultsDCEcxm = S.fitresultsDCEcxm;
                simulatedCXM = S.simulatedCXM;
            else
                fprintf("Processing simulation %s in %s...\n", name, output)
                [fitresultsDCEcxm, simulatedCXM] = CXMB.SOLVER_CXM_BOUND_ITERATIONS(aif, myo_curves, mask, mask, timepoints);
            end
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
            x = linspace(1, size(CMRmap, 1), 256);
    
            % Interpolate the colormap to create a smooth gradient
            smoothCMRmap = interp1(1:size(CMRmap, 1), CMRmap, x, 'linear');
            
            figure('visible','off');
            t = tiledlayout(4,3);
            width = 900;    % Replace with your desired figure width
            height = 1200;   % Replace with your desired figure height
            %fig = gcf;      % Get the current figure handle
            set(gcf, 'Position', [100, 100, width, height]);
            title(t, name)
            labels = {'Flow', 'PS', 'Vp', 'Ve'};
            vmin = [0,0,0,0];
            vmax = [DRO.flow(2)*1.25,DRO.ps(2)*1.25,DRO.vp(2)*1.25,DRO.ve(2)*1.25];
            vals = [DRO.flow; DRO.ps; DRO.vp; DRO.ve];
            for i = 1:4
                pred = fitresultsDCEcxm(:,:,i)';
                gt = myo_pars(:, :, i)';
                MSE = immse(gt, pred);
                SSIM = ssim(pred, gt);
                results(e_idx, i, :) = [MSE, SSIM];
                nexttile
                h1 = imshow(gt, [vmin(i), vmax(i)]); colormap(smoothCMRmap);
                colorbar;
                title(labels{i} + sprintf(" lower: %.3f upper: %.3f", vals(i, 1), vals(i, 2)))
                nexttile
                h2 = imshow(pred, [vmin(i), vmax(i)]); colormap(smoothCMRmap);
                title(sprintf("mse: %.3f\nssim: %.3f", MSE, SSIM))
                colorbar;
                nexttile
                h3 = imshow(abs(gt-pred), [vmin(i), vmax(i)]); colormap(smoothCMRmap);
                title("Absolute Diff")
                colorbar;
                t.Padding = 'none';
                t.TileSpacing = 'tight';
            end
            saveas(gcf, fullfile(output, [filename, '_gt_pred.png']))
            save(fullfile(output, [filename, '_data.mat']), 'results', 'fitresultsDCEcxm', 'simulatedCXM');

            % plot the aif and myo_curves
            figure('visible','off');
            h1 = plot(timepoints, aif, '.'); hold on
            h2 = plot(timepoints, myo_curves(:, 1, 1), '-.'); hold on
            h3 = plot(timepoints, simulatedCXM(:, 1, 1), '-.'); hold on
            leg_all = legend([h1 h2 h3], 'AIF', 'Myo curves', 'Simulated myo curves', 'fontsize',12, 'FontName', 'Arial', ...
                'Location', 'northeast', 'FontWeight', 'bold');
            title(sprintf("Flow:%.1f, PS:%.1f, Ve:%.2f, Vp:%.2f", myo_pars(1, 1, 1), myo_pars(1, 1, 2), myo_pars(1, 1, 3), myo_pars(1, 1, 4)));
            
            saveas(gcf, fullfile(output, [filename, '_aif_myo_curve.png']))
            % Save the simulated DRO and predicted DRO
            % mkdir(fullfile(output, name))
            % CXMB.mat2Nifti(permute(myo_curves, [2, 3, 1]), fullfile(output, name, 'SimulationDRO.nii'), [1, 1, 1]);
            % CXMB.mat2Nifti(permute(simulatedCXM, [2, 3, 1]), fullfile(output, name, 'PredictedDRO.nii'), [1, 1, 1]);
            % save(fullfile(output, name, 'workspace.mat'))
            % close all 
        end
    end
end
    
clear all
addpath(genpath('.'))
dir = "results/comparison#4";
dir = "results/#13 - pinnpaper";
mkdir(dir)
output = "./"+dir;

flow = [0.01, 5];
ps = [0.01, 3];
vp = [0.01, .3];
ve = [0.01, 0.5];
vecs = [flow; ps; vp; ve];
vars = ["flow", "ps", "vp", "ve"];
n = 20;
elem_per_iter = 15;
fontSize = 50;

DRO = SimulatedDRO();
simulation(DRO, output, elem_per_iter);
quit;

for param = 1:4
    vec = vecs(param, :);
    gr = (vec(2)/vec(1))^(1/(n-1));
    iter = vec(1) * gr.^(0:n-1);
    results = zeros(n-1, elem_per_iter, 4, 2);
    for variation = 1:length(iter)-1
        currpath = "/"+vars(param)+"/"+sprintf("lower%.3f_upper%.3f", iter(variation), iter(variation+1));
        [status, msg] = mkdir(dir+currpath);
        DRO = SimulatedDRO();
        if param == 1
            DRO.flow = [iter(variation), iter(variation+1)];
        elseif param == 2
            DRO.ps = [iter(variation), iter(variation+1)];    
        elseif param == 3
            DRO.vp = [iter(variation), iter(variation+1)];    
        elseif param == 4
            DRO.ve = [iter(variation), iter(variation+1)];
        else
            error("Values not changed")
        end
        results(variation, :, :, :) = simulation(DRO, output + currpath, elem_per_iter);
    end
    for time = 1:elem_per_iter
        fig = gcf;
        t = tiledlayout(4,4);
        width = 1200;
        height = 1200;
        set(fig, 'Position', [100, 100, width, height]);
        % t.Padding = 'none';
        t.TileSpacing = 'tight';
        title(t, sprintf(vars(param)+" variation between %.3f and %.3f at time %s", vecs(param, 1), vecs(param, 2)), time); 
        for allparam = 1:4
            nexttile;
            histogram('BinEdges',iter,'BinCounts',results(:, time, allparam, 1)');
            set(gca,'xscale','log')
            title(vars(allparam)+" mse log")
            nexttile;
            histogram('BinEdges',iter,'BinCounts',max(results(:, time, allparam, 2)',0));
            set(gca,'xscale','log')
            title(vars(allparam)+" ssim log")
            ylim([min(results(:, time, allparam, 2)'*.9), 1]);
            nexttile;
            histogram('BinEdges',iter,'BinCounts',results(:, time, allparam, 1)');
            title(vars(allparam)+" mse")
            nexttile;
            histogram('BinEdges',iter,'BinCounts',max(results(:, time, allparam, 2)',0));
            title(vars(allparam)+" ssim")
            ylim([min(results(:, time, allparam, 2)'*.9), 1]);
        end
        % exportgraphics(t,fullfile(output+"/"+vars(param)+"/compare_end"+time+".png"),'ContentType', 'vector')
        saveas(t,fullfile(output+"/"+vars(param)+"/compare_end"+time+".png"));
    end
    
end