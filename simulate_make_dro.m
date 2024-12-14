% Simulate the four parameters, vp, ve, ps, flow
% depend on Classdef SimulatedDRO and CXM_bound
% Mona, April 2024

function results = simulate_make_dro(DRO, output, t_start, end_val, start_end_val)
    results = zeros(end_val-fix(start_end_val/60), 4, 2);
    resume = true;
    reimage = false;

    DRO.snr = 0;
    DRO.t_start = 0;
    DRO.t_end = 21;
    DRO.ps_slice = 1;

    if any(mod(t_start, DRO.t_intval))
        error("starting value is not on a timepoint")
    end
    % Load DRO Params
    if exist(fullfile(output, 'dro_data.mat'), "file") && resume
        fprintf("Loading dro\n");
        S = load(fullfile(output, 'dro_data.mat'));
        dims = S.dims; timepoints_all = S.timepoints_all; aif_all = S.aif_all; myo_curves_all = S.myo_curves_all; myo_pars = S.myo_pars;
        DRO = S.DRO;
    else
        [dims, timepoints_all, aif_all, myo_curves_all, myo_pars] = DRO.SIMULATE_2CXM();
        save(fullfile(output, 'dro_data.mat'), 'DRO', 'dims', 'timepoints_all', 'aif_all', 'myo_curves_all', 'myo_pars');
    end

    % Plot AIF and one myo curve
    figure('visible','off');
    h1 = plot(timepoints_all, aif_all, '-o'); hold on
    h2 = plot(timepoints_all, myo_curves_all(:, 1, 1), '-o'); hold on
    leg_all = legend([h1 h2], 'AIF', 'Myo curves', 'fontsize', 12, 'FontName', 'Arial', ...
        'Location', 'northeast', 'FontWeight', 'bold');
    saveas(gcf, fullfile(output, 'alltimepoints_aif_myo_curve.png'));
    
    % Plot all myo curves
    fig = figure('visible','off');
    t = tiledlayout(4,4);
    tiles = arrayfun(@(x) nexttile, 1:16);
    set(fig, 'Position', [100, 100, 1200, 1200]);
    title(t, "All Myo Curves")
    plotMyos(timepoints_all, myo_curves_all, DRO, tiles, "gt")
    saveas(fig, fullfile(output, 'allmyocurves.png'));

    for s_idx = 1:length(t_start)
        t_end = (t_start(s_idx)+start_end_val:60:1200);
        t_end = t_end(1:end_val);

        % One iteration of time
        for e_idx = 1:length(t_end)
            start_index = t_start(s_idx)/60/DRO.t_intval;
            end_index = t_end(e_idx)/60/DRO.t_intval;

            % timepoints and aif start from 0 and are sliced after fed
            % through 2cxm 
            timepoints = timepoints_all(1:end_index);
            aif = aif_all(1:end_index);
            myo_curves = myo_curves_all(start_index+1:end_index,:,:);
            % end
            name = sprintf('start%.1f end%.1f interval%.2f snr%.2f', t_start(s_idx)/60, t_end(e_idx)/60, DRO.t_intval, DRO.snr);
            filename = sprintf('start%.1f_end%.1f_interval%.2f_snr%.2f', t_start(s_idx)/60, t_end(e_idx)/60, DRO.t_intval, DRO.snr);

            CXMB = CXM_bound();
            mask = ones(dims(1), dims(2));
            if exist(fullfile(output, [filename, '_gt_pred.png']), "file") && resume
                fprintf("Already calculated %s in %s\n", name, output)
                S = load(fullfile(output, [filename, '_data.mat']));
                results = S.results;
                if reimage == false
                    continue
                end
                fitresultsDCEcxm = S.fitresultsDCEcxm;
                simulatedCXM = S.simulatedCXM;
                fullsimulatedCXM = S.fullsimulatedCXM;
            else
                fprintf("Processing simulation %s in %s...\n", name, output)
                [fitresultsDCEcxm, simulatedCXM, fullsimulatedCXM] = CXMB.SOLVER_CXM_BOUND_ITERATIONS(aif, aif_all, myo_curves, myo_curves_all, mask, mask, timepoints, timepoints_all, start_index+1);
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

            % Save predicted 
            saveas(gcf, fullfile(output, [filename, '_gt_pred.png']))
            save(fullfile(output, [filename, '_data.mat']), 'results', 'fitresultsDCEcxm', 'simulatedCXM', 'fullsimulatedCXM');
            
            % Plot calculated myo curves 
            figure('visible','off');
            t = tiledlayout(4,4);
            tiles = arrayfun(@(x) nexttile, 1:16);
            set(gcf, 'Position', [100, 100, 1200, 1200]);
            title(t, "All Myo Curves for "+name)
            plotMyos(timepoints_all, myo_curves_all, DRO, tiles, "gt")
            plotMyos(timepoints_all(start_index+1:end_index), myo_curves, DRO, tiles, "current")
            plotMyos(timepoints_all, fullsimulatedCXM, DRO, tiles, "calculated")
            saveas(gcf, fullfile(output, [filename, '_allmyocurves.png']));
        end
    end
end

function plotMyos(time, myo_curves, DRO, tiles, displayname)
    hold(tiles(1), 'on'); plot(tiles(1), time, myo_curves(:, 1, 1), '-', 'DisplayName', displayname); title(tiles(1), "min flow, min ps, min vp, min ve"); legend(tiles(1), 'show'); 
    hold(tiles(2), 'on'); plot(tiles(2), time, myo_curves(:, 1, DRO.vpstep), '-', 'DisplayName', displayname); title(tiles(2), "min flow, min ps, min vp, max ve"); legend(tiles(2), 'show'); 
    hold(tiles(3), 'on'); plot(tiles(3), time, myo_curves(:, 1, end-DRO.vpstep), '-', 'DisplayName', displayname); title(tiles(3), "min flow, max ps, min vp, min ve"); legend(tiles(3), 'show'); 
    hold(tiles(4), 'on'); plot(tiles(4), time, myo_curves(:, 1, end), '-', 'DisplayName', displayname); title(tiles(4), "min flow, max ps, min vp, max ve"); legend(tiles(4), 'show'); 

    hold(tiles(5), 'on'); plot(tiles(5), time, myo_curves(:, DRO.psstep, 1), '-', 'DisplayName', displayname); title(tiles(5), "max flow, min ps, min vp, min ve"); legend(tiles(5), 'show'); 
    hold(tiles(6), 'on'); plot(tiles(6), time, myo_curves(:, DRO.psstep, DRO.vpstep), '-', 'DisplayName', displayname); title(tiles(6), "max flow, min ps, min vp, max ve"); legend(tiles(6), 'show'); 
    hold(tiles(7), 'on'); plot(tiles(7), time, myo_curves(:, DRO.psstep, end-DRO.vpstep), '-', 'DisplayName', displayname); title(tiles(7), "max flow, max ps, min vp, min ve"); legend(tiles(7), 'show'); 
    hold(tiles(8), 'on'); plot(tiles(8), time, myo_curves(:, DRO.psstep, end), '-', 'DisplayName', displayname); title(tiles(8), "max flow, max ps, min vp, max ve"); legend(tiles(8), 'show'); 

    hold(tiles(9), 'on'); plot(tiles(9), time, myo_curves(:, end-DRO.psstep, 1), '-', 'DisplayName', displayname); title(tiles(9), "min flow, min ps, max vp, min ve"); legend(tiles(9), 'show'); 
    hold(tiles(10), 'on'); plot(tiles(10), time, myo_curves(:, end-DRO.psstep, DRO.vpstep), '-', 'DisplayName', displayname); title(tiles(10), "min flow, min ps, max vp, max ve"); legend(tiles(10), 'show'); 
    hold(tiles(11), 'on'); plot(tiles(11), time, myo_curves(:, end-DRO.psstep, end-DRO.vpstep), '-', 'DisplayName', displayname); title(tiles(11), "min flow, max ps, max vp, min ve"); legend(tiles(11), 'show'); 
    hold(tiles(12), 'on'); plot(tiles(12), time, myo_curves(:, end-DRO.psstep, end), '-', 'DisplayName', displayname); title(tiles(12), "min flow, max ps, max vp, max ve"); legend(tiles(12), 'show'); 

    hold(tiles(13), 'on'); plot(tiles(13), time, myo_curves(:, end, 1), '-', 'DisplayName', displayname); title(tiles(13), "max flow, min ps, max vp, min ve"); legend(tiles(13), 'show'); 
    hold(tiles(14), 'on'); plot(tiles(14), time, myo_curves(:, end, DRO.vpstep), '-', 'DisplayName', displayname); title(tiles(14), "max flow, min ps, max vp, max ve"); legend(tiles(14), 'show'); 
    hold(tiles(15), 'on'); plot(tiles(15), time, myo_curves(:, end, end-DRO.vpstep), '-', 'DisplayName', displayname); title(tiles(15), "max flow, max ps, max vp, min ve"); legend(tiles(15), 'show'); 
    hold(tiles(16), 'on'); plot(tiles(16), time, myo_curves(:, end, end), '-', 'DisplayName', displayname); title(tiles(16), "max flow, max ps, max vp, max ve"); legend(tiles(16), 'show'); 
    
end
    
