addpath(genpath('.'))
dir = "./results/comparison#5";
mkdir(dir)

flow = [0.01, 5];
ps = [0.01, 3];
vp = [0.01, .3];
ve = [0.01, 0.5];
vecs = [flow; ps; vp; ve];
vars = ["flow", "ps", "vp", "ve"];
n = 20;
elem_per_iter = 8;

for param = 1:4
    vec = vecs(param, :);
    gr = (vec(2)/vec(1))^(1/(n-1));
    iter = vec(1) * gr.^(0:n-1);
    results = zeros(n-1, elem_per_iter, 4, 2);
    for variation = 1:length(iter)-1
        currpath = "/"+vars(param)+"/"+sprintf("lower%.3f_upper%.3f", iter(variation), iter(variation+1));
        [status, msg] = mkdir(dir+currpath);
        DRO = SimulatedDRO();
        
        % for the comparison
        DRO.flow = [0.15, .5];
        DRO.ps = [0.05, 0.12];
        DRO.vp = [0.06, .17];
        DRO.ve = [0.07, 0.2];

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
        results(variation, :, :, :) = simulate_make_dro(DRO, dir + currpath, [0], elem_per_iter);
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
        % exportgraphics(t,fullfile(dir+"/"+vars(param)+"/compare_end"+time+".png"),'ContentType', 'vector')
        saveas(t,fullfile(dir+"/"+vars(param)+"/compare_end"+time+".png"));
    end
    
end
