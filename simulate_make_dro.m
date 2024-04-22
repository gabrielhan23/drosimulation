% Simulate the four parameters, vp, ve, ps, flow

% Mona, April 2024
dims = [40, 120, 3];
t_start = 0;
t_end = 2;
t_intval = 0.02;
ps_slice = 2;
snr = 0;
[timepoints, aif, myo_curves, myo_pars] = simulate_2CXM(dims, t_start, t_end, t_intval, ps_slice, snr);

figure, plot(timepoints, aif)
x0=[0.5 2 0.1 0.5];
lb=[0 0 0 0];
ub=[3.0, 4.0, 0.24, 0.60];

fitresultsDCEcxm = Double_cxm(aif, myo_curves, timepoints, x0, lb, ub, 0.85, 20);
%%
F_cxm = fitresultsDCEcxm(:,:,1);
PS_cxm = fitresultsDCEcxm(:,:,2);
vp_cxm = fitresultsDCEcxm(:,:,3);
ve_cxm = fitresultsDCEcxm(:,:,4);
Rsq_cxm = fitresultsDCEcxm(:,:,end);

pred_params = fitresultsDCEcxm(:, :, [1, 3, 4, 2, 5]);

myo_pars = squeeze(myo_pars);

% demonstrate the results
CMRmap=[0 0 0;.15 .15 .5;.3 .15 .75;.6 .2 .50;1 .25 .15;.9 .5 0;.9 .75 .1;.9 .9 .5;1 1 1];

t = tiledlayout(4,2);
width = 300;    % Replace with your desired figure width
height = 800;   % Replace with your desired figure height
fig = gcf;      % Get the current figure handle
set(fig, 'Position', [100, 100, width, height]);
labels = {'Flow', 'Vp', 'Ve', 'PS'};
vmax = [3.0, 0.24, 0.60, 4.0];
for i = 1:4
    pred = pred_params(:,:,i)';
    gt = myo_pars(:, :, i)';
    MSE = immse(gt, pred);
    SSIM = ssim(pred, gt);
    nexttile
    h1 = imshow(gt, [0, vmax(i)]); colormap(CMRmap);
    title(labels{i})
    nexttile
    h2 = imshow(pred, [0, vmax(i)]); colormap(CMRmap);
    title(sprintf("mse: %.3f\nssim: %.3f", MSE, SSIM))
    colorbar
    t.Padding = 'none';
    t.TileSpacing = 'tight';
end

function [fitresultsDCEcxm] = Double_cxm(aif, myo_curves, time, x0, lb, ub, threshold, iters)
NityFivpercent=0;

[len, rows, cols] = size(myo_curves);

Rsqc=-1*ones([rows*cols 1]);
x0c=repmat(x0,[rows*cols 1]);
% cp: blood pool concentration mean
% sigin: ctoi, concentration within the heartmask
iter = 0;
while(NityFivpercent<threshold & iter<iters) % initial guess optimization
    for i = 1:rows
        for j = 1:cols
            c = (i-1) * rows + j;
            if Rsqc(c)<0.9   
                sigin = myo_curves(:,i,j);
                cp = aif';
                tempfitDCE = fitdcemri_etk(sigin, cp,time,x0,lb,ub,'cxm_bound');
                fitresultsDCEcxm(i, j, :) = tempfitDCE;
                Rsqc(c)=tempfitDCE(end);
                if tempfitDCE(end) > Rsqc(c)%initial guess otimization
                    x0c(c,:)=x0;
                    Rsqc(c)=tempfitDCE(end);
                end
                % adjustable
                x0 = x0c(c,:).*((1-(rand(1,length(x0))-0.5)*1.5));
            end
        end
    end
    fitQ = fitresultsDCEcxm(:,:,end);
    maskadj= (fitQ>0);
    f = sort(fitQ(maskadj==1));
    NityFivpercent = f(floor(sum(maskadj(:))*0.05));% get the 95% R2
    fprintf("Iter %d ---- 95 percentage R2 in iter: %.3f\n", iter, NityFivpercent)
    iter=iter+1;
end

end


function [timepoints, aif, myo_curves_slice, myo_pars_slice] = simulate_2CXM(dims, t_start, t_end, t_intval, ps_slice, snr)
% Simulate the DRO (digital reference object) data
timepoints = t_start:t_intval:t_end;
timepoints = timepoints(1:end-1);

% aif, arterial input function
aif = gamma_function(10^5/3, 2.0, 0.1, timepoints);

% set parameters
flows = zeros(dims);
vps = zeros(dims);
ves = zeros(dims);
pss = zeros(dims);
flow = [0.5, 1.0, 1.5, 2.0];
for i = 0:3
    flows(i*10+1:(i+1)*10, :, :) = flow(i+1);
end
vp = [0.02, 0.05, 0.1, 0.2];
for i = 0:3
    for j = 0:2
        vps(:,10*(i+4*j)+1:10*(i+1+4*j), :) = vp(i+1);
    end
end
ve = [0.1, 0.2, 0.5];
for i = 0:2
    ves(:, i*40+1:(i+1)*40, :) = ve(i+1);
end
ps = [0.5, 1.5, 2.5];
for i = 0:2
    pss(:, :, i+1) = ps(i+1);
end

% simulate the signal within the myo
for i = 1:dims(1)
    for j = 1:dims(2)
        for k = 1:dims(3)
            % fprintf("i: %d, j: %d, k: %d\n", i, j, k)
            pars = [flows(i, j, k), vps(i,j,k), ves(i,j,k), pss(i,j,k)];
            signal = two_comp_ex_model(num2cell(pars), timepoints, aif);

            % v approximate v_e + v_p
            myo_curves(:, i, j, k, 1) = signal;
            myo_pars(i, j, k, 1, :) = pars;
        end
    end
end

if snr ~= 0
    signal_db_mean = 10 * log10(mean(myo_curves(:)));
    noise_db_mean = signal_db_mean - snr;
    noise_pw_mean = 10 ^ (noise_db_mean / 10);
    noise_pw = normrnd(0, sqrt(noise_pw_mean), size(myo_curves));
    myo_curves = myo_curves + noise_pw;
end

% remove impossible values
myo_curves(1,:,:,:,:) = 0;

% input and output normalization
timepoints = normalize(timepoints);
maxx = max(aif(:));
aif = aif / maxx;
myo_curves = myo_curves / maxx;

% select slices
myo_pars_slice = myo_pars(:, :, ps_slice, :, :);
myo_curves_slice = myo_curves(:, :, :, ps_slice);
end

function [gamma] = gamma_function(K, alpha, beta, time)
gamma = (K*(time.^alpha)) .* double(exp(-time/beta));
end

function [signal] = two_comp_ex_model(parameters, time, AIF)
[flow, vol_p, vol_e, PS] = deal(parameters{:});
[tb, tp, te, km, k_p, E, delta] = phys_params_to_model_params(flow, vol_p, vol_e, PS);

signal = (time(2) - time(1)) * flow * ...
    (E * conv(exp(-km * time), AIF) + (1 - E) * conv(exp(-((km + delta) * time)), AIF));
signal = signal(1:length(AIF));
end

function [tb, tp, te, k_m, k_p, e, d] = phys_params_to_model_params(f, vp, ve, perm_surf)
tb = vp / f;
tp = vp / (perm_surf + f);
te = ve / perm_surf;

k_m = 0.5 * (1.0 / tp + 1.0 / te - sqrt((1.0 / tp + 1.0 / te) ^ 2 - 4.0 * (1.0 / tb) * (1.0 / te)));
k_p = 0.5 * (1.0 / tp + 1.0 / te + sqrt((1.0 / tp + 1.0 / te) ^ 2 - 4.0 * (1.0 / tb) * (1.0 / te)));
e = (k_p - 1.0 / tb) / (k_p - k_m);
d = k_p - k_m;
end
