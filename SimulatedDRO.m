classdef SimulatedDRO
    properties

        alpha (1, 1) {mustBeNumeric} = 2.0
        beta (1, 1) {mustBeNumeric} = 0.1
        snr (1, 1) {mustBeNumeric} = 0 % gaussian noise level
        ps_slice (1, 1) {mustBeNumeric} = 1 % referring to a PS full-slice value of (0.5, 1.5, 2.5) respectively
        t_start (1, 1) {mustBeNumeric} = 0 % started timepoints [min]
        t_intval (1, 1) {mustBeNumeric} = 0.02 % temporal resolution [min]
        t_end (1, 1) {mustBeNumeric} = 2 % end timepoints [min]

        flow (1, 4) {mustBeNumeric} = [0.5, 1.0, 1.5, 2.0]
        vp (1, 4) {mustBeNumeric} = [0.02, 0.05, 0.1, 0.2]
        ve (1, 3) {mustBeNumeric} = [0.1, 0.2, 0.5]
        ps (1, 3) {mustBeNumeric} = [0.5, 1.5, 2.5]

    end
    methods

        function [timepoints, aif, myo_curves_slice, myo_pars_slice] = SIMULATE_2CXM(obj, dims)
            % Simulate the DRO (digital reference object) data
            timepoints = obj.t_start:obj.t_intval:obj.t_end;
            timepoints = timepoints(1:end-1);

            % aif, arterial input function
            aif = obj.f_GAMMA(10^5/3, timepoints);

            % set parameters
            flows = zeros(dims);
            vps = zeros(dims);
            ves = zeros(dims);
            pss = zeros(dims);

            for i = 0:3
                flows(i*10+1:(i + 1)*10, :, :) = obj.flow(i+1);
            end

            for i = 0:3
                for j = 0:2
                    vps(:, 10*(i + 4 * j)+1:10*(i + 1 + 4 * j), :) = obj.vp(i+1);
                end
            end
            
            for i = 0:2
                ves(:, i*40+1:(i + 1)*40, :) = obj.ve(i+1);
            end
            
            for i = 0:2
                pss(:, :, i+1) = obj.ps(i+1);
            end


            % simulate the signal within the myo
            for i = 1:dims(1)
                for j = 1:dims(2)
                    for k = 1:dims(3)
                        % fprintf("i: %d, j: %d, k: %d\n", i, j, k)
                        pars = [flows(i, j, k), pss(i, j, k), vps(i, j, k), ves(i, j, k)];
                        % signal = obj.TWO_COMPONENT_EX_MODEL(num2cell(pars), timepoints, aif);
                        signal = obj.CXM_BOUND(pars, [timepoints', aif']);

                        % v approximate v_e + v_p
                        myo_curves(:, i, j, k, 1) = signal;
                        myo_pars(i, j, k, 1, :) = pars;
                    end
                end
            end

            if obj.snr ~= 0
                signal_db_mean = 10 * log10(mean(myo_curves(:)));
                noise_db_mean = signal_db_mean - obj.snr;
                noise_pw_mean = 10^(noise_db_mean / 10);
                noise_pw = normrnd(0, sqrt(noise_pw_mean), size(myo_curves));
                myo_curves = myo_curves + noise_pw;
            end

            % remove impossible values
            myo_curves(1, :, :, :, :) = 0;

            % input and output normalization
            % timepoints = normalize(timepoints);
            maxx = max(aif(:));
            % aif = aif / maxx;
            % myo_curves = myo_curves / maxx;

            % select slices
            myo_pars_slice = myo_pars(:, :, obj.ps_slice, :, :);
            myo_curves_slice = myo_curves(:, :, :, obj.ps_slice);
        end

        function [gamma] = f_GAMMA(obj, K, time)
            gamma = (K * (time.^obj.alpha)) .* double(exp(-time/obj.beta));
        end

        function [signal] = TWO_COMPONENT_EX_MODEL(obj, parameters, time, AIF)
            [f, vol_p, vol_e, PS] = deal(parameters{:});
            [tb, tp, te, km, k_p, E, delta] = obj.PHY2MODEL_PARAMS(f, vol_p, vol_e, PS);

            signal = (time(2) - time(1)) * f * ...
                (E * conv(exp(-km*time), AIF) + (1 - E) * conv(exp(-((km + delta) * time)), AIF));
            signal = signal(1:length(AIF));
        end
    end

    methods (Static)
        function [tb, tp, te, k_m, k_p, e, d] = PHY2MODEL_PARAMS(f, vp, ve, perm_surf)
            tb = vp / f;
            tp = vp / (perm_surf + f);
            te = ve / perm_surf;

            k_m = 0.5 * (1.0 / tp + 1.0 / te - sqrt((1.0 / tp + 1.0 / te)^2-4.0*(1.0 / tb)*(1.0 / te)));
            k_p = 0.5 * (1.0 / tp + 1.0 / te + sqrt((1.0 / tp + 1.0 / te)^2-4.0*(1.0 / tb)*(1.0 / te)));
            e = (k_p - 1.0 / tb) / (k_p - k_m);
            d = k_p - k_m;
        end

        function Ctoi = CXM_BOUND(beta, params)
            % ---------------------------------------------------------------
            % Calculate the simulated Ctoi according to the given parameters
            % and AIF in corresponding time, using the CXM_bound model.
            % /::INPUTS::\
            %     beta: {Flow, PS, Vp, Ve}, given parameters
            %     params: (time, Cp), e.g. size(100, 2), a column vector of time in minutes corresponding to Cp,
            % an Arterial Input Function (AIF) for the tissue in ROI.
            % /::OUTPUT::\
            %     Ctoi: Contrast Agent(CA) Concentration over time
            %     for the Tissue of Interest
            % ---------------------------------------------------------------
            F = beta(1);
            PS = beta(2);
            vp = beta(3);
            ve = beta(4);

            time = params(:, 1);
            Cp = params(:, 2);

            a = (F + PS) / vp;
            b = PS / ve;
            c = F / vp;
            M1 = 0.5 * (a + b + sqrt((a + b).^2-4*b*c));
            M2 = 0.5 * (a + b - sqrt((a + b).^2-4*b*c));
            B = (M2 - c) / (M2 - M1);
            H = B * exp(-M1*time) + (1 - B) * exp(-M2*time);
            Ctoi = F * convolution(H, Cp);
        end
    end
end