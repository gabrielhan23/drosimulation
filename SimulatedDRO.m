classdef SimulatedDRO
    properties

        alpha (1, 1) {mustBeNumeric} = 2.0
        beta (1, 1) {mustBeNumeric} = 0.1
        snr (1, 1) {mustBeNumeric} = 0 % gaussian noise level
        ps_slice (1, 1) {mustBeNumeric} = 1 % referring to a PS full-slice value of (0.5, 1.5, 2.5) respectively
        t_start (1, 1) {mustBeNumeric} = 0 % started timepoints [min]
        t_intval (1, 1) {mustBeNumeric} = 0.02 % temporal resolution [min]
        t_end (1, 1) {mustBeNumeric} = 2 % end timepoints [min]
        
        % start, stop
        flow (1, 2) {mustBeNumeric} = [0.15, .5]
        ps (1, 2) {mustBeNumeric} = [0.05, 0.12]
        vp (1, 2) {mustBeNumeric} = [0.06, .17]
        ve (1, 2) {mustBeNumeric} = [0.07, 0.2]

        % flow (1, 2) {mustBeNumeric} = [0.05, 1]
        % ps (1, 2) {mustBeNumeric} = [0.01, 0.2]
        % vp (1, 2) {mustBeNumeric} = [0.01, .4]
        % ve (1, 2) {mustBeNumeric} = [0.01, 0.5]
        
        flowstep (1, 1) {mustBeNumeric} = 10
        psstep (1, 1) {mustBeNumeric} = 10
        vpstep (1, 1) {mustBeNumeric} = 10
        vestep (1, 1) {mustBeNumeric} = 10
        
        % parameters for gaussian kernal curve
        B (1, 1) {mustBeNumeric} =  1.050 % mM
        a1 (1, 1) {mustBeNumeric} = 0.809 % mM min
        a2 (1, 1) {mustBeNumeric} = 0.330 % mM min
        m1 (1, 1) {mustBeNumeric} = 0.1682 % min-1
        m2 (1, 1) {mustBeNumeric} = 38.078 % min-1
        sigma1 (1, 1) {mustBeNumeric} = 0.0563 % gaussian distribution
        sigma2 (1, 1) {mustBeNumeric} = 0.132 % gaussian distribution
        mu1 (1, 1) {mustBeNumeric} = 0.170 % gaussian distribution
        mu2 (1, 1) {mustBeNumeric} = 0.365 % gaussian distribution
        tc (1, 1) {mustBeNumeric} = 0.1483 % sigmod width and center

        % parameters for dual exponential curve
        D (1, 1) {mustBeNumeric} =  0.1 %0.1
        dual_a1 (1, 1) {mustBeNumeric} = 8
        dual_a2 (1, 1) {mustBeNumeric} = 3
        dual_m1 (1, 1) {mustBeNumeric} = 0.5
        dual_m2 (1, 1) {mustBeNumeric} = 0.01
        

    end
    methods
        % function obj = SimulatedDRO()
        % 
        % end 
        function [dims, timepoints, aif, myo_curves_slice, myo_pars_slice] = SIMULATE_2CXM(obj)
            % dims = (40, 120, 3)
            flowvec = linspace(obj.flow(1), obj.flow(2), obj.flowstep);
            psvec = linspace(obj.ps(1), obj.ps(2), obj.psstep).';
            vpvec = linspace(obj.vp(1), obj.vp(2), obj.vpstep);
            vevec = linspace(obj.ve(1), obj.ve(2), obj.vestep).';

            dims = [obj.flowstep * obj.vpstep, obj.vestep * obj.psstep];
            % Simulate the DRO (digital reference object) data
            timepoints = obj.t_start:obj.t_intval:obj.t_end;
            timepoints = timepoints(1:end-1);

            % aif, arterial input function
            % aif = obj.f_GAMMA(10^5/3, timepoints);
            aif = obj.DUAL_EXP(timepoints);
            % aif = obj.GAUSSIAN_EXP(timepoints);

            % used for 100 by 100
            flows = repmat(flowvec, obj.vestep * obj.psstep, obj.vpstep);
            ves = repmat(vevec, obj.psstep, obj.vpstep * obj.flowstep);
            vps = repelem(vpvec, obj.vestep * obj.psstep, obj.flowstep);
            pss = repelem(psvec, obj.vestep, obj.flowstep * obj.vpstep);

            % used for a 40 by 120 matrix 
            % for i = 0:(length(obj.flow)-1)
            %     flows(i*10+1:(i + 1)*10, :, :) = obj.flow(i+1);
            % end
            % 
            % for i = 0:(length(obj.vp)-1)
            %     for j = 0:2
            %         vps(:, 10*(i + 4 * j)+1:10*(i + 1 + 4 * j), :) = obj.vp(i+1);
            %     end
            % end
            % 
            % for i = 0:(length(obj.ve)-1)
            %     ves(:, i*(dims(2)/length(obj.ve))+1:(i + 1)*(dims(2)/length(obj.ve)), :) = obj.ve(i+1);
            % end
            % 
            % for i = 0:(length(obj.ps)-1)
            %     pss(:, :, i+1) = obj.ps(i+1);
            % end


            % simulate the signal within the myo
            for i = 1:dims(1)
                for j = 1:dims(2)
                    % fprintf("i: %d, j: %d, k: %d\n", i, j, k)
                    pars = [flows(i, j), pss(i, j), vps(i, j), ves(i, j)];
                    signal = obj.TWO_COMPONENT_EX_MODEL(pars, timepoints, aif);
                    % signal = obj.CXM_BOUND(pars, [timepoints', aif']);

                    % v approximate v_e + v_p
                    myo_curves(:, i, j) = signal;
                    myo_pars(i, j, :) = pars;
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
            myo_curves(1, :, :, :) = 0;

            % input and output normalization
            % timepoints = normalize(timepoints);
            maxx = max(aif(:));
            % aif = aif / maxx;
            % myo_curves = myo_curves / maxx;

            % select slices
            myo_pars_slice = myo_pars;
            myo_curves_slice = myo_curves;
        end

        function [gamma] = f_GAMMA(obj, K, time)
            gamma = (K * (time.^obj.alpha)) .* double(exp(-time/obj.beta));
        end

        function [Cp] = GAUSSIAN_EXP(obj, time)
            % G. J. M. Parker, etc "Experimentally-derived functional form for a population-averaged 
            % hightemporal-resolution arterial input function for dynamic contrast-enhanced MRI,"
            p1 = (obj.B .* exp(-obj.m1 .* time)) ./ (1 + exp(-obj.m2 .* (time - obj.tc)));
            p2 = (obj.a1 / obj.sigma1 / sqrt(2 * pi)) .* exp(-((time - obj.mu1)./obj.sigma1.*sqrt(2)).^2);
            p3 = (obj.a2 / obj.sigma2 / sqrt(2 * pi)) .* exp(-((time - obj.mu2)./obj.sigma2.*sqrt(2)).^2);
            Cp = p1 + p2 + p3;
        end

        function [Cp] = DUAL_EXP(obj, time)
            % "Cp(t) = D(a1e^(−m1t) + a2e^(−m2t))" ([Khalifa et al., 2014, p. 1243016]
            Cp = obj.D .* (obj.dual_a1 .* exp(-obj.dual_m1 * time) + obj.dual_a2 .* exp(-obj.dual_m2 * time));
            
        end

        function [signal] = TWO_COMPONENT_EX_MODEL(obj, parameters, time, AIF)
            F = parameters(1);
            PS = parameters(2);
            vp = parameters(3);
            ve = parameters(4);
            
            [tb, tp, te, km, k_p, E, delta] = obj.PHY2MODEL_PARAMS(F, vp, ve, PS);

            signal = (time(2) - time(1)) * F * ...
                (E * conv(exp(-km*time), AIF) + (1 - E) * conv(exp(-((km + delta) * time)), AIF));
            signal = signal(1:length(AIF))';
               
            % same model different variables 
            % [f, PS, vol_p, vol_e] = deal(parameters{:});
            % [tb, tp, te, km, k_p, E, delta] = obj.PHY2MODEL_PARAMS(f, vol_p, vol_e, PS);
            % 
            % signal = (time(2) - time(1)) * f * ...
            %     (E * conv(exp(-km*time), AIF) + (1 - E) * conv(exp(-((km + delta) * time)), AIF));
            % signal = signal(1:length(AIF));
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
            % and AIF in corresponding time, using the CXM_bound model(2CXM).
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