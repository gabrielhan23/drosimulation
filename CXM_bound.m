classdef CXM_bound
    properties
        debug = false
        threshold (1, 1) {mustBeNumeric} = 0.85
        iters (1, 1) {mustBeNumeric} = 1
        x0(1, 4) {mustBeNumeric} = [0.01, 0.001, 0.1, 0.1] % F, PS, Vp, Ve
        lb(1, 4) {mustBeNumeric} = [0, 0, 0, 0]
        ub(1, 4) {mustBeNumeric} = [0.5, 0.5, 0.5, 0.5]
        CMRmap = [0 0 0;.15 .15 .5;.3 .15 .75;.6 .2 .50;1 .25 .15;.9 .5 0;.9 .75 .1;.9 .9 .5;1 1 1];
    end
    methods



        function [fitresultsDCEcxm, simulatedCXM, fullsimulatedCXM] = SOLVER_CXM_BOUND_ITERATIONS(obj, aif, aif_all, myo_curves, myo_curves_all, mask, myo_mask, time, time_all)
            % ---------------------------------------------------------------
            % solve the cxm bound model in an iterative optimization way
            % /::INPUTS::\
            %   aif: (1, length(time)), Arterial Input Function (AIF)
            %   myo_curves: (length(time), rows, cols), Contrast Agent(CA) Concentration over time for
            %   simulated region
            %   mask: ROI mask for myo_curves
            %   myo_mask: myocardium region for measuring accuracy of the
            %   curve fitting
            %   time: a column vector of time in minutes corresponding
            %   to aif
            % /::OUTPUT::\
            %   fitresultsDCEcxm: (rows, cols, 6), F, PS, Vp, Ve, rmse, R2
            % ---------------------------------------------------------------

            [len, rows, cols] = size(myo_curves);
            % initial guess optimization
            NityFivpercent = 0;
            iter = 0;
            Rsqc = -1 * ones([rows * cols, 1]);
            x0c = repmat(obj.x0, [rows * cols, 1]);
            simulatedCXM = zeros(size(myo_curves));
            fullsimulatedCXM = zeros(size(myo_curves_all));
            fitresultsDCEcxm = zeros(rows, cols, 6);
            while (NityFivpercent < obj.threshold && iter < obj.iters)
            % while (iter < obj.iters)
                for i = 1:rows
                    for j = 1:cols
                        c = (i - 1) * cols + j;

                        if mask(i, j) == 1
                            if Rsqc(c) < 0.9
                                sigin = myo_curves(:, i, j);
                                cp = aif';
                                % cp: blood pool concentration mean
                                % sigin: ctoi, concentration within the heartmask
                                % return: F, PS, Vp, Ve, rmse, R2
                                % x0c, F, PS, Vp, Ve
                                % tempfitDCE = obj.curve_fit(sigin, cp, time, x0c(c,:), obj.lb, obj.ub, 'cxm_bound');
                                tempfitDCE = obj.curve_fit(sigin, aif, time, x0c(c,:), obj.lb, obj.ub,'2cxm');
                                if tempfitDCE(end) > Rsqc(c)
                                    % fit success, update
                                    % simulatedCXM(:, i, j) = obj.CXM_BOUND(tempfitDCE(1:4), [time, aif]);
                                    simulatedCXM(:, i, j) = obj.TWO_COMPONENT_EX_MODEL(tempfitDCE(1:4), [time', aif']);
                                    fullsimulatedCXM(:, i, j) = obj.TWO_COMPONENT_EX_MODEL(tempfitDCE(1:4), [time_all', aif_all']);
                                    fitresultsDCEcxm(i, j, :) = tempfitDCE;
                                    Rsqc(c) = tempfitDCE(end);
                                end
                                % initial guess optimization
                                tmp = x0c(c,:);
                                x0c(c, :) =  tmp .* ((1 - (rand(1, length(tmp)) - 0.5) * 1.5));

                                % adjustable
                            end
                        end
                    end
                end

                fitQ = fitresultsDCEcxm(:, :, end);
                maskadj = myo_mask.*(fitQ > 0);
                f = sort(fitQ(maskadj == 1));
                try
                    NityFivpercent = f(max([(floor(sum(maskadj(:))*0.05)),1])); % get the 95% R2
                catch
                    NityFivpercent = 0;
                end
                fprintf("Iter %d ---- 95 percentage R2 in iter: %.3f\n", iter, NityFivpercent)
                iter = iter + 1;
            end

        end

        function [pars] = curve_fit(obj, varargin)

            toi = varargin{1};
            Cp = varargin{2};
            time = varargin{3};
            x0 = varargin{4};
            lb = varargin{5};
            ub = varargin{6};
            type = varargin{7};

            % Make sure all vectors are column, not row
            if size(toi,1) == 1
                toi = toi';
            end
            if size(Cp,1) == 1
                Cp = Cp';
            end
            if size(time,1) == 1
                time = time';
            end

            % Make sure data vectors are same size
            if size(toi,1) ~= size(time,1) || size(toi,1) ~= size(Cp,1) || size(Cp,1) ~= size(time,1)
                error('dce_mri_fit:toi,Cp,and time vectors must all be same size');
            end

            % Set options for non-linear LSQ fitting

            S=optimset;
            % S.Algorithm='interior-point';

            S.Algorithm='trust-region-reflective';
            % S.TolFun=1e-3; S.TolX=1e-3;% RY change fitting threshold
            S.TolFun=1e-6; S.TolX=1e-6;% RY change fitting threshold
            S.MaxIter=1000;
            % step size


            if obj.debug
                S.Display='iter-detailed';
                S.PlotFcns=["optimplotstepsize"];
                S.OutputFcn = @outfun;
                history.x = [];
            else
                S.Display='off';
            end

            function stop = outfun(x, optimValues, state)
                stop = false;

                switch state
                    case 'init'
                        hold on
                    case 'iter'
                        % Concatenate current point and objective function
                        % value with history. x must be a row vector.
                        if size(x, 1) == 1
                            history.x = [history.x; x];
                        else
                            history.x = [history.x; x'];
                        end
                    case 'done'
                        hold off
                    otherwise
                end
            end

            % Perform the fitting (the function "rrm" is in the nested functions below)
            if strcmp(type, "cxm_bound")
                [B,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@obj.CXM_BOUND, x0, [time, Cp], toi, lb, ub, S);
                times = linspace(time(1),time(end))';
                calculated = obj.CXM_BOUND(B, [time, Cp])
            else
                AIF = Cp;
                [B,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@obj.TWO_COMPONENT_EX_MODEL, x0, [time, AIF], toi, lb, ub, S);
                times = linspace(time(1),time(end))';
                calculated = obj.TWO_COMPONENT_EX_MODEL(B, [time, AIF]);
            end
            % for debug use
            if obj.debug
                disp(history.x)
                figure('Position', [100, 300, 1200, 300]), plot(time, toi,'ko',time,calculated,'g-', time, Cp, 'ro')
                legend('toi','Fitted exponential', 'aif')
                text = sprintf('FS: %.3f\nPs: %.3f\nVp: %.3f\nVe: %.3f', B(1),B(2),B(3),B(4));
                annotation('textbox',[.91 .4 .1 .3],'String',text,'FitBoxToText','on')
            end

            % Store each parameter (see reference)
            pars(1,1)=B(1);    %Fs
            pars(2,1)=B(2);    %PS
            pars(3,1)=B(3);    %Vp
            pars(4,1)=B(4);    %Ve

            pred = obj.CXM_BOUND(B, [time, Cp]);
            %pars(5,1)=rsquare(toi,pred);%ry
            [pars(6,1),pars(5,1)]=obj.rsquare(toi,pred);%ry
        end


    end

    methods (Static)
        function [] = mat2Nifti(volume, savepath, voxelSize)
            % save Nifti
            % reference https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image
            temp_nii = make_nii(volume);
            temp_nii.hdr.dime.pixdim(2:4) = voxelSize;

            save_nii(temp_nii, savepath);
            disp(["Suceess to save the nifti files" + savepath])
        end
        function [tb, tp, te, k_m, k_p, e, d] = PHY2MODEL_PARAMS(f, vp, ve, perm_surf)
            tb = vp / f;
            tp = vp / (perm_surf + f);
            te = ve / perm_surf;

            k_m = 0.5 * (1.0 / tp + 1.0 / te - sqrt((1.0 / tp + 1.0 / te)^2-4.0*(1.0 / tb)*(1.0 / te)));
            k_p = 0.5 * (1.0 / tp + 1.0 / te + sqrt((1.0 / tp + 1.0 / te)^2-4.0*(1.0 / tb)*(1.0 / te)));
            e = (k_p - 1.0 / tb) / (k_p - k_m);
            d = k_p - k_m;
        end
        function Ctoi = TWO_COMPONENT_EX_MODEL(beta, params)
            F = beta(1);
            PS = beta(2);
            vp = beta(3);
            ve = beta(4);

            time = params(:, 1);
            AIF = params(:, 2);
            time = time';
            AIF = AIF';
            
            [tb, tp, te, km, k_p, E, delta] = CXM_bound.PHY2MODEL_PARAMS(F, vp, ve, PS);

            signal = (time(2) - time(1)) * F * ...
                (E * conv(exp(-km*time), AIF) + (1 - E) * conv(exp(-((km + delta) * time)), AIF));
            Ctoi = signal(1:length(AIF))';
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
            % Ctoi = F * conv(H, Cp, 'same');
        end

        function [r2, rmse] = rsquare(y,f,varargin)
            % Compute coefficient of determination of data fit model and RMSE
            %
            % [r2 rmse] = rsquare(y,f)
            % [r2 rmse] = rsquare(y,f,c)
            %
            % RSQUARE computes the coefficient of determination (R-square) value from
            % actual data Y and model data F. The code uses a general version of
            % R-square, based on comparing the variability of the estimation errors
            % with the variability of the original values. RSQUARE also outputs the
            % root mean squared error (RMSE) for the user's convenience.
            %
            % Note: RSQUARE ignores comparisons involving NaN values.
            %
            % INPUTS
            %   Y       : Actual data
            %   F       : Model fit
            %
            % OPTION
            %   C       : Constant term in model
            %             R-square may be a questionable measure of fit when no
            %             constant term is included in the model.
            %   [DEFAULT] TRUE : Use traditional R-square computation
            %            FALSE : Uses alternate R-square computation for model
            %                    without constant term [R2 = 1 - NORM(Y-F)/NORM(Y)]
            %
            % OUTPUT
            %   R2      : Coefficient of determination
            %   RMSE    : Root mean squared error
            %
            % EXAMPLE
            %   x = 0:0.1:10;
            %   y = 2.*x + 1 + randn(size(x));
            %   p = polyfit(x,y,1);
            %   f = polyval(p,x);
            %   [r2 rmse] = rsquare(y,f);
            %   figure; plot(x,y,'b-');
            %   hold on; plot(x,f,'r-');
            %   title(strcat(['R2 = ' num2str(r2) '; RMSE = ' num2str(rmse)]))
            %
            % Jered R Wells
            % 11/17/11
            % jered [dot] wells [at] duke [dot] edu
            %
            % v1.2 (02/14/2012)
            %
            % Thanks to John D'Errico for useful comments and insight which has helped
            % to improve this code. His code POLYFITN was consulted in the inclusion of
            % the C-option (REF. File ID: #34765).

            if isempty(varargin); c = true;
            elseif length(varargin)>1; error 'Too many input arguments';
            elseif ~islogical(varargin{1}); error 'C must be logical (TRUE||FALSE)'
            else c = varargin{1};
            end

            % Compare inputs
            if ~all(size(y)==size(f)); error 'Y and F must be the same size'; end

            % Check for NaN
            tmp = ~or(isnan(y),isnan(f));
            y = y(tmp);
            f = f(tmp);

            if c; r2 = max(0,1 - sum((y(:)-f(:)).^2)/sum((y(:)-mean(y(:))).^2));
            else r2 = 1 - sum((y(:)-f(:)).^2)/sum((y(:)).^2);
                if r2<0
                    % http://web.maths.unsw.edu.au/~adelle/Garvan/Assays/GoodnessOfFit.html
                    warning('Consider adding a constant term to your model') %#ok<WNTAG>
                    r2 = 0;
                end
            end

            rmse = sqrt(mean((y(:) - f(:)).^2));

        end


    end

end