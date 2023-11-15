%fit Cp
count=0
clear pars
%T0=-600:10:600;
for T0=0:10:600
    count=count+1;
x0=[T0 Cp(2) 0 0 0 0 0];
lb=[-1000 -10 -10 -10 -10 -10];
ub=[1000 10 10 10 10 10]

  % Set options for non-linear LSQ fitting
            S=optimset; S.Algorithm='trust-region-reflective'; S.Display='off';
            %S.TolFun=1e-3; S.TolX=1e-3;% RY change fitting threshold
            S.TolFun=1e-10; S.TolX=1e-10;% RY change fitting threshold
            S.MaxIter=1000;

            % Perform the fitting (the function "rrm" is in the nested functions below)
            [B,~,~] = lsqcurvefit(@TKCp,x0,timein',Cp(2:end)',lb,ub,S);

            % Store each parameter (see reference)
            pars(1,count)=B(1);    %Ktrans
            pars(2,count)=B(2);    %kepTOI
            pars(3,count)=max(B(3),B(4));    %kepTOI
            pars(4,count)=min(B(3),B(4));    %Ktrans
            pars(5,count)=B(find(B==max(B(3),B(4)))+2);    %kepTOI
            pars(6,count)=B(find(B==min(B(3),B(4)))+2);
            pred=TKCp(B,timein');
            %pars(4,1)=rsquare(toi,pred);
                      [pars(7,count),pars(8,count)]=rsquare(Cp(2:end)',pred);%ry

                      
end

function Cp = TKCp(beta,X)
    % beta is 
    % [T_first D a1 a2 m1 m2]
    time=X(:,1);%sec
    
    T_f = beta(1); %sec
    D = beta(2); %1/t=sec^-1
    a1=beta(5); % 1/100
    a2=beta(6); %1/100
    m1=beta(3); % 1/100
    m2=beta(4); %1/100
    Cp=D*(a1*exp(-m1*(T_f+time))+a2*exp(-m2*(T_f+time)));
    %C_toi = F*convolution(H,Cp);

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
