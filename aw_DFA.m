function [dfa_exponent, dfa_intercept, scales, rmses] = aw_DFA(signal, nScales, to_plot)
%
% [dfa_exponent, dfa_intercept, scales, rmses] = aw_DFA(signal, nScales)
% 
% Detrended fluctuation analysis according to the methods and code of MX
% Cohen.
% 
%   Inputs:
%       signal -                    signal to be measured in a nx1 vector
%       nScales -               integer (ex: 20) of number of timescales to be
%                                           used in DFA calculation
%                                           (default = 20).
%       to_plot -                   Boolean for plotting of figure (default
%                                           = 1). Set as 0 if do not want
%                                           to plot results.
% 
%   Outputs:
%       dfa_exponent -  DFA exponent (DFAe) which is the slope of the
%                               line of best fit for the log-log plot between scales and RMS
%                               values. Also known as the Hurst exponent.
%       dfa_intercept -   y-intercept of the best fit line for the log-log
%                               plot between scales and RMS values.
%       scales -            log-spaced scales between upper and lower
%                                   bounds, length of nScales
%       rmses -         computed root mean squares
% 
% written by annemarie wolff, awolf037@uottawa.ca
%
%% setup parameters for DFA
if nargin == 1
    nScales = 20;
    to_plot = 1;
end

Num_pts = length(signal);
ranges  = round(Num_pts*[.01 .2]); % in fraction of the total length, so 1% of length of total signal to 20% of length of total signal
scales  = ceil(logspace(log10(ranges(1)), log10(ranges(2)), nScales)); % log-spaced scales between lower and upper bound
rmses   = zeros(1, nScales);

%% COMPUTE DFA
% integrate and mean-center the signals
integ_sig = cumsum(signal - nanmean(signal));

% compute RMS over different time scales
for scalei = 1:nScales
    % number of epochs for this scale
    num_epochs = floor(Num_pts/scales(scalei)); 
    % compute RMS for the signal
    epochs = reshape(integ_sig(1:num_epochs*scales(scalei)), scales(scalei), num_epochs);
    % detrend
    depochs = detrend(epochs);
    % the root mean square computation
    rmses(1,scalei) = nanmean(sqrt(nanmean(depochs.^2, 1)));
end

% fit a linear model to quantify scaling Hurst exponent
A = [ ones(nScales,1) log10(scales)' ];  % linear model
dfa = (A'*A) \ (A'*log10(rmses(1,:))'); % fit to signal; first number corresponds to the y-intercept, the second number to the slope of the line

dfa_intercept = dfa(1);
dfa_exponent = dfa(2);

%% plot results
if to_plot == 1
    hold on;
    plot(log10(scales), log10(rmses), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
    plot(log10(scales), (dfa_intercept + dfa_exponent * log10(scales)), 'k--', 'linew', 2);
    xlabel('Data scale (log)');
    ylabel('RMS (log)');
    title('DFA/Hurst exponent results', 'FontWeight', 'normal');
    legend({'Data';['Fit (DFA = ' num2str(round(dfa_exponent)) ')']}, 'Location', 'northwest');
end

