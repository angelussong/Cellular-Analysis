function [fitresult, gof] = createFit(Neuron_name,t, v,cur)
%CREATEFIT(T,V)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : t
%      Y Output: v
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 09-May-2018 17:08:36


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( t, v );

% Set up fittype and options.
ft = fittype( 'a+b*exp(-x/tau)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-100 0.5 5];
opts.Robust = 'Bisquare';
opts.StartPoint = [0.695013095702136 0.439117405741358 0.913375856139019];
opts.Upper = [0 5 100];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure(2);
h = plot( fitresult, xData, yData );
legend( h, 'v vs. t', 'exponential fit', 'Location', 'NorthEast' );
% Label axes
ttl = sprintf('exponential fit to %dpA: baseline%g mV\t slope%g\t tau%g ms',cur,fitresult.a,fitresult.b,fitresult.tau);
title(ttl)
xlabel t
ylabel v
grid on
xlabel('Time (ms)');
ylabel('Vm (mV)');
saveas(gcf,[Neuron_name,'_tau_fitting.png']);

close



