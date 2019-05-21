function [inputs] = predruns_inputs

% general
inputs.tozvar = 'toz';
inputs.var = 'ICEFRAC'; %TS %SNOWHICE %ICEFRAC
inputs.detrend = 1;
inputs.detrend_ozone = 1;
inputs.lats = [63,90];
inputs.ClLevel = {'highCl'}; %{'highCl','lowCl','lowGHG'}
inputs.timeperiodtoz = [1995,2024];%[1955,1979];
inputs.timeperiodvar = [1995,2024];%[1955,1979];
% inputs.timeperiodtoz = [1955,1979];
% inputs.timeperiodvar = [1955,1979];
inputs.percentile = 20;
inputs.tozmonth = 3;
inputs.varmonth = [3,4,5,6,7,8,9,10,11];%[3,4,5,6,7];
inputs.varmonthtomean = [9];
inputs.removeENSO = 1;
inputs.lastfiveyears = 0;
%compare to obs
inputs.compareERA = 1;
inputs.takediff = 1;
inputs.obsperc = 0;
inputs.includemarkers = 1;
inputs.lastfiveyears_obs = 0;
inputs.justextractobs = 0;

% plotting inputs
inputs.plotcorr = 0;
inputs.plotdiff = 0;
inputs.plotind = 0;
inputs.plotindmonths = 0;
inputs.plotcombine = 0;
end
