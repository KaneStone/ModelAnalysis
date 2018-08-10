function [inputs] = predruns_inputs

% general
inputs.tozvar = 'toz';
inputs.var = 'TS';
inputs.detrend = 0;
inputs.detrend_ozone = 0;
inputs.lats = [63,90];
inputs.ClLevel = {'highCl'}; %{'highCl','lowCl','lowGHG'}
inputs.timeperiodtoz = [1995,2024];
inputs.timeperiodvar = [1995,2024];
inputs.percentile = 20;
inputs.tozmonth = 3;
inputs.varmonth = [3,4,5,6,7];
inputs.varmonthtomean = [4];
inputs.removeENSO = 0;
inputs.lastfiveyears = 1;
%compare to obs
inputs.compareERA = 1;
inputs.takediff = 1;
inputs.obsperc = 0;
inputs.includemarkers = 1;
inputs.lastfiveyears_obs = 1;
inputs.justextractobs = 1;

% plotting inputs
inputs.plotcorr = 0;
inputs.plotdiff = 0;
inputs.plotind = 0;
inputs.plotindmonths = 0;
inputs.plotcombine = 0;
end