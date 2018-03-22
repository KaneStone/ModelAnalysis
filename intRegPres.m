function [dataRegPres,regpres] = intRegPres(data,pressure,varargin)

% integrate regular pressure over 1 dimension (pressure, dim1);

if nargin == 3
    regpres = varargin{1};
else 
    regpres = [1000;975;950;925;900;875;850;825;800;775;750;700;650;600;550;500;450;400;350;...
        300;250;225;200;175;150;125;100;70;50;30;20;10;7;5;3;2;1;0.7;0.5;0.3;0.2;0.1;.07;.05;...
        .03;.02;.01;.007;.005;.003;.002;.001];
end

for j = 1:size(data,2)    
    dataRegPres(:,j)  = interp1(log(pressure(:,j)),data(:,j),log(regpres),'linear','extrap');
end

end