function conc=vmr2conc(vmr,t,p,trace,varargin)
%converts volume mixing ratios to concentrations in molec.cm-3
%t - temperature in K
%p - pressure in Pa
%for keyword 'ice' conc is returned as gm-3

k=1.38066e-23;  %JK-1molec-1
Na=6.02214179e23; %mol-1
molmasswater=18.015286; %Molecular mass of water (g mol-1)
R=k*Na; %JK-1mol-1
air = 28.97; %(g mol-1)

if (strcmp(trace, 'O3'))
    tr = 48; %(g mol-1)
elseif (strcmp(trace, 'NO2'))
    tr = 46.0055; %(g mol-1)
elseif (strcmp(trace, 'N2O'))
    tr = 44.0130; %(g mol-1)
elseif (strcmp(trace, 'NONO2'))
    tr = 76.0116; %(g mol-1)
elseif (strcmp(trace, 'BrO'))
    tr = 95.904; %(g mol-1)
elseif (strcmp(trace, 'CLO'))
    tr = 51.45; %(g mol-1)
elseif (strcmp(trace, 'HCL'))
    tr = 36.46; %(g mol-1)
elseif (strcmp(trace, 'CL'))
    tr = 35.453; %(g mol-1)
elseif (strcmp(trace, 'CH4'))
    tr = 16.04; %(g mol-1)
end

if (~isempty(varargin))
    if (strcmp(varargin(1), 'ppm'))
        conv=1e-6;
        conc=1/k*1e-6*conv*vmr.*p./t; %element wise arithmetic.
    elseif (strcmp(varargin(1), 'ppb')) 
        conv=1e-9;
        conc=1/k*1e-6*conv*vmr.*p./t;
    elseif (strcmp(varargin(1), 'ppt')) 
        conv=1e-12;
        conc=1/k*1e-6*conv*vmr.*p./t;
    elseif (strcmp(varargin(1), 'mmr'))
        conv=1;
        conc=((air/(tr))*(1e-6*conv*((vmr.*p)./t)))/k;%converts from mmr (kg/kg)
        %conc=1/k*1e-6*conv*vmr.*p./t;
        %=((1e-6*conv*((vmr.*p)./t)))/k;
    elseif (strcmp(varargin(1), 'conc')) 
        conv=1;
        conc=1/k*1e-6*conv*vmr.*p./t;
    elseif (strcmp(varargin(1), 'ice'))
        conc=1/R*vmr*p/t; %mol.m-3
        conc=molmasswater*conc; %g.m-3
    end
else
    conv=1;
    conc=1/k*1e-6*conv*vmr*p/t;
end
