function [areaout] = gridcellarea(latin,lonin)

[LON,LAT] = meshgrid(lonin,latin);
areaout = cdtarea(LAT,LON,'km2');

end