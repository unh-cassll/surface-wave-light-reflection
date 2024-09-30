% Given slope variances in two dimensions and number of standard deviations
% over which to constrain distribution, produces slope PDF
function [Pxy,sx,sy] = produce_slope_distribution(slope_var_x,slope_var_y,num_standard_dev)

sx = linspace(-sqrt(slope_var_x)/2*num_standard_dev,sqrt(slope_var_x)/2*num_standard_dev,100)';
sy = linspace(-sqrt(slope_var_y)/2*num_standard_dev,sqrt(slope_var_y)/2*num_standard_dev,100);

Pxy = 1/(2*pi*sqrt(slope_var_x)*sqrt(slope_var_y))*exp(-0.5*(sx.^2/slope_var_x+sy.^2/slope_var_y));

sx = repmat(sx,1,length(sy));
sy = repmat(sy,length(sy),1);

