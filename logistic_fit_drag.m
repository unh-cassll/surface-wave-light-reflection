% Logistic fit to drag coefficient Cd data from:
% Edson et al. (2013) [U10 < 25]
% Curcic & Haus (2020) [U10 > 25]
%
% Nathan Laxague, 2020
%
function [Cd_out,forcing_out] = logistic_fit_drag(input_forcing,forcing_type)

input_forcing = input_forcing(:)';

switch forcing_type
    
    case 'U10'
        
        U10 = input_forcing;
        
    case 'ustar'
        
        U10 = (0.01:0.01:60);
        
end

% Logistic fit parameters, U10 -> Cd
U10_offset = 17;
s_curve_strength = 0.2;
Cd_min = 0.8e-3;
Cd_max = 3.2e-3;
Cd_span = Cd_max - Cd_min;

Cd_out = Cd_span*(1+exp(-s_curve_strength*(U10-U10_offset))).^-1+Cd_min;

ustar_out = sqrt(Cd_out).*U10;

if strcmp(forcing_type,'U10')
    
    forcing_out = ustar_out;
    
else
    
    forcing_out = NaN*input_forcing;
    
    for i = 1:length(input_forcing)
        
        ustar_diff = abs(input_forcing(i)-ustar_out);
        ind = find(ustar_diff==min(ustar_diff),1,'first');
        forcing_out(i) = U10(ind);
        
    end
    
    Cd_out = (input_forcing./forcing_out).^2;
    
end
