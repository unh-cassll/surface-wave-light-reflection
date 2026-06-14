% Computes unified omnidirectional wavenumber spectrum
% Elfouhaily et al. (1997)
%
% Coded by N. Laxague between 2013-2024
%
function [out_K,out_S] = Elfouhaily_omni(U10,min_k,max_k,num_k,fetch_m)

if nargin < 5
    fetch_m = 1e6;
end

% Define constants
g = 9.81; %m/s/s
sigma = 0.072;
rho = 1030;

% Compute drag/friction velocity
%CD = 1.2*10^-3;
%CD = 1e-3*(0.8+0.065*U10);
%u_star = sqrt(CD).*U10;
[~,u_star] = logistic_fit_drag(U10,'U10');

% Compute wavenumber, frequency, and phase speed
k = logspace(log10(min_k),log10(max_k),num_k)';     % rad/m
ang_f = sqrt(g*k+sigma/rho*k.^3);                   % rad/s
c = ang_f./k;                                       % phase speed

% Capillarity effects
cm = 0.23;                                          % minimum phase speed
km = 2*g/(cm.*cm);                                  % wavenumber at which cm occurs

% Retrieve wavenumber peak
k0 = g./(U10.*U10);                                  % rad/m
X0 = 2.2e4;                                         % "empirical constant"
omega_c = 0.84*tanh((fetch_m*k0/X0).^0.4).^-0.75;     % critical inverse wave age
kp = k0.*omega_c.^2;                                  % peak wavenumber
ang_fp = sqrt(g*kp);                                % peak angular frequency
cp = ang_fp./kp;                                     % phase speed at peak wavenumber
omega = U10./cp;                                     % inverse wave age
gamma = 1.7;                                        % JONSWAP parameter
% gamma(omega_c>1) = 1.7 + 6*log10(omega_c(omega_c>1));
alpha_p = 0.006*omega_c.^0.55;
sigma = 0.08*(1+4*omega_c.^-3);

% Long wave portion
Lpm = exp(-5/4*(kp./k).^2);                         % Pierson-Moskowitz shape spectrum
big_gamma = exp(-(sqrt(k./kp)-1).^2./(2*sigma.^2));
Jp = gamma.^big_gamma;                              % JONSWAP peak enhancement
Fp = Lpm.*Jp.*exp(-omega/sqrt(10).*(sqrt(k./kp)-1));  % Long wave side effect function
Bl = 0.5*alpha_p.*cp./c.*Fp;

% Short wave portion
alpha_m = 1 + log(u_star/cm);
alpha_m(u_star>cm) = 1 + 3*log(u_star(u_star>cm)/cm);
alpha_m = alpha_m*10^-2;
Fm = exp(-1/4*(k/km-1).^2);
Bh = 0.5*alpha_m.*cm./c.*Fm;

% Trim off short wave contribution for low wavenumbers
B_rat = Bl./Bh;
B_rat(k>100,:) = NaN;
Bl(k>2e3,:) = 0;
for n = 1:length(U10)
    ind = find(B_rat(:,n) > 0.9*max(B_rat(:,n),[],'all','omitnan'),1,'first');
    Bh_val = Bh(ind,n);
    Bh(1:ind,n) = exp(-0.5*(k(ind)./k(1:ind)).^2).*Bh(1:ind,n);
    Bh(1:ind,n) = Bh(1:ind,n)*Bh_val/Bh(ind,n);
end

% Combine portions of spectra
S = k.^-3.*(Bl+Bh);

out_K = k;
out_S = S;
