%
% Random phase approach to simulation of realistic sea surface
% Methodology rests on the shoulders of many, especially J. Tessendorf
%
% Script executing water surface wave simulation
% v1 - 2015/06/02 - Hoki
% MIT licence
%
% The original full function is shared at the following link:
% https://www.dropbox.com/s/qow7cdf5z95t7hx/ocean_simulator.m?dl=0
%
% v2 - 2023 Nathan Laxague
% *incorporation of Elfouhaily spectrum
% *enforcement to keep Hs of original spectrum

function [x,y,eta] = ocean_simulator_compact_elfouhaily(U10,grid_size,dx,nframes,fps)

warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId')

x = (0:grid_size-1)*dx;
y = x;

rng(13);

n = grid_size * [1 1] ;
[x_grid,y_grid] = meshgrid(1:grid_size,1:grid_size) ;
Grid_Sign = ones( n ) ;
Grid_Sign(mod(x_grid+y_grid,2)==0) = -1 ;

[H0,omega,Hs] = wave_height_spectrum(U10,grid_size,dx);

s = struct();

for i = 1:nframes

    time = i/fps;

    Z = calc_wave(H0,omega,time,Grid_Sign);

    s(i).eta = Z;

end

eta = reshape([s.eta],[grid_size grid_size nframes]);

Hm0 = 4*mean(mean(std(eta,0,3)));
eta = eta*Hs/Hm0/2;

    function [H0,omega,Hs] = wave_height_spectrum(U10,grid_size,dx)

        num_wave_dir = 360;
        g = 9.81;
        sigma = 0.072;
        rho_w = 1030;

        [~,u_star] = logistic_fit_drag(U10,'U10');

        num_k = 1000;
        kmin = 2*pi/(grid_size*dx);
        kmax = pi/dx;

        [k,Fk] = Elfouhaily_omni(U10,kmin,kmax,num_k);

        Fk(Fk<0) = 0;

        Hs = 4*sqrt(trapz(k,Fk));

        % Repeat blocks that vary in wavenumber and wave direction
        wave_dir_vec = linspace(-180,180,num_wave_dir);
        k_block = repmat(k,[1 num_wave_dir]);
        Fk_block = repmat(Fk,[1 num_wave_dir]);
        wave_dir_block = repmat(wave_dir_vec,[num_k 1]);

        % Compute wave frequency and celerity given wavenumber and current profile
        %Doppler_block = k_block.*U_current_block.*cosd(dir_current_rel_wind_block);
        omega_block = sqrt(g*k_block+sigma/rho_w*k_block.^3);%+Doppler_block;
        cp_block = omega_block./k_block;
        waveage_block = cp_block./U10;

        % Compute directional spreading function
        [deltak] = Elfouhaily_spread(k_block,U10,waveage_block,u_star);
        s_block = atanh(deltak)*2/log(2);
        spreading_block = real(cosd(wave_dir_block).^(2*s_block)/(2*pi));

        % Produce 'F_block', the directional wavenumber elevation variance spectrum
        F_block = k_block.^-1.*Fk_block.*spreading_block;

        [kx,ky] = pol2cart(mod(wave_dir_block+90,360)*pi/180,k_block);

        kx_vec = reshape(kx,[],1);
        ky_vec = reshape(ky,[],1);
        F_vec = reshape(F_block,[],1);

        [kx_new,ky_new] = meshgrid(-kmax:kmin:kmax,-kmax:kmin:kmax);
        kx_new(:,1) = [];
        ky_new(:,1) = [];
        kx_new(end,:) = [];
        ky_new(end,:) = [];
        ky_new = flipud(ky_new);
        k_new = sqrt(kx_new.^2+ky_new.^2);

        omega = sqrt(g*k_new+sigma/rho_w*k_new.^3);

        S = scatteredInterpolant(kx_vec,ky_vec,F_vec);
        F_block_new = abs(S(kx_new,ky_new));

        F_block_filt = 10.^medfilt2(log10(F_block_new),[5 5]);

        k1 = 20;

        w1 = 1*(1+exp(-k1*(ky_new-kmin/2))).^-1;
        w2 = 0*w1 + 1;
        w2 = circular_tukey(w2,0.5);
        w = w1.*w2;
        w = w*sum(0*w+1,'all')/sum(w,'all');
        
        F_block_filt = F_block_filt.*w;

        F_block_filt(isnan(F_block_filt)) = min(F_block_filt,[],'all','omitnan');

        randmat1 = randn(grid_size);
        randmat2 = randn(grid_size);

        H0 = 1/sqrt(2) .* (randmat1 + 1i .* randmat2) .* sqrt(F_block_filt); % height field at time t = 0

    end

    function Z = calc_wave( H0,omega,time,Grid_Sign )
        wt = exp(1i .* omega .* time ) ;
        Ht = H0 .* wt + conj(flipud(H0)) .* conj(wt) ;
        Z = real( ifft2(Ht) .* Grid_Sign ) ;
    end

end