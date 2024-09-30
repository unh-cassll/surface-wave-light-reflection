% Simulates a sea surface and models the polarization state of light
% reflected from that sea surface
%
% Inputs:
% * wind speed
% * grid size
% * spatial resolution
% * temporal sampling interval
% * number of frames desired
%
% Outputs:
% * horizontal displacements (x,y)
% * vertical surface displacement field (wse_m_stack)
% * Stokes parameters (S0_stack, S1_stack, S2_stack)
%
% Polarized light reflection model is a simplified (and evolving)
% adaptation of the technique described by
% Mobley, C. D. (2015). "Polarized reflectance and transmittance properties
% of windblown sea surfaces". Applied optics, 54(15), 4828-4849.
%
% Code by N. Laxague 2024
%
function [x,y,wse_m_stack,S0_stack,S1_stack,S2_stack] = produce_simulated_sea_surface_modeled_reflection(U10_m_s,grid_size,spatial_resolution_m,temporal_interval_s,nframes);

t_on = tic;

% Model sea surface
[x,y,wse_m_stack] = ocean_simulator_compact_elfouhaily(U10_m_s,grid_size,spatial_resolution_m,nframes,1/temporal_interval_s);

t_off = toc(t_on);

dstr_now = datestr(now,'mm/dd/yyyy HH:MM:SS');
tstr = datestr(t_off/86400,'HH:MM:SS');

disp([dstr_now '... WAVE SURFACE SIMULATED AFTER ' tstr])

t_on = tic;

% Incident Stokes vector
S_inc1 = 1;
S_inc2 = 0;
S_inc3 = 0;
S_inc4 = 0;

% Solar angles (for future use)
% sun_polar_ang = 0;
% sun_az_ang = 0;

% Real index of refraction
n = 1.34;

% Virtual camera parameters
camera_incidence_deg = 30;
camera_azimuth_deg = 0;
camera_height_m = 100;
camera_location = [sind(camera_azimuth_deg+180)*tand(camera_incidence_deg) cosd(camera_azimuth_deg+180)*tand(camera_incidence_deg) 1]'*camera_height_m;

% Compute unresolved slope variance (below the subpixel)
min_k = 2e-3;
max_k = 2e4;
num_k = 1000;
fetch_m = 1e6;
[k,Fk] = Elfouhaily_omni(U10_m_s,min_k,max_k,num_k,fetch_m);
k_high = pi/spatial_resolution_m;
slope_var = trapz(k(k>k_high),k(k>k_high).^2.*Fk(k>k_high));
slope_var_x = slope_var/2;
slope_var_y = slope_var/2;

% Produce subpixel slope distributions
num_standard_dev = 3;
[Pxy,sx,sy] = produce_slope_distribution(slope_var_x,slope_var_y,num_standard_dev);
Pxy_copy = floor(Pxy);
Pxy_copy = Pxy_copy/min(Pxy_copy,[],'all');
slope_PDF_struc = struct();
increment = 0;
for i = 1:size(sx,1)
    for j = 1:size(sx,2)
        Pxy_count = Pxy_copy(i,j);
        for counter = 1:Pxy_count
            increment = increment + 1;
            slope_PDF_struc(increment).sx = sx(i,j);
            slope_PDF_struc(increment).sy = sy(i,j);
        end
    end
end
sx_vec = [slope_PDF_struc.sx];
sy_vec = [slope_PDF_struc.sy];
num_rays_per_ensemble = 100;

% Preallocate Stokes vector components
S0_stack = NaN*wse_m_stack;
S1_stack = S0_stack;
S2_stack = S0_stack;

t_off = toc(t_on);

dstr_now = datestr(now,'mm/dd/yyyy HH:MM:SS');
tstr = datestr(t_off/86400,'HH:MM:SS');

disp([dstr_now '... SUBPIXEL SLOPE DISTRIBUTION COMPUTED AFTER ' tstr])

for frame_ind = 1:nframes

    t_on = tic;

    % Pluck out a particular time step's surface wave field
    wse_m = squeeze(wse_m_stack(:,:,frame_ind));

    % Compute the surface normal field
    [gx,gy] = gradient(wse_m);
    N_mat = NaN*repmat(wse_m,[1 1 3]);
    N_mat(:,:,1) = gx/spatial_resolution_m;
    N_mat(:,:,2) = gy/spatial_resolution_m;
    N_mat(:,:,3) = 0*gy + 1;
    N_mat = N_mat./sqrt(sum(N_mat.^2,3));
    nx = N_mat(:,:,1);
    ny = N_mat(:,:,2);
    nz = N_mat(:,:,3);

    % Set up and run loop over all surface facets
    nx_vec = reshape(nx,[],1);
    ny_vec = reshape(ny,[],1);
    nz_vec = reshape(nz,[],1);
    X_vec = reshape(repmat(x,length(y),1),[],1);
    Y_vec = reshape(repmat(y',1,length(x)),[],1);
    Z_vec = reshape(wse_m,[],1);
    [total_rows,total_cols] = size(wse_m);
    total_inds = total_rows*total_cols;

    image_holder_struc = struct();

    parfor ind = 1:total_inds

        % Grab surface normal vector and facet position
        N0 = [nx_vec(ind) ny_vec(ind) nz_vec(ind)]';
        X = X_vec(ind);
        Y = Y_vec(ind);
        Z = Z_vec(ind);

        % Compute reflected ray
        ray_refl = camera_location - [X Y Z]';
        ray_refl = ray_refl/norm(ray_refl);

        % Add some Gaussian noise in slope
        rand_ind = rand(num_rays_per_ensemble,1);
        rand_ind = floor((rand_ind-min(rand_ind))/(max(rand_ind)-min(rand_ind))*(length(sx_vec)-1))+1;
        wiggle_x = sx_vec(rand_ind);
        wiggle_y = sy_vec(rand_ind); %#ok<PFBNS>
        wiggle_z = 0*wiggle_x;

        N = N0 + [wiggle_x;wiggle_y;wiggle_z];
        N = N./sqrt(sum(N.^2));

        % Compute incident ray
        ray_inc = ray_refl - (2*ray_refl'*N).*N;

        % Compute incidence angle with respect to surface facet
        inc_dot_N = sum(ray_inc.*N,1);
        theta_i = abs(acosd(inc_dot_N*-1));
        theta_t = asind(sind(theta_i)./n);

        % Enforce sky-leaving radiance conditions (todo)
        % polar_ang = acosd([0 0 1]*N);
        % az_ang = atan2(N(1,:),N(2,:))*180/pi;

        % NEED TO REPLACE THIS WITH A SKY RADIANCE MODEL
        %weight = cosd(polar_ang-sun_polar_ang).*cosd(az_ang-sun_az_ang);
        weight = 0*N(1,:) + 1;

        % Ensure that ray is coming from above horizon
        % NEED TO REPLACE THIS WITH MULTIPLE REFLECTION CALCULATION
        pass_cond = ray_inc(3,:) < 0;
        N = N(:,pass_cond);
        ray_inc = ray_inc(:,pass_cond);
        theta_i = theta_i(pass_cond);
        theta_t = theta_t(pass_cond);
        weight = weight(pass_cond);

        if ~isempty(N)

            % Compute h and s vectors (Mobley 2015)
            h_inc = cross(repmat([0;0;1],1,length(theta_i)),ray_inc);
            h_inc = h_inc./sqrt(sum(h_inc.^2));
            s = cross(ray_inc,N);
            s = s./sqrt(sum(s.^2));
            h_cross_s = cross(h_inc,s);

            % Compute Stokes vector rotation angle alpha (Mobley 2015)
            h_inc_dot_s = sum(h_inc.*s,1);
            alpha_ang = abs(acosd(h_inc_dot_s));
            antiparallel_inds = h_cross_s(3,:).*ray_inc(3,:) < 0;
            alpha_ang(antiparallel_inds) = 360 - alpha_ang(antiparallel_inds);
            c2a = cosd(2*alpha_ang);
            s2a = sind(2*alpha_ang);

            % Rotation matrix, Stokes vector
            % R = [1,0,            0,             0;
            %      0,cosd(2*alpha),-sind(2*alpha),0;
            %      0,sind(2*alpha),cosd(2*alpha), 0;
            %      0,0,            0,             1];

            % Mueller matrix for reflection off surface, KA89
            alpha = 1/2*(tand(theta_i-theta_t)./tand(theta_i+theta_t)).^2;
            eta = 1/2*(sind(theta_i-theta_t)./sind(theta_i+theta_t)).^2;
            gamma_Re = (tand(theta_i-theta_t).*sind(theta_i-theta_t))./(tand(theta_i+theta_t).*sind(theta_i+theta_t));

            % R_AM = [alpha+eta,alpha-eta,0,       0;
            %         alpha-eta,alpha+eta,0,       0;
            %         0,        0,        gamma_Re,0;
            %         0,        0,        0,       gamma_Re];

            % Compute reflected Stokes vector
            %S_refl = R*R_AM*S_inc;
            S_refl_0 = S_inc1*(alpha+eta)+S_inc2*(alpha-eta);
            S_refl_1 = S_inc1*c2a.*(alpha-eta)+S_inc2*c2a.*(alpha+eta)-S_inc3*s2a.*gamma_Re;
            S_refl_2 = S_inc1*s2a.*(alpha-eta)+S_inc2*s2a.*(alpha+eta)+S_inc3*c2a.*gamma_Re;
            S_refl_3 = S_inc4*gamma_Re;

            S_refl = median([S_refl_0; S_refl_1; S_refl_2; S_refl_3].*weight,2,'omitnan')./median(weight,'omitnan');

            % Save Stokes vector in running array
            image_holder_struc(ind).S_refl = S_refl;
            image_holder_struc(ind).theta_i = median(theta_i.*weight)./median(weight);

        else

            % Save Stokes vector in running array
            image_holder_struc(ind).S_refl = NaN*ones(4,1);
            image_holder_struc(ind).theta_i = NaN;

        end

    end

    Stokes_array = permute(reshape([image_holder_struc.S_refl],4,total_rows,total_cols),[2 3 1]);

    % Extract computed Stokes vector field
    S0_stack(:,:,frame_ind) = real(squeeze(Stokes_array(:,:,1)));
    S1_stack(:,:,frame_ind) = real(squeeze(Stokes_array(:,:,2)));
    S2_stack(:,:,frame_ind) = real(squeeze(Stokes_array(:,:,3)));

    t_off = toc(t_on);

    dstr_now = datestr(now,'mm/dd/yyyy HH:MM:SS');
    tstr = datestr(t_off/86400,'HH:MM:SS');

    disp([dstr_now '... DONE WITH FRAME # ' sprintf(['%0' num2str(ceil(log10(nframes))) 'u'],frame_ind) '/' num2str(nframes) ' AFTER ' tstr])

end
