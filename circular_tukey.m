%
% Circular Tukey (tapered cosine) window
% Works for 3D arrays such that [x,y,z] array will be tapered in x and y,
% repeated along z axis in cylindrical coordinates.
% Additional cosine taper in z coordinate
%
% Nathan J. M. Laxague, 4/2019
%
function out_array = circular_tukey(in_array,taper_width)

% Obtain size of input array, find minimum horizontal dimension
[s1,s2,s3] = size(in_array);
min_dim = min([s1 s2]);

% Build radial array
[x,y] = meshgrid(1:s2,1:s1);
x = x - mean(mean(x));
y = y - mean(mean(y));
r = sqrt(x.^2+y.^2);
r = r/(min_dim/2);

% Compute cosine taper portion
cos_period = taper_width*2;
rshift = 1 - taper_width;
cos_mat = (cos(2*pi/cos_period*(r-rshift))+1)/2;

% Flatten window elsewhere
w2D = cos_mat;
w2D(r<(1-taper_width)) = 1;
w2D(r>1) = 0;

% If input is three-dimensional, build 3D window
if s3 > 1
   
    w = repmat(w2D,[1 1 s3]);
    
else
    
    w = w2D;
    
end

% Create temporal window
% t = (1:s3*taper_width)'/(s3*taper_width)*2*pi - pi;
% t = t(1:floor(length(t)/2));
% w_t = [(cos(t)+1)/2; ones(s3-2*length(t),1); flip((cos(t)+1)/2)];
% w_t_rep = permute(repmat(w_t,[1 s1 s2]),[2 3 1]);
% w = w.*w_t_rep;

% Compute norm of window
test_array = ones(s1,s2);
out_array = mean(w,3);
C = (norm(test_array)/norm(out_array));

% Output
out_array = C*w.*in_array;