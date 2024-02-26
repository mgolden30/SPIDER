
%% Load data
clear;
restoredefaultpath();

addpath("library_generation/");

load('trajectory.mat');

% Add noise to u
%%
amplitude = 0*1e-7;
noise = gaussian_noise( u );

u = u + amplitude*noise;


seeds = [1, 543212345];

for seed = seeds

%Integrate library
addpath("../SPIDER_functions/");

number_of_library_terms = 3;   %under-estimate this
number_of_windows       = 1024; %number of domains we integrate over 
degrees_of_freedom      = 1;   %vectors have one degree of freedom
dimension               = 2;   %how many dimensions does our data have?
envelope_power          = 8;   %weight is (1-x^2)^power
size_vec                = [16, 512]; %how many gridpoints should we use per integration?
buffer                  = 0; %Don't use points this close to boundary

%define shorthand notation
nl = number_of_library_terms;
nw = number_of_windows;
dof= degrees_of_freedom;

%Make important objects for integration
pol      = envelope_pol( envelope_power, dimension );
G        = zeros( dof*nw, nl );
labels   = cell(nl, 1);
scales   = zeros(1,nl);

%we also need pol to be odd in time. multiply it by t
pol0 = pol; %save for normalization

size_of_data = size(u, 1:dimension);
corners = pick_subdomains_manual_seed( size_of_data, size_vec, buffer, nw, seed );

grid = { (0:N-1)/N*L, linspace(0,T,M+1)};

x = grid{1};
t = grid{2};
dxs = [ 1, 1];

%Nondimensionalize and compute mean and variation of data
%TODO;

%%%%%%%%%%%%%%%%%%%%%%
% LIBRARY INTEGRATION
%%%%%%%%%%%%%%%%%%%%%%
a = 1;

dt = T/(M-1);

nd = 10;%num derivs
du = cell(nd,1);
for l = 1:nd
  k = 0:N-1;
  k(k>N/2) = k(k>N/2) - N;
  k = k';
  k = 2*pi/L*k;
  du{l} = real(ifft( (1i * k).^l .* fft(u) ));
end

u_t = (circshift(u,-1,2) - u)/dt; %finite difference in time
u_t(:,1) = 0;
u_t(:,end) = 0;
time_scale = mean(abs(u),"all") / mean(abs(u_t), "all");

labels{a} =  "\partial_t u";
G(:,a)    = SPIDER_integrate( u, [2], grid, corners, size_vec, pol0 );
scales(a) = std(u,0,"all")/time_scale;
a = a+1;

addpath("library_generation/");
%All other library terms are computed with the auto library generation
for i = 1:3^(nd)
  [canonical, derivs, v] = check_library_term(i);
  if(canonical)
    %Add to library
    [term, str, scale] = generate_library_term(u,du,derivs);
    G(:,a)    = SPIDER_integrate( term, [], grid, corners, size_vec, pol0 );
    scales(a) = scale;
    labels{a} = str;
    a = a+1;
  end
end

%normalize with non-odd polynomial.
norm_vec = SPIDER_integrate( 0*u + 1, [], grid, corners, size_vec, pol0 );      

G = G./norm_vec;

%scales = vecnorm(G);
G = G./scales;

save("seed_" + seed + ".mat", "G", "labels", "scales");

end


function noise = gaussian_noise( u )
  seed = 1;
  rng(seed);

  cr = 2*rand( size(u) ) - 1;
  ci = 2*rand( size(u) ) - 1;
  noise = real(ifft2( cr + 1i*ci ));

  noise = noise / max(max(abs(noise)));
end