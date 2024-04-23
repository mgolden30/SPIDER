%{
PURPOSE:
This script will generate a feature matrix G from a 6th order numerical
solution of moidified KS dynamics.
%}



%% Load data
clear;
restoredefaultpath();

addpath("library_generation/");

load('data/trajectory.mat');

%Add noise to u
u = add_noise( u, 1e-7 );

%Seeds for two different feature matrices
seeds = [1, 543212345];

for seed = seeds

%Integrate library
addpath("../SPIDER_functions/");

number_of_library_terms = 3;    %under-estimate this
number_of_windows       = 1024; %number of domains we integrate over 
degrees_of_freedom      = 1;    %vectors have one degree of freedom
dimension               = 2;    %how many dimensions does our data have?
envelope_power          = 8;    %weight is (1-x^2)^power
size_vec                = [16, 512]; %how many gridpoints should we use per integration?
buffer                  = 0;    %Don't use points this close to boundary

%define shorthand notation
nl = number_of_library_terms;
nw = number_of_windows;
dof= degrees_of_freedom;

%Make important objects for integration
pol      = envelope_pol( envelope_power, dimension );
G        = zeros( dof*nw, nl );
labels   = cell(nl, 1);
scales   = zeros(1,nl);

size_of_data = size(u, 1:dimension);
corners = pick_subdomains_manual_seed( size_of_data, size_vec, buffer, nw, seed );

%Create a numerical x and t grid
x = (0:N-1)/N*L;
t = linspace(0,T,M+1);
grid = {x,t};
dt = T/(M-1); %needed to estimate timescale


%%%%%%%%%%%%%%%%%%%%%%
% LIBRARY INTEGRATION
%%%%%%%%%%%%%%%%%%%%%%

addpath("library_generation/");

%This controls the complexity of terms in the library
symbol_length = 10; 
n_f = 1; %number of fields: u
n_d = 1; %number of dimensions: x (to differentiate)
base = n_f + n_d + 1; %used in counting library terms


%Precompute spatial derivatives
du = compute_spatial_derivatives(u, symbol_length, N, L);

a = 1; %running index over library

%First library term is (\partial_t u)
labels{a} = "\partial_t u";
G(:,a)    = SPIDER_integrate( u, [2], grid, corners, size_vec, pol );
scales(a) = std(u,0,"all")/estimate_timescale(u, dt);
a = a+1;

for i = 1:base^(symbol_length)
  [valid, fields, derivs, digits] = check_library_term( i, n_f, n_d, symbol_length );
  if( valid )
    %Add to library
    [term, str, scale] = generate_library_term( u, du, derivs);
    G(:,a)    = SPIDER_integrate( term, [], grid, corners, size_vec, pol );
    scales(a) = scale;
    fprintf(str + "\n");
    labels{a} = str;
    a = a+1;
  end
end

%normalize with integral of unity
norm_vec = SPIDER_integrate( 0*u + 1, [], grid, corners, size_vec, pol );      
G        = G./norm_vec;

%Nondimensionalize
G = G./scales;

save("seed_" + seed + ".mat", "G", "labels", "scales");

end

function du = compute_spatial_derivatives(u, nd, N, L)
  %{
  PURPOSE:
  Precompute a few spatial derivatives of u so that we can quickly
  construct library terms as needed. Store in a cell du such that du{i} 
  is the ith spatial derivative of u.
  %}

  %spatial wavevector
  k = 0:N-1;
  k(k>N/2) = k(k>N/2) - N;
  k = k';
  k = 2*pi/L*k;

  du = cell(nd,1);
  for l = 1:nd
    du{l} = real(ifft( (1i * k).^l .* fft(u) ));
  end
end

function time_scale = estimate_timescale(u, dt)
  %For nondimensionalization of the time derivative
  u_t = (circshift(u,-1,2) - u)/dt; %finite difference in time
  %ignore endpoints of u_t since not periodic in time
  time_scale = mean(abs(u),"all") / mean(abs(u_t(:,2:end-1)), "all");
end

function u = add_noise( u, amplitude )
  %Add white Gaussian noise to our numerical solution

  seed = 1;
  rng(seed);

  cr = 2*rand( size(u) ) - 1;
  ci = 2*rand( size(u) ) - 1;
  noise = real(ifft2( cr + 1i*ci ));

  %affine transformation to zero mean and unity std
  noise = noise - mean(noise, "all");
  noise = noise / std(  noise,0,"all");

  %{
  std(noise, 0,"all")  
  mean( noise,  "all")
  histogram( noise, 100 );
  %}

  u = u + amplitude*noise;
end