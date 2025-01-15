%{

A tutorial on using model discovery on 1D data.

Press ctrl + enter to run individual sections

%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 1: generate some time depenendent data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;

% Simulate the Thomas attractor
% https://en.wikipedia.org/wiki/Thomas%27_cyclically_symmetric_attractor

dt = 0.005;  
skip = 16; %record every skip timesteps
nu = 1e-2;
n  = 128;  %grid resolution
timesteps = 128;

[x,y,t,u,v,p] = navier_stokes( n, dt, timesteps, nu, skip );




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2: compute a library matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath("../SPIDER_functions/");

number_of_library_terms = 2;    %under-estimate this
number_of_windows       = 128;  %number of domains we integrate over 
degrees_of_freedom      = 1;    %scalars have one degree of freedom
dimension               = 3;    %how many dimensions does our data have?
envelope_power          = 4;    %weight is (1-x^2)^power
size_vec                = [64,64,32]; %how many gridpoints should we use per integration?
buffer                  = 0;    %Don't use points this close to boundary

% BEGIN BOILERPLATE
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
seed = 1;
corners = pick_subdomains_manual_seed( size_of_data, size_vec, buffer, nw, seed );

grid = { y,x,t };
a = 1; %running index over library
% END BOILERPLATE



labels{a} = "d_x u";
G(:,a)    = SPIDER_integrate( u, [2], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} = "d_y v";
G(:,a)    = SPIDER_integrate( v, [1], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} = "d_t u";
G(:,a)    = SPIDER_integrate( u, [3], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} = "d_x(uu)";
G(:,a)    = SPIDER_integrate( u.^2, [2], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} = "d_y(uv)";
G(:,a)    = SPIDER_integrate( u.*v, [1], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} = "d_xx u";
G(:,a)    = SPIDER_integrate( u, [2,2], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} = "d_yy u";
G(:,a)    = SPIDER_integrate( u, [1,1], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} = "d_x p";
G(:,a)    = SPIDER_integrate( p, [2], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

%HOMEWORK: find the time derivative of v

%normalize the feature matrix.
norm_vec = SPIDER_integrate( 0*u + 1, [], grid, corners, size_vec, pol );      
G        = G./norm_vec;
%Nondimensionalize
G = G./scales;









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 3: sparse regression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gamma = 1.25;

%search for three models
for i = 1:2
  %Do sparse regression
  [cs, residuals] = greedy_regression_pure_matlab( G );

  %Print the discovered model
  k = report_identified_model(cs, residuals, scales, labels, gamma);

  %remove the most important term
  [~, kill] = max( vecnorm(G*diag(cs(:,k))) );
  G(:, kill) = [];
  labels(kill) = [];
  scales(kill) = [];
end
