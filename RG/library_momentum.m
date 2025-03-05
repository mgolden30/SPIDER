%{
PURPOSE:
This script will generate a feature matrix G from a 6th order numerical
solution of moidified KS dynamics.
%}


tic
%% Load data
clear;
restoredefaultpath();


load('data.mat');
rho = rho_rg;

p = estimate_momentum( rho );
p = imgaussfilt(p, 1);

%Visualize the data
tiledlayout(1,2);

buffer = 3;
rho = rho(buffer:end-buffer,buffer:end-buffer);
p = p(buffer:end-buffer,buffer:end-buffer);

nexttile
imagesc(rho);
title("\rho")

nexttile
imagesc(p);
title("p")



%% Integrate library
addpath("../SPIDER_functions/");

number_of_library_terms = 3;    %under-estimate this
number_of_windows       = 1024; %number of domains we integrate over 
degrees_of_freedom      = 1;    %vectors have one degree of freedom
dimension               = 2;    %how many dimensions does our data have?
envelope_power          = 8;    %weight is (1-x^2)^power
size_vec                = [32, 32]; %how many gridpoints should we use per integration?
buffer                  = 3;    %Don't use points this close to boundary

%define shorthand notation
nl = number_of_library_terms;
nw = number_of_windows;
dof= degrees_of_freedom;

%Make important objects for integration
pol      = envelope_pol( envelope_power, dimension );
G        = zeros( dof*nw, nl );
labels   = cell(nl, 1);
scales   = zeros(1,nl);

size_of_data = size(rho, 1:dimension);
seed = 1;
corners = pick_subdomains_manual_seed( size_of_data, size_vec, buffer, nw, seed );

%Create a numerical x and t grid
nx = size(p,2);
nt = size(p,1);
x = linspace(0,1,nx);
t = linspace(0,1,nt);
grid = {t,x};


%%%%%%%%%%%%%%%%%%%%%%
% LIBRARY INTEGRATION
%%%%%%%%%%%%%%%%%%%%%%


a = 1; %running index over library

labels{a} = "\partial_t p";
G(:,a)    = SPIDER_integrate( p, [1], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} = "\partial_x( p^2)";
G(:,a)    = SPIDER_integrate( p.^2, [2], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} = "\partial_x( p^2 / rho)";
G(:,a)    = SPIDER_integrate( p.^2./rho, [2], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;


labels{a} = "\partial_x( \rho^2)";
G(:,a)    = SPIDER_integrate( rho.^2, [2], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} = "\partial_x \rho";
G(:,a)    = SPIDER_integrate( rho, [2], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} = "\partial_x^2 \rho";
G(:,a)    = SPIDER_integrate( rho, [2,2], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} = "\partial_x^3 \rho";
G(:,a)    = SPIDER_integrate( rho, [2,2,2], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} = "\partial_x p";
G(:,a)    = SPIDER_integrate( p, [2], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} = "\partial_x^2 p";
G(:,a)    = SPIDER_integrate( p, [2,2], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} = "\partial_x^3 p";
G(:,a)    = SPIDER_integrate( p, [2,2,2], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} = "\rho";
G(:,a)    = SPIDER_integrate( rho, [], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} = "\rho^2";
G(:,a)    = SPIDER_integrate( rho.^2, [], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} = "p";
G(:,a)    = SPIDER_integrate( p, [], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} = "p.^2";
G(:,a)    = SPIDER_integrate( p.^2, [], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

%normalize with integral of unity
norm_vec = SPIDER_integrate( 0*p + 1, [], grid, corners, size_vec, pol );      
G        = G./norm_vec;

%Nondimensionalize
G = G./scales;

%save("G_momentum.mat", "G", "labels", "scales");

%%
%load("G_momentum.mat");


%Split up the library into a target term b
b = G(:, 1);
target = labels{1};
target_scale = scales(1);

%And the rest of the library
A = G(:, 2:end);
labels = labels(2:end);
scales(1) = [];

tic
[cs, residuals] = inhomogeneous_regression( A, b );
toc

scatter( 1:numel(residuals), residuals, 'o', 'filled' );

toc


%normalization = 1;
for interest=1:numel(residuals)
  c = cs(:,interest); 
  str = interest + ":  residual = " + residuals(interest) +  ", \quad   " + target + " = ";
  for i = 1:numel(c)
    if( c(i) == 0 )
      continue; 
    end

   str = str + " + " + c(i)  + labels{i};

  end

  %For printing, replace \ with \\
  str = strrep(str, "\", "\\");

  fprintf(str + "\n");
end
