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
b = 0.15;
v = @(x) sin(circshift(x,-1)) - b*x;

dt =  0.05;  %timetep
n  = 1024*32; %steps

%Do RK4 to fill out a trajectory
data = zeros( n, 3 );
data(1,:) = [1;0;0]; %initial data

for i = 2:n
  k1 = dt*v(data(i-1,:) );
  k2 = dt*v(data(i-1,:) +k1/2);
  k3 = dt*v(data(i-1,:) +k2/2);
  k4 = dt*v(data(i-1,:) +k3);

  data(i,:) = data(i-1,:) + (k1 + 2*k2 + 2*k3 + k4)/6;
end

%Name each coordinate
x = data(:,1);
y = data(:,2);
z = data(:,3);

ms = 1; %marker size
scatter3( x, y, z, ms, z, 'o');

xlim([-6,6]);
ylim(xlim);
zlim(ylim);
pbaspect([1 1 1]);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2: compute a library matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath("../SPIDER_functions/");

number_of_library_terms = 2;    %under-estimate this
number_of_windows       = 128;  %number of domains we integrate over 
degrees_of_freedom      = 1;    %scalars have one degree of freedom
dimension               = 1;    %how many dimensions does our data have?
envelope_power          = 4;    %weight is (1-x^2)^power
size_vec                = [128]; %how many gridpoints should we use per integration?
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

size_of_data = size(data, 1:dimension);
seed = 1;
corners = pick_subdomains_manual_seed( size_of_data, size_vec, buffer, nw, seed );

t = dt*(0:n-1);
grid = { t };
a = 1; %running index over library
% END BOILERPLATE



labels{a} = "dx/dt";
G(:,a)    = SPIDER_integrate( x, [1], grid, corners, size_vec, pol );
scales(a) = mean(abs(x));
a = a+1;

labels{a} = "dy/dt";
G(:,a)    = SPIDER_integrate( y, [1], grid, corners, size_vec, pol );
scales(a) = mean(abs(y));
a = a+1;

labels{a} = "dz/dt";
G(:,a)    = SPIDER_integrate( z, [1], grid, corners, size_vec, pol );
scales(a) = mean(abs(z));
a = a+1;


labels{a} = "x";
G(:,a)    = SPIDER_integrate( x, [], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} = "y";
G(:,a)    = SPIDER_integrate( y, [], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} = "z";
G(:,a)    = SPIDER_integrate( z, [], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;



labels{a} = "sin(x)";
G(:,a)    = SPIDER_integrate( sin(x), [], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} = "sin(y)";
G(:,a)    = SPIDER_integrate( sin(y), [], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} = "sin(z)";
G(:,a)    = SPIDER_integrate( sin(z), [], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

%normalize the feature matrix.
norm_vec = SPIDER_integrate( 0*x + 1, [], grid, corners, size_vec, pol );      
G        = G./norm_vec;
%Nondimensionalize
G = G./scales;









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 3: sparse regression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gamma = 1.25;

%search for three models
for i = 1:3
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
