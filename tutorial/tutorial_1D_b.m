%{

A tutorial on using model discovery on 1D data.

Press ctrl + enter to run individual sections

This example is 1D in time, but assume the vector nature of the governing
equations.
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 1: generate some time depenendent data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;

% Simulate the 3-body problem
%Try 2, 5, 6,...
initial_conditions = readmatrix("three_body_solutions/25.txt");

dt =  0.0025;  %timetep
n  = 1024*8; %steps

%Do RK4 to fill out a trajectory
data = zeros( n, 12 );
data(1,:) = initial_conditions(1:end-1); %initial data

for i = 2:n
  k1 = dt*v(data(i-1,:) );
  k2 = dt*v(data(i-1,:) +k1/2);
  k3 = dt*v(data(i-1,:) +k2/2);
  k4 = dt*v(data(i-1,:) +k3);

  data(i,:) = data(i-1,:) + (k1 + 2*k2 + 2*k3 + k4)/6;
end

%Name each coordinate
r1 = data(:,1:3);
r2 = data(:,4:6);
r3 = -r1-r2;

r12 = r1-r2;
r13 = r1-r3;
r23 = r2-r3;

ms = 1; %marker size
scatter3( r1(:,1), r1(:,2), r1(:,3), ms, 'o', 'filled', 'markerfacecolor', 'red' );
hold on
scatter3( r2(:,1), r2(:,2), r2(:,3), ms, 'o', 'filled', 'markerfacecolor', 'blue' );
scatter3( r3(:,1), r3(:,2), r3(:,3), ms, 'o', 'filled', 'markerfacecolor', 'black' );
hold off

xlim([-2,2]);
ylim(xlim);
zlim(ylim);
pbaspect([1 1 1]);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2: compute a library matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath("../SPIDER_functions/");

number_of_library_terms = 2;    %under-estimate this
number_of_windows       = 128;  %number of domains we integrate over 
degrees_of_freedom      = 3;    %vectors have three degrees of freedom
dimension               = 1;    %how many dimensions does our data have?
envelope_power          = 4;    %weight is (1-x^2)^power
size_vec                = [1024]; %how many gridpoints should we use per integration?
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



labels{a} = "ddot{r}_1";
G(:,a)    = [SPIDER_integrate( r1(:,1), [1,1], grid, corners, size_vec, pol );
             SPIDER_integrate( r1(:,2), [1,1], grid, corners, size_vec, pol );
             SPIDER_integrate( r1(:,3), [1,1], grid, corners, size_vec, pol );];
scales(a) = 1;
a = a+1;

labels{a} = "ddot{r}_2";
G(:,a)    = [SPIDER_integrate( r2(:,1), [1,1], grid, corners, size_vec, pol );
             SPIDER_integrate( r2(:,2), [1,1], grid, corners, size_vec, pol );
             SPIDER_integrate( r2(:,3), [1,1], grid, corners, size_vec, pol );];
scales(a) = 1;
a = a+1;

labels{a} = "ddot{r}_3";
G(:,a)    = [SPIDER_integrate( r3(:,1), [1,1], grid, corners, size_vec, pol );
             SPIDER_integrate( r3(:,2), [1,1], grid, corners, size_vec, pol );
             SPIDER_integrate( r3(:,3), [1,1], grid, corners, size_vec, pol );];
scales(a) = 1;
a = a+1;

%try several different power laws
for p = 2:7
  F12 = r12./vecnorm(r12,2,2).^(p+1);
  F13 = r13./vecnorm(r13,2,2).^(p+1);
  F23 = r23./vecnorm(r23,2,2).^(p+1);

  labels{a} = "hat{r}_{12} / r_{12}^" + p;
  G(:,a)    = [SPIDER_integrate( F12(:,1), [], grid, corners, size_vec, pol );
               SPIDER_integrate( F12(:,2), [], grid, corners, size_vec, pol );
               SPIDER_integrate( F12(:,3), [], grid, corners, size_vec, pol );];
  scales(a) = 1;
  a = a+1;
 
  labels{a} = "hat{r}_{13} / r_{13}^" + p;
  G(:,a)    = [SPIDER_integrate( F13(:,1), [], grid, corners, size_vec, pol );
               SPIDER_integrate( F13(:,2), [], grid, corners, size_vec, pol );
               SPIDER_integrate( F13(:,3), [], grid, corners, size_vec, pol );];
  scales(a) = 1;
  a = a+1;
 
  labels{a} = "hat{r}_{23} / r_{23}^" + p;
  G(:,a)    = [SPIDER_integrate( F23(:,1), [], grid, corners, size_vec, pol );
               SPIDER_integrate( F23(:,2), [], grid, corners, size_vec, pol );
               SPIDER_integrate( F23(:,3), [], grid, corners, size_vec, pol );];
  scales(a) = 1;
  a = a+1;
end



%normalize the feature matrix.
norm_vec = SPIDER_integrate( 0*r1(:,1) + 1, [], grid, corners, size_vec, pol );      
norm_vec = repmat(norm_vec, [3,1]);
G        = G./norm_vec;
%Nondimensionalize
G = G./scales;









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 3: sparse regression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gamma = 10;

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

function vel = v( x )
  x = x.'; %take transpose

  r1 = x(1:3);
  r2 = x(4:6);

  p1 = x(7:9);
  p2 = x(10:12);

  r12 = r1-r2;
  r13 = 2*r1 + r2;
  r23 = r1 + 2*r2;
  
  F12 = r12 / norm(r12)^3;
  F13 = r13 / norm(r13)^3;
  F23 = r23 / norm(r23)^3;
  
  vel = [p1; p2; -F12 - F13; F12 - F23];
  vel = vel.'; %take transpose
end
