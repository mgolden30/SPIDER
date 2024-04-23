clear;
restoredefaultpath();
addpath("../SPIDER_functions/");


load("data1e10.mat");
U = data(:,:,:,1);
W = data(:,:,:,2);
P = data(:,:,:,3);
T = data(:,:,:,4);


number_of_library_terms = 3;   %under-estimate this
number_of_windows       = 1024; %number of domains we integrate over 
degrees_of_freedom      = 1;   %vectors have one degree of freedom
dimension               = 3;   %how many dimensions does our data have?
envelope_power          = 8;   %weight is (1-x^2)^power
size_vec                = [32,32,32,32]; %how many gridpoints should we use per integration?
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

size_of_data = size(U, 1:dimension);

seed = 1;
corners = pick_subdomains_manual_seed( size_of_data, size_vec, buffer, nw, seed );

grid = {x,z,t};

%Nondimensionalize and compute mean and variation of data
%TODO;

%%%%%%%%%%%%%%%%%%%%%%
% LIBRARY INTEGRATION
%%%%%%%%%%%%%%%%%%%%%%
a = 1;

labels{a} =  "dWdt";
G(:,a)    = SPIDER_integrate( W, [3], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} =  "d(UW)dx";
G(:,a)    = SPIDER_integrate( U.*W, [1], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} =  "d(WW)dz";
G(:,a)    = SPIDER_integrate( W.*W, [2], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} =  "Px";
G(:,a)    = SPIDER_integrate( P, [1], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} =  "Pz";
G(:,a)    = SPIDER_integrate( P, [2], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} =  "W_xx";
G(:,a)    = SPIDER_integrate( W, [1,1], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} =  "W_zz";
G(:,a)    = SPIDER_integrate( W, [2,2], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} =  "T";
G(:,a)    = SPIDER_integrate( T, [], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} =  "T_x";
G(:,a)    = SPIDER_integrate( T, [1], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} =  "T_y";
G(:,a)    = SPIDER_integrate( T, [2], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;

labels{a} =  "1";
G(:,a)    = SPIDER_integrate( 0*T+1, [], grid, corners, size_vec, pol );
scales(a) = 1;
a = a+1;


%normalize with non-odd polynomial.
norm_vec = SPIDER_integrate( 0*U + 1, [], grid, corners, size_vec, pol );      
%norm_vec = repmat( norm_vec, [dof*numel(u0s),1] );

G = G./norm_vec;
G = G./scales;

