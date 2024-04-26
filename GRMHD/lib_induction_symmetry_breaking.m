%{
Vector library for GRMHD
%}


%% Load data
%clear;

name = "sim1";

%cartesian MHD data
tic
[rho, v, B, P, t, x1, x2, x3] = load_MHD_data_stitch();
toc

%Integrate library
addpath("../SPIDER_functions/");

number_of_library_terms = 2;  %under-estimate this
number_of_windows       = 128;%number of domains we integrate over 
degrees_of_freedom      = 1;  %vectors have one degree of freedom
dimension               = 4;  %how many dimensions does our data have?
envelope_power          = 4;  %weight is (1-x^2)^power
size_vec                = [48, 48, 48, 32]; %how many gridpoints should we use per integration?
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
pol = multiply_pol( pol, legendre_pol( [1,0,0,0], dimension ) );

size_of_data = size(B, 1:dimension);
corners = pick_subdomains( size_of_data, size_vec, buffer, nw );

%Nondimensionalize and compute mean and variation of data
%TODO;

%Create grid variable
grid = {t,x1,x2,x3};

%%%%%%%%%%%%%%%%%%%%%%
% LIBRARY INTEGRATION
%%%%%%%%%%%%%%%%%%%%%%
a = 1;

%Specify the component you want
labels{a} = "\partial_t B_3";
G(:,a)    =  SPIDER_integrate( B(:,:,:,:,3), [1], grid, corners, size_vec, pol );  
scales(a) = 1;
a         = a+1;

for k = 1:3
for j = 2:4
  labels{a} = "\partial_" + j + "\partial_" + j + "B_" + k;
  G(:,a)    =  SPIDER_integrate( B(:,:,:,:,k), [j,j], grid, corners, size_vec, pol );  
  scales(a) = 1;
  a         = a+1;
end
labels{a} = "B_" + k;
G(:,a)    =  SPIDER_integrate( B(:,:,:,:,k), [], grid, corners, size_vec, pol );  
scales(a) = 1;
a         = a+1;
end

%v times B for magnetodynamics

source_x = v(:,:,:,:,2).*B(:,:,:,:,3) - v(:,:,:,:,3).*B(:,:,:,:,2);
source_y = v(:,:,:,:,3).*B(:,:,:,:,1) - v(:,:,:,:,1).*B(:,:,:,:,3);
source_z = v(:,:,:,:,1).*B(:,:,:,:,2) - v(:,:,:,:,2).*B(:,:,:,:,1);

labels{a} = "curl( v x B )_z";
G(:,a)    = SPIDER_integrate( source_y, [2], grid, corners, size_vec, pol ) - SPIDER_integrate( source_x, [3], grid, corners, size_vec, pol );
scales(a) = 1;
a         = a+1;

labels{a} = "curl( v x B )_z cubic correction 1";
G(:,a)    = SPIDER_integrate( source_y, [2,2,2], grid, corners, size_vec, pol );
scales(a) = 1;
a         = a+1;

labels{a} = "curl( v x B )_z cubic correction 2";
G(:,a)    = SPIDER_integrate( source_x, [3,3,3], grid, corners, size_vec, pol );
scales(a) = 1;
a         = a+1;

%normalize with non-odd polynomial.
norm_vec = SPIDER_integrate( 0*B(:,:,:,:,1) + 1, [], grid, corners, size_vec, pol0 );      
norm_vec = repmat( norm_vec, [dof,1] );


G = G./norm_vec;
G = G./scales;