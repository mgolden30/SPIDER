%{
This script provides some methods for finding sparse relations from feature
matrices G
%}



N = 128;
M = 1024;

x = (0:(N-1))/N*(2*pi);

u0 = cos(3*x) + cos(x-1);
u0 = u0';

T = 100;
L = 22;
U = generate_KS_timeseries(u0, T, L ,M);


imagesc(U)

%Integrate library
addpath("../SPIDER_functions/");


number_of_library_terms = 3;   %under-estimate this
number_of_windows       = 1024; %number of domains we integrate over 
degrees_of_freedom      = 1;   %vectors have one degree of freedom
dimension               = 2;   %how many dimensions does our data have?
envelope_power          = 4;   %weight is (1-x^2)^power
size_vec                = [32,32]; %how many gridpoints should we use per integration?
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
%pol = multiply_pol( pol, legendre_pol( [0,0,1], dimension ) );

size_of_data = size(U, 1:dimension);
seed = 1;
corners = pick_subdomains_manual_seed( size_of_data, size_vec, buffer, nw, seed );

%Nondimensionalize and compute mean and variation of data
%TODO;

%Create grid variable
grid = {x*L/(2*pi),  T*(1:M)/M};

%%%%%%%%%%%%%%%%%%%%%%
% LIBRARY INTEGRATION
%%%%%%%%%%%%%%%%%%%%%%
a = 1;




labels{a} = "u";
G(:,a)    = SPIDER_integrate( U, [], grid, corners, size_vec, pol );
scales(a) = 1;
a         = a+1;

labels{a} = "u_t";
G(:,a)    = SPIDER_integrate( U, [2], grid, corners, size_vec, pol );
scales(a) = 1;
a         = a+1;

labels{a} = "(u^2)_x";
G(:,a)    = SPIDER_integrate( U.^2, [1], grid, corners, size_vec, pol );
scales(a) = 1;
a         = a+1;

labels{a} = "(u)_xx";
G(:,a)    = SPIDER_integrate( U, [1,1], grid, corners, size_vec, pol );
scales(a) = 1;
a         = a+1;

labels{a} = "(u)_xxxx";
G(:,a)    = SPIDER_integrate( U, [1,1,1,1], grid, corners, size_vec, pol );
scales(a) = 1;
a         = a+1;

%normalize with non-odd polynomial.
norm_vec = SPIDER_integrate( 0*U + 1, [], grid, corners, size_vec, pol0 );      
norm_vec = repmat( norm_vec, [dof,1] );


G = G./norm_vec;
G = G./scales;

