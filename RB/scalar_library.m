clear;
restoredefaultpath();
addpath("../SPIDER_functions/");


number_of_library_terms = 3;   %under-estimate this
number_of_windows       = 1024; %number of domains we integrate over 
degrees_of_freedom      = 1;   %vectors have one degree of freedom
dimension               = 2;   %how many dimensions does our data have?
envelope_power          = 8;   %weight is (1-x^2)^power
size_vec                = [round(N/3), round(M/8)]; %how many gridpoints should we use per integration?
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

x = grid{1};
t = grid{2};
dxs = [ 1, 1];

%Nondimensionalize and compute mean and variation of data
%TODO;

%%%%%%%%%%%%%%%%%%%%%%
% LIBRARY INTEGRATION
%%%%%%%%%%%%%%%%%%%%%%
a = 1;

labels{a} =  "dudt";
G(:,a)    = SPIDER_integrate( u, [2], grid, corners, size_vec, pol0 );
scales(a) = 1;
a = a+1;


%normalize with non-odd polynomial.
norm_vec = SPIDER_integrate( 0*u + 1, [], grid, corners, size_vec, pol0 );      
norm_vec = repmat( norm_vec, [dof*numel(u0s),1] );

G = G./norm_vec;
G = G./scales;


save("seed_" + seed + ".mat", "G", "labels", "scales");

end