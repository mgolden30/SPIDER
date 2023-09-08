%{
Scalar library for GRMHD
%}


%% Load data
clear;

%cartesian MHD data
tic
[rho, v, B, P, t, x1, x2, x3] = load_MHD_data_stitch();
toc

%%%Compute correlations
%tau = correlation_length( rho )
%plot(-1./tau);
%return

%Integrate library
addpath("../SPIDER_functions/");

number_of_library_terms = 3;   %under-estimate this
number_of_windows       = 276*2; %number of domains we integrate over 
degrees_of_freedom      = 1;   %vectors have one degree of freedom
dimension               = 4;   %how many dimensions does our data have?
envelope_power          = 5;   %weight is (1-x^2)^power
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
%pol = multiply_pol( pol, legendre_pol( [0,0,1], dimension ) );

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

%macro for a string partial derivative
%subtract 1 to go from MATLAB indexing to what theorists usually take for
%spacetime x0= t, x1 = x, x2 = y, etc.
str  = @(j) "\partial_" + (j-1);

%{
for i = 1:3
  for j = 1:4 %since 2:4 are the spatial coordinates
    labels{a} = "\partial_" + j + " B_" + i;
    G(:,a)    = SPIDER_integrate_meshes( B(:,:,:,:,i), [j], grid, corners, size_vec, pol );
    scales(a) = 1;
    a         = a+1; 
  end
end
%}


fields = cell(4,1); %4-momentum
fields{1} = rho;
fields{2} = rho.*v(:,:,:,:,1);
fields{3} = rho.*v(:,:,:,:,2);
fields{4} = rho.*v(:,:,:,:,3);

names = { "\rho", "(\rho u_1)", "(\rho u_2)", "(\rho u_3)" };

for f = 1:4 %loop over fields

%first derivatives of rho
for j = 1:4 %since 1:4 are the spatial coordinates
  labels{a} = "\partial_" + (j-1) + names(f);
  G(:,a)    = SPIDER_integrate( fields{f}, [j], grid, corners, size_vec, pol );
  scales(a) = 1;
  a         = a+1;
end

%second derivatives of rho
for j = 1:4 %since 1:4 are the spatial coordinates
  for k = j:4
    labels{a} = str(j) + str(k) + names(f);
    G (:,a)    = SPIDER_integrate( fields{f}, [j,k], grid, corners, size_vec, pol );
    scales(a) = 1;
    a         = a+1;
  end
end

%third derivatives of rho
for j = 1:4 %since 1:4 are the spatial coordinates
for k = j:4
for l = k:4
    labels{a} = str(j) + str(k) + str(l) + names(f);
    G (:,a)    = SPIDER_integrate(fields{f} , [j,k,l], grid, corners, size_vec, pol );
    scales(a) = 1;
    a         = a+1;
end
end
end

%fourth derivatives of rho
for j = 1:4 %since 1:4 are the spatial coordinates
for k = j:4
for l = k:4
for m = l:4
    labels{a} = str(j) + str(k) + str(l) + str(m) + names(f);
    G (:,a)    = SPIDER_integrate( fields{f}, [j,k,l,m], grid, corners, size_vec, pol );
    scales(a) = 1;
    a         = a+1;
end
end
end
end

end %end loop over fields


%normalize with non-odd polynomial.
norm_vec = SPIDER_integrate( 0*B(:,:,:,:,1) + 1, [], grid, corners, size_vec, pol0 );      
norm_vec = repmat( norm_vec, [dof,1] );


G = G./norm_vec;
G = G./scales;