%{
Vector library for GRMHD
%}


%% Load data
%clear;


%cartesian MHD data
tic
[rho, v, B, P, t, x1, x2, x3] = load_MHD_data_stitch();
toc



%Integrate library
addpath("../SPIDER_functions/");

number_of_library_terms = 2;  %under-estimate this
number_of_windows       = 128;%number of domains we integrate over 
degrees_of_freedom      = 3;  %vectors have one degree of freedom
dimension               = 4;  %how many dimensions does our data have?
envelope_power          = 5;  %weight is (1-x^2)^power
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
pol  = multiply_pol( pol, legendre_pol( [1,0,0,0], dimension ) );

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

%{
labels{a} = "\partial_t B_i";
G(:,a)    = [SPIDER_integrate( B(:,:,:,:,1), [1], grid, corners, size_vec, pol );
             SPIDER_integrate( B(:,:,:,:,2), [1], grid, corners, size_vec, pol );
             SPIDER_integrate( B(:,:,:,:,3), [1], grid, corners, size_vec, pol );];  
scales(a) = 1;
a         = a+1;
%}

for j = 2:4
  labels{a} = "\partial_" + j + "\partial_" + j + "B_i";
  G(:,a)    = [SPIDER_integrate( B(:,:,:,:,1), [j,j], grid, corners, size_vec, pol );
               SPIDER_integrate( B(:,:,:,:,2), [j,j], grid, corners, size_vec, pol );
               SPIDER_integrate( B(:,:,:,:,3), [j,j], grid, corners, size_vec, pol );]; 
  scales(a) = 1;
  a         = a+1;
end


for j = 2:4
  labels{a} = "\partial_" + j + "\partial_" + j + "(rho u_i)";
  G(:,a)    = [SPIDER_integrate( rho.*v(:,:,:,:,1), [j,j], grid, corners, size_vec, pol );
               SPIDER_integrate( rho.*v(:,:,:,:,2), [j,j], grid, corners, size_vec, pol );
               SPIDER_integrate( rho.*v(:,:,:,:,3), [j,j], grid, corners, size_vec, pol );]; 
  scales(a) = 1;
  a         = a+1;
end

%{
for j = 2:4
  for jj = j:4
    labels{a} = "\partial_" + j + "\partial_" + j + "\partial_" + jj + "\partial_" + jj + "(rho u_i)";
    G(:,a)    = [SPIDER_integrate( rho.*v(:,:,:,:,1), [j,j,jj,jj], grid, corners, size_vec, pol );
                 SPIDER_integrate( rho.*v(:,:,:,:,2), [j,j,jj,jj], grid, corners, size_vec, pol );
                 SPIDER_integrate( rho.*v(:,:,:,:,3), [j,j,jj,jj], grid, corners, size_vec, pol );]; 
    scales(a) = 1;
    a         = a+1;
  end
end
%}


labels{a} = "B_i";
G(:,a)    = [SPIDER_integrate( B(:,:,:,:,1), [], grid, corners, size_vec, pol );
             SPIDER_integrate( B(:,:,:,:,2), [], grid, corners, size_vec, pol );
             SPIDER_integrate( B(:,:,:,:,3), [], grid, corners, size_vec, pol );];  
scales(a) = 1;
a         = a+1;

%v times B for magnetodynamics
source_x = v(:,:,:,:,2).*B(:,:,:,:,3) - v(:,:,:,:,3).*B(:,:,:,:,2);
source_y = v(:,:,:,:,3).*B(:,:,:,:,1) - v(:,:,:,:,1).*B(:,:,:,:,3);
source_z = v(:,:,:,:,1).*B(:,:,:,:,2) - v(:,:,:,:,2).*B(:,:,:,:,1);


labels{a} = "curl( v x B )";
G(:,a)    = [SPIDER_integrate( source_z, [3], grid, corners, size_vec, pol ) - SPIDER_integrate( source_y, [4], grid, corners, size_vec, pol );
             SPIDER_integrate( source_x, [4], grid, corners, size_vec, pol ) - SPIDER_integrate( source_z, [2], grid, corners, size_vec, pol );
             SPIDER_integrate( source_y, [2], grid, corners, size_vec, pol ) - SPIDER_integrate( source_x, [3], grid, corners, size_vec, pol );];  
scales(a) = 1;
a         = a+1;



labels{a} = "\partial_t( rho v_i)";
G(:,a)    = [SPIDER_integrate( rho.*v(:,:,:,:,1), [1], grid, corners, size_vec, pol );
             SPIDER_integrate( rho.*v(:,:,:,:,2), [1], grid, corners, size_vec, pol );
             SPIDER_integrate( rho.*v(:,:,:,:,3), [1], grid, corners, size_vec, pol );];  
scales(a) = 1;
a         = a+1;

v1 = v(:,:,:,:,1);
v2 = v(:,:,:,:,2);
v3 = v(:,:,:,:,3);

labels{a} = "\nabla_j( \rho v_j v_i)";
G(:,a)    = [SPIDER_integrate( rho.*v1.*v1, [2], grid, corners, size_vec, pol )+SPIDER_integrate( rho.*v1.*v2, [3], grid, corners, size_vec, pol )+SPIDER_integrate( rho.*v1.*v3, [4], grid, corners, size_vec, pol );
             SPIDER_integrate( rho.*v2.*v1, [2], grid, corners, size_vec, pol )+SPIDER_integrate( rho.*v2.*v2, [3], grid, corners, size_vec, pol )+SPIDER_integrate( rho.*v2.*v3, [4], grid, corners, size_vec, pol );
             SPIDER_integrate( rho.*v3.*v1, [2], grid, corners, size_vec, pol )+SPIDER_integrate( rho.*v3.*v2, [3], grid, corners, size_vec, pol )+SPIDER_integrate( rho.*v3.*v3, [4], grid, corners, size_vec, pol );]; 
scales(a) = 1;
a         = a+1;

labels{a} = "\nabla_i P";
G(:,a)    = [SPIDER_integrate( P, [2], grid, corners, size_vec, pol );
             SPIDER_integrate( P, [3], grid, corners, size_vec, pol );
             SPIDER_integrate( P, [4], grid, corners, size_vec, pol );];  
scales(a) = 1;
a         = a+1;

labels{a} = "\nabla_i B^2";
B_sq = sum( B.^2, 5 );
G(:,a)    = [SPIDER_integrate( B_sq, [2], grid, corners, size_vec, pol );
             SPIDER_integrate( B_sq, [3], grid, corners, size_vec, pol );
             SPIDER_integrate( B_sq, [4], grid, corners, size_vec, pol );];  
scales(a) = 1;
a         = a+1;

v1 = B(:,:,:,:,1);
v2 = B(:,:,:,:,2);
v3 = B(:,:,:,:,3);
labels{a} = "\nabla_j( B_j B_i)";
G(:,a)    = [SPIDER_integrate( v1.*v1, [2], grid, corners, size_vec, pol )+SPIDER_integrate( v1.*v2, [3], grid, corners, size_vec, pol )+SPIDER_integrate( v1.*v3, [4], grid, corners, size_vec, pol );
             SPIDER_integrate( v2.*v1, [2], grid, corners, size_vec, pol )+SPIDER_integrate( v2.*v2, [3], grid, corners, size_vec, pol )+SPIDER_integrate( v2.*v3, [4], grid, corners, size_vec, pol );
             SPIDER_integrate( v3.*v1, [2], grid, corners, size_vec, pol )+SPIDER_integrate( v3.*v2, [3], grid, corners, size_vec, pol )+SPIDER_integrate( v3.*v3, [4], grid, corners, size_vec, pol );];
scales(a) = 1;
a         = a+1;


v1 = v(:,:,:,:,1);
v2 = v(:,:,:,:,2);
v3 = v(:,:,:,:,3);

labels{a} = "\hat{z} \times \rho \vec{v}";
G(:,a)    = [SPIDER_integrate(  rho.*v2,  [], grid, corners, size_vec, pol );
             SPIDER_integrate( -rho.*v1,  [], grid, corners, size_vec, pol );
             SPIDER_integrate( 0*v1, [], grid, corners, size_vec, pol );];  
scales(a) = 1;
a         = a+1;


labels{a} = "\hat{y} \times \rho \vec{v}";
G(:,a)    = [SPIDER_integrate(  rho.*v3,  [], grid, corners, size_vec, pol );
             SPIDER_integrate( 0*v1,  [], grid, corners, size_vec, pol );
             SPIDER_integrate( -rho.*v1, [], grid, corners, size_vec, pol );];  
scales(a) = 1;
a         = a+1;


labels{a} = "\hat{x} \times \rho \vec{v}";
G(:,a)    = [SPIDER_integrate(   0.*v2,  [], grid, corners, size_vec, pol );
             SPIDER_integrate(  rho.*v3,  [], grid, corners, size_vec, pol );
             SPIDER_integrate( -rho.*v2, [], grid, corners, size_vec, pol );];  
scales(a) = 1;
a         = a+1;

[T,X,Y,Z] = ndgrid( t, x1, x2, x3 );

labels{a} = "rho x \hat{x}";
G(:,a)    = [SPIDER_integrate(  rho.*X,  [], grid, corners, size_vec, pol );
             SPIDER_integrate( 0*X,  [], grid, corners, size_vec, pol );
             SPIDER_integrate( 0*X, [], grid, corners, size_vec, pol );];  
scales(a) = 1;
a         = a+1;


labels{a} = "rho y \hat{y}";
G(:,a)    = [SPIDER_integrate( 0*X,  [], grid, corners, size_vec, pol );
             SPIDER_integrate( rho.*Y,  [], grid, corners, size_vec, pol );
             SPIDER_integrate( rho.*0,  [], grid, corners, size_vec, pol );];  
scales(a) = 1;
a         = a+1;

labels{a} = "rho z \hat{z}";
G(:,a)    = [SPIDER_integrate( 0*X,  [], grid, corners, size_vec, pol );
             SPIDER_integrate( 0*X,  [], grid, corners, size_vec, pol );
             SPIDER_integrate( rho.*Z,  [], grid, corners, size_vec, pol );];  
scales(a) = 1;
a         = a+1;

%normalize with non-odd polynomial.
norm_vec = SPIDER_integrate( 0*B(:,:,:,:,1) + 1, [], grid, corners, size_vec, pol0 );      
norm_vec = repmat( norm_vec, [dof,1] );

G = G./norm_vec;
G = G./scales;