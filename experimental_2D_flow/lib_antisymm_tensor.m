%{
Generate a feature matrix G for vectors using simulated QED data
%}


%% Load data
if ~exist("U", "var")
  fprintf("Loading in data...\n");
  load("UVexp_Re36.34_dim_1.mat"); %experimental data

  %Create a vector object U
  U = U_t;
  V = V_t;

  U = pagetranspose(U);
  V = pagetranspose(V);

  U( ~isfinite(U) ) = 0;
  V( ~isfinite(V) ) = 0;
end


%Integrate library
addpath("../SPIDER_functions/");

number_of_library_terms = 3;  %under-estimate this
number_of_windows       = 128;%number of domains we integrate over 
degrees_of_freedom      = 1;  %vectors have one degree of freedom
dimension               = 3;  %how many dimensions does our data have?
envelope_power          = 5;  %weight is (1-x^2)^power
size_vec                = [ 64, 64, 64 ]; %how many gridpoints should we use per integration?


buffer                  = 3; %Don't use points this close to boundary

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
pol = multiply_pol( pol, legendre_pol( [0,0,1], dimension ) );

size_of_data = size(U, 1:dimension);
corners = pick_subdomains( size_of_data, size_vec, buffer, nw );


%Nondimensionalize!
mean_U = mean( sqrt(U.^2 + V.^2), 'all' ); %mean flow velocity
mean_W = mean( abs(dxV_t - dyU_t), 'all' ); %mean vorticity
U = U./mean_U;
V = V./mean_U;
dx = dx/mean_U*mean_W;
dy = dy/mean_U*mean_W;
dt = dt*mean_W;



%Define characteristic scales
mean_U  = mean( abs(U), 'all' );
std_U   = std( U, 0, 'all' );

%Create grid variable
x = (1:size(U,1))*dx;
y = (1:size(U,2))*dx;
t = (1:size(U,3))*dt;
grid = {x,y,t};

%length and time scales
L = x(end) - x(1); %Use size of domain
T = t(end) - t(1); %Use size of domain

a = 1; %running index over library space

labels{a} = "u_i";
G(:,a)    = SPIDER_integrate_curl( U, V, [], grid, corners, size_vec, pol );
scales(a) = mean_U;
a         = a+1;

labels{a} = "u^2 u_i";
G(:,a)    = SPIDER_integrate_curl( (U.^2 + V.^2).*U, (U.^2 + V.^2).*V, [], grid, corners, size_vec, pol );
scales(a) = mean_U^3;
a         = a+1;

labels{a} = "\partial_t u_i";
G(:,a)    = SPIDER_integrate_curl( U, V, [3], grid, corners, size_vec, pol );
scales(a) = std_U/T;
a         = a+1;

labels{a} = "\partial_t^2 u_i";
G(:,a)    = SPIDER_integrate_curl( U, V, [3,3], grid, corners, size_vec, pol );
scales(a) = std_U/T/T;
a         = a+1;

labels{a} = "\nabla_j(u_i u_j)";
G(:,a)    = SPIDER_integrate_curl( U.*U, U.*V, [1], grid, corners, size_vec, pol ) ...
          + SPIDER_integrate_curl( U.*V, V.*V, [2], grid, corners, size_vec, pol );
scales(a) = mean_U*std_U/L;
a         = a+1;

labels{a} = "\nabla^2 u_i";
G(:,a)    = SPIDER_integrate_curl( U, V, [1,1], grid, corners, size_vec, pol ) + ...
          + SPIDER_integrate_curl( U, V, [2,2], grid, corners, size_vec, pol );
scales(a) = std_U/L/L;
a         = a+1;



%normalize with non-odd polynomial.
norm_vec = SPIDER_integrate( 0*U(:,:,:,1) + 1, [], grid, corners, size_vec, pol0 );      
norm_vec = repmat( norm_vec, [dof,1] );

scales = scales/L; %since we took another spatial derivative of eevery term

G = G./norm_vec;
G = G./scales;



function vals = SPIDER_integrate_curl( U, V, derivs, grid, corners, size_vec, pol )
  % U is x component of vector, V is y component
  vals = SPIDER_integrate( V, [derivs,1], grid, corners, size_vec, pol ) ...
       - SPIDER_integrate( U, [derivs,2], grid, corners, size_vec, pol );
end