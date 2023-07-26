%{
Generate a feature matrix G for vectors using simulated QED data
%}


%% Load data
clear;

%load in vorticity, flow velocity, ....
load("vorticity.mat")

N = size(o, 1);
dx = 2*pi/N;

omega = o;
U = 0*omega;
V = U;

omega_x = U;
omega_y = U;

for t = 1:size(omega,3)
  k = 0:N-1;
  k(k>N/2) = k(k>N/2) - N;

  k_sq = k.^2 + k'.^2;
  k_sq(1,1) = 1;

  o = fft2( omega(:,:,t) );

  V(:,:,t) = -real(ifft2( o.*1i.*k./k_sq ));
  U(:,:,t) =  real(ifft2( k'.*o.*1i./k_sq ));

  omega_x(:,:,t) =  real(ifft2( o.*1i.*k ));
  omega_y(:,:,t) =  real(ifft2( k'.*o.*1i ));
  
  %{
  imagesc( omega(:,:,t) );
  hold on
    quiver( V(:,:,t), U(:,:,t));
  hold off
  drawnow;
  %}
end


%Integrate library
addpath("../SPIDER_functions/");

number_of_library_terms = 3;   %under-estimate this
number_of_windows       = 128; %number of domains we integrate over 
degrees_of_freedom      = 1;   %vectors have one degree of freedom
dimension               = 3;   %how many dimensions does our data have?
envelope_power          = 4;   %weight is (1-x^2)^power
size_vec                = [ 32, 32, 32 ]; %how many gridpoints should we use per integration?

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
%pol = multiply_pol( pol, legendre_pol( [0,0,1], dimension ) );

size_of_data = size(U, 1:dimension);
corners = pick_subdomains( size_of_data, size_vec, buffer, nw );

% Nondimensionalize!
% TODO

% Create grid variable
x = (1:size(U,1))*dx;
y = (1:size(U,2))*dx;
t = (1:size(U,3))*dt;
grid = {x,y,t};

%length and time scales
L = x(end) - x(1); %Use size of domain
T = t(end) - t(1); %Use size of domain








%% END OF BOILERPLATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO TOUCH PAST HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = 1; %running index over library space

labels{a} = "w";
G(:,a)    = SPIDER_integrate( omega, [], grid, corners, size_vec, pol );
scales(a) = 1;
a         = a+1;

labels{a} = "\partial_t w";
G(:,a)    = SPIDER_integrate( omega, [3], grid, corners, size_vec, pol );
scales(a) = 1;
a         = a+1;

labels{a} = "u \partial_x w";
G(:,a)    = SPIDER_integrate( U.*omega_x, [], grid, corners, size_vec, pol );
scales(a) = 1;
a         = a+1;

labels{a} = "v \partial_x w";
G(:,a)    = SPIDER_integrate( V.*omega_x, [], grid, corners, size_vec, pol );
scales(a) = 1;
a         = a+1;

labels{a} = "u \partial_y w";
G(:,a)    = SPIDER_integrate( U.*omega_y, [], grid, corners, size_vec, pol );
scales(a) = 1;
a         = a+1;

labels{a} = "v \partial_y w";
G(:,a)    = SPIDER_integrate( V.*omega_y, [], grid, corners, size_vec, pol );
scales(a) = 1;
a         = a+1;

labels{a} = "\nabla^2 w";
G(:,a)    = SPIDER_integrate( omega, [1,1], grid, corners, size_vec, pol ) + ...
            SPIDER_integrate( omega, [2,2], grid, corners, size_vec, pol );
scales(a) = 1;
a         = a+1;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DON't TOUCH PAST HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% normalize
norm_vec = SPIDER_integrate( 0*U(:,:,:,1) + 1, [], grid, corners, size_vec, pol );      
norm_vec = repmat( norm_vec, [dof,1] );

G = G./norm_vec;
G = G./scales;