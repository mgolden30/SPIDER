%{
Generate a feature matrix G for scalars using simulated QED data
%}

%% Load data
if ~exist("Es", "var")
  fprintf("Loading in data...\n");
  load("QED_simulation_data.mat");
end

% Integrate library
addpath("../SPIDER_functions/");

number_of_library_terms = 5;  %under-estimate this
number_of_windows       = 128;%number of domains we integrate over 
degrees_of_freedom      = 1;  %scalars have one degree of freedom
dimension               = 3;  %how many dimensions does our data have?
envelope_power          = 4;  %weight is (1-x^2)^power
size_vec                = [ 32, 32, 32]; %how many gridpoints should we use per integration?
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

size_of_data = size(As, 1:dimension);
corners = pick_subdomains( size_of_data, size_vec, buffer, nw );

%Compute secondary fields we might want (charge density and current)

rho = sum( abs(Qs).^2 - abs(Ps).^2, 4 );

J = 0*Es; %J will be current density
J(:,:,:,1) = conj( Qs(:,:,:,1) ).*Qs(:,:,:,4) + ...
             conj( Qs(:,:,:,2) ).*Qs(:,:,:,3) + ...
             conj( Qs(:,:,:,3) ).*Qs(:,:,:,2) + ...
             conj( Qs(:,:,:,4) ).*Qs(:,:,:,1) - ...
             conj( Ps(:,:,:,1) ).*Ps(:,:,:,4) - ...
             conj( Ps(:,:,:,2) ).*Ps(:,:,:,3) - ...
             conj( Ps(:,:,:,3) ).*Ps(:,:,:,2) - ...
             conj( Ps(:,:,:,4) ).*Ps(:,:,:,1);
 
J(:,:,:,2) = -1i*conj( Qs(:,:,:,1) ).*Qs(:,:,:,4) + ...
              1i*conj( Qs(:,:,:,2) ).*Qs(:,:,:,3) + ...
             -1i*conj( Qs(:,:,:,3) ).*Qs(:,:,:,2) + ...
              1i*conj( Qs(:,:,:,4) ).*Qs(:,:,:,1) - ...
             -1i*conj( Ps(:,:,:,1) ).*Ps(:,:,:,4) - ...
              1i*conj( Ps(:,:,:,2) ).*Ps(:,:,:,3) - ...
             -1i*conj( Ps(:,:,:,3) ).*Ps(:,:,:,2) - ...
              1i*conj( Ps(:,:,:,4) ).*Ps(:,:,:,1);
J = real(J);

phi = get_scalar_potential(Ps,Qs);

%Define characteristic scales
mean_E  = mean( abs(Es), 'all' );
mean_A  = mean( abs(As), 'all' );
mean_phi= mean(abs(phi), 'all' );
mean_rho= mean(abs(rho), 'all' );
mean_J  = mean( abs(J),  'all' );

std_E  = std( Es, 0, 'all' );
std_A  = std( As, 0, 'all' );
std_rho= std(rho, 0, 'all' );
std_phi= std(phi, 0, 'all' );
std_J  = std(  J, 0, 'all' );

%length and time scales
x = grid{1};
t = grid{3};
L = x(end) - x(1); %Use size of domain
T = t(end) - t(1); %Use size of domain




a = 1; %running index over library space


labels{a} = "\nabla_i E_i";
G(:,a)    = SPIDER_integrate( Es(:,:,:,1), [1], grid, corners, size_vec, pol ) + ...
          + SPIDER_integrate( Es(:,:,:,2), [2], grid, corners, size_vec, pol );
scales(a) = std_E/L;
a         = a+1;

%{
labels{a} = "\nabla_i A_i";
G(:,a)    = SPIDER_integrate( As(:,:,:,1), [1], grid, corners, size_vec, pol ) + ...
          + SPIDER_integrate( As(:,:,:,2), [2], grid, corners, size_vec, pol );
scales(a) = std_A/L;
a         = a+1;
%}

labels{a} = "E^2";
G(:,a)    = SPIDER_integrate( sum(Es.^2,4), [], grid, corners, size_vec, pol );
scales(a) = mean_E^2;
a         = a+1;

labels{a} = "A^2";
G(:,a)    = SPIDER_integrate( sum(As.^2,4), [], grid, corners, size_vec, pol );
scales(a) = mean_E^2;
a         = a+1;

labels{a} = "\rho";
G(:,a)    = SPIDER_integrate( rho, [], grid, corners, size_vec, pol );
scales(a) = mean_rho;
a         = a+1;

labels{a} = "\partial_t \rho";
G(:,a)    = SPIDER_integrate( rho, [3], grid, corners, size_vec, pol );
scales(a) = std_rho/T;
a         = a+1;


labels{a} = "\nabla_i j_i";
G(:,a)    = SPIDER_integrate( J(:,:,:,1), [1], grid, corners, size_vec, pol ) + ...
          + SPIDER_integrate( J(:,:,:,2), [2], grid, corners, size_vec, pol );
scales(a) = std_J/L;
a         = a+1;

labels{a} = "\nabla^2 \phi";
G(:,a)    = SPIDER_integrate( phi, [1,1], grid, corners, size_vec, pol ) + ...
          + SPIDER_integrate( phi, [2,2], grid, corners, size_vec, pol );
scales(a) = std_J/L;
a         = a+1;

norm_vec = SPIDER_integrate( 0*rho + 1, [], grid, corners, size_vec, pol );      
norm_vec = repmat( norm_vec, [dof,1] );

G = G./norm_vec;
G = G./scales;



function phi = get_scalar_potential(Ps,Qs)
  %P is electron, Q is positron
  rho = sum( conj(Qs).*Qs, 4 ) - sum( conj(Ps).*Ps, 4 );
  rho = squeeze(rho);
  
  phi = 0*rho;

  N = size(phi,1);
  k = 0:N-1; k( k>N/2 ) = k( k>N/2 ) - N;

  for t = 1:size(phi,3)
    %solve nabla^2 \phi = - \rho
    V_fft = fft2( squeeze(rho(:,:,t)) )./(k.^2 + k'.^2);
    V_fft(1,1) = 0; %fix the mean voltage
  
    phi(:,:,t) = real( ifft2( V_fft ) );

    %V_x = real( ifft2( 1i*k'.*V_fft    ) );
    %V_y = real( ifft2( 1i*    V_fft.*k ) );
  end
end