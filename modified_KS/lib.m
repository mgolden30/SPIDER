
%% Load data
clear;
restoredefaultpath();

addpath("library_generation/");

%grid variables
%{
N = 128;
M = 128;

% physical parameters
L = 22;
T = 20;
nu = 0.1;

x = (0:(N-1))/N*2*pi;
x = x';

%define many initial conditions to train on
u0s = { 
        cos(x);
        exp( -(x-pi).^4 );
        sin(2*x) + exp(cos(x));
        sin(x) - 1./(1 + 0.2*cos(x));
        cos(x)./(1 + 0.2*sin(x));
        cos(2*x);
        sin(3*x-1);
        cos(x-1)./(2+cos(x));
        cos(sin(x));
        cos(2*sin(x));
      };

t = (0:(M-1))*T/M;
[xx,tt] = meshgrid( x, t );
force = cos(2*xx).*cos(tt*2*pi);
%}

load('trajectory.mat');

seeds = [1, 543212345];

for seed = seeds

G_total = [];

for state_index = 1:numel(u0s)

u0 = u0s{state_index};
%[u] = analytic_solution( u0, nu, L, T, M );
force = 0*force;

addpath("../../KS/")
[u, dx_vec] = generate_KS_timeseries( u0, T, L, M-1 );
%[u] = second_order_with_forcing( u0, M, N, L, T, nu, force );

imagesc((1:M)/M*T, x/(2*pi)*L, u); xlabel('t'); ylabel('x'); 
colormap jet
set(gca, 'ydir', 'normal');
colorbar();
clim([-1 1]*5);
title("u(x,t)")
drawnow;
saveas(gcf, "training_"+state_index+".png");

grid = { x/(2*pi)*L, linspace(0,T,M) };


%%%Compute correlations
%tau = correlation_length( rho )
%plot(-1./tau);
%return

%Integrate library
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

dt = T/(M-1);
u_t = (circshift(u,-1,2) - u)/dt; %finite difference in time

nd = 10;%num derivs
du = cell(nd,1);
for l = 1:nd
  k = 0:N-1;
  k(k>N/2) = k(k>N/2) - N;
  k = k';
  k = 2*pi/L*k;
  du{l} = real(ifft( (1i * k).^l .* fft(u) ));
end

labels{a} =  "dudt";
G(:,a)    = SPIDER_integrate( u, [2], grid, corners, size_vec, pol0 );
scales(a) = 1;
a = a+1;

addpath("library_generation/");
%All other library terms are computed with the auto library generation
for i = 1:3^(nd)
  [canonical, derivs, v] = check_library_term(i);
  if(canonical)
    %Add to library
    [term, str] = generate_library_term(u,du,derivs);
    %{
    str
    v
    i
    derivs
    %}
    G(:,a)    = SPIDER_integrate( term, [], grid, corners, size_vec, pol0 );
    scales(a) = 1;
    labels{a} = str;
    a = a+1;
  end
end

G_total = [G_total; G];

end

G = G_total;


%normalize with non-odd polynomial.
norm_vec = SPIDER_integrate( 0*u + 1, [], grid, corners, size_vec, pol0 );      
norm_vec = repmat( norm_vec, [dof*numel(u0s),1] );

G = G./norm_vec;

%scales = vecnorm(G);
G = G./scales;

%{
addpath("../../GISR-main/GISR/");
[cs, residuals] = greedy_regression_pure_matlab( G );


semilogy(1:numel(residuals), residuals, 'o');
return
%}

save("seed_" + seed + ".mat", "G", "labels", "scales");

end