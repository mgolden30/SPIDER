
%% Load data
clear;


Hts = [0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4]; %dimensionless
t_per_subdomain = 0.8; %hold this constant
x_per_subdomain = 1; %hold this constant

for sweep_idx = 1:numel( Hts )
  Ht = Hts(sweep_idx)


%grid variables
N = 256;
M = 2048;

% physical parameters
L = 2*pi;
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

t = (0:(M-1))*T/(M-1);
[xx,tt] = meshgrid( t, x );

sp_forcing = 5./(1+cos(xx).^2);
sp_forcing = sp_forcing - mean(sp_forcing, 2);

force = sp_forcing.*(cos(tt*2*pi) + cos(tt*4*pi));

seeds = [1, 543212345];

for seed = seeds

G_total = [];

for state_index = 1:numel(u0s)

u0 = u0s{state_index};
%[u] = analytic_solution( u0, nu, L, T, M );
[u] = second_order_with_forcing( u0, M, N, L, T, nu, force );

imagesc((1:M)/M*T, x/(2*pi)*L, u); xlabel('t'); ylabel('x'); 
%return
colormap jet
set(gca, 'ydir', 'normal');
colorbar();
clim([-1 1]);
title("u(x,t)")
pause(0.1);

grid = { x/(2*pi)*L, linspace(0,T,M) };


%%%Compute correlations
%tau = correlation_length( rho )
%plot(-1./tau);
%return

%Integrate library
addpath("../SPIDER_functions/");


number_of_library_terms = 3;   %under-estimate this
number_of_windows       = 64; %number of domains we integrate over 
degrees_of_freedom      = 1;   %vectors have one degree of freedom
dimension               = 2;   %how many dimensions does our data have?
envelope_power          = 8;   %weight is (1-x^2)^power
size_vec                = [round(x_per_subdomain/2/pi*N), round(t_per_subdomain/T*M*Ht) ]; %how many gridpoints should we use per integration?
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
u_t = (circshift(u,-1,2) - u)/dt;

names = {"u", "u^2"};
fields = {u, u.^2};


labels{a} =  "f";
G(:,a)    = SPIDER_integrate( force, [], grid, corners, size_vec, pol );
scales(a) = 1;
a         = a+1;

for f = 1:numel(fields)

labels{a} =  names{f};
G(:,a)    = SPIDER_integrate( fields{f}, [], grid, corners, size_vec, pol );
scales(a) = 1;
a         = a+1;

for i = 1:2
labels{a} = "\partial_" + i + names{f};
G(:,a)    = SPIDER_integrate( fields{f}, [i], grid, corners, size_vec, pol );
scales(a) = dxs(i);
a         = a+1;
end


for i = 1:2
for j = i:2
labels{a} = "\partial_" + i + "\partial_" + j + names{f};
G(:,a)    = SPIDER_integrate( fields{f}, [i,j], grid, corners, size_vec, pol );
scales(a) = dxs(i)*dxs(j);
a         = a+1;
end
end


for i = 1:2
for j = i:2
for k = j:2
labels{a} = "\partial_" + i + "\partial_" + j + "\partial_" + k + names{f};
G(:,a)    = SPIDER_integrate( fields{f}, [i,j,k], grid, corners, size_vec, pol );
scales(a) = dxs(i)*dxs(j)*dxs(k);
a         = a+1;
end
end
end


for i = 1:2
for j = i:2
for k = j:2
for l = k:2
labels{a} = "partial_" + i + "\partial_" + j + "\partial_" + k + "\partial_" + l + names{f};
G(:,a)    = SPIDER_integrate( fields{f}, [i,j,k,l], grid, corners, size_vec, pol );
scales(a) = dxs(i)*dxs(j)*dxs(k)*dxs(l);
a         = a+1;
end
end
end
end


for i = 1:2
for j = i:2
for k = j:2
for l = k:2
for m = l:2
labels{a} = "partial_" + i + "\partial_" + j + "\partial_" + k + "\partial_" + l + "\partial_" + m + names{f};
G(:,a)    = SPIDER_integrate( fields{f}, [i,j,k,l,m], grid, corners, size_vec, pol );
scales(a) = dxs(i)*dxs(j)*dxs(k)*dxs(l)*dxs(m);
a         = a+1;
end
end
end
end
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

save("seed_" + seed + "_Ht_" + Ht + ".mat", "G", "labels", "scales");
end
end