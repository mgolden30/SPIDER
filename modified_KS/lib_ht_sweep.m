
%% Load data
clear;

%grid variables
N = 256;
M = 512;

% physical parameters
L = 2*pi;
T = 4; %doubled from 2 when M=256
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

us = {};

for i = 1:numel(u0s)
  us{i} = second_order( u0s{i}, M, N, L, T, nu );
end

seeds = [1, 543212345];

for ht = 32:8:256
ht

for seed = seeds

G_total = [];

for state_index = 1:numel(u0s)

%u0 = u0s{state_index};
%[u] = analytic_solution( u0, nu, L, T, M );
%[u] = second_order( u0, M, N, L, T, nu );

u = us{state_index};

%{
imagesc((1:M)/M*T, x/(2*pi)*L, u); xlabel('t'); ylabel('x'); 
colormap jet
set(gca, 'ydir', 'normal');
colorbar();
clim([-1 1]);
title("u(x,t)")
pause(0.1);
%}

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
size_vec                = [64, ht]; %how many gridpoints should we use per integration?
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

save("Ht_sweep/seed2_" + seed + "_ht_" + ht + ".mat", "G", "labels", "scales");

end
end