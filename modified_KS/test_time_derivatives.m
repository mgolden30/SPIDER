%{
Generate a library for Kuramoto-Sivashinsky 
%}


%% Load data
clear;
restoredefaultpath();

addpath("library_generation/");

load('trajectory.mat');

%%

imagesc(u)

tiledlayout(2,2);
nexttile %tile 1 is interpolated data at a spatial slice

%I need time derivatives of this data. Let's pick a spatial slice and look
%at the time series
t = linspace(0,T,M+1);
data = u(1,:);
plot(t, data);

%This data is nonperiodic, so I cannot just take a FFT to compute
%derivatives. Let's try a Chebyshev transform. I'll use this function
%provided here: https://www.mathworks.com/matlabcentral/fileexchange/44034-fast-chebyshev-differentiation

addpath("fchd\");

num_nodes = 512; %number of Chebyshev-Gauss nodes to project onto (actually uses one more than this)
nodes = cos(pi*(0:num_nodes)/num_nodes);

%shift nodes onto [0,T]
nodes = (nodes+1)/2*T;

%interpolate onto the Chebyshev grid

data2 = interpn( t,data,nodes );

hold on
plot( nodes, data2, 'o' );
hold off
title("u(t)")

%Take the chebyshev derivative
d_data2 = fchd( data2 ) * 2 / T; %Don;t forget to rescale

%Take the Finite Difference derivative
dt = t(2) - t(1);
d_data = (circshift( data,-1 ) - circshift(data,1))/(2*dt);
%Don't trust end points
d_data(1)   = inf;
d_data(end) = inf;

nexttile
plot( t, d_data );
hold on
  plot(nodes, d_data2, 'o');
hold off
title("first derivative");


dd_data = (circshift( data,-1 ) -2*data + circshift(data,1))/(dt*dt);
%Don't trust end points
dd_data(1)   = inf;
dd_data(end) = inf;
dd_data2 = fchd( d_data2 ) * 2 / T; %Don't forget to rescale

nexttile
plot( t, dd_data );
hold on
  plot(nodes, dd_data2, 'o');
hold off
title("second derivative");

nexttile
ddd_data = (-1/2*circshift(data,-2) + circshift( data,-1 ) - circshift(data,1) + 1/2*circshift(data,2))/(dt*dt*dt);
%Don't trust end points
ddd_data(1:2)   = inf;
ddd_data(end-1:end) = inf;
ddd_data2 = fchd( dd_data2 ) * 2 / T; %Don't forget to rescale
plot( t, -ddd_data );
hold on
  plot(nodes, ddd_data2, 'o');
hold off
ylim([-1 1]);
title("third derivative");


%% Since this seems to work, let's generate mixed partial derivative of u up to tenth order

nodes = fliplr(nodes);


max_d = 10;
deriv_matrix = cell( max_d +1, max_d + 1 );

[X,  TT ] = ndgrid( x, t );
[X2, TT2] = ndgrid( x, nodes );


for m = 1:max_d+1 %derivative in time
  u2 = interpn( X, TT, u, X2, TT2 ); %interpolate onto Chebyshev points in time 
  
  %compute the (m-1)th derivative in time
  for i = 1:(m-1)
    for j = 1:size(u2,1)
      u2(j,:) = fchd(u2(j,:)) *2 / T;
    end
  end
  
  %interpolate back to uniform grid
  u2= interpn( X2, TT2, u2, X, TT );

  for n = 1:max_d+1
    k = 0:N-1;
    k(k>N/2) = k(k>N/2) - N;
    k = k*2*pi/L;
    k = k';

    %Take the spatial derivative with a Fourier transform
    deriv_matrix{m,n} = real(ifft( (1i*k).^(n-1).*fft(u2) ));
  end
end

save("data/derivative_matrix.mat", "deriv_matrix", "num_nodes", "max_d");
