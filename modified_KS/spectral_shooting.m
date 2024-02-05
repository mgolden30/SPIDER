function [u] = spectral_shooting( u0, nu, L, T, M )
  %Use spectral in space, first order in time
  
  N = numel(u0);

  k = 0:(N-1);
  k(k>N/2) = k(k>N/2)-N;
  k = k';
  k0= k;
  k = (2*pi)*k/L;


  dt = T/(M-1);
  for i = 1:M
    u(:,i) = u0;

    rhs = real(ifft( k.^2 .* fft(u0.^2) ))/2;

    u0 = u0 + dt*rhs;
    
  t = (0:(M-1))/(M-1) * T;
  x0= (0:N-1)/N*2*pi;
  x0 = x0';

  %turn into spacetime matrix

  phase = exp( 1i*k0.*x0 );
  decay = exp( - nu*k.^2.*t );

  size(decay);
  size(phase);
  size(phi2);

  %decay = 0*decay + 1;

  phi2 = phi2 .* decay ;

  phi = real(ifft( phi2 ));

  %imagesc(phi); pause(5);

  u = -2*nu*real(ifft( 1i*k.*fft(phi) ))./phi;
end