function p = estimate_momentum( rho )
  
  m = size(rho,1);
  n = size(rho,2);

  kx = 0:n-1;
  kx(kx>n/2) = kx(kx>n/2) - n;

  kt = 0:m-1;
  kt(kt>m/2) = kt(kt>m/2) - m;
  kt = kt.';

  rho = fft2(rho);
  rho = -kt.*rho./kx;
  rho(~isfinite(rho)) = 0;

  p = real(ifft2(rho));
end