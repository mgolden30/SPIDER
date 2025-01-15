function [x,y,t,u,v,p] = navier_stokes( n, dt, timesteps, nu, skip )
  grid1d = (0:(n-1))/n*2*pi;
  [x,y] = meshgrid( grid1d );

  omega = 3*sin(3*x).*sin(3*y) + sin(x).*sin(y);

  u = zeros(n,n,timesteps);
  v = zeros(n,n,timesteps);
  p = zeros(n,n,timesteps);
  
  k = 0:n-1;
  k(k>n/2) = k(k>n/2) - n;

  %dealiasing mask
  mask = abs(k) <= n/3;
  mask = mask & mask.';

  kx = k;
  ky = k.';

  to_u =  1i*ky ./(kx.^2 + ky.^2);
  to_v = -1i*kx./(kx.^2 + ky.^2);
  to_p = 1./(kx.^2 + ky.^2);

  to_p(1,1) = 0;
  to_u(1,1) = 0;
  to_v(1,1) = 0;
  
  to_u = to_u.*mask;
  to_v = to_v.*mask;
  to_p = to_p.*mask;
  


  kx = kx .*mask;
  ky = ky .*mask;

  omega = fft2(omega);

  e = exp(-dt/2 * nu * (kx.^2 + ky.^2) );

  for t = 1:timesteps
    u(:,:,t) = real(ifft2( to_u.*omega ));
    v(:,:,t) = real(ifft2( to_v.*omega ));

    ux = real(ifft2(1i*kx.*to_u.*omega));
    uy = real(ifft2(1i*ky.*to_u.*omega));
    vx = real(ifft2(1i*kx.*to_v.*omega));
    vy = real(ifft2(1i*ky.*to_v.*omega));
    
    p(:,:,t) = real(ifft2( to_p.*fft2( ux.^2 + vy.^2 + 2*uy.*vx ) ));
  
  
    for i = 1:skip
    k1 = dt*advection(omega, kx, ky, to_u, to_v, mask);
    omega = e.*omega;
    k1 = e.*k1;
    k2 = dt*advection(omega + k1/2, kx, ky, to_u, to_v, mask);
    k3 = dt*advection(omega + k2/2, kx, ky, to_u, to_v, mask);
    omega = e.*omega;
    k1 = e.*k1;
    k2 = e.*k2;
    k3 = e.*k3;
    k4 = dt*advection(omega + k3, kx, ky, to_u, to_v, mask);
    
    omega = omega + (k1 + 2*k2 + 2*k3 + k4)/6;
    end

    imagesc( real(ifft2(omega)) );
    axis square;
    colorbar();
    drawnow;
  end

  x = grid1d;
  y = grid1d;
  t = (0:timesteps-1)*dt*skip;
end



function advec = advection(omega, kx, ky, to_u, to_v, mask)
  u = real(ifft2( to_u.*omega ));
  v = real(ifft2( to_v.*omega ));
  
  wx = real(ifft2( 1i*kx.*omega ));
  wy = real(ifft2( 1i*ky.*omega ));

  advec = mask .* fft2(-u.*wx - v.*wy);
  
end