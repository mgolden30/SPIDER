function [x] = TSGLRK2_iterative( x0, f, dt, N, threshold )
  %{
  PURPOSE:
  This is a time-symetric implicit integrator based on Gauss-Legendre
  quadrature

  INPUT:
  x0 - initial state
  f  - velocity
  j  - jacobian of velocity
  dt - timestep
  N - number of timesteps
  %}

  if( size(x0,2) ~= 1 )
    error('x0 needs to be a column vector');
  end

  sq3 = sqrt(3); %needed for quadrature
  n = numel(x0);
  x = zeros( n, N );
  x(:,1) = x0;
  for i = 2:N
    xl = x(:,i-1); %x last
    kl = f(xl); %f(x_last). Use this to generate an initial guess
    k1 = f(xl + sq3/6*dt*kl);
    k2 = f(xl + (0.5 + sq3/6)*dt*kl);
    xn = xl + dt*kl; %x_next
    
    %Make a giant vector
    v = [k1; k2; xn];

    it = 1;
    while( true )
      N = numel(x);
      
      %macros for 
      K1 = @(v) v(0*N + (1:N));
      K2 = @(v) v(1*N + (1:N));
      xf = @(v) v(2*N + (1:N));

      F = @(v) [ K1(v) - f( (xl+xf(v))/2 - sq3/6*dt*K2(v) ); K2(v) - f( (xl+xf(v))/2 + sq3/6*dt*K1(v) ); (xf(v) - xl)/dt- (K1(v) + K2(v))/2;];
      F0 = F(v);

      %Check for convergence
      if( norm(F) < threshold )
        [it norm(F)]
        break;
      end

      h = 1e-5;
      J = @(s) (F(v + h*s) - F0)/h;
      
      inner = 64;
      tol = 1e-6;
      outer = 1;
      step = gmres( J, F0, inner, tol, outer );

      v = v - 0.1 * step;
      xn = xf(v); it = it + 1;
    end
    x(:,i) = xn;
  end
end