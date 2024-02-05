N = 128;
x = (0:(N-1))/N * 2 * pi;

u0 = cos(3*x) - 0.5*sin(x);
L = 22;
M = 256 * 20;
T = 10 * 20;

maxit = 128;
threshold = 1e-9; %acceptable error in state-velocity

epsilon = 1e-6;

[u, iterations, flag] = KS_implicit_6th( u0, T, M, L, maxit, threshold, epsilon );


tiledlayout(1,2);

nexttile
imagesc(u);

nexttile
scatter( 1:M+1, iterations )

save("trajectory.mat");

function [u, iterations, flag] = KS_implicit_midpoint( u0, T, M, L, maxit, threshold )
  %{
  PURPOSE:
  To avoid small timesteps, use a fully implicit method to solve the KS equation

  Use Newton-Raphson-GMRES to converge the implicit step
  %}

  N = numel(u0);
  u = zeros(N,M+1);
  u(:,1) = u0;

  dt = T/M;
  
  iterations = zeros(M+1,1);

  %Make a fourier wavevector k
  k = 0:N-1; k(k>N/2) = k(k>N/2) - N;
  k = k*2*pi/L;
  k = k'; %make it a column vector

  deriv = @(u,n) real(ifft( (1i*k).^n .* fft(u) ));
  f     = @(u) - deriv(u.^2,1)/2 - deriv(u,2) - deriv(u,4); 
  df    = @(u,du) -u.*deriv(du,1) - du.*deriv(u,1) - deriv(du,2) - deriv(du,4);

  for i = 2:M+1
    % use a trivial implicit guess
    % I need k1 = f( (x + xf)/2) and xf = x + dt * k1. 
    % Substituting, I get the implicit equation k1 = f( x + dt*k1/2 ) 
    
    %Make a guess
    %u_mid = u(:,i-1) + dt * deriv(u(:,i-1).^2, 2)/2; %explicit advective step
    %u_mid = real(ifft( fft(u_mid)./(1 - dt/2*k.^2 + dt/2*k.^4) ));
    k1 = f( u(:,i-1) );
 
    size(k1)

    %anonymous functions for implicit condition and its Jacobian
    F  = @( k1)  k1 -   f( u(:,i-1) + dt*k1/2 );

    for it = 1:maxit
      dF = @(dk1) dk1 -  dt/2*df( u(:,i-1) + dt*k1/2, dk1 );
      
      F0 = F(k1);
      plot(k1);
      drawnow
      fprintf( "it %d: |F| = %.6e\n", it, norm(F0) );

      if( norm(F0) < threshold )
        iterations(i) = it;
        break;
      end
      if( it == maxit)
        %we didn't converge
        flag = false;
        return;
      end

      %otherwise we need to Newton step
      %h = 1e-5;
      %dF = @(v) (F(k1 + v) - F0)/h;
      
      inner = 64;
      outer = 2;
      tol = 1e-6;
      [step, ~] = gmres( dF, F0, inner, tol, outer );

      k1 = k1 - step;
    end
    u(:,i) = u(:,i-1) + k1*dt;
  end
  flag = true; %we converged!
end


function [u, iterations, flag] = KS_implicit_4th( u0, T, M, L, maxit, threshold )
  %{
  PURPOSE:
  To avoid small timesteps, use a fully implicit method to solve the KS
  equation. Fourth order GLRK

  Use Newton-Raphson-GMRES to converge the implicit step
  %}

  %Runge-Kutta coefficients
  a11 = 1/4;
  a12 = 1/4 - sqrt(3)/6;
  a21 = 1/4 + sqrt(3)/6;
  a22 = 1/4;

  N = numel(u0);
  u = zeros(N,M+1);
  u(:,1) = u0;

  dt = T/(M);
  
  %Make a fourier wavevector k
  k = 0:N-1; k(k>N/2) = k(k>N/2) - N;
  k = k*2*pi/L;
  k = k'; %make it a column vector

  deriv = @(u,n) real(ifft( (1i*k).^n .* fft(u) ));
  
  f     = @(u) - deriv(u.^2,1)/2 - deriv(u,2) - deriv(u,4); 
  df    = @(u,du) -u.*deriv(du,1) - du.*deriv(u,1) - deriv(du,2) - deriv(du,4);

  iterations= zeros(M+1,1);

  %keep track of converged vs
  vs = zeros(2*N,M+1);

  for i = 2:M+1
    % use a trivial implicit guess
    % I need k1 = f( (x + xf)/2) and xf = x + dt * k1. 
    % Substituting, I get the implicit equation k1 = f( x + dt*k1/2 ) 
    
    % Make a guess
    if( i < 6 )
      %Not enough data to extrapolate
      c1 = a11 + a12;
      c2 = a21 + a22;
    
      u1 = u(:,i-1) - c1*dt*deriv(u(:,i-1).^2, 1)/2; %explicit advection
      u1 = real(ifft( exp(c1*dt*(k.^2 - k.^4)) .* fft(u1) ));
      k1 = f(u1);
    
      u2 = u(:,i-1) - c2*dt*deriv(u(:,i-1).^2, 1)/2; %explicit advection
      u2 = real(ifft( exp(c2*dt*(k.^2 - k.^4)) .* fft(u2) ));
      k2 = f(u2);
      
      v = [k1;k2];
    else
      %We have past values of v to work with
      v = 2*vs(:,i-1) - vs(:,i-2);
     
      %Use a third order polynomial
      %v = 4*vs(:,i-1) - 6*vs(:,i-2) + 4*vs(i-3) - vs(:,i-4);
    end

    %define accessors
    K1 = @(v) v(      1:N );
    K2 = @(v) v( N + (1:N));

    %anonymous functions for implicit condition and its Jacobian
    F1     = @(k1,k2)  k1 - f( u(:,i-1) + dt*(a11*k1 + a12*k2)      );
    F2     = @(k1,k2)  k2 - f( u(:,i-1) + dt*(a21*k1 + a22*k2)      );
    F      = @(v) [ F1(K1(v), K2(v)); F2(K1(v), K2(v)) ];


    for it = 1:maxit
      F0 = F(v);

      k1 = K1(v);
      k2 = K2(v);

      %define the Jacobian since the base-point will be changing.
      dF1_d1 = @(dk1)   dk1 -  df( u(:,i-1) + dt*(a11*k1 + a12*k2), a11*dt*dk1 );
      dF1_d2 = @(dk2)       -  df( u(:,i-1) + dt*(a11*k1 + a12*k2), a12*dt*dk2 );
      dF2_d1 = @(dk1)       -  df( u(:,i-1) + dt*(a21*k1 + a22*k2), a21*dt*dk1 );
      dF2_d2 = @(dk2)   dk2 -  df( u(:,i-1) + dt*(a21*k1 + a22*k2), a22*dt*dk2 );
      
      dF = @(dv) [ dF1_d1(K1(dv)) + dF1_d2(K2(dv)); 
                   dF2_d1(K1(dv)) + dF2_d2(K2(dv)) ];

      %plot(v);
      %drawnow
      fprintf( "timestep %04d: it %03d: |F| = %.6e\n", i, it, norm(F0) );
      
      if( norm(F0) < threshold )
        iterations(i) = it;
        break;
      end

      if( it == maxit)
        %we didn't converge
        flag = false;
        return;
      end
      
      inner = 64;
      outer = 1;
      tol = 1e-6;
      [step, ~] = gmres( dF, F0, inner, tol, outer );

      v = v - step;
    end
    vs(:,i) = v;
    fprintf("\n\n");

    k1 = K1(v);
    k2 = K2(v);
    %Take the fourth order step
    u(:,i) = u(:,i-1) + dt*(k1 + k2)/2;
  end
  flag = true; %we converged!
end


function [u, iterations, flag] = KS_implicit_6th( u0, T, M, L, maxit, threshold, epsilon )
  %{
  PURPOSE:
  To avoid small timesteps, use a fully implicit method to solve the KS
  equation. Sixth order GLRK

  Use Newton-Raphson-GMRES to converge the implicit step
  %}

  %Runge-Kutta coefficients
  a11 = 5/36;
  a12 = 2/9  - sqrt(15)/15;
  a13 = 5/36 - sqrt(15)/30;

  a21 = 5/36 + sqrt(15)/24;
  a22 = 2/9;
  a23 = 5/36 - sqrt(15)/24;

  a31 = 5/36 + sqrt(15)/30;
  a32 = 2/9  + sqrt(15)/15;
  a33 = 5/36;

  N = numel(u0);
  u = zeros(N,M+1);
  u(:,1) = u0;

  dt = T/(M);
  
  %Make a fourier wavevector k
  k = 0:N-1; k(k>N/2) = k(k>N/2) - N;
  k = k*2*pi/L;
  k = k'; %make it a column vector

  deriv = @(u,n) real(ifft( (1i*k).^n .* fft(u) ));
  
  % define a modification term and its Jacobian
  coeffs  = [1,1,1,1];

  %add a strange nonlinearity
  modif = @(u)    -epsilon*( 3*coeffs(1)*u.^2 + 4*coeffs(2)*u.^3 + 5*coeffs(3)*u.^4 + 6*coeffs(4)*u.^5) .* deriv( u, 1 );
  dmodif= @(u,du) -epsilon*( 3*coeffs(1)*u.^2 + 4*coeffs(2)*u.^3 + 5*coeffs(3)*u.^4 + 6*coeffs(4)*u.^5) .* deriv( du, 1 ) -epsilon*du.*( 6*coeffs(1)*u + 12*coeffs(2)*u.^2 + 20*coeffs(3)*u.^3 + 30*coeffs(4)*u.^4) .* deriv( u, 1 );

  f     = @(u) - deriv(u.^2,1)/2 - deriv(u,2) - deriv(u,4) + modif(u); 
  df    = @(u,du) -u.*deriv(du,1) - du.*deriv(u,1) - deriv(du,2) - deriv(du,4) + dmodif(u,du);

  iterations= zeros(M+1,1);

  %keep track of converged vs
  vs = zeros(3*N,M+1);

  for i = 2:M+1
    % use a trivial implicit guess
    % I need k1 = f( (x + xf)/2) and xf = x + dt * k1. 
    % Substituting, I get the implicit equation k1 = f( x + dt*k1/2 ) 
    
    % Make a guess
    if( i < 6 )
      %Not enough data to extrapolate
      c1 = a11 + a12 + a13;
      c2 = a21 + a22 + a23;
      c3 = a31 + a32 + a33;
    
      u1 = u(:,i-1) - c1*dt*deriv(u(:,i-1).^2, 1)/2; %explicit advection
      u1 = real(ifft( exp(c1*dt*(k.^2 - k.^4)) .* fft(u1) ));
      k1 = f(u1);
    
      u2 = u(:,i-1) - c2*dt*deriv(u(:,i-1).^2, 1)/2; %explicit advection
      u2 = real(ifft( exp(c2*dt*(k.^2 - k.^4)) .* fft(u2) ));
      k2 = f(u2);

      u3 = u(:,i-1) - c3*dt*deriv(u(:,i-1).^2, 1)/2; %explicit advection
      u3 = real(ifft( exp(c3*dt*(k.^2 - k.^4)) .* fft(u3) ));
      k3 = f(u3);

      v = [k1;k2;k3];
    else
      %We have past values of v to work with
      v = 2*vs(:,i-1) - vs(:,i-2);
    end

    %define accessors
    K1 = @(v) v(       (1:N));
    K2 = @(v) v(   N + (1:N));
    K3 = @(v) v( 2*N + (1:N));

    %anonymous functions for implicit condition and its Jacobian
    F1     = @(k1,k2,k3)  k1 - f( u(:,i-1) + dt*(a11*k1 + a12*k2 + a13*k3)      );
    F2     = @(k1,k2,k3)  k2 - f( u(:,i-1) + dt*(a21*k1 + a22*k2 + a23*k3)      );
    F3     = @(k1,k2,k3)  k3 - f( u(:,i-1) + dt*(a31*k1 + a32*k2 + a33*k3)      );
    F      = @(v) [ F1(K1(v), K2(v), K3(v)); 
                    F2(K1(v), K2(v), K3(v));
                    F3(K1(v), K2(v), K3(v))];


    for it = 1:maxit
      F0 = F(v);

      k1 = K1(v);
      k2 = K2(v);
      k3 = K3(v);

      %define the Jacobian since the base-point will be changing.
      dF1_d1 = @(dk1)   dk1 -  df( u(:,i-1) + dt*(a11*k1 + a12*k2 + a13*k3), a11*dt*dk1 );
      dF1_d2 = @(dk2)       -  df( u(:,i-1) + dt*(a11*k1 + a12*k2 + a13*k3), a12*dt*dk2 );
      dF1_d3 = @(dk3)       -  df( u(:,i-1) + dt*(a11*k1 + a12*k2 + a13*k3), a13*dt*dk3 );
      
      dF2_d1 = @(dk1)       -  df( u(:,i-1) + dt*(a21*k1 + a22*k2 + a23*k3), a21*dt*dk1 );
      dF2_d2 = @(dk2)   dk2 -  df( u(:,i-1) + dt*(a21*k1 + a22*k2 + a23*k3), a22*dt*dk2 );
      dF2_d3 = @(dk3)       -  df( u(:,i-1) + dt*(a21*k1 + a22*k2 + a23*k3), a23*dt*dk3 );
      
      dF3_d1 = @(dk1)       -  df( u(:,i-1) + dt*(a31*k1 + a32*k2 + a33*k3), a31*dt*dk1 );
      dF3_d2 = @(dk2)       -  df( u(:,i-1) + dt*(a31*k1 + a32*k2 + a33*k3), a32*dt*dk2 );
      dF3_d3 = @(dk3)   dk3 -  df( u(:,i-1) + dt*(a31*k1 + a32*k2 + a33*k3), a33*dt*dk3 );

      dF = @(dv) [ dF1_d1(K1(dv)) + dF1_d2(K2(dv)) + dF1_d3(K3(dv)); 
                   dF2_d1(K1(dv)) + dF2_d2(K2(dv)) + dF2_d3(K3(dv));
                   dF3_d1(K1(dv)) + dF3_d2(K2(dv)) + dF3_d3(K3(dv));
                   ];

      %plot(v);
      %drawnow
      fprintf( "timestep %04d: it %03d: |F| = %.6e\n", i, it, norm(F0) );
      
      if( norm(F0) < threshold )
        iterations(i) = it;
        break;
      end

      if( it == maxit)
        %we didn't converge
        flag = false;
        return;
      end
      
      inner = 64;
      outer = 1;
      tol = 1e-6;
      [step, ~] = gmres( dF, F0, inner, tol, outer );

      v = v - step;
    end
    vs(:,i) = v;
    fprintf("\n\n");

    k1 = K1(v);
    k2 = K2(v);
    k3 = K3(v);

    %Take the sixth order step
    b1 = 5/18;
    b2 = 4/9;
    b3 = 5/18;
    u(:,i) = u(:,i-1) + dt*(b1*k1 + b2*k2 + b3*k3);
  end
  flag = true; %we converged!
end