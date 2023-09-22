function [U, dx_vec] = generate_KS_timeseries( u, T, L, M )
  %{
  This solves the KS equation with periodic boundary conditions using a
  time-symmetric scheme inspired by Crank-Nicolson (but not exactly
  Crank-Nicolson)
  
  INPUT:
  u - initial state
  T - integration time
  L - physical size of domain (chaos is L~40)
  M - number of timesteps
 
  OUTPUT:
  U - timeseries matrix. Contains the numerical trajectory.
  dx_vec - [dx, dt] for model discovery
  %}

  N = numel(u);     %number of gridpoints
  U = zeros(N,M+1); %Matrix which will contain the timeseries 
  h  = T/M;         %timestep for integration
  I  = eye(N);      %identity matrix
    
  x = (0:N-1)/N*2*pi; % spatial coordinates
  F = exp( 1i*(x').*(0:N-1) );
  k = 0:(N-1); k(k>N/2) = k(k>N/2) - N;
  k = k*2*pi/L; %Don't forget to rescale!!!
  
  D  =  real(1i*F*diag(k   )*F'/N);
  D2 = -real(   F*diag(k.^2)*F'/N);
  D4 =  real(   F*diag(k.^4)*F'/N);
      
  B  = I - h/2*(D2 + D4); %constant matrix used in timestepping
  A0 = I + h/2*(D2 + D4); %constant part of other matrix used in timestepping
  
  U(:,1) = u;
  for i = 1:M
    A  = A0 + h/2*D*diag(u);
    u = linsolve(A, B*u);

    %restart = 32;
    %outer   = 32;
    %tol     = 1e-9;
    %[u, ~] = gmres(A, B*u, restart, tol, outer);

    U(:,i+1) = u;
  end
  
  dx_vec = [L/N, h]; %[dx, dt]
end
