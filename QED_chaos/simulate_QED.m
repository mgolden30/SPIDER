%{
This script solves the Dirac-Maxwell equations in 2+1D as specified in 
"Homogeneous Quantum Electrodynamic Turbulence" by John Shebalin

We use the Coulomb gauge (div(A) = 0) in numerics
%}

clear; %clear memory

N =  128; %number of points per spatial axis
M = 1024*4; %number of timesteps

save_every = 32; %for movie/ simulation output

A = zeros(N,N,2); % vector potential
E = zeros(N,N,2); % electric field
P = zeros(N,N,4); % bispinor electron
Q = zeros(N,N,4); % bispinor positron (I realize this is a bad name...)

%Timeseries matrices to save for SPIDER
As = zeros(N,N,M/save_every,2);
Es = zeros(N,N,M/save_every,2);
Ps = zeros(N,N,M/save_every,4);
Qs = zeros(N,N,M/save_every,4);



gp = (0:N-1)/N*2*pi;
[x,y] = ndgrid(gp,gp);


% Initial data! 
% Remember div(A)=0 is a constraint
% Also net charge must be zero for periodic boundary conditions to work
% |P|^2 - |Q|^2 = 0


P(:,:,1) = exp(1i*x);
P(:,:,2) = exp(1i*2*y);
P(:,:,3) = exp(1i*(x-0.5));
P(:,:,4) = exp(1i*y);

Q(:,:,1) = exp( 1i*(x+y)  );
Q(:,:,2) = exp( 1i*x );
Q(:,:,3) = exp( 1i*y );
Q(:,:,4) = exp( 1i*(2*x+y) );

normalization = 2;
P = P/normalization;
Q = Q/normalization;

A(:,:,1) = sin(y);
A(:,:,2) = sin(x);

E = coulomb_E(P,Q); %Solve Gauss's law for E

state = pack_state(A,E,P,Q);

h = 1e-3;

vidObj = VideoWriter('QED.avi');
open(vidObj);
fprintf("Starting...\n");
div_E = [];
div_A = [];
counter = 1;
for step = 1:M
  step
  [c1, c2] = check_constraints(state);
  
  div_A(step) = c2;
  div_E(step) = c1;
  
  state= rk4_step( state, @state_velocity, h);
  
  %Decide if it is time to save/plot
  if mod(step, save_every) ~= 0
    continue 
  end


  [A,E,P,Q] = unpack_state( state );
  As(:,:,counter,:) = A;
  Es(:,:,counter,:) = E;
  Ps(:,:,counter,:) = P;
  Qs(:,:,counter,:) = Q;
  counter = counter + 1;


  plot_state(state);
  drawnow
  
  % Write each frame to the file.
  currFrame = getframe(gcf);
  writeVideo(vidObj,currFrame);
end
% Close the file.
close(vidObj);

%% Plot constraints
clf
semilogy(div_E);
hold on
  semilogy(div_A);
hold off
xlim([0, M]);
xlabel("timestep");
ylabel("constraint value");
legend({"div(E)", "div(A)"});
saveas(gcf, "constraints.png");

%% Save data for SPIDER
grid = {gp, gp, (1:save_every:M)*h}; %{x,y,t} arrays
save('QED_simulation_data.mat', "As", "Es", "Ps", "Qs", "dx", "dt", "grid");







function state = pack_state( A,E,P,Q )
  dim   = 3;
  state = cat( dim, A, E, P, Q);
end

function [A,E,P,Q] = unpack_state( state )
  A = state(:,:,1:2);
  E = state(:,:,3:4);
  P = state(:,:,5:8);
  Q = state(:,:,9:12);
end

function dstate = state_velocity(state)
  [A,E,P,Q] = unpack_state( state );
  
  N = size(A,1);
  k = 0:N-1; k( k>N/2 ) = k( k>N/2 ) - N;
  
  %Step 1: solve for the electric potential via Poisson
  %P is electron, Q is positron
  rho = sum( conj(Q).*Q, 3 ) - sum( conj(P).*P, 3 );
  rho = squeeze(rho);
  
  %solve nabla^2 \phi = - \rho
  V_fft = fft2(rho)./(k.^2 + k'.^2);
  V_fft(1,1) = 0; %fix the mean voltage
  
  V   = real( ifft2(        V_fft    ) );
  V_x = real( ifft2( 1i*k'.*V_fft    ) );
  V_y = real( ifft2( 1i*    V_fft.*k ) );
  
  dA = zeros(N,N,2);
  dA(:,:,1) = -E(:,:,1) - V_x;
  dA(:,:,2) = -E(:,:,2) - V_y;
  
  %Step 2: Compute time derivative of E 
  lap_A = zeros(N,N,2); %need laplacian of A
  for i=1:2
    lap_A(:,:,i) = -real(ifft2( (k.^2 + k'.^2) .* fft2(squeeze(A(:,:,i))))); 
  end
  
  %Need current from spinors
  % P is electorn, Q is positron
  J = zeros(N,N,2);
  J(:,:,1) = conj( Q(:,:,1) ).*Q(:,:,4) + ...
             conj( Q(:,:,2) ).*Q(:,:,3) + ...
             conj( Q(:,:,3) ).*Q(:,:,2) + ...
             conj( Q(:,:,4) ).*Q(:,:,1) - ...
             conj( P(:,:,1) ).*P(:,:,4) - ...
             conj( P(:,:,2) ).*P(:,:,3) - ...
             conj( P(:,:,3) ).*P(:,:,2) - ...
             conj( P(:,:,4) ).*P(:,:,1);
  
         
  J(:,:,2) = -1i*conj( Q(:,:,1) ).*Q(:,:,4) + ...
              1i*conj( Q(:,:,2) ).*Q(:,:,3) + ...
             -1i*conj( Q(:,:,3) ).*Q(:,:,2) + ...
              1i*conj( Q(:,:,4) ).*Q(:,:,1) - ...
             -1i*conj( P(:,:,1) ).*P(:,:,4) - ...
              1i*conj( P(:,:,2) ).*P(:,:,3) - ...
             -1i*conj( P(:,:,3) ).*P(:,:,2) - ...
              1i*conj( P(:,:,4) ).*P(:,:,1);
          
  dE = - lap_A - real(J);
  
  
  %Step 3: Dirac time
  g0 = [1, 0, 0, 0;
        0, 1, 0, 0;
        0, 0,-1, 0;
        0, 0, 0,-1];
    
  a1 = [0,0,0,1;
        0,0,1,0;
        0,1,0,0;
        1,0,0,0];
  
  a2 = [0,   0,  0,-1i;
        0,   0, 1i,  0;
        0, -1i,  0,  0;
       1i,   0,  0,  0];
  
  P_x = ifft2( 1i*k'.*fft2(P)    );
  P_y = ifft2( 1i*    fft2(P).*k );
  
  Q_x = ifft2( 1i*k'.*fft2(Q)    );
  Q_y = ifft2( 1i*    fft2(Q).*k );
  
  P   = reshape(P,  [N*N,4]);
  P_x = reshape(P_x,[N*N,4]);
  P_y = reshape(P_y,[N*N,4]);
  
  Q   = reshape(Q,  [N*N,4]);
  Q_x = reshape(Q_x,[N*N,4]);
  Q_y = reshape(Q_y,[N*N,4]);

  V = reshape(V, [N*N,1]);
  A = reshape(A, [N*N,2]);
  
  % P is electron, Q is positron
  dP =  1i*V.*P + (-1i*A(:,1).*P - P_x)*a1.' + (-1i*A(:,2).*P - P_y)*a2.' - 1i*P*g0.';
  dQ = -1i*V.*Q + ( 1i*A(:,1).*Q - Q_x)*a1.' + ( 1i*A(:,2).*Q - Q_y)*a2.' - 1i*Q*g0.'; %change sign of charge coupling terms
  
  dP = reshape(dP, [N,N,4]);
  dQ = reshape(dQ, [N,N,4]);
  
  dstate = pack_state( dA, dE, dP, dQ );
end

function E = coulomb_E(P,Q)
  N = size(P,1);
  k = 0:N-1; k( k>N/2 ) = k( k>N/2 ) - N;
  
  %Step 1: solve for the electric potential via Poisson
  rho = sum( conj(Q).*Q, 3 ) - sum( conj(P).*P, 3 );
  rho = squeeze(rho);
  
  V_fft = fft2(rho)./(k.^2 + k'.^2);
  V_fft(1,1) = 0; %fix the mean voltage
  
  V_x = real( ifft2( 1i*k'.*V_fft    ) );
  V_y = real( ifft2( 1i*    V_fft.*k ) );
  
  dim = 3;
  E   = -cat(dim, V_x, V_y);
end

function xp = rk4_step( x, v, h )
  k1 = h*v(x);
  k2 = h*v(x+k1/2);
  k3 = h*v(x+k2/2);
  k4 = h*v(x+k3);
  
  xp = x + (k1 + 2*k2 + 2*k3 + k4)/6;
end

function plot_state(state)
  tiledlayout(2,2);
  
  [A,E,P,Q] = unpack_state(state);
  
  rho_e = sum( conj(P).*P, 3 );
  rho_p = sum( conj(Q).*Q, 3 );
  
  rho_e = squeeze(rho_e)';
  rho_p = squeeze(rho_p)';
  
  nexttile
  imagesc(rho_e);
  title("electron density")
  axis square
  colorbar(); 
  set(gca, 'ydir', 'normal' );
  caxis([0 1]);
  
  nexttile
  imagesc(rho_p);
  title("positron density");
  axis square
  colorbar(); 
  set(gca, 'ydir', 'normal' );
  caxis([0 1]);
  
  nexttile
  imagesc( rho_p - rho_e );
  axis square
  title('charge density');
  colorbar(); 
  set(gca, 'ydir', 'normal' );
  caxis([-0.5 0.5])
  
  nexttile
  imagesc( squeeze(sqrt(sum(E.^2, 3)))' )
  axis square
  title('|E|');
  colorbar(); 
  set(gca, 'ydir', 'normal' );
  caxis([0 2])
  drawnow
end

function [c1,c2] = check_constraints(state)
  %See if A is div free, and if E satisfies Gauss
  
  [A,E,P,Q] = unpack_state(state);
  
  rho_e = sum( conj(P).*P, 3 );
  rho_p = sum( conj(Q).*Q, 3 );
  
  rho_e = squeeze(rho_e);
  rho_p = squeeze(rho_p);
  
  N = size(A,1);
  
  k = 0:N-1; k( k>N/2 ) = k( k>N/2 ) - N;

  Gauss = real(ifft2( 1i*k'.*fft2(squeeze(E(:,:,1)))    )) + ...
          real(ifft2( 1i*    fft2(squeeze(E(:,:,2))).*k )) - rho_p + rho_e;
  
  Gauge = real(ifft2( 1i*k'.*fft2(squeeze(A(:,:,1)))    )) + ...
          real(ifft2( 1i*    fft2(squeeze(A(:,:,2))).*k ));
  
  c1 = max(max(abs(Gauss)));
  c2 = max(max(abs(Gauge)));
end