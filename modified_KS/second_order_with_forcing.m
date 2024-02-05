function U = second_order_with_forcing( u0, M, N, L, T, nu, f )
  %Solve inviscid Burgers' equation

  U = zeros(N,M);
  U(:,1) = u0;

  I = eye(N);

  h  = L/N;
  D  = ( circshift(I,-1) - circshift(I,1) ) / (2*h);
  D2 = ( circshift(I,-1) - 2*I + circshift(I,1) )/h/h;

  dt = T/(M-1);
  for i = 2:M
    A = I + dt/2*D*diag(U(:,i-1)) - nu*dt/2*D2;
    U(:,i) = linsolve( A, (I + nu*dt/2*D2)*U(:,i-1) + dt*(f(:,i) + f(:,i-1))/2 );
  end
end