function [cs, residuals] = greedy_regression_pure_matlab_naive( A )
  %{
  PURPOSE:
  The minimum of L = |A*c| in nested sparse subspaces, such that L
  increases minimally at each stage.
  %}

  m = size(A,1);
  n = size(A,2);

  cs = 0;

  if m<n
    fprintf("error: matrix is underdetermined.");
    return;
  end

  %first rotate A so it is square and upper triangular
  [~, A] = qr(A);  A = A(1:n, 1:n);
  
  %keep a copy
  A0 = A; 

  cs = zeros(n,n);
  I  = ones(n,1); I = (I == 1); %logical vector indicating sparsity
  residuals = zeros(n,1);

  while( n > 0 )
    [U, S, V] = svd(A, 'econ'); 
    cs(I,n) = V(:,n);   %save out the smallest singular vector
    residuals(n) = S(n,n);

    if( n == 1 )
      break;
    end

    candidates = zeros(n,1);
    for i = 1:n
      a = A(:,i);
      %{
      alpha = 1/norm(a);
      w = alpha*U'*a;

      s = diag(S); %turn singular vectors into array
      bounds = [s(end), s(end-1)];
      
      f0  = @(sigma)  1 - 1/alpha^2 * sum( w.^2 ./ ( s.^2 - sigma.^2 ) );
      reg = @(sigma)  ( s(end)^2 - sigma.^2 ) .* ( s(end-1)^2 - sigma.^2 ) * alpha^2 / (s(end)^2 - s(end-1)^2);
      f   = @(sigma)  f0(sigma) .* reg(sigma);

      maxit = 128;
      threshold = 1e-130;
      g = 0;
      for j = 1:maxit
        g  = sum(bounds)/2; %bisection guess
        fg = f(g);
        if(abs(fg) < threshold)
          break;
        end

        if fg > 0
          bounds(1) = g;
        else
          bounds(2) = g;
        end
      end
      candidates(i) = g;
      %}
      A2 = A;
      A2(:,i) = [];
      s2 = svd(A2);
      candidates(i) = min(s2);
    end

    [~, i_min] = min( candidates );

    j = find(I);
    I(j(i_min)) = 0;
    A = A0(:, I);
    size(A)
    n = n-1;
  end

  %rescale the residual
  residuals = residuals / sqrt(m);
end