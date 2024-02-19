function [cs, residuals, fs] = greedy_regression_pure_matlab_reg( A, eps )
  %{
  PURPOSE:
  The minimum of L = |A*c| in nested sparse subspaces, such that L
  increases minimally at each stage.

  INPUT:
  A - a matrix to look for sparse null vectors of

  OUTPUT:
  cs - columns of this matrix are the increasingly sparse approximate null
       vectors.
  reisudals - vecnorm( A*cs );
  %}

  m = size(A,1);
  n = size(A,2);

  fbounds = zeros(n,2);

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
      alpha = 1/norm(a);
      w = alpha*U'*a;

      s = diag(S); %turn singular vectors into array
      
      bounds = [1, s(end-1)/s(end)];
      
      f = @(sigma)  1 - 1/alpha^2 * sum( w.^2 ./ ( s.^2 - sigma.^2 - eps) );
      
      maxit = 128;
      g = 0;

      fbounds(i,1) = f(bounds(1));
      fbounds(i,2) = f(bounds(2));
      for j = 1:maxit
        g  = sum(bounds)/2; %bisection guess
        fg = f(g);
        
        if fg < 0
          bounds(1) = g;
        else
          bounds(2) = g;
        end
      end
      candidates(i) = g;
    end

    [~, i_min] = min( candidates );

    fs(n) = f(candidates(i_min));

    j = find(I);
    I(j(i_min)) = 0;
    A = A0(:, I);
    size(A)
    n = n-1;
  end
end