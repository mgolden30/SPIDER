function [cs, residuals] = exhaustive_search( A )
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
    [~, S, V] = svd(A, 'econ'); 
    cs(I,n) = V(:,n);   %save out the smallest singular vector
    residuals(n) = S(n,n);

    if( n == 1 )
      break;
    end

    candidates = zeros(n,1);
    for i = 1:n
      a = A(:,i);
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