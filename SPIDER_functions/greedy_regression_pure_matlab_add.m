function [cs, residuals] = greedy_regression_pure_matlab_add( G, c0, n_max )
  %{
  PURPOSE:
  The minimum of L = |G*c|_2 / |c|_2 in nested sparse subspaces, such that L
  decreases maximally at each stage.

  INPUT:
  G - a matrix to look for sparse null vectors of

  OUTPUT:
  cs - columns of this matrix are the increasingly sparse approximate null
       vectors.
  resiudals - vecnorm( G*cs );
  %}

  m = size(G,1);
  n = size(G,2);
  if m<n
    fprintf("error: matrix is underdetermined.");
    return;
  end

  %first rotate G so it is square and upper triangular
  [~, G] = qr(G);  
  G = G(1:n, 1:n);
  
  cs = zeros(n,n);
  I  = (c0~=0); %logical vector indicating sparsity
  residuals = zeros(n,1);

  while( sum(I) < n_max )
    ns = sum(I);

    [U, S, V] = svd( G(:,I), 'econ'); %take SVD of current submatrix 
    cs(I,ns) = V(:,ns);   %save out the smallest singular vector
    residuals(ns) = S(ns,ns);

    candidates = zeros(n,1);
    %Loop over terms we haven't used
    for i = 1:n
      if I(i) == 1
        candidates(i) = inf;
        continue;
      end

      %otherwise, we should try adding this term to the relation
      a = G(:,i);
      alpha = 1/norm(a);
      w = alpha*U'*a;

      s = diag(S); %turn singular vectors into array
      bounds = [0, s(end)]; 
      %must be lower than current singular value
      %must be greater than or equal to zero

      tau = sqrt( 1 - sum(w.^2) );
      f0  = @(sigma)  1 + 1/alpha^2 * sum( w.^2 ./ ( s.^2 - sigma.^2 ) ) - tau^2/alpha^2/sigma^2;
      reg = @(sigma)  sigma^2 * ( s(end)^2 - sigma.^2 ) * alpha^2;
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

        if (fg < 0  )
          bounds(1) = g;
        else
          bounds(2) = g;
        end

        %fprintf("iteration %d: g = %e, f(g) = %e\n", j, g, fg);
      end
      candidates(i) = g;
    end

    [~, i_min] = min( candidates );
    I(i_min) = 1; %include this term
  end

  %rescale the residual
  residuals = residuals / sqrt(m);
end