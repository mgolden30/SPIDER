function [res, cs] = GISR(G,c0)
  %{
  PURPOSE:
  Do sparse SVD stuff.
  %}

  %Boolean vector of terms which are turned on
  I = (c0 ~= 0);

  n0 = sum(I); %initial size of model
  N  = size(G,2); %library size

  res = zeros(N,1); %residuals of each sparse subspace
  cs  = zeros(N,N); %coefficients of sparse approximate null vectors

  %threshold is NOT uncertainty in singular values
  %rather it is the value of f(x) upon which root finding exits
  threshold = 1e-20;

  %Add terms
  for n = n0:N
    Gr = G( :, I );
    [~,S,V] = svd(Gr, 'econ');
    
    %read coefficients and residuals
    cs(I,n) = V(:,n);
    res(n)  = S(n,n);

    %Call C code to dermine modified singular values
    s_min = add_column( Gr, G(:,~I), threshold  );

    [~,j] = min(s_min);
    off = find( I == 0);
    I( off(j) )=1;
  end

  %Remove terms
  I = (c0 ~= 0);
  for n = n0:-1:1
    Gr = G( :, I );
    [~,S,V] = svd(Gr);
    %read coefficients and residuals
    cs(I,n) = V(:,n);
    res(n)  = S(n,n);

    if(n==1); break; end

    %Call the C code to compute modified singular values
    s_min = remove_column( Gr, threshold  );

    [~,j] = min( s_min );

    on = find(I);
    I(on(j)) = 0;
  end
end