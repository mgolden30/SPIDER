function [c, residual] = inhomogeneous_regression( A, b )
  %QR trick to remove scaling with number of windows
  G = [b, A];
  [~, G] = qr(G, "econ");
  b = G(:,1);
  A = G(:,2:end);

  n = size(A,2);
  c = zeros(n,n);
  residual = zeros(n,1);

  size(A)
  size(b)
  c(:,n) = A\b;
  residual(n) = norm(A*c(:,n) - b);
  I = ones(n,1);

  for i = (n-1):-1:1
    res = [];
    for j = 1:sum(I)
      Ar = A(:,I == 1);
      Ar(:,j) = []; %try eliminating
      
      cr = Ar\b;
      res = [res, norm(Ar*cr - b) ];
    end
    [~, j] = min(res);

    idx = find(I);
    I(idx(j)) = 0;

    Ar = A(:,I == 1);
    cr = Ar\b;
    residual(i) = norm(Ar*cr - b);
    c(I==1,i) = cr;
  end

  residual = residual / norm(b);
end