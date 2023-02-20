function pol = legendre_pol( l, dimension )
  
  pols = cell(dimension,1);
  for i = 1:dimension
    syms x;
    pols{i} = flip( coeffs( legendreP(l(i), x), 'all' ) ); 
  end
  
  lens = zeros(dimension,1);
  for i = 1:dimension
    lens(i) = numel( pols{i} ); 
  end
  
  pol = zeros( max(lens), dimension );
  for i=1:dimension
    pol(1:numel(pols{i}), i) = pols{i}; 
  end
  
  pol = double(pol);
end