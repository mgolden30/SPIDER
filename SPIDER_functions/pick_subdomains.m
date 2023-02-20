function corners = pick_subdomains( size_of_data, size_vec, buffer, nw )
  %For reproducibility
  seed = 1;
  rng(seed);

  dim = numel(size_of_data); %dimension of space we integrate in

  corners = zeros( dim, nw );
  for d = 1:dim
    corners(d,:) = randi( size_of_data(d) - size_vec(d) - 2*buffer, [1,nw] ) + buffer;
  end
end