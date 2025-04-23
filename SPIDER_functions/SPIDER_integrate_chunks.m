function vals = SPIDER_integrate_safe_fast( data, derivs, grid, betas)
  %{
  PURPOSE:
  This function is meant to integrate space-time fields with nice
  polynomial weight functions to be used for model discovery. Here are some
  key points.

  1. Integration by parts is handled ANALYTICALLY
  2. Non-uniform grid spacings are allowed. Integration is approximated
     with the trapezoid rule.
  3. The "data" argument should only have dimensions corresponding to space
     and time. If you have a "velocity" matrix of size [3,128,128,...]
     where the first index is the components, you are going to have to feed
     components in one by one. This is a pain, but you can write a wrapper
     that uses this function as a building block.
  
  
  INPUT:
  data    - the data you need weighted integrals of - pre-subsampled, (data_size, n_windows)

  derivs  - [1,1] for instance uses integration by parts to integrate two
            derivatives along the first dimension (could be x or y
            depending on your convention)

  grid    - a cell containing [---points in x---;
             ---points in y---;
             ...] - this is already subdamples, size - (data_size, n_windows)

  
  
  pol - a [l, ndim] matrix describing the polynomial weight function

 
  OUTPUT:
  vals - a vector of integrated values.
  %}
  s = size(data);
  n           = size(s',1) - 1;        %dimensionality of the data
  num_windows = s(end);        %number of spacetime windows we want to integrate over

  %vals        = zeros(num_windows,1);   %vals are the integrated values
  
  %data = squeeze(data);        % only matters if data is one dimensional

  %{
  %Do a sanity check and make sure the number of coordinates given for the
  %spacetime cubes matches the dimensionality of the data
  data_size = size(data);
  if numel(data_size) == 2
    %If this is the case, then your data is stored as a simple matrix.
    %If data is 1D, it still is stored as a matrix since MATLAB is built on
    %matrices. Check for one of these dimensions being trivial.
    data_size( data_size == 1 ) = []; %kill trivial dimension
  end
  assert( numel(grid)      == n );
  %}
  
  % subsample_indices is a cell that will be used to subsample the data for 
  % spacetime domains.


  % Do integration by parts to figure out the new polynomial weight.
    %{
  poly_sign = 1;
  l = size(pol,1);  %upper bound of degree of the polynomial weight
  derivs0 = derivs; %save a copy since we will destroy the derivs object
  while( numel(derivs) ~= 0 )
    d = derivs(1); %direction to differentiate
      
    poly_sign = -poly_sign; %each integration by parts accumulates a -1
    for ll = 2:l
      pol( ll-1, d ) = (ll-1) * pol( ll, d );
    end
    pol(ll, d) = 0; %Fixing a bug! Don't forget to kill the last term manually!
    
    derivs(1) = []; %Delete this derivative as we have taken it analytically
  end
    %}

  %Now we are ready to start the integration loop
  

     %Construct literal weight "functions"

  Fs = cell(n,1);
  
  for i = 1:n
      Fs{i} = weight_function( betas(i), sum(derivs==i) );
  end

    
    conversion_factor = 1; 
    %^this will account for coordinate changes and the sign change from
    %integration by parts.
    
    


    for j = 1:n  % loop is over spacetime dimensions
      x_spec = grid{j};


      % Let x be the grid passed in (assumed to have units), and let y be 
      % Some subset of x translated and rescaled to the interval [-1,1]
      % Our weight functions are naturally functions of y instead of x
      % y = mx+b
      m = 2./(x_spec(end,:) - x_spec(1,:));
      b = -(x_spec(end,:) + x_spec(1,:))./(x_spec(end,:) - x_spec(1,:));
      
      y = m.*x_spec + b; %[-1,1] coordinates

      conversion_factor = conversion_factor .* m .^ sum(derivs == j);
      %^this accounts for changing derivatives from d/dx = m d/dy
      

      %poly_weight_1d = polyval( fliplr( pol(:,j)' ), y ); %evaluate the polynomial weight for this dimension

      desired_shape = ones(n-j+2,1);
      desired_shape(1) = size(x_spec,1);
      desired_shape(end) = size(x_spec,2);

      %phi = weight_function( y, betas(j), sum(derivs==j) );

      phi = Fs{j}(y);
      
      
      phi = reshape(phi, desired_shape'); %enforce that phi is a column vector


      %size(poly_weight_1d), desired_shape, size(data)

      %poly_weight_1d = reshape( poly_weight_1d, desired_shape'); %make sure its a column vector
      %poly_weight_1d([1,end]) = 0;
      
      %compute grid spacings to left and right
      hr      = circshift(x_spec, -1) - x_spec;
      hr(end) = 0;
      hl      = x_spec - circshift(x_spec, 1);
      hl(1)   = 0;

      %size(hl), size(hr), "mark"
      hl      = reshape( hl, desired_shape' );
      hr      = reshape( hr, desired_shape' );
      
      %hl      = reshape( hl, [size(y,1), 1] );
      %hr      = reshape( hr, [size(y,1), 1] );
      %size(hl), size(hr),size(data),'sizes'

      %"mark", size(data), size(poly_weight_1d), desired_shape
      %size(data), size(phi), size(y)


      data = phi.* data; %add polynomial weight
      
      data = (hr+hl)/2     .* data; %add trapezoid rule weights
      data = sum(data);             %sum along first dimension
      data = kill_first_dimension( data ); %kill first dimension
      
      %size(data)
    end
    %{toc,"integration"}
    %5{size(data),"data"}
    
    
    %size(conversion_factor)
    vals = data .* conversion_factor';
    
    %{size(vals),"vals"}

  end



function data = kill_first_dimension( data )
  %Squeeze almost always does what you want, until you have something of
  %size [1,10]. Squeeze does nothing in this case, so we can just take the
  %transpose.
  data = squeeze(data);
  if(size(data,1) == 1)
    data = data';
  end
end


function F = weight_function( beta, derivs )
  syms x;

  phi = (1-x^2)^beta;
  
  for i = 1:derivs
    phi = diff(phi,x);
  end

  phi = phi * (-1)^derivs;
 
  F = matlabFunction(phi, 'vars', {[x]});
end