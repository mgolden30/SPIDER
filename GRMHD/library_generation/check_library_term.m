function [valid, fields, derivs, digits] = check_library_term( i, num_fields, dim )
  %{
  PURPOSE:
  I need words made up of the letters u, d/dx
  To accomplish this, I will take the number i, and I will expand it in
  binary. 
  0's    will correspond to nothing
  1-8's  will correspond to ux, uy, uz, Bx, By, Bz, rho, P, respectively
  9-13's will correspond to d/dt, d/dx, d/dy, d/dz, respectively

  INPUT:
  i - a number that will uniquely identify a word

  OUTPUT:
  valid  - a Boolean flag for if you should consider this library term
  fields - [1,num_f] a vector communicating the fields of interest
  derivs - [4,num_f] a matrix telling you the number of derivatives on each
           term
  %}

  f_idx = 1:num_fields;         %field indices
  d_idx = num_fields + (1:dim); %derivative indices

  base  = num_fields + dim + 1; %expand in this base

  %Compute the length of the binary sequence
  n = ceil( max(log(i),1)/log(base) ) + 1;

  digits = zeros(1,n);
  for j = 1:n
    digits(n+1-j) = mod(i,base);
    i = (i - mod(i,base))/base;
  end

  %Now we should check if the term is "valid". Set this to False by default
  valid  = false;
  num_f  = sum(digits <= num_fields) - sum(digits == 0); %count the number of fields IN THIS TERM
  derivs = zeros(4, num_f);
  fields = zeros(1, num_f);

  %First note that if we do not end with a field, the term is invalid
  if( ~any( digits(end) == 1:num_fields)  )
    return;
  end

  j = 1; %index over digits
  %Read through and ignore leading zeros
  while( digits(j) == 0 )
    j = j + 1;
  end



  %Read the rest
  k = 1; %k is our counter over fields
  last_deriv = 0;
  while( j <= n)
    if( digits(j) == 0 )
      return; 
      %found a useless space, which invalidates this term
    end

    %Check if we landed on a field
    if( any(digits(j) == f_idx) )
      fields(k) = digits(j);
      k = k+1; %Go to next u
      last_deriv = 0; %reset the last deriv we read
    end

    %Check if we landed on a derivative
    if( any(digits(j) == d_idx) )
      d = digits(j) - num_fields; %add a derivative in this dimension
      derivs(d,k) = derivs(d,k) + 1;
      if( last_deriv > d )
        %invalid ordering of partial derivatives.
        return;
      end
      last_deriv = d;
    end

    j = j+1;
  end

  %Check that the fields are sorted
  if( any(fields ~= sort(fields)) )
    return;
  end

  %Check that among like fields, derivatives are sorted
  for i = 1:num_fields
    I = (fields == i);
    if( sum(I) == 0)
      continue; %nothing to check
    end
 
    der = derivs(:, I); %subsample to size [dim, sum(I)];
    der = sum( der.*(dim.^(0:dim-1)') ); %map to a [1, sum(I)];
    if( any(der ~= sort(der)) )
      return; %derivatives not sorted between like fields!
    end
  end
  valid = true;
end