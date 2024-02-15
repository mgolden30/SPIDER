function [valid, derivs, digits] = check_library_term( i )
  %{
  PURPOSE:
  I need words made up of the letters u, d/dx
  To accomplish this, I will take the number i, and I will expand it in
  binary. 
  0's will correspond to nothing
  1's will correspond to u
  2's will correspond to d/dx

  INPUT:
  i - a number that will uniquely identify a word

  OUTPUT:
  valid - a Boolean flag for if you should consider this library term or
  not
  %}

  base = 3; %expand in this base

  %Compute the length of the binary sequence
  n = ceil( max(log(i),1)/log(base) ) + 1;

  digits = zeros(1,n);
  for j = 1:n
    digits(n+1-j) = mod(i,base);
    i = (i - mod(i,base))/base;
  end

  %Now we should check if the term is "valid". Set this to False by default
  valid  = false;
  num_u  = sum(digits == 1); %count the number of u's
  derivs = zeros( 1, num_u);

  %First note that if we do not end with a u, the term is invalid
  if( digits(end) ~= 1 )
    return;
  end

  j = 1; %index over digits
  %Read through and ignore leading zeros
  while( digits(j) == 0 )
    j = j + 1;
  end
  
  %Read the rest
  k = 1; %k is our counter over u's
  while( j <= n)
    if( digits(j) == 0 )
      return; %found a useless space
    end

    if( digits(j) == 1 )
      k = k+1; %Go to next u
    end

    if( digits(j) == 2 )
      derivs(k) = derivs(k) + 1;
    end

    j = j+1;
  end

  %lastly, sort derivs
  if( any(derivs ~= sort(derivs)))
    return;
  end
  
  valid = true;
end