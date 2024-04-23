function [term, str, scales] = generate_library_term( u, deriv_matrix, derivs )
  %{
  PURPOSE:
  Generate the term of interest and a string to describe it based on
  "derivs". us is a cell, where us{k,l} is the kth x derivative and lth time derivative of u.

  OUTPUT:
  term - a numerical representation of this term to integrate
  str  -  string representation of the term
  scales - the scale S_j as described in the arxiv preprint.
  %}


  u_mean = mean(abs(u),"all");
  u_std  = std(abs(u),0,"all");
  T      = u_std/mean( abs(deriv_matrix{2,1}),"all" ); %time scale
  L      = u_std/mean( abs(deriv_matrix{1,2}),"all" ); %length scale
  
  %Initialize blank slates
  term = 0*u+1;
  str = "";
  scales = 1;

  if numel(derivs) == 0
    return;
  end

  for i = 1:size(derivs, 2)
    m = derivs(1,i); %time  derivatives
    n = 0;%derivs(2,i); %space derivatives
    if( m + n  == 0 )
      term = term .* u;
      str = str + "u";
      scales = scales * u_mean;
    else
      %use a derivative
      term = term .* deriv_matrix{m+1,n+1};
      str  = str + "(dx^" + (m) + " u)";
      scales = scales * u_std / L^(m) / T^(n);
    end
  end 
end