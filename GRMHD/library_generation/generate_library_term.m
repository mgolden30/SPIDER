function [term, str, scales] = generate_library_term( u, us, derivs )
  %{
  PURPOSE:
  Generate the term of interest and a string to describe it based on
  "derivs". us is a cell, where us{k} is the kth derivative of u.

  OUTPUT:
  term - a numerical representation of this term to integrate
  str  -  string representation of the term
  scales - the scale S_j as described in the arxiv preprint.
  %}


  u_mean = mean(abs(u),"all");
  u_std  = std(abs(u),0,"all");
  L      = u_mean/mean( abs(us{1}),"all" );

  %Initialize blank slates
  term = 0*u+1;
  str = "";
  scales = 1;

  if numel(derivs) == 0
    return;
  end

  for i = 1:numel(derivs)
    k = derivs(i);
    if( k == 0 )
      term = term .* u;
      str = str + "u";
      scales = scales * u_mean;
    else
      %use a derivative
      term = term .* us{k};
      str  = str + "(dx^" + k + " u)";
      scales = scales * u_std / L^k;
    end
  end 
end