function [term, str] = generate_library_term( u, us, derivs )
  %{
  PURPOSE:
  Generate the term of interest and a string to describe it based on
  "derivs". us is a cell, where us{k} is the kth derivative of u.
  %}

  term = 0*u+1;
  str = "";
  if numel(derivs) == 0
    return;
  end


  for i = 1:numel(derivs)
    k = derivs(i);
    if( k == 0 )
      term = term .* u;
      str = str + "u";
    else
      %use a derivative
      term = term .* us{k};
      str  = str + "(dx^" + k + " u)";
    end
  end 
end