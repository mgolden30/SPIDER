%{
%}

clear;


num_fields = 8; %\vec{u}, \vec{B}, rho, P
dim = 4; %(t,x,y,z)

base = num_fields+dim+1;
n    = 4; %number of letters we allow

nl = 0;

for i = 1:base^(n)
  i;
  [valid, fields, derivs, digits] = check_library_term( i, num_fields, dim );
  if valid
    digits;
    [fields;derivs];
    nl = nl + 1;
    %return
  end
end

nl