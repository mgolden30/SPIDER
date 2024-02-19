%{
Write the 
%}

clear;
load("trajectory.mat");

file = fopen("library.tex", "w");

fwrite(file, "\mathcal{L} = \{ & ");


per_line = 3;
for i = 1:numel(labels)
  fwrite( file, labels{i} );
  fwrite( file, ", \, " );
  
  if mod( i, per_line ) == 0 
    fwrite(file, "\\");
    fwrite(file, "&");
  end
 
end
fwrite(file, "\}")
fclose(file);

