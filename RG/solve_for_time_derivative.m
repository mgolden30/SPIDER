
clear;

load("Gcomp3.mat");

%Split up the library into a target term b
b = G(:, 1);
target = labels{1};
target_scale = scales(1);

%And the rest of the library
A = G(:, 2:end);
labels = labels(2:end);
scales(1) = [];

tic
[cs, residuals] = inhomogeneous_regression( A, b );
toc

scatter( 1:numel(residuals), residuals, 'o', 'filled' );


%normalization = 1;
for interest=1:numel(residuals)
  c = cs(:,interest); 
  str = interest + ":  residual = " + residuals(interest) +  ", \quad   " + target + " = ";
  for i = 1:numel(c)
    if( c(i) == 0 )
      continue; 
    end

   str = str + " + " + c(i)  + labels{i};

  end

  str = str + " = 0";

  %For printing, replace \ with \\
  str = strrep(str, "\", "\\");

  fprintf(str + "\n");
 end
