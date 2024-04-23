%{
%}

clear;

num_fields = 2; %\vec{u}, \vec{B}, rho, P
dim = 2; %(t,x,y,z)

ns = [];
ws = 4;
for w = ws
  tic
  [nl] = compute_library_size( w, num_fields, dim );
  fprintf("walltime %e\n", toc);

  ns = [ns, nl]
end

ns
return


%%
clf;

load("lib_size.mat");

ws = ws(1:numel(ns)); %since I killed this counting early

ms = 100;
fs = 16;
scatter( ws, ns, ms, 's', 'filled', "markerfacecolor", "black" );
hold on
  plot( ws, ns, "linewidth", 2, "color", "black" );
hold off

set(gca, "yscale", "log");
xlim([0.5, max(ws)+0.5]);
xticks( ws )
xlabel("number of letters $n$", "interpreter", "latex");
ylabel("$|\mathcal{L}_n|$", "interpreter", "latex", "rotation", 0 );

ylim([3, 3e4]);
yticks();
set(gcf, "color", "w");
set(gca, "fontsize", fs);

function [nl] = compute_library_size( word_limit, num_fields, dim )
  base = num_fields+dim+1;
  nl = 0;

  file = fopen("representions.tex", "w");
  fprintf(file, "\\documentclass{article}\n \\usepackage{color} \\usepackage{amsmath} \\begin{document} \\noindent");

  for i = 1:base^(word_limit)
    [valid, fields, derivs, digits] = check_library_term( i, num_fields, dim );
    if valid
      %digits;
      representation = [derivs;fields];
      nl = nl + 1;

      %Add as many preceeding columns of zeros as needed.
      rep = zeros(dim+1, word_limit);
      rep(:, end-size(representation,2)+1:end) = representation;

      representation = rep;

      fprintf(file, "\nterm %d:\n$$\\begin{pmatrix}", nl);
      for i = 1:size(representation,1)
        for j = 1:size(representation,2)
          if representation(i,j) == 0
            fprintf( file, "{\\color{red} %d}", representation(i,j) );
          else
            fprintf( file, "%d", representation(i,j) );
          end
          if j ~= size(representation,2)
            fprintf( file, " & " );
          end
        end
        fprintf(file, "\\\\\n");
        if (i == size(representation,1) - 1)
          fprintf(file, "\\hline");
        end
      end
      fprintf(file, "\\end{pmatrix}$$");
      %return
    end
  end
  fprintf(file, "\\end{document}");

  fclose(file);
end