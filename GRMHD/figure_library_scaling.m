%{
%}

clear;

addpath("library_generation");

num_fields = 8; %\vec{u}, \vec{B}, rho, P
dim = 4; %(t,x,y,z)

ns = [];
ws = 1:5;
for w = ws
  [nl] = compute_library_size( w, num_fields, dim );
  ns = [ns, nl];
end

%%
clf;
ms = 100;
fs = 16;
scatter( ws, ns, ms, 's', 'filled', "markerfacecolor", "black" );
hold on
  plot( ws, ns, "linewidth", 2, "color", "black" );
  plot( ws, 12.^ws, "linewidth", 2, "color", "red" );
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

x = [0.45 0.45];
y = [0.7 0.6];
annotation('textarrow',x,y,'String','$12\,{}^n$', ...
      'interpreter', 'latex', 'fontsize', fs, "color", "red" );

drawnow;
exportgraphics(gcf, "../modified_KS/figs/library_scaling.pdf");

function [nl] = compute_library_size( word_limit, num_fields, dim )
  base = num_fields+dim+1;
  nl = 0;  
  for i = 1:base^(word_limit)
    [valid, fields, derivs, digits] = check_library_term( i, num_fields, dim );
    if valid
      %digits;
      [fields;derivs];
      nl = nl + 1;
      %return
    end
  end
end