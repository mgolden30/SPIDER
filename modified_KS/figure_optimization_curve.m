clear;
restoredefaultpath();
addpath("../SPIDER_functions/")

%remove time derivative from library?
remove_td = false;
seed = 1;
[r_e, r_p, r_m ] = optimization_curves(seed, remove_td);



%%

figure(1);
clf
tiledlayout(2,2);
nexttile
%imagesc(cs ~= 0);
title("naive");
nexttile
imagesc(cs2~= 0);
title("fast");
nexttile
   % loglog(residuals, 'o');
nexttile
loglog(residuals2, 'o');

%% Make a nice figure
clf;

color2 = [62, 218, 130]/256;
color1 = [0 ,0, 0];

fs = 16;
ms = 100;
scatter( 1:numel(residuals2), residuals2, ms, 'o', 'filled', 'markerfacecolor', color1 );
hold on
scatter( 1:numel(residuals), residuals, ms/3, 's', 'filled', 'markerfacecolor', color2 );
hold off
set(gca, "xscale", "log");
set(gca, "yscale", "log");

xlabel("terms in relation $k$", "interpreter", "latex");
ylabel("residual $r_k$", "interpreter", "latex");
set(gca, "fontsize", fs );
set(gcf, "color", "w");

pbaspect([1,1,1]);


if(~remove_td)
%Add annotations
x = [0.49 0.42];
y = [0.65 0.48];
annotation('textarrow',x,y,'String','$\partial_t u + u \partial_x u + \partial_x^2 u + \partial_x^4u$', ...
    'interpreter', 'latex', 'fontsize', fs );

x = [0.61 0.48];
y = [0.55 0.21];
annotation('textarrow',x,y,'String','$\cdots + \sum_{k=3}^6 c_k \partial_x u^k$', ...
    'interpreter', 'latex', 'fontsize', fs );
%{
x = [0.63 0.17];
y = [0.5  0.35];
dim = [x,y];
annotation('textbox',dim,'String','$+c_5 \partial_x u^5 + c_6 \partial_x u^6 = 0$', ...
    'interpreter', 'latex', 'fontsize', fs, 'EdgeColor', [1,1,1,0] );
%}
else
  x = [0.44 0.44];
  y = [0.55 0.67];
  annotation('textarrow',x,y,'String','$C( u ,\partial_x^n u) = 0$', ...
      'interpreter', 'latex', 'fontsize', fs );
end

if(remove_td)
  ks = [13];
else
  ks = [4,8];
end

ylim([1e-7, 1e0]);
yticks([1e-6, 1e-4, 1e-2, 1e0]);


hold on
scatter( ks, residuals(ks), ms, "x", "linewidth", 3, "MarkerEdgeColor", "red" );
hold off

legend( {'ANUBIS-', 'Thorough', ''}, 'location', 'NorthEast' );


f = gcf;
set(gcf, "Position", [91.4000 365.8000 716.8000 442.4000] );

if(remove_td)
  exportgraphics( gcf, "figs/optimization2.pdf" );
else
  exportgraphics( gcf, "figs/optimization.pdf" );
end


function [r_e, r_p, r_m ] = optimization_curves(seed, remove_td)
  %{
  PURPOSE:
  Do sparse regression three ways and record the resulting curves
  
  r_e - exhaustive search
  r_p/m - ANUBIS+-
  %}
  
  load("seed_"+seed+".mat");

  if(remove_td)
    G(:,1) = [];
    labels(1) = [];
  end

  %Exhaustive
  tic
  [cs_e,  r_e ] = exhaustive_search( G );
  t_e = toc;
  fprintf("Exhasutive Search: %.2f seconds...\n", t_e);

  %ANUBIS-minus
  tic
  [cs_m, r_m] = greedy_regression_pure_matlab( G );
  t_m = toc
  fprintf("ANUBIS-: %.2f seconds...\n", t_m);

  %Need an initial guess
  c0 = zeros( size(G,2),1 ); 
  desired_terms = {"\partial_t u", "u(dx^1 u)"};%, "(dx^2dt^0 u)", "(dx^4dt^0 u)"};
  for i = 1:numel(labels)
    for j = 1:numel(desired_terms)
      if( labels{i} == desired_terms{j} )
        c0(i) = 1;
        i
      end
    end
  end
  n_max = 50;
  tic
  [cs_p, r_p] = greedy_regression_pure_matlab_add( G, c0, n_max );
  t_p = toc
  fprintf("ANUBIS+: %.2f seconds...\n", t_p);
end