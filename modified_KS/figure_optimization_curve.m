clear;
restoredefaultpath();
addpath("../SPIDER_functions/")


remove_td = false; %remove time derivative from library?
seed = 1;
[r_e, r_p, r_m ] = optimization_curves(seed, remove_td);


%% Make a nice figure
clf;

%color3 = [62, 218, 130]/256;
%color2 = [62, 218, 130]/256/2;
color1 = [0 ,0, 0];
color2 = [0.5, 0.5, 1];
color3 = [0.8, 0, 0];

fs = 16;
ms = 200;
scatter( 1:numel(r_e), r_e, ms, 'o', 'filled', 'markerfacecolor', color1 );
hold on
scatter( 1:numel(r_m), r_m, ms/2, 's', 'filled', 'markerfacecolor', color2 );
if(~remove_td)
scatter( 1:numel(r_p), r_p, ms/4, 'd', 'filled', 'markerfacecolor', color3 );
end
hold off
set(gca, "xscale", "log");
set(gca, "yscale", "log");
set(gca, "box", "on");

if(~remove_td)
  ylabel("residual $r_k$", "interpreter", "latex");
end
xlabel("terms in relation $k$", "interpreter", "latex");

set(gca, "fontsize", fs );
set(gcf, "color", "w");

pbaspect([1,1,1]);


if(~remove_td)
%Add annotations
x = [0.45 0.42];
y = [0.6 0.5];
annotation('textarrow',x,y,'String','$\partial_t u + u \partial_x u + \partial_x^2 u + \partial_x^4u$', ...
    'interpreter', 'latex', 'fontsize', fs );

x = [0.54 0.48];
y = [0.4 0.26];
annotation('textarrow',x,y,'String','$+ \sum_{k=3}^6 c_k \partial_x u^k$', ...
    'interpreter', 'latex', 'fontsize', fs );
%{
x = [0.63 0.17];
y = [0.5  0.35];
dim = [x,y];
annotation('textbox',dim,'String','$+c_5 \partial_x u^5 + c_6 \partial_x u^6 = 0$', ...
    'interpreter', 'latex', 'fontsize', fs, 'EdgeColor', [1,1,1,0] );
%}
else
  x = [0.46 0.46];
  y = [0.55 0.67] - 0.04;
  annotation('textarrow',x,y,'String','$C( u ,\partial_x^n u) = 0$', ...
      'interpreter', 'latex', 'fontsize', fs );
end

if(remove_td)
  ks = [13];
else
  ks = [4,8];
end

xlim([0.8, 200])
ylim([1e-8/2, 5e0]);
yticks([1e-8, 1e-4, 1e0]);
if remove_td
  %yticks([]);
end


%hold on
%scatter( ks, r_e(ks), ms, "x", "linewidth", 3, "MarkerEdgeColor", "red" );
%hold off

if(~remove_td)
legend( {'exhaustive', 'ANUBIS--', 'ANUBIS+'}, 'location', 'NorthEast', 'interpreter', 'latex' );
end

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
  t_m = toc;
  fprintf("ANUBIS-: %.2f seconds...\n", t_m);

  %find out the last term to use
  labels( cs_m(:,1) ~=0 )

  %Need an initial guess
  %{
  c0 = zeros( size(G,2),1 ); 
  desired_terms = {"\partial_t u"};%, "u(dx^1 u)"};%, "(dx^2dt^0 u)", "(dx^4dt^0 u)"};
  for i = 1:numel(labels)
    for j = 1:numel(desired_terms)
      if( labels{i} == desired_terms{j} )
        c0(i) = 1;
        i
      end
    end
  end
  %}
  c0 = cs_m(:,2);
  n_max = numel(c0);
  tic
  [cs_p, r_p] = greedy_regression_pure_matlab_add( G, c0, n_max );
  t_p = toc;
  fprintf("ANUBIS+: %.2f seconds...\n", t_p);
end