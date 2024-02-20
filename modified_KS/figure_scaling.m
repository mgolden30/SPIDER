%{
Scaling figures
%}

clear;
addpath("../SPIDER_functions")

seed = 1;
rng(seed);

ns = round( 32*1.1.^(0:30) )
%return

times_slow = [];
times_fast = [];

for n = ns
  A = 2*rand(n,n) - 1;
  
  tic
  [cs, residuals] = greedy_regression_pure_matlab_naive( A );
  times_slow = [times_slow, toc];

  tic
  eps = 0;
  [cs, residuals] = greedy_regression_pure_matlab( A, eps );
  times_fast = [times_fast, toc];
end

save("times");

%%
load("times.mat");

fs = 12;
lw = 2;
ms = 100;
clf
%plot(ns, times_slow, "linewidth", lw, "color", "black");
%hold on
plot(ns, (ns/90 ).^4, "linewidth", lw, "color","red", "linestyle", "-");
hold on
plot(ns, (ns/160).^3, "linewidth", lw, "color","red", "linestyle", "-");

scatter( ns, times_slow, ms, 'd', 'filled', 'markerfacecolor', 'black')
hold on
%plot(ns, times_fast,  "linewidth", lw, "color", "black");
scatter( ns, times_fast, ms, 's', 'filled','markerfacecolor', 'black')
hold off
set( gca, "xscale", "log" );
set( gca, "yscale", "log" );
ylim([1e-2, 1.3e2]);
xticks([32, 64, 128, 256]);
xlim([60, 270]);
yticks([1e-2, 1e0, 1e2]);
ylabel("walltime (s)", "interpreter", "latex");
xlabel("library size $n$", "interpreter", "latex");
legend({ '','','naive', 'accelerated'}, "location", "SouthEast", "interpreter", "latex");
set(gca,"fontsize", fs);

%Now decide how big I want this figure in inches
figure = gcf;
im_size = round([500,300]*1);
figure.Position = [400,400,im_size];
set(gcf, "color", "white");

%Add annotations for the power law scaling
annotation('textarrow',[.4,.4],[.7,.6],'String','$n^4$', 'Interpreter', 'latex', 'fontsize', fs, 'color', 'red');
annotation('textarrow',[.4,.4],[.25,.35]+0.03,'String','$n^3$', 'Interpreter', 'latex', 'fontsize', fs, 'color', 'red');
exportgraphics( gcf, "figs/scaling.pdf" );