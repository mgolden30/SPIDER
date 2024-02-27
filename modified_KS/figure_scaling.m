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
times_add = [];
for n = ns
  n
  A = 2*rand(n,n) - 1;
  %{
  tic
  [cs, residuals] = greedy_regression_pure_matlab_naive( A );
  times_slow = [times_slow, toc];

  tic
  eps = 0;
  [cs, residuals] = greedy_regression_pure_matlab( A, eps );
  times_fast = [times_fast, toc];
  %}

  c0 = zeros(n,1);
  c0([1,3]) = 1;
  tic
  max_n = min( [64, n] );
  [cs, residuals] = greedy_regression_pure_matlab_add( A, c0, max_n );
  times_add = [times_add, toc];
end

save("times2");

%%
load("times2.mat");
load("times.mat");

fs = 12;
lw = 2;
ms = 100;
clf
%plot(ns, times_slow, "linewidth", lw, "color", "black");
%hold on


p1 = polyfit( log(ns), log(times_slow), 1 );
p2 = polyfit( log(ns), log(times_fast), 1 );
p3 = polyfit( log(ns), log(times_add),  1 );

plot( ns, exp( p1(1)*log(ns) + p1(2) ), "linewidth", lw, "color","red", "linestyle", "-");
hold on
plot( ns, exp( p2(1)*log(ns) + p2(2) ), "linewidth", lw, "color","red", "linestyle", "-");
plot( ns, exp( p3(1)*log(ns) + p3(2) ), "linewidth", lw, "color","red", "linestyle", "-");

scatter( ns, times_slow, ms, 'd', 'filled', 'markerfacecolor', 'black')
scatter( ns, times_fast, ms, 's', 'filled', 'markerfacecolor', 'black')
scatter( ns, times_add,  ms, 'o', 'filled', 'markerfacecolor', 'black')
hold off
set( gca, "xscale", "log" );
set( gca, "yscale", "log" );
ylim([1e-2, 2e3]);
xticks([32, 64, 128, 256, 512]);
xlim([60, 530]);
yticks([1e-2, 1e0, 1e2]);
ylabel("walltime (s)", "interpreter", "latex");
xlabel("library size $n$", "interpreter", "latex");
legend({ '','','','naive', 'ANUBIS--','ANUBIS+'}, "location", "SouthEast", "interpreter", "latex");
set(gca,"fontsize", fs);

%Now decide how big I want this figure in inches
figure = gcf;
im_size = round([500,300]*1);
figure.Position = [400,400,im_size];
set(gcf, "color", "white");

%Add annotations for the power law scaling
annotation('textarrow',[.4,.4],[.7,.6],'String','$n^4$', 'Interpreter', 'latex', 'fontsize', fs, 'color', 'red');
annotation('textarrow',[.5,.5],[.28,.48],'String','$n^3$', 'Interpreter', 'latex', 'fontsize', fs, 'color', 'red');
annotation('textarrow',[.3,.3],[.24,.35],'String','$n^2$', 'Interpreter', 'latex', 'fontsize', fs, 'color', 'red');

exportgraphics( gcf, "figs/scaling.pdf" );


polyfit( log(ns), log(times_slow), 1 )
polyfit( log(ns), log(times_fast), 1 )
polyfit( log(ns), log(times_add),  1 )
polyfit( log(ns(1:end-5)), log(times_add(1:end-5)), 1 )