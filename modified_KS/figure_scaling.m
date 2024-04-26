%{
Scaling figures
%}

clear;
restoredefaultpath();

addpath("../SPIDER_functions")
addpath("D:/GISR-main/GISR/SINDy/SINDy-PI-master/Comparison/PDE_Comparison/Implicit_SINDy/Functions/");
%addpath("SINDy/SINDy-PI-master/Comparison/PDE_Comparison/Implicit_SINDy/TempFunctions/");

seed = 1;
rng(seed);

ns = round( 32*1.1.^(0:30) )
%return

num_trials = 8;

times_add = zeros( numel(ns), num_trials );
times_pi  = zeros( numel(ns), num_trials );

times_slow = [];
times_fast = [];

i = 1;
for n = ns
  n

  for trial = 1:num_trials
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

   %{
    c0 = zeros(n,1);
    c0([1,3]) = 1;
    max_n = min( [64, n] );
    tic
    [cs, residuals] = greedy_regression_pure_matlab_add( A, c0, max_n );
    times_add(i, trial) =  toc;
  %}

    tic
    Theta = A;
    tol = 1e-3;
    pflag = 0; %no plotting
    [Xi, indTheta, lambdavec, numterms, errorv] = ADMpareto( Theta, tol, pflag );
    times_pi(i, trial) = toc;
  end
  i = i + 1;  
end

save("times3");

%%
load("times3.mat");
load("times2.mat");
load("times.mat");

fs = 12;
lw = 2;
ms = 100;
clf
%plot(ns, times_slow, "linewidth", lw, "color", "black");
%hold on

%average over trials
times_add = mean(times_add, 2);
times_pi  = mean(times_pi,  2);

%restrict data to x-axis we are plotting so the power law fit looks decent
rel = (ns<520) & (ns>120);
ns         = ns(rel);
times_slow = times_slow(rel);
times_add  = times_add(rel);
times_fast = times_fast(rel);
times_pi   = times_pi(rel);

p1 = polyfit( log(ns), log(times_slow), 1 );
p2 = polyfit( log(ns), log(times_fast), 1 );
p3 = polyfit( log(ns), log(times_add),  1 );
p4 = polyfit( log(ns), log(times_pi),  1 );



%use these for scatter and line plots
my_colors = cell(4,1);
my_colors{1} = [0.0, 0.7, 0.4];
my_colors{2} = [0.6, 0.0, 0.5];
my_colors{3} = [0.5, 1.0, 1.0]/2;
my_colors{4} = [0.4, 0.6, 0.8];

plot( ns, exp( p1(1)*log(ns) + p1(2) ), "linewidth", lw, "color", my_colors{1}, "linestyle", "-");
hold on
plot( ns, exp( p2(1)*log(ns) + p2(2) ), "linewidth", lw, "color",my_colors{2}, "linestyle", "-");
plot( ns, exp( p3(1)*log(ns) + p3(2) ), "linewidth", lw, "color",my_colors{3}, "linestyle", "-");
plot( ns, exp( p4(1)*log(ns) + p4(2) ), "linewidth", lw, "color",my_colors{4}, "linestyle", "-");
pbaspect([1,2,1]);

scatter( ns, times_slow, ms, 'd', 'filled', 'markerfacecolor', my_colors{1} )
scatter( ns, times_fast, ms, 's', 'filled', 'markerfacecolor', my_colors{2} )
scatter( ns, times_add,  ms, '^', 'filled', 'markerfacecolor', my_colors{3} )
scatter( ns, times_pi,   ms, 'o', 'filled', 'markerfacecolor', my_colors{4} )

hold off
set( gca, "xscale", "log" );
set( gca, "yscale", "log" );

ylim([1e-1, 2e3]);

xticks([32, 64, 128, 256, 512]);
xlim([120, 530]);
yticks([1e-1, 1e1, 1e3]);

ylabel("walltime (s)", "interpreter", "latex");
xlabel("library size $n$", "interpreter", "latex");
%legend({ '','','','','ExhaustiveSearch', 'ANUBIS--','ANUBIS+','SINDy'}, "location", "SouthEast", "interpreter", "latex");
set(gca,"fontsize", fs);

%Now decide how big I want this figure in inches
figure = gcf;
im_size = round([300,500]*2);
figure.Position = [200,200,im_size];
set(gcf, "color", "white");

%Add annotations for the power law scaling
annot_col = "black";
annotation('textarrow',[.34,.34],[.78,.7 ],'String', sprintf("implicit-SINDy ($n^{%0.2f}$)", p4(1)), 'Interpreter', 'latex', 'fontsize', fs, 'color', annot_col);
annotation('textarrow',[.7,.7],[.67,.785 ], 'String', sprintf("exhaustive search ($n^{%0.2f}$)", p1(1)), 'Interpreter', 'latex', 'fontsize', fs, 'color',  annot_col);
annotation('textarrow',[.5,.5],[.48,.42],  'String', sprintf("ANUBIS-- ($n^{%0.2f}$)", p2(1)), 'Interpreter', 'latex', 'fontsize', fs, 'color',  annot_col);
annotation('textarrow',[.6,.6],[.28,.34],  'String', sprintf("ANUBIS+ $(n^{%0.2f}$)", p3(1)), 'Interpreter', 'latex', 'fontsize', fs, 'color',  annot_col);

exportgraphics( gcf, "figs/scaling.pdf" ); 

polyfit( log(ns), log(times_slow), 1 )
polyfit( log(ns), log(times_fast), 1 )
polyfit( log(ns), log(times_add),  1 )
polyfit( log(ns), log(times_pi),  1 )


%% Estimate time to evaluate a library of size 10,000

load("../GRMHD/library_generation/lib_size.mat");
ws = ws(1:numel(ns));
%ns = 10.^(3:8);

%use these for scatter and line plots
my_colors = cell(4,1);
my_colors{1} = [0.0, 0.7, 0.4];
my_colors{2} = [0.6, 0.0, 0.5];
my_colors{3} = [0.5, 1.0, 1.0]/2;
my_colors{4} = [0.4, 0.6, 0.8];
          

sindypi_time = exp( p4(1)*log(ns) + p4(2) )
exhaust_time = exp( p1(1)*log(ns) + p1(2) )
anubisp_time = exp( p2(1)*log(ns) + p2(2) )
anubism_time = exp( p3(1)*log(ns) + p3(2) )

clf;
scatter([], []);

hold on
plot( ws, exp( p1(1)*log(ns) + p1(2) ), "linewidth", lw, "color", my_colors{1}, "linestyle", "-");
plot( ws, exp( p2(1)*log(ns) + p2(2) ), "linewidth", lw, "color", my_colors{2}, "linestyle", "-");
plot( ws, exp( p3(1)*log(ns) + p3(2) ), "linewidth", lw, "color", my_colors{3}, "linestyle", "-");
plot( ws, exp( p4(1)*log(ns) + p4(2) ), "linewidth", lw, "color", my_colors{4}, "linestyle", "-");

ms = nan;
scatter( ws, exhaust_time,  ms, 'd', 'filled', 'markerfacecolor', my_colors{1})
scatter( ws, anubisp_time,  ms, 's', 'filled', 'markerfacecolor', my_colors{2})
scatter( ws, anubism_time,  ms, '^', 'filled', 'markerfacecolor', my_colors{3})
scatter( ws, sindypi_time,  ms, 'o', 'filled', 'markerfacecolor', my_colors{4})

hold off

day  = 60*60*24;
year = 60*60*24*365;
hubble = 14.4 * 1e9 * year;

ylim( [1, hubble*8] )

yline( day,  "linewidth", 2, "color", "black", "alpha", 1 );
yline( year, "linewidth", 2, "color", "black", "alpha", 1 );
yline(hubble, "linewidth", 2, "color", "black", "alpha", 1);
%set( gca, "xscale", "log" );
set( gca, "yscale", "log" );

xlabel("word size", "interpreter", "latex" );
ylabel("predicted walltime (s)", "interpreter", "latex" );

yticks( [1, 1e4, 1e8, 1e12, 1e16] );

%BOX 1: Day
px = 0.15 ; %position
py = 0.33;
wx = 0.16 ; %width
wy = 0.023;
dims = [px,py,wx,wy];
annotation('textbox', dims, 'String', 'Day', 'Interpreter', 'latex', 'fontsize', fs, 'color', 'black', 'edgecolor', 1*[1,1,1] );

%BOX 2: Year
px = 0.15 ; %position
py = 0.445;
wx = 0.16 ; %width
wy = 0.023;
dims = [px,py,wx,wy];
annotation('textbox', dims, 'String', 'Year', 'Interpreter', 'latex', 'fontsize', fs, 'color', 'black', 'edgecolor', 1*[1,1,1] );

%BOX 3: Hubble
px = 0.15 ; %position
py = 0.89;
wx = 0.32 ; %width
wy = 0.023;
dims = [px,py,wx,wy];
annotation('textbox', dims, 'String', 'Hubble time', 'Interpreter', 'latex', 'fontsize', fs, 'color', 'black', 'edgecolor', 1*[1,1,1] );


% Add labels to all the curves as well


figure;
ax = gca;
ax.Box = 'on';  % This sets the axis box to be visible


%ANUBIS+
px = 0.5; %position
py = 0.2;
wx = 0.32 ; %width
wy = 0.023;
dims = [px,py,wx,wy];
annotation('textbox', dims, 'rotation', 45, 'String', 'ANUBIS+', 'Interpreter', 'latex', 'fontsize', fs, 'color', 'black', 'edgecolor', 1*[1,1,1] );


%ANUBIS-
px = 0.5; %position
py = 0.34;
wx = 0.32 ; %width
wy = 0.023;
dims = [px,py,wx,wy];
annotation('textbox', dims, 'rotation', 65, 'String', 'ANUBIS--', 'Interpreter', 'latex', 'fontsize', fs, 'color', 'black', 'edgecolor', 1*[1,1,1] );


%SINDy
px = 0.7; %position
py = 0.53;
wx = 0.32 ; %width
wy = 0.023;
dims = [px,py,wx,wy];
annotation('textbox', dims, 'rotation', 55, 'String', 'implicit-SINDy', 'Interpreter', 'latex', 'fontsize', fs, 'color', 'black', 'edgecolor', 1*[1,1,1] );

%exhaustive search
px = 0.6; %position
py = 0.7;
wx = 0.32 ; %width
wy = 0.023;
dims = [px,py,wx,wy];
annotation('textbox', dims, 'rotation', 65, 'String', 'exhasutive search', 'Interpreter', 'latex', 'fontsize', fs, 'color', 'black', 'edgecolor', 1*[1,1,1] );



exportgraphics( gcf, "figs/predicted_scaling.pdf" );