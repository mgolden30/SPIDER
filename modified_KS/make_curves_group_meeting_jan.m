%{
Make a "publication worthy" optimization curve
%}

clear;
restoredefaultpath();
addpath("../SPIDER_functions/")

load("seed_1.mat");
[cs, residuals] = greedy_regression_pure_matlab_naive( G );

load("seed_543212345.mat");
residuals2 = vecnorm( G*cs );

%%
clf
ms1 = 50; %marker size
ms2 = 80;

fs = 16; %fontsize

color1 = [0.5 0 0];
color2 = [0 0 0.5];
scatter(1:numel(residuals), residuals,  ms1, 'd', 'filled', 'MarkerFaceColor', color1, 'MarkerEdgeColor', color1 );
hold on
scatter(1:numel(residuals), residuals2, ms2, 's', 'filled', 'MarkerFaceColor', color2, 'MarkerEdgeColor', color2);
hold off
set(gca, 'yscale', 'log' );
set(gca, 'xscale', 'log' );
xlabel("nonzero terms (model complexity)", "interpreter", "latex", "FontSize", fs);
ylabel("residual (model error)", "interpreter", "latex", "FontSize", fs);
xticks([1, 10, 100]);

set(gca,"fontsize", fs);
set(gcf, 'color', 'w');

%Save figure
saveas( gcf, "figs/modified_KS_curve.pdf" );