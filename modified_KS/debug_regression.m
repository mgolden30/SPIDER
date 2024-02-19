%{
My optimization curves are currently ever so slightly, but importantly
different. I wnt to check for numerical instability in my bisection method.
%}

clear;
restoredefaultpath();
addpath("../SPIDER_functions/")
load("seed_1.mat");

tic
[cs,  residuals ] = greedy_regression_pure_matlab_naive( G );
t2 = toc

%%
eps = 0e-7;
tic
[cs2, residuals2, fs] = greedy_regression_pure_matlab( G, eps );
t3 = toc

[t2, t3, t2/t3]

figure(1);
tiledlayout(2,2);
nexttile
imagesc(cs ~= 0);
title("naive");
nexttile
imagesc(cs2~= 0);
title("fast");
nexttile
loglog(residuals, 'o');
nexttile
loglog(residuals2, 'o');
return

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
    legend( {'fast', 'slow'}, 'location', 'southwest' );

xlabel("terms in relation $k$", "interpreter", "latex");
ylabel("residual $r_k$", "interpreter", "latex");
set(gca, "fontsize", fs );
set(gcf, "color", "w");

pbaspect([2,1,1]);

%Add annotations
x = [0.5 0.35];
y = [0.8 0.44];
%annotation('textarrow',x,y,'String','$\partial_t u + u \partial_x u + \partial_x^2 u + \partial_x^4u = 0$', ...
%    'interpreter', 'latex', 'fontsize', fs );

x = [0.6 0.46];
y = [0.6 0.24];
%annotation('textarrow',x,y,'String','$\cdots + c_3 \partial_x u^3 + c_4 \partial_x u^4$', ...
%    'interpreter', 'latex', 'fontsize', fs );

x = [0.63 0.2];
y = [0.6 0.4];
dim = [x,y];
%annotation('textbox',dim,'String','$+c_5 \partial_x u^5 + c_6 \partial_x u^6 = 0$', ...
%    'interpreter', 'latex', 'fontsize', fs, 'EdgeColor', [1,1,1,0] );

f = gcf;
set(gcf, "Position", [91.4000 365.8000 716.8000 442.4000] );

saveas(gcf, "figs/optimiation2.pdf");


%% Make a figure showing the alignment of the two algorithms
nl = size(G,2);
cos_mat = zeros( [nl,nl] );

for i = 1:nl
  for j = 1:nl
    ci = cs(:,i) ~= 0;
    cj = cs2(:,j) ~= 0;
    
    cos_mat(i,j) = dot(ci,cj)/sqrt(sum(ci)*sum(cj));
  end
end
clf
imagesc( cos_mat );

loglog( 1:nl, diag(cos_mat) );