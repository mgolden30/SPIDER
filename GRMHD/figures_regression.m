clear;


%%Vector Library First

addpath("../SPIDER_functions/");
load("vector_G.mat");

%starting model is \partial_t B_i
c0 = zeros( size(G,2),1  );
c0(1) = 1;

tic
n_max = numel(c0);
[cs, res] = greedy_regression_pure_matlab_add(G, c0, n_max);
toc

xt = [1, 4,16]; %xticks
optimization_curve( res, xt );
drawnow;
exportgraphics( gcf, "figs/induction.pdf" ); 

k=2;
print_model( k, labels, cs)





%% Now momentum equation

%Discard first element
G(:,1)    = [];
labels(1) = [];
scales(1) = [];

c0 = zeros( size(G,2),1  );
c0(15) = 1;

tic
%[cs, res] = greedy_regression_pure_matlab_add(G, c0, n_max);
[cs, res] = greedy_regression_pure_matlab(G );
toc

optimization_curve( res, xt );
drawnow;

exportgraphics( gcf, "figs/momentum.pdf" ); 

k=5;
print_model( k, labels, cs)


%% Now scalar 

clear;
load("scalar_G.mat");
%starting model is \partial_t B_i
c0 = zeros( size(G,2),1  );
c0(1) = 1;

tic
n_max = numel(c0);
[cs, res] = greedy_regression_pure_matlab(G);
toc

xt = [1, 16, 256]; %xticks
optimization_curve( res, xt );
drawnow;
exportgraphics( gcf, "figs/gauss.pdf" ); 

%%
k=3;
print_model( k, labels, cs)

%% Now 

bad = find(cs(:,3));
labels(bad) = [];
G(:,bad) = [];
scales(bad) = [];

tic
n_max = numel(c0);
[cs, res] = greedy_regression_pure_matlab(G);
toc

xt = [1, 16, 256]; %xticks
optimization_curve( res, xt );
drawnow;
exportgraphics( gcf, "figs/continuity.pdf" ); 

%%
k=4;
print_model( k, labels, cs)


%% Lastly, save the library nicely in latex

load("scalar_G.mat");

% Specify the file name
filename = 'latex_equations.txt';

% Open the file for writing
fileID = fopen(filename, 'w');

% Write each equation to the file
for i = 1:numel(labels)
    fprintf(fileID, '%s\n', labels{i});
end

% Close the file
fclose(fileID);

disp('Latex equations written to file successfully.');
function optimization_curve( res, xt )
  ms  = 100;
  fs  = 24;
  fsa = fs/2;

  col = [96 224 151]/256;
  scatter(1:numel(res), res, ms, 'o', "MarkerFaceColor", col, "MarkerEdgeColor", "black", "LineWidth", 2 );
  set(gca, 'xscale', 'log');
  set(gca, 'yscale', 'log');
  xticks(xt);

  set(gca,"fontsize", fsa);

  xlabel("terms in model", "fontsize", fs, "interpreter", "latex");
  ylabel("residual", "fontsize", fs, "interpreter", "latex");
  set(gcf, "color", "w");
end

function print_model( k, labels, cs)
  c = cs(:,k);

  idx = find(c ~= 0)
  
  labels = labels(idx);
  c      = c(idx);
  norm = max(abs(c));
  str = "";
  for i = 1:k
    str = str + " + " + c(i)/norm + " " + labels{i};
  end
  str
end