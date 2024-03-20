%{
This script provides some methods for finding sparse relations from feature
matrices G
%}

%% Greedy_regression

addpath('SPIDER_functions/');

%Pick a residual to minimize. I have played with three choices, and choice
%1 and 3 are the most sensible. Choice #2 will likely be deprecated
%eventually.
%1: |Gc|
%2: |Gc|/|G_restricted|
%3: |Gc|/max|G.*c| <- relative residual
index = 3;


%This choice tells sparse regression to start with the full library and
%work our way backwards, discarding terms in a greedy fashion.
starting_model = ones( size(G,2), 1 );

number_of_subsamples = 100; %Do parameter estimation 100 times to estimate uncertainty
fraction_of_windows_in_subsample = 1/2; %use half of windows randomly each time we estimate parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start Boilerplate. Probably don't touch this.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%We can pick any of these residuals to minimize
res1 = @(G,c) norm(G*c)/norm(c);
res2 = @(G,c) norm(G*c)/norm(c)/norm(G);
res3 = @(G,c) norm(G*c)/max( vecnorm(G*diag(c)) );
residual_labels = { "$|Gc|$", "$|Gc|/|G_{c\neq 0}|$", "$|Gc|/\textrm{max}\{|G c_n|\}$" };
residual_functions = {res1, res2, res3}; %residuals to compute during sparsification
res_func_for_sparsification = residual_functions{index};
[res_ave, res_std, cs_ave, cs_std] = greedy_regression( G, starting_model, ...
                                                                   number_of_subsamples, ...
                                                                   fraction_of_windows_in_subsample, ...
                                                                   residual_functions, ...
                                                                   res_func_for_sparsification );
          

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End Boilerplate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Plot results                                                  
fs = 16; %Font Size
lim = size(res_ave,2);
h = scatter( 1:lim, res_ave(index,1:lim), 's', 'MarkerEdgeColor', 'black' );
set(get(h,'Parent'), 'YScale', 'log');

xlabel('N', 'interpreter', 'latex', 'FontSize', fs);
ylabel( residual_labels{index}, 'interpreter', 'latex', 'FontSize', fs );
h1 = get(gca,'YLabel');
xticks(1:numel(res_ave))
xlim( [0.5, size(res_ave,2)+0.5] );
set(gcf,'color','w');
pbaspect([2 1 1])

%Print some models to screen
for interest = 2:min(8, size(G,2))
  print_model( interest, cs_ave, res_ave, index, labels, scales )
end

return;

%% Combinatoric Search

addpath('SPIDER_functions/');

number_of_terms = 3;
%residual_type = "|Gc|"; %Use this residual or the next
residual_type = "|Gc|/max";
[models, residuals, cs] = combinatoric_search( G, labels, number_of_terms, scales, residual_type );

%Print to screen best 10 models
models(1:min(10, numel(models)))





function str = print_model( interest, cs_ave, res_ave, index, labels, scales )
  c = cs_ave(:,interest);
  c = c./scales';

  normalization = max(abs(c)); %normalize by the maximal term
  
  str = interest + ":  residual = " + res_ave(index,interest) + ", \quad   ";
  for i = 1:numel(c)
    if( c(i) == 0 )
      continue; 
    end
  
    str = str + " + " + c(i)/normalization + "  " + labels{i};
  end
  str = str + " = 0";
end