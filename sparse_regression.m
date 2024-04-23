%{
This script provides some methods for finding sparse relations from feature
matrices G
%}


%% Combinatoric Search
addpath('SPIDER_functions/');

number_of_terms = 3;
%residual_type = "|Gc|"; %Use this residual or the next
residual_type = "|Gc|/max";
[models, residuals, cs] = combinatoric_search( G, labels, number_of_terms, scales, residual_type );

%Print to screen best 10 models
models(1:min(10, numel(models)))





%% Greedy_regression
% Add/remove a term at a time to see how sparsity effects the residual

index = 3; %use 1 2 or 3 depending on your residual choice.
%1: |Gc|
%2: |Gc|/|G_restricted|
%3: |Gc|/max|G.*c| <- relative residual

clf
%{
if exist('cs_ave')
  %You've run greedy regression before
  starting_K = 5;
  starting_model = cs_ave(:, starting_K);
else
  %You have not run greedy regression yet. Start with combinatoric result
  starting_model = cs(1,:);
  starting_K = sum( starting_model ~= 0 );
end
%}

starting_model = ones( size(G,2),1 );
number_of_subsamples = 100; %Do parameter estimation 100 times to estimate uncertainty
fraction_of_windows_in_subsample = 1/2; %use half of windows randomly each time we estimate parameters

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

%xline( starting_K );

%Print some models 
for interest = 2:7
  print_model( interest, cs_ave, res_ave, index, labels, scales )
end






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