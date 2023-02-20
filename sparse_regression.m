%{
This script provides some methods for finding sparse relations from feature
matrices G
%}


%% Combinatoric Search
addpath('SPIDER_functions/');

number_of_terms = 3;
residual_type = "|Gc|"; %Use this residual or the next
%residual_type = "|Gc|/max";
[models, residuals, cs] = combinatoric_search( G, labels, number_of_terms, scales, residual_type );

%Print to screen best 10 models
models(1:min(10, numel(models)))


