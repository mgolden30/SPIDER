%{
This script provides some methods for finding sparse relations from feature
matrices G
%}


%% Combinatoric Search

clear;

%load training data
load('seed2_1.mat')
test = load('seed2_543212345.mat');


addpath('../SPIDER_functions/');
addpath('../SPIDER_functions/GISR/');

nl = size(G,2);

c0 = ones(size(G,2),1); %start with full library
%[res, cs] = GISR(G, c0);


%Alternatively, use the old regression code
number_of_subsamples = 1; %Do parameter estimation 100 times to estimate uncertainty
fraction_of_windows_in_subsample = 1; %use half of windows randomly each time we estimate parameters

%We can pick any of these residuals to minimize
index = 1;
res1 = @(G,c) norm(G*c)/norm(c);
res2 = @(G,c) norm(G*c)/norm(c)/norm(G);
res3 = @(G,c) norm(G*c)/max( vecnorm(G*diag(c)) );
residual_labels = { "$|Gc|$", "$|Gc|/|G_{c\neq 0}|$", "$|Gc|/\textrm{max}\{|G c_n|\}$" };
residual_functions = {res1, res2, res3}; %residuals to compute during sparsification
res_func_for_sparsification = residual_functions{index};
[res_ave, res_std, cs, cs_std] = greedy_regression( G, c0, ...
                                                                   number_of_subsamples, ...
                                                                   fraction_of_windows_in_subsample, ...
                                                                   residual_functions, ...
                                                                   res_func_for_sparsification );

res = res_ave(1,:);


scatter( 1:nl, res, 's', 'markeredgecolor', 'black' )

hold on

rescale = scales./(test.scales);

res_test = vecnorm( test.G * (cs) );
scatter( 1:nl, res_test, 's', 'markeredgecolor', 'blue' )
hold off

set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');
ylabel('absolute residual');
xlabel('number of terms');
%{
interest = 6;
I = (cs(:,interest) ~= 0);
save('reference_indices.mat', "I");
%}
%xline(interest)
print_model( interest, cs, res, labels, scales )
legend({'training', 'testing', ''});

%% See scaling with Ht

cs = [];
for i =32:8:128
  load( 'reference_indices.mat' );
  load( "Ht_sweep/seed2_1_ht_" + i + ".mat" );

  G = G(:,I);
  labels = labels(I);

  [~,~,V] = svd(G);

  c = V(:,end);

  [~,idx] = max(abs(c));
  c = c/c(idx);

  cs = [cs c];
end

imagesc(cs)

tiledlayout(3,2);
for i = 1:6
  nexttile;
  plot( 32:8:128, cs(i,:));
  title( labels{i} );
  xlabel("H_t");

  
  m = mean(cs(i,:));
  ylim( sort([0.8, 1.2]*m) );
  drawnow;
end



function str = print_model( interest, cs_ave, res, labels, scales )
  c = cs_ave(:,interest);
  c = c./scales';

  normalization = max(abs(c)); %normalize by the maximal term
  
  str = interest + ":  residual = " + res(interest) + ", \quad   ";
  for i = 1:numel(c)
    if( c(i) == 0 )
      continue; 
    end
  
    str = str + " + " + c(i)/normalization + "  " + labels{i};
  end
  str = str + " = 0";
end