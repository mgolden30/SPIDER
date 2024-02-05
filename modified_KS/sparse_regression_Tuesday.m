%{
This script provides some methods for finding sparse relations from feature
matrices G
%}


%% Combinatoric Search

clear;

Ms = [128, 256, 512, 1024, 2048, 4096, 2^13, 2^14, 2^15, 2^16, 2^17, 2^18, 2^19];
t_per_subdomain = 0.8; %hold this constant

res_4  = zeros(numel(Ms),1);
res_as = zeros(numel(Ms),1);

for sweep_idx = 1:numel( Ms )
  M = Ms(sweep_idx)
%  M = 4096;
%load training data
load("seed_1_M_" + M + ".mat")
test = load("seed_543212345_M_" + M + ".mat");


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

interest = 6;
%xline(interest)
print_model( interest, cs, res, labels, scales )
legend({'training', 'testing', ''});
%return;

res_4(sweep_idx ) = res_ave(1,4);
res_as(sweep_idx ) = res_ave(1,end);

end

%%
clf();

tiledlayout(1,2);

nexttile

T = 2;
points_per_integral = round( t_per_subdomain/T*Ms );

scatter( points_per_integral, res_4, 100, 'filled' );
title("four term residual");
xlabel("$n_t$", 'Interpreter', 'latex' );
ylabel("residual", 'Interpreter', 'latex' );

set(gca, 'fontsize', 32 );
set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');



nexttile
scatter( points_per_integral, res_as, 100, 'filled' );
set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');

title("full library residual");
xlabel("$n_t$", 'Interpreter', 'latex' );
ylabel("residual", 'Interpreter', 'latex' );
set(gca, 'fontsize', 32 );

saveas(gcf, 'figs/4_term_and_full.png')

%%

Ms = [4096];
t_per_subdomain = 0.8; %hold this constant

res_4  = zeros(numel(Ms),1);
res_as = zeros(numel(Ms),1);

for sweep_idx = 1:numel( Ms )
  M = Ms(sweep_idx)
%  M = 4096;
%load training data
load("seed_1_M_" + M + ".mat")
test = load("seed_543212345_M_" + M + ".mat");


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
ylabel('residual');
xlabel('number of terms');

interest = 6;
%xline(interest)
print_model( interest, cs, res, labels, scales )
end

yticks([1e-5, 1])
set(gca, 'fontsize', 32 );
xline(6);
legend({'training', 'testing'});

saveas(gcf, 'figs/4096_opt_curve.png')


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