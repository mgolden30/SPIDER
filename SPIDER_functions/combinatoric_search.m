function [models, residuals, cs] = combinatoric_search( G, labels, desired_terms, scales, residual_type )
  %{
  PURPOSE:
  Do a combinatoric search for relations 

  INPUT:
  G      - integrated matrix
  labels - cell of strings labeling each term
  desired_terms - integer > 1 
  scales - scale of terms 
  residual_type - "|Gc|" or "|Gc|/max"
  %}

  [~, nl] = size(G); %num windows and num library
  big = nchoosek(nl, desired_terms); %Combinations to check

  residual_function = 0;
  if residual_type == "|Gc|"
    residual_function = @(G,c) norm(G*c);
  end
  if residual_type == "|Gc|/max"
    residual_function = @(G,c) norm(G*c)/max(vecnorm(G*diag(c)));
  end

  cs        = zeros(big, nl);
  residuals = zeros(big, 1 );
  models    = cell( big, 1 );
  combinations = nchoosek( 1:nl, desired_terms );
  
  for i = 1:big
    combination = combinations(i,:);
    G_restricted = G(:, combination);
    
    %compute smallest singular value
    [~, ~, V] = svd(G_restricted, 'econ');
    c = V(:,end);
    
    cs(i,combination) = c;
    residuals(i)      = residual_function(G_restricted,c);
    model_str = "residual = " + residuals(i) + ", ";
    for j=1:desired_terms
      %rescale coefficients!
      scales_restricted = scales( combination );
      c_true = c./scales_restricted';
      c_true = c_true/max(abs(c_true)); %normalize by max coefficient
      model_str =  model_str + c_true(j) + " " + labels( combination(j) );
      if j < desired_terms
        model_str = model_str + "  +  ";
      else
        model_str = model_str + "  = 0";
      end
    end
    models{i} = model_str;
  
    %If any term in the model is negligible, artificially inflate eta so it
    %gets thrown out.
    %mags = vecnorm(G_restricted*diag(c));
    %if( min(mags)/max(mags) < 1e-2 ) %Check if there is a negligible term
    %  etas(i) = 1; 
    %end
  end
  
  %bad = (etas < 0.001);     %Throw out identities
  %bad = 0*bad;
  %bad = bad | (etas > 0.5); %Throw out if eta is too high. Nonsense relations
  
  %models( bad ) = [];
  %etas( bad )   = [];
  %cs( bad, : )  = [];
  
  [residuals, I] = sort( residuals );
  cs = cs(I, :);
  models = models(I);
end