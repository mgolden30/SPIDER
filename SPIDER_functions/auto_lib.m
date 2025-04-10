function [G, labels, scales] = auto_lib( alphabet, fields, deriv_order, dimension, word_size, grid, corners, size_vec, pol )
  %{
  PURPOSE:
  Automatically generate a library and integrate it.
  This function generates AND integrates the library in this form.
  
  INPUT:
  alphabet - a cell {"u", "v", ..., "cos(u)", ...}. Does not contain
             derivatives. Derivatives will be added as this function runs.

  fields - a cell {u, v, ..., cos(u), ...} that contains the corresponding data matrices of alphabet. 
  
  OUTPUT:
  G - feature matrix
  %}

  num_fields = numel(fields);

  alphabet0 = alphabet;
  fields0 = fields;

  tic
  [alphabet, fields] = add_derivatives_to_alphabet( alphabet, fields, deriv_order, dimension, grid );
  walltime = toc;

  fprintf("Took %.3e seconds to precompute derivatives.\n", walltime);

  G = [];
  labels = {};
  scales = [];

  construction_time = [];
  integration_time  = [];

  A = 1; %incrmemental index over library
  for w = 1:word_size %word size
    %Generate all combinations with replacement
    P = combinator( numel(alphabet), w, 'c', 'r');

    %Loop over all combinations.
    %This is the loop over library terms
    for j = 1:size(P,1)

      %Construct our term letter by letter
      tic;
      term = 1; %start with 1
      label = "";

      %Decide while we construct this term if we can instead replace it
      %with a weak-form term
      can_replace = true; %innocent until proven guilty
      derivs = []; %by default we will do no integration by parts

      for k = 1:size(P,2)
        letter = alphabet{ P(j,k) };
      
        %Check if letters (besides the final letter) have derivatives
        if contains(letter, "_") && k < size(P,2)
          %Not canonically replacable by a weak-form term
          can_replace = false;
        end

        if k == size(P,2) && can_replace && contains(letter,"_")
          %we can replace this term with a weak-form term!
          c = char(letter); %turn our string into a character array

          %find the index of the base field
          for index = 1:num_fields
            if alphabet0{index} == c(1)
              break;
            end
          end

          %Use this base field INSTEAD
          term = term .* fields0{index};
          
          label = "(" + label + " " + c(1) + ")_";
          
          % This line is MATLAB spooky magic. 
          % Apparently '123' - '0' = [1,2,3]
          derivs = c(3:end) - '0';
          
          %append these derivatives to the label
          for d = derivs
            label = label + d;
          end

          %label = label + " " + letter;
        else
          %Proceed normally
          %We precomputed each term
          term  = term.*fields{P(j,k)};
          label = label + " " + letter;
        end
      end
      construction_time(A) = toc;

      tic;
      %Do the weak integral instead
      labels{A} = label;
      G(:,A) = SPIDER_integrate( term, derivs, grid, corners, size_vec, pol );
      scales(A) = 1;
      integration_time(A) = toc;
      A         = A+1;

      if mod(A, 128) ==0
        plot(integration_time, 'o');
        hold on
          plot(construction_time, 'd')
        hold off

        legend({'integration', 'construction'})
        ylim([0, 0.05]);
        drawnow;
      end
    end
  end
end


function [alphabet, fields] = add_derivatives_to_alphabet( alphabet0, fields0, deriv_order, dimension, grid )
  %{
  Helper function for auto_lib that adds derivatives to the library.

  INPUT:
  alphabet - {"u", "v", ....}
  fields - {u, v, ...}
  deriv_order - max number of derivatives to apply (inclusive)
  dimension - spacetime dimension of the data. How many different partial
              derivatives can we apply?
  grid - cell of grid info
  %}

  %To compute combinations of derivatives, we will use combinator
  addpath("combinator/");

  n = numel(alphabet0); %number of letters

  alphabet = {};
  fields = {};

  for i=1:n %loop over each letter to take derivatives
    alphabet{end+1} = alphabet0{i};
    fields{end+1} = fields0{i};
    
    for num_derivs = 1:deriv_order
      %Combinations with replacement
      P = combinator( dimension, num_derivs, 'c', 'r');
      
      %Add each combination of derivatives
      for j = 1:size(P,1)
        new_letter = alphabet0{i} + "_"; %underscore denotes derivatives
        data = fields0{i}; %data of this letter
        
        for k = 1:size(P,2)
          %add the derivative to the string
          new_letter = new_letter + P(j,k);

          %Take the derivative of the data
          direction = P(j,k); %axis we differentiate
          x  = grid{direction};
          dx = x(2) - x(1);

          %Do centered differencing for each derivative          
          data = (circshift( data, -1, direction ) - circshift( data, 1, direction ))/(2*dx);
        end %end derivative loop
    
        %append this differentiated field to the alphabet
        alphabet{end+1} = new_letter;
        fields{end+1}   = data;
      end %end combination
    end %end number of derivatives
  end %end field loop
end