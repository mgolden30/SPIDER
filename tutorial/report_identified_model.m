function k = report_identified_model(cs, residuals, scales, labels, gamma)
    
    figure(); %make a new figure
    scatter( 1:numel(residuals), residuals, 'o', 'filled');
    set(gca, 'yscale', 'log');
    xlabel("# of terms");
    ylabel('residual');
    
    %Select number of terms
    jumps = find(residuals(1:end-1)./residuals(2:end) > gamma);
    k = jumps(end) + 1; %identified model
    
    %Highlight the identified model
    hold on 
    scatter( k, residuals(k), 200, 's', 'MarkerEdgeColor', 'red', 'LineWidth', 2 );
    hold off
    
    %Print identified model
    c = cs(:,k)./scales.'; %rescale to physical units
    c = c(c~= 0);
    c = c/max(abs(c));
    labels(cs(:,k) ~= 0);

    I = find(cs(:,k));
    for i = 1:k
      fprintf("%.5f " + labels{I(i)}, c(i));
      if i == k
        continue;
      end 
      if c(i+1) >0
        fprintf(" +");
      end
    end
    fprintf(" = 0\n");
end