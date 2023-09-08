function tau = correlation_length( rho )
  T = size(rho,1);

  buffer = 4; %timesteps used in estimating correlation time

  %subtract the mean over space
  rho = rho - mean(rho, [2,3,4]);

  corr  = @(x,y) sum(x.*y, 'all')/sqrt( sum(x.*x, 'all') * sum(y.*y, 'all') );
  tau   = zeros(T-buffer, 1);
  corrs = zeros(buffer,1);

  for t = 1:T-buffer
    for i = 1:buffer
      corrs(i) = corr( rho(t,:,:,:), rho(t+i,:,:,:) );
    end

    %plot(corrs);
    %pause(1)

    %fit corr to a decyaing exponential
    P = polyfit( 1:buffer, log(corrs), 1 )
    tau(t) = P(1);  
  end


end