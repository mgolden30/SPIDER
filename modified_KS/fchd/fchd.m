function w = fchd(v)
%FCHD	Fast Chebyshev derivative
% FCHD(V) computes the first derivative of the data in V located along
% the N+1 Chebyshev points cos(pi*(0:N)/N).
% 
% 
% Example 1:
% Use FCHT to differentiate the function f(x) = tan(x) over [-1,1], and
% compare with the exact derivate f'(x) = sec(x)^2.
% 
% x = cos(pi*(0:10)/10); % create sparse Chebyshev-spaced grid of 11 points
% xx = linspace(-1,1); % create dense, linearly spaced grid
% plot(xx,sec(xx).^2,x,fchd(tan(x))); % compare Chebyshev derivative to exact
%
% 
% Example 2:
% To show the spectral convergence property of the Chebyshev derivative,
% compute the error between the Chebyshev derivative and the exact
% derivative of f(x) = tan(x) for several N.
% 
% N = 1:30;
% err = zeros(1,length(N));
% 
% for n = N
%     x = cos(pi*(0:n)/n)'; % establish grid
%     err(n) = max(sec(x).^2 - fchd(tan(x))); % compute error
% end
% 
% loglog(N,err); %display

v = v(:)';
M = length(v);
N = M - 1;

% fast Chebyshev transforms
b = real(fft([v fliplr(v(2:N))]));
c = real(ifft(1i*[0:M-2 0 2-M:-1].*b));
n2b = (0:M-2).^2.*b(1:N);

% derivative construction
w = zeros(M,1);
w(1) = sum(n2b)/N + 0.5*N*b(M);
w(2:N) = -csc(pi/N*(1:M-2)).*c(2:N);
w(M) = sum((-1).^(1:N).*n2b)/N + 0.5*(-1)^M*N*b(M);

end