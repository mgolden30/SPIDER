


N = 64;
M = 128;

L = 2*pi;
T = 4;

x = (0:(N-1))/N*2*pi;
x = x';

u0 = exp(-(x-pi).^4);
u0 = cos(x);

plot(u0);

nu = 0.1;

[u] = analytic_solution( u0, nu, L, T, M );

imagesc(u)
