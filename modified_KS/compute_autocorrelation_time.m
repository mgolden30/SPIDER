%{
Get correlation time
%}



%grid variables
N = 256;
M = 2048;

% physical parameters
L = 2*pi;
T = 20;
nu = 0.1;

x = (0:(N-1))/N*2*pi;
x = x';

%define many initial conditions to train on
u0s = { 
        cos(x);
        exp( -(x-pi).^4 );
        sin(2*x) + exp(cos(x));
        sin(x) - 1./(1 + 0.2*cos(x));
        cos(x)./(1 + 0.2*sin(x));
        cos(2*x);
        sin(3*x-1);
        cos(x-1)./(2+cos(x));
        cos(sin(x));
        cos(2*sin(x));
      };

t = (0:(M-1))*T/(M-1);
[xx,tt] = meshgrid( t, x );

sp_forcing = 5./(1+cos(xx).^2);
sp_forcing = sp_forcing - mean(sp_forcing, 2);

force = sp_forcing.*(cos(tt*2*pi) + cos(tt*4*pi));

%imagesc(force);
%return;

state_index = 2;
u0 = u0s{state_index};
[u] = second_order_with_forcing( u0, M, N, L, T, nu, force );

imagesc((1:M)/M*T, x/(2*pi)*L, u); xlabel('t'); ylabel('x'); 
colormap jet
set(gca, 'ydir', 'normal');
colorbar();
%clim([-1 1]);
title("u(x,t)")
set(gca, 'fontsize', 32);
set(gcf, 'color', 'w');
pause(0.1);

with_square = true

if with_square
    hold on
      dx = 2*pi/N;
      dt = T/M;
      rectangle('position', [10.7, 3, 0.8, 1], 'linewidth', 3);
    hold off
    saveas(gcf, "figs/trajectory_with_domain.png");    
else
  saveas(gcf, "figs/trajectory.png");
end

grid = { x/(2*pi)*L, linspace(0,T,M) };


%% compute correlation in time and space now.

corr = @(x,y) sum(x.*y)./sqrt( sum(x.^2) * sum(y.^2) );
%corr = @(x,y) norm(x-y);

sub_mean = @(x) x - mean(x);
%sub_mean = @(x) x;


max_future = 256;
time_corr = zeros( max_future, 1);
for i = 1:M-max_future
  for j = 1:max_future
    u1 = sub_mean( u(:,i) );
    u2 = sub_mean( u(:,i+j) );
    
    time_corr(j) = time_corr(j) + corr( u1, u2 );
  end
end
time_corr = time_corr / (M - max_future);

dt = T/M;
plot( dt*(0:max_future-1), time_corr, 'linewidth', 3 );
ylabel('correlation');
xlabel('t');

yline(exp(-1));

title('$\langle \textrm{corr}[ u(s), u(s+t) ] \rangle_s$', 'Interpreter', 'latex');

xlim()
ylim([1e-1 1]);
yticks([0.1, exp(-1) 1]);
xline(0.8);
xlim([0 1]);
xticks([0, 0.8, 1]);
set(gca, 'fontsize', 32);

saveas(gcf, "figs/correlation_time.png");

%% Find correlation distance


max_future = N/2;
time_corr = zeros( max_future, 1);
for i = 1:N-max_future
  for j = 1:max_future
    u1 = sub_mean( u(:,i) );
    u2 = sub_mean( circshift(u(:,i), j ) );
    
    time_corr(j) = time_corr(j) + corr( u1, u2 );
  end
end
time_corr = time_corr / (N - max_future);

dt = 2*pi/N;
plot( dt*(0:max_future-1), time_corr, 'linewidth', 3 );
ylabel('correlation');
xlabel('x');

yline(exp(-1));

title('$\langle \textrm{corr}[ u(y), u(y+x) ] \rangle_y$', 'Interpreter', 'latex');

xlim()
ylim([1e-1 1]);
yticks([0.1, exp(-1) 1]);
x0 = 1
xline( x0 );
%xlim([0 1]);
xticks([0, x0, pi]);
set(gca, 'fontsize', 32);

saveas(gcf, "figs/correlation_space.png");