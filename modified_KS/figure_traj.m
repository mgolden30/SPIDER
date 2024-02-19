%{
Make a figure of modified KS trajectory
%}

clear;
load("trajectory.mat");

x = (0:N)/N*L;
t = (0:M)/M*(T/2);

imagesc(t,x,u(:,1:M/2));
pbaspect([2,1,1]);
%cb = colorbar();
%set( cb, "Location","north" );
clim([-1 1]*2);
xticks([0 100 200]);
yticks([0 22]);
colormap bone
set( gca, 'ydir', 'normal' );

%Draw a weak-form subdomain
hold on
  seed = 1;
  rng(seed);

  %[x,y,w,h]
  dx = x(2) - x(1);
  dt = t(2) - t(1);
  size_vec = [512*dt, 64*dx];
 
  num_boxes = 10;
  for i = 1:num_boxes
    center = [rand()*(T - size_vec(1)), rand()*(L - size_vec(2))];
    pos = [center, size_vec];
    rectangle( "position", pos, "EdgeColor", "red", "LineWidth", 2, "FaceColor", [1, 0, 0, 0.3]);
  end
hold off

fs = 12;
xlabel("$t$", "interpreter", "latex", "fontsize", fs);
ylabel("$x$", "interpreter", "latex", "fontsize", fs, "rotation", 0);
set(gca, "fontsize", fs);
set(gcf, "Color", "w");

saveas( gcf, "figs/traj.png" )