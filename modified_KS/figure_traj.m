%{
Make a figure of modified KS trajectory
%}

clear;
load("trajectory.mat");

x = (0:N)/N*L;
t = (0:M)/M*(T);

imagesc(t,x,u(:,1:M));
pbaspect([4,1,1]);
%cb = colorbar();
%set( cb, "Location","north" );
clim([-1 1]*2);
xticks([0 100 200]);
yticks([0 22]);
colormap bone
set( gca, 'ydir', 'normal' );

%Draw a weak-form subdomain
hold on
  seed = 11;
  rng(seed);

  %[x,y,w,h]
  dx = x(2) - x(1);
  dt = t(2) - t(1);
  size_vec = [512*dt, 64*dx];
 
  num_boxes = 5;
  for i = 1:num_boxes
    center = [rand()*(T - size_vec(1)), rand()*(L - size_vec(2))];
    pos = [center, size_vec];
    rectangle( "position", pos, "EdgeColor", "red", "LineWidth", 1, "FaceColor", [1, 0, 0, 0.3]);
  end
hold off

%Now decide how big I want this figure in inches
figure = gcf;
im_size = round([500*2,300]*0.6);
figure.Position = [400,400,im_size];
set(gcf, "color", "white");

fs = 12;
xlabel("$t$", "interpreter", "latex", "fontsize", fs);
ylabel("$x$", "interpreter", "latex", "fontsize", fs, "rotation", 0);
set(gca, "fontsize", fs);
set(gcf, "Color", "w");

exportgraphics( gcf, "figs/traj.pdf" )