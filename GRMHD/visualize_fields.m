%{
visualize density
%}



[rho, v, B, t, x1, x2, x3] = load_MHD_data_stitch( );


%%
c = 2;
t = 16;
clf
plot(  squeeze(v(t,:,1,1,c) ))
hold on
plot(  squeeze(v(t,1,:,1,c) ))
plot(  squeeze(v(t,1,1,:,c) ))
legend

%%
for t = 1:40
  clf
  imagesc( squeeze(v(t,:,:,10,2)) );
  axis square
  colorbar
  %clim([0.5 1.5]);
  title(t + "")
  colormap jet
  drawnow
end