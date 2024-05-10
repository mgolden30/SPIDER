%{
Visualize/clean data
%}

[eta0, x, y, t] = load_data();

sigma = 1;
eta = imgaussfilt(eta0, sigma);

%%
ny = size(eta,1);
nx = size(eta,2);
nt = size(eta,3);


vidObj = VideoWriter('peaks.avi');
    vidObj.FrameRate=10;
open(vidObj);


for t = 1:600

  %imagesc(eta(:,:,t));
  surf(eta(:,:,t));
  shading interp;

  colormap parula;

  %colorbar;

  pbaspect([nx*nt, ny*nt, 10*nx*ny])

  zlim([-8, 8]);
  xlim([0, 300]);
  ylim([0, 300]);

  drawnow

% Write each frame to the file.
       currFrame = getframe(gcf);
       writeVideo(vidObj,currFrame);
    end
  
    % Close the file.
    close(vidObj);