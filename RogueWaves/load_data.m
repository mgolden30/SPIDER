function [eta, x, y, t] = load_data()

% ncdisp('xygrid_50cm_20200211_1040_plane_sub.nc')

% Dimensions:
%            time = 5679
%            x    = 280
%            y    = 300
% Variables:
%     time 
%            Size:       5679x1
%            Dimensions: time
%            Datatype:   double
%            Attributes:
%                        units = 'microseconds since 2020-02-11 10:40:00:383874'
%     x    
%            Size:       280x1
%            Dimensions: x
%            Datatype:   single
%            Attributes:
%                        units = 'm'
%     y    
%            Size:       300x1
%            Dimensions: y
%            Datatype:   single
%            Attributes:
%                        units = 'm'
%     xgrid
%            Size:       300x280
%            Dimensions: y,x
%            Datatype:   single
%     ygrid
%            Size:       300x280
%            Dimensions: y,x
%            Datatype:   single
%     eta  
%            Size:       300x280x5679
%            Dimensions: y,x,time
%            Datatype:   single

fsample=5 ; % 5 Hz 

%Xg=ncread('xygrid_50cm_20200211_1040_plane_sub.nc','xgrid');
%Yg=ncread('xygrid_50cm_20200211_1040_plane_sub.nc','ygrid');
%Zg=ncread('xygrid_50cm_20200211_1040_plane_sub.nc','eta');



x =ncread('xygrid_50cm_20231124_1300_plane_sub.nc','xgrid');
y =ncread('xygrid_50cm_20231124_1300_plane_sub.nc','ygrid');
eta=ncread('xygrid_50cm_20231124_1300_plane_sub.nc','eta');

t = 0:size(eta,3)-1/fsample;

sigma = 1;
eta  =imgaussfilt3( eta, sigma);
%{
for j=1:size(Zg, 3)
   eta2=squeeze( eta(:,:,j)); 
   imagesc(eta2)


  xlabel('X [meters]')
  ylabel('Y [meters]')
  title(['Hmax=',num2str(max(eta(:))),' meters'])
  colorbar()
  clim( [-2, 6]);
  pause(0.05)

end
%}



end