%  
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
clear all

fsample=5 ; % 5 Hz 

  %=1 storm Nov 24 2023 
  %=2 storm Jan 4 2020
  %=3 storm Feb 11 2020
  %=4 storm  Nov 19 2020
  %=5 storm April 6 2021
  %=6 storm Aug 17 2023  ** VERY NOISY DATA
  %=7 storm March 28 2020

%1 3 4 5

  STORM=5

  % theta ; angle between the Y axis and the mean wave direction 

if STORM==1
filename='xygrid_50cm_20231124_1300_plane_sub.nc';
theta=0; 
end

if STORM==2
filename='xygrid_50cm_20200104_1040_plane_sub.nc';
theta=0; 
end

if STORM==3
filename='xygrid_50cm_20200211_1040_plane_sub.nc';
theta=pi/2; 
end

if STORM==4
    filename='xygrid_50cm_20201119_1100_plane_sub.nc';
    theta=0; 
end

if STORM==5
    filename='xygrid_50cm_20210406_1700_plane_sub.nc';
    theta=0; 
end

if STORM==6
    filename='xygrid_50cm_20230817_1000_plane_sub.nc';
    theta=0; 
end


if STORM==7
    filename='xygrid_50cm_20200328_1620_plane_sub.nc';
    theta=-pi/10; 
end



Xg=ncread(filename,'xgrid');
Yg=ncread(filename,'ygrid');
Zg=ncread(filename,'eta');
Zg(isnan(Zg))=0; 

Zgn=Zg; % full data set on the trapeziodal domain



[MY MX MT]=size(Zg);


%% rotate axis (counterclock rotation from X axis) 



Xr=Xg*cos(theta) - Yg*sin(theta);
Yr=Xg*sin(theta) + Yg*cos(theta);

figure(1)
for j=1:200
   eta=squeeze( Zg(:,:,j)); 

   subplot(1,2,1)
   pcolor(Xg,Yg,eta)
   shading flat

xlabel('X [meters]')
ylabel('Y [meters]')
title(['Hmax=',num2str(max(eta(:))),' meters'])


 subplot(1,2,2)
   pcolor(Xr,Yr,eta)
   shading flat

xlabel('X [meters]')
ylabel('Y [meters]')
title(['Hmax=',num2str(max(eta(:))),' meters'])

pause(0.05)
end


%%  data in a rectangular domain

% fat rectangle 
i1=2:151; 
i2=50:211; 

% elongated rectangle 
i1=2:201; 
i2=66:197; 
Zgs=Zg(i1,i2,:); % data set on a rectangular domain

Hs=4*std(Zgs,0,3);  % average wave height 

figure
pcolor(Xg,Yg, 4*std(Zg,0,3))
shading flat
hold on
plot(Xg(i1,i2),Yg(i1,i2),'.')
shading flat
colorbar
title('average wave height ')

%%
for j=1:MT
   eta=squeeze( Zg(:,:,j)); 

   pcolor(Xg,Yg,eta)
   shading flat



xlabel('X [meters]')
ylabel('Y [meters]')
title(['Hmax=',num2str(max(eta(:))),' meters'])
pause(0.05)
end
%%