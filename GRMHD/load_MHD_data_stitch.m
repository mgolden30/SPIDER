function [rho, v, B, P, t, x1, x2, x3] = load_MHD_data_stitch( )
  %{
  PURPOSE:
  Load in fields for GRMHD, as well as numerical grids
   
  INPUT:

  OUTPUT:
  %}
  
  % B is a 3 component field and prim is 4 component
  % I have to assume prim is density and velocity

  %filename = sprintf("sim1/HGB.out2.%05d.athdf", 0);
  %filename = sprintf("C:/Users/wowne/Downloads/sim4/HGB.out2.%05d.athdf", 0);
  %filename = sprintf("E:/sim3/turb.out2.%05d.athdf", 200);
  %filename = sprintf("C:/Users/wowne/Downloads/sim11/Turb.out2.%05d.athdf", 0);
  filename = sprintf("C:/Users/wowne/Downloads/sim13/Turb.out2.%05d.athdf", 0);

  log_loc = h5read( filename, '/LogicalLocations' );
  nm = max( log_loc, [], 2 ) + 1; %number of meshes in each direction
  
  gs = h5readatt( filename, '/', 'RootGridSize' ); %true grid size


  timesteps = 300; %number of timesteps I want to read

  B    = zeros(timesteps, gs(1), gs(2), gs(3), 3);
  prim = zeros(timesteps, gs(2), gs(2), gs(3), 5); %4 or 5 depending on if pressure is needed
  
  gs
  nm
  n = int64(gs)./nm %grid sizes of each mesh

  t = zeros(1, timesteps); %we can read the time directly from each data file
  
  for slice = 0:timesteps-1
    
    filename = sprintf("C:/Users/wowne/Downloads/sim13/Turb.out2.%05d.athdf", slice + 0);
    %filename = sprintf("C:/Users/wowne/Downloads/sim4/HGB.out2.%05d.athdf", t + 500 );
    %filename = sprintf("sim1/HGB.out2.%05d.athdf", t);
    %filename = sprintf("E:/sim3/turb.out2.%05d.athdf", 200 + t);

    %read out this file's time
    t(slice+1) = h5readatt( filename, '/', 'Time' );

    B_temp    = h5read( filename, '/B' );
    prim_temp = h5read( filename, '/prim' );
    
    for i = 1:prod(nm)
      i1 = n(1)*log_loc(1,i) + (1:n(1));
      i2 = n(2)*log_loc(2,i) + (1:n(2));
      i3 = n(3)*log_loc(3,i) + (1:n(3));

      %go to t+1 since we start with zero
      prim( slice+1, i1, i2, i3, :) = prim_temp(:,:,:,i,:);
      B(    slice+1, i1, i2, i3, :) = B_temp(:,:,:,i,:);
    end
  end

  rho = squeeze( prim(:,:,:,:,1) );
  P   = squeeze( prim(:,:,:,:,2) );
  v   = squeeze( prim(:,:,:,:,3:5) );
  
  %now figure out the spatial grid

  L1 = h5readatt( filename, '/', 'RootGridX1' );
  L2 = h5readatt( filename, '/', 'RootGridX2' );
  L3 = h5readatt( filename, '/', 'RootGridX3' );

  f = @(D,N) (0:N-1)*D/N;
  x1 = f( L1(3), size(B,2) );
  x2 = f( L2(3), size(B,3) );
  x3 = f( L3(3), size(B,4) );
end