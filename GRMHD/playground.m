t = 1;
filename = sprintf("sim1/HGB.out2.%05d.athdf", t);
filename = sprintf("E:/sim3/turb.out2.%05d.athdf", 200);

hf = h5info( filename );

hf.Datasets(1)
x1f = h5read( filename, '/x1f' );
x1v = h5read( filename, '/x1v' );

%read in the submesh structure
%log_loc = h5read( filename, '/LogicalLocations' );

hf.Attributes(13)
%h5readatt( filename, '/', 'RootGridSize' )
h5readatt( filename, '/', 'VariableNames' )