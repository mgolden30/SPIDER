function s_min = add_column( A, a, threshold )   
  %A is current matrix
  %a are the possible columns to add

  [U,S,V] = svd(A);

  %compute information about the column modification
  alpha = 1./vecnorm(a);  %compute norm
  w     = U' * a .* alpha; %rotate and normalize

  wbar  = w(1:size(A,2), :);

  dim = 1;
  tau = 1 - sum(wbar.^2, dim);

  al = 1./alpha.^2;
  
  s = diag(S);
   
  n = numel(s);
  p = size(wbar,2);


  filename = mfilename('fullpath');
  [filepath,~,~] = fileparts(filename);

  sfile = filepath + "/temp/S.bin";
  afile = filepath + "/temp/al.bin";
  wfile = filepath + "/temp/w.bin";
  tfile = filepath + "/temp/tau.bin";
  ofile = filepath + "/temp/out.bin";

  %write data to disk
  write_data( s,    sfile );
  write_data( al,   afile );
  write_data( wbar, wfile );
  write_data( tau,  tfile );

  %call binary
  binary = filepath + "/add_column";
  system(binary + " " + n + " " + p + sprintf(" %.16e ", threshold) + " " + sfile + " " + wfile + " " + afile + " " + tfile + " " + ofile );

  %Read output
  s_min = read_data( ofile );
end