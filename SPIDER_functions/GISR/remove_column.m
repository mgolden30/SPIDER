function s_min = remove_column( A, threshold  )
  %{
  PURPOSE:
  Given a tall matrix A, find the smallest singular values of A_k, where
  A_k is A without the kth column.

  INPUT:
  A: the matrix of interest
  threhsold: 
  %}

  [U,S_mat,~] = svd(A, 'econ');

  alpha = 1./vecnorm(A);

  w  = U'*A.*alpha; %normalized vectors

  al = 1./alpha.^2;

  S = diag(S_mat); %get the diagonal as a vector

  n = numel(S);

  filename = mfilename('fullpath');
  [filepath,~,~] = fileparts(filename);

  sfile = filepath + "/temp/S.bin";
  afile = filepath + "/temp/al.bin";
  wfile = filepath + "/temp/w.bin";
  ofile = filepath + "/temp/out2.bin";

  %write data to disk
  write_data( S,    sfile );
  write_data( al,   afile );
  write_data( w,    wfile );

  %call binary
  binary = filepath + "/remove_column";
  system(binary + " " + n + " " + sprintf(" %.16f ", threshold) + " " + sfile + " " + wfile + " " + afile + " " + ofile );

  %Read output
  s_min = read_data( ofile );
end