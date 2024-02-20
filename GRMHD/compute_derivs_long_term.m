%{
PURPOSE:
Compute derivatives of our data that we can load as we need.
%}

clear;

[rho, v, B, P, t, x1, x2, x3] = load_MHD_data_stitch();

fields = {"rho", "v", "B", "P"};
dx = x1(2) - x1(1);
dy = x2(2) - x2(1);
dz = x3(2) - x3(1);

%%
for i = 1:numel(fields)
  %Do x deriv
  f = fields{i};

  tic
  eval( sprintf("df = (circshift( %s, -1, 2 ) - circshift( %s, 1, 2 ))/dx;", f, f ) );
  toc
  save( "E:\Matts_Data\fields\" + f + "_dx", "df" );
  clearvars df;

  tic
  eval( sprintf("df = (circshift( %s, -1, 3 ) - circshift( %s, 1, 3 ))/dx;", f, f ) );
  toc
  save( "E:\Matts_Data\fields\" + f + "_dy", "df" );
  clearvars df;

  tic
  eval( sprintf("df = (circshift( %s, -1, 4 ) - circshift( %s, 1, 4 ))/dx;", f, f ) );
  toc
  save( "E:\Matts_Data\fields\" + f + "_dz", "df" );
  clearvars df;
end