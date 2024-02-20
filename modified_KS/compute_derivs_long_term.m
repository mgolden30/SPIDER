%{
PURPOSE:
Compute derivatives of our data that we can load as we need.
%}

clear;

[rho, v, B, P, t, x1, x2, x3] = load_MHD_data_stitch();
