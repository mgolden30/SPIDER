function w = read_data( filename )
  f = fopen(filename, "r");

  N = fread(f, 1, "int");

  w = fread(f, N, "double");

  fclose(f);
end