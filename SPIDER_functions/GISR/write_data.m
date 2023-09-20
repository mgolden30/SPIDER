
function write_data(data, filename)
  N    = numel(data);

  f = fopen(filename, "w");
  
  fwrite(f, N,    "int");
  fwrite(f, data, "double");
  
  fclose(f);
end