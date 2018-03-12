function classname = build_variant(formulation, solver)

  classname = sprintf('pdco_%s_%s', formulation, solver);
  fname = ['Variants/' , classname, '.m'];
  fid = fopen(fname, 'w');
  code = format_template(formulation, solver);
  fprintf(fid, code);
  fclose(fid);
end
