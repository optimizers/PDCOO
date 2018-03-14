function classname = build_variant(pdcoo_home, formulation, solver)

  classname = sprintf('pdco_%s_%s', formulation, solver);
  variants_dir = fullfile(pdcoo_home, 'Variants')
  fname = fullfile(variants_dir, sprintf('%s.m', classname));
  if exist(variants_dir, 'dir') ~= 7
    mkdir(variants_dir);
  end
  fid = fopen(fname, 'w');
  code = format_template(formulation, solver);
  fprintf(fid, code);
  fclose(fid);
end
