function classname = build_variant(formulation, solver)

  classname = sprintf('pdco_%s_%s', formulation, solver);
  fname = ['Variants/' , classname, '.m'];
  if exist('Variants','dir')~= 7
    % Variants directory does not exist, make it.
    mkdir('Variants');
  end
  fid = fopen(fname, 'w');
  code = format_template(formulation, solver);
  fprintf(fid, code);
  fclose(fid);
end