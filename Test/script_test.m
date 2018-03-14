import model.lpmodel;
import model.slackmodel;

options_pdco.file_id = 1;

formulation = 'K1x';
solver = 'Cholesky';
classname = build_variant(pdcoo_home, formulation, solver);

list_problem = {'afiro.mps'};

n_problem = length(list_problem);

options_pdco.d1 = 1.0e-2;
options_pdco.d2 = 1.0e-2;
options_pdco.OptTol = 1.0e-9;
options_solv.atol1 = 1.0e-4;
options_solv.atol2 = 1.0e-10;
options_solv.itnlim = 10;
options_pdco.Print = 1;

fprintf(options_pdco.file_id, ...
        '\n    Name    Objectif   Presid   Dresid   Cresid   PDitns   Inner     Time      D2 * r\n\n');

for i = 1:n_problem

  mps_name = list_problem{i};
  fprintf('%s\n', mps_name);

  % Read .mps file
  mps_stru = readmps(mps_name);
  lp = mpstolp(mps_stru);
  slack = slackmodel(lp);
  Anorm = normest(slack.gcon(slack.x0), 1.0e-3);

  options_pdco.x0 = slack.x0;
  options_pdco.x0(slack.jLow) = slack.bL(slack.jLow) + 1;
  options_pdco.x0(slack.jUpp) = slack.bU(slack.jUpp) - 1;
  options_pdco.x0(slack.jTwo) = (slack.bL(slack.jTwo) + slack.bU(slack.jTwo)) / 2;
  options_pdco.xsize = max(norm(options_pdco.x0, inf), 1);
  options_pdco.zsize = max(norm(slack.gobj(slack.x0), inf) + sqrt(slack.n) * Anorm, 1);
  options_pdco.z0 = options_pdco.zsize * ones(slack.n, 1);
  options_pdco.y0 = zeros(slack.m, 1);
  options_pdco.mu0 = options_pdco.zsize;
  options_pdco.Maxiter = min(max(30, slack.n), 100);
  options_form = struct();

  Problem = eval([classname, '(slack, options_pdco,options_form,options_solv)']);
  Problem.solve;
  fprintf(Problem.file_id, ...
          '\n%12s   %11.4e   %6.0f   %6.0f   %6.0f   %6d   %6d   %7.2f s   %11.4e\n', ...
          mps_name, slack.fobj(Problem.x),                                            ...
          log10(Problem.Pinf), log10(Problem.Dinf), log10(Problem.Cinf0),             ...
          Problem.PDitns, Problem.inner_total, Problem.time, options_pdco.d2^2 * norm(Problem.y));
end

fclose('all');
