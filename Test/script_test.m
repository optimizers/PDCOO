import model.qpmodel;
import model.slackmodel;

options_pdco.file_id = 1;

formulation = 'K25';
solver = 'LDL';
classname = build_variant(pdcoo_home, formulation, solver);

% list_problem = {'afiro.mps'};
% list_problem = {'CVXQP1_S.QPS'};
% list_problem = {'CVXQP2_S.QPS'};
list_problem = {'CVXQP3_S.QPS'};

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
  lp = qpstoqp(mps_stru);
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

  Problem = eval([classname, '(slack, options_pdco, options_form, options_solv)']);
  Problem.solve;
  fprintf(Problem.file_id, ...
          '\n%12s   %11.4e   %6.0f   %6.0f   %6.0f   %6d   %6d   %7.2f s   %11.4e\n', ...
          mps_name, slack.fobj(Problem.x),                                            ...
          log10(Problem.Pinf), log10(Problem.Dinf), log10(Problem.Cinf0),             ...
          Problem.PDitns, Problem.inner_total, Problem.time, options_pdco.d2^2 * norm(Problem.y));

   figure;
   semilogy(Problem.conds, 'b-', 'LineWidth', 2); hold on
   semilogy(Problem.conds_pred, 'r-', 'LineWidth', 2);
   legend('condition number', 'estimate');
   title('K2.5: Condition number');
   xlabel('Iterations');
   figure;
   eigmin_neg = semilogy(-Problem.eigmin_neg, 'b.', 'MarkerSize', 15); hold on;
   eigmax_neg = semilogy(-Problem.eigmax_neg, 'r.', 'MarkerSize', 15);
   neg_lbnds = semilogy(-Problem.neg_lbnds, 'b-');
   neg_ubnds = semilogy(-Problem.neg_ubnds, 'r-');
   legend([eigmin_neg, eigmax_neg, neg_lbnds, neg_ubnds], ...
          'smallest negative', 'largest negative', 'lower bound', 'upper bound');
   title('K2.5: Negative eigenvalues');
   xlabel('Iterations');
   figure;
   eigmin_pos = semilogy(Problem.eigmin_pos, 'b.', 'MarkerSize', 15); hold on;
   eigmax_pos = semilogy(Problem.eigmax_pos, 'r.', 'MarkerSize', 15);
   pos_lbnds = semilogy(Problem.pos_lbnds, 'b-');
   pos_ubnds = semilogy(Problem.pos_ubnds, 'r-');
   legend([eigmin_pos, eigmax_pos, pos_lbnds, pos_ubnds], ...
          'smallest positive', 'largest positive', 'lower bound', 'upper bound');
   title('K2.5: Positive eigenvalues');
   xlabel('Iterations');
end

fclose('all');
