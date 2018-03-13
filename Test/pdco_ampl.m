% Benchmark variants of PDCO on a list of problems in AMPL format.
% D. Orban, June 2014.

import model.amplmodel;
import model.slackmodel;

% Lists of problems.
brunel_qps = {'q25fv47',  'qadlittl', 'qafiro',   'qbandm',   'qbeaconf', 'qbore3d',  ...
              'qbrandy',  'qcapri',   'qe226',    'qetamacr', 'qfffff80', 'qforplan', ...
              'qgfrdxpn', 'qgrow15',  'qgrow22',  'qgrow7',   'qisrael',  'qpilotno', ...
              'qrecipe',  'qsc205',   'qscagr25', 'qscagr7',  'qscfxm1',  'qscfxm2',  ...
              'qscfxm3',  'qscorpio', 'qscrs8',   'qscsd1',   'qscsd6',   'qscsd8',   ...
              'qsctap1',  'qsctap2',  'qsctap3',  'qseba',    'qshare1b', 'qshare2b', ...
              'qshell',   'qship04l', 'qship04s', 'qship08l', 'qship08s', 'qship12l', ...
              'qship12s', 'qsierra',  'qstair',   'qstandat'};

cute_qps = {'aug2d', 'aug2dc', 'aug2dcqp', 'aug2dqp', 'aug3d', 'aug3dc', 'aug3dcqp', ...
            'aug3dqp', 'cvxqp1_l', 'cvxqp1_m', 'cvxqp1_s', 'cvxqp2_l', 'cvxqp2_m', ...
            'cvxqp2_s', 'cvxqp3_l', 'cvxqp3_m', 'cvxqp3_s', 'dtoc3', 'dual1', 'dual2', ...
            'dual3', 'dual4', 'dualc1', 'dualc2', 'dualc5', 'dualc8', 'genhs28', ...
            'gouldqp2', 'gouldqp3', 'hs118', 'hs21', 'hs21mod', 'hs268', 'hs35', ...
            'hs35mod', 'hs51', 'hs52', 'hs53', 'hs76', 'hues-mod', 'huestis', 'ksip', ...
            'liswet1', 'liswet10', 'liswet11', 'liswet12', 'liswet2', 'liswet3', ...
            'liswet4', 'liswet5', 'liswet6', 'liswet7', 'liswet8', 'liswet9', ...
            'lotschd', 'mosarqp1', 'mosarqp2', 'powell20', 'primal1', 'primal2', ...
            'primal3', 'primal4', 'primalc1', 'primalc2', 'primalc5', 'primalc8', ...
            'qpcblend', 'qpcboei1', 'qpcboei2', 'qpcstair', 's268', 'stcqp1', ...
            'stcqp2', 'tame', 'ubh1', 'yao', 'zecevic2'};

mm_qps = {'cont-050', 'cont-100', 'cont-101', 'cont-200', 'cont-201', ...
          'dpklo1', 'exdata', 'laser', 'qptest', 'stadat1', 'stadat2', ...
          'stadat3', 'values'};
mm_large_qps = {'boyd1', 'boyd2', 'cont-300'};

% qp_set = brunel_qps;
%qp_set = cute_qps;
%qp_set = mm_qps;
%qp_set = mm_large_qps;
%qp_set = {'liswet1'};
qp_set = {'q25fv47'};

% PDCO options.
options_pdco = pdcoSet;
options_pdco.Print = true;
%options_pdco.Print = false;
options_solv.atol1 = 1.0e-6;  % Used as relative stopping tolerance by MINRES.
options_solv.atol2 = 1.0e-9; % Acts as atol_max.
options_solv.itnlim = 2;     % Proportion of the number of variables.

% Other options.
options_pdco.d1 = 1.0e-2;
options_pdco.d2 = 1.0e-2;

stats = zeros(length(qp_set), 3);  % We keep PDitns, CGitns, time.

% Print header.
fprintf('%15s', 'name');
fprintf('  %4s  %5s  %6s  %7s %3s', 'meth', 'outer', 'inner', 'time', 'err');
fprintf('\n');

for k = 1 : length(qp_set)
  qp_name = qp_set{k};
  fprintf('%15s', qp_name);

  % Convert problem so the constraints have the form: Ax = b, l ≤ x ≤ u.
  % Use sparse interface.
  qp = slackmodel(amplmodel(qp_name, true));

  % Setup problem data.
  zn = zeros(qp.n, 1);
  q0 = qp.fobj(zn); g0 = qp.gobj(zn); H  = qp.hobj(qp.x0);
  A  = qp.gcon(qp.x0); cL = qp.cL;
  bL = qp.bL; bU = qp.bU;

  % Set initial guess.
  x0 = zeros(qp.n, 1);
  x0(qp.jLow) = bL(qp.jLow) + 1;
  x0(qp.jTwo) = (bL(qp.jTwo) + bU(qp.jTwo)) / 2;
  x0(qp.jUpp) = bU(qp.jUpp) - 1;

  Anorm = normest(A, 1.0e-3);
  xsize = max(norm(x0, inf), 1);
  zsize = max(norm(qp.gobj(x0), inf) + sqrt(qp.n) * Anorm, 1);

  if options.Print
    fprintf('\nCalling PDCO with xsize = %7.1e and zsize = %7.1e\n', xsize, zsize);
  end
  z0 = zsize * ones(qp.n, 1);
  y0 = zeros(size(cL));

  options_pdco.mu0 = zsize;
  options_pdco.Maxiter = min(max(30, qp.n), 100);
  
  build_variant('K2,LDL');
  
  options_form = struct();

  solver = pdco_K2_LDL(qp, options_pdco,options_form,options_solv);
  solver.solve()

  PDitns = solver.PDitns;
  CGitns = solver.CGitns;
  time = solver.time;
  if solver.inform ~= 0
    PDitns = -PDitns;
    CGitns = -CGitns;
    time   = -time;
  end
  fprintf('  %5d  %6d  %7.2f', solver.PDitns, solver.CGitns, solver.time);
  if solver.inform == 1
    fprintf(' %3s', 'i');
  elseif solver.inform == 2
    fprintf(' %3s', 'l');
  elseif solver.inform == 3
    fprintf(' %3s', 's');
  else
    fprintf(' %3s', '');
  end
  stats(k, m, :) = [solver.PDitns, solver.CGitns, solver.time];
  fprintf('\n');
  qp.delete();
end
