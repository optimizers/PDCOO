function [outProblem] = mpstolp(inProblem)
  % Function "mpstolp" build the lpmodel equivalent to the problem
  % described in the structure input "inProblem" (for maximization problem,
  % the objective function will be turned in it opposite)
  %
  % Input args   :
  %   inProblem : a structure built by readmps.m, this structure should
  %               represent a linear optimization problem.
  %
  % Outputs args :
  %   outProblem : a lpmodel object

  % Extract constraints
  rowtypes = cell2mat(inProblem.rowtypes');

  A_E = find(rowtypes == 'E');
  A_L = find(rowtypes == 'L');
  A_G = find(rowtypes == 'G');
  A_N = find(rowtypes == 'N');

  % Compute dimensions
  [m, n] = size(inProblem.A);
  n_L = size(A_L, 1);
  n_G = size(A_G, 1);

  % Extract c
  c = inProblem.A(A_N(1), :)';
  if strcmp(inProblem.objsense, 'MAX') || strcmp(inProblem.objsense, 'MAXIMIZE')
    c = - c;
  end

  % Build A
  A = inProblem.A([A_E ; A_L ; A_G], :);

  % Build cL and cU
  if size(inProblem.rhs, 1) == 0
    inProblem.rhs = sparse(m, 1);
  end

  if size(inProblem.ranges, 1) ~= 0
    % Compute ranges
    ranges_E = inProblem.ranges(A_E);
    ranges_L = inProblem.ranges(A_L);
    ranges_G = inProblem.ranges(A_G);

    cL = [inProblem.rhs(A_E) - (ranges_E < 0) .* abs(ranges_E) ;    ...
          inProblem.rhs(A_L) - abs(ranges_L) ; inProblem.rhs(A_G)];
    cU = [inProblem.rhs(A_E) + (ranges_E > 0) .* abs(ranges_E) ;    ...
          inProblem.rhs(A_L) ; inProblem.rhs(A_G) + abs(ranges_G)];
  else
    cL = [inProblem.rhs(A_E) ; - Inf * ones(n_L, 1) ; inProblem.rhs(A_G)];
    cU = [inProblem.rhs(A_E) ; inProblem.rhs(A_L) ; Inf * ones(n_G, 1)];
  end

  % Build bL
  bL = sparse(n, 1);
  if size(inProblem.lbnds, 1) ~= 0
     bL = inProblem.lbnds';
  end

  % Build bu
  bU = Inf * ones(n, 1);
  if size(inProblem.ubnds, 1) ~= 0
     bU = inProblem.ubnds';
  end

  % Set initial guess.
  x0 = zeros(n, 1);

  outProblem = model.lpmodel(inProblem.name, x0, cL, cU, bL, bU, A, c);

end
