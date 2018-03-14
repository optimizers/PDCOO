classdef pdcoO < handle
    %-----------------------------------------------------------------------
    % pdcoO.m: Primal-Dual Barrier Method for Convex Objectives
    %-----------------------------------------------------------------------
    % This abstract class deals with solving optimization problems of the form
    %
    %    minimize    phi(x) + 1/2 norm(D1*x)^2 + 1/2 norm(r)^2
    %      x, r
    %    subject to  A*x + D2*r = b,   bl <= x <= bu,   r unconstrained,
    %
    % where
    %    phi(x) is a smooth convex function  defined by function pdObj;
    %    A      is an m x n matrix defined by opMatrix or opFunction pdMat;
    %    b      is a given m-vector;
    %    D1, D2 are positive-definite diagonal matrices defined from d1, d2.
    %           In particular, d2 indicates the accuracy required for
    %           satisfying each row of Ax = b.
    %
    % D1 and D2 (via d1 and d2) provide primal and dual regularization
    % respectively.  They ensure that the primal and dual solutions
    % (x, r) and (y, z) are unique and bounded.
    %
    % A scalar d1 is equivalent to d1 = ones(n, 1), D1 = diag(d1).
    % A scalar d2 is equivalent to d2 = ones(m, 1), D2 = diag(d2).
    % Typically, d1 = d2 = 1.0e-4.
    % These values perturb phi(x) only slightly  (by about 1.0e-8) and request
    % that A*x = b be satisfied quite accurately (to about 1.0e-8).
    % Set d1 = 1.0e-4, d2 = 1 for least-squares problems with bound constraints.
    % The problem is then equivalent to
    %
    %    minimize    phi(x) + 1/2 norm(d1*x)^2 + 1/2 norm(A*x - b)^2
    %    subject to  bl <= x <= bu.
    %
    % More generally, d1 and d2 may be n and m vectors containing any positive
    % values (preferably not too small, and typically no larger than 1).
    % Bigger elements of d1 and d2 improve the stability of the solver.
    %
    % At an optimal solution, if x(j) is on its lower or upper bound,
    % the corresponding z(j) is positive or negative respectively.
    % If x(j) is between its bounds, z(j) = 0.
    % If bl(j) = bu(j), x(j) is fixed at that value and z(j) may have
    % either sign.
    %
    % Also, r and y satisfy r = D2 y, so that Ax + D2^2 y = b.
    % Thus if d2(i) = 1.0e-4, the i-th row of Ax = b will be satisfied to
    % approximately 1.0e-8.  This determines how large d2(i) can safely be.

    properties

        slack       % slackmodel or slackmodel_spot object

        % Input extracted from "slack"
        A           % Matrix or opSpot
        b           % vector of RHS
        bl          % vector of lowerbounds
        bu          % vector of upper bounds
        n           % number of variables
        m           % number of constraints

        % Input provided by the user via "options"
        d1          % positive scalar or positive vector (see above)
        d2          % positive scalar or positive vector (see above)
        x0          % initial primal variables x
        y0          % initial dual   variables y
        z0          % initial dual   variables z
        xsize       % estimates of the largest x at the solution (for scaling)
        zsize       % estimates of the largest z at the solution (for scaling)
        maxitn      % max iterations
        featol      % required primal feasibility accuracy
        opttol      % required dual feasibility and complementarity accuracy
        Prilev      % 1 for verbose, else 0
        steptol
        stepSame    % true if stepx and stepz should be the same
        x0min       % min distance between x0 and bl or bu AFTER SCALING
        z0min       % min distance between z0 and zero AFTER SCALING
        mu0         % initial barrier parameter
        backtrack   % specifies whether to do a backtracking linesearch
        file_id

        % Output
        x           % primal solution
        y           % dual solution associated with Ax + D2 r = b
        z           % dual solution associated with bl <= x <= bu
        inform      % final solver status
                    % inform = 0 if a solution is found
                    % = 1 if too many iterations were required
                    % = 2 if the linesearch failed too often
                    % = 3 if the step lengths became too small
                    % = 4 if other errors occured.
        PDitns      % number of primal-dual barrier iterations required
        inner_itns  % number of inner iterations required in a given outer iteration
        inner_total % total number of inner iterations required
        time        % the cpu time used (via cputime)
                    % (we also use tic/toc to allow for multicore systems.)

        explicitA
        zn

        % Constraints type
        low
        upp
        one
        fix
        two
        zlo
        zup
        nfix

        % Used to establish convergence
        eta
        maxf
        maxfail
        bigcenter
        gamma
        delta
        beta
        zeta
        theta
        rL
        rU
        cL
        cU
        Pinf
        Dinf
        Cinf
        Cinf0
        fmerit

        % Used to compute mu
        mufirst
        mulast
        mu
        center

        % Variables
        x1
        x2
        z1
        z2
        r1
        r2
        dx1
        dx2
        dz1
        dz2
        dx
        dy

        obj
        grad
        hess
        H

        step
        nf

        % Some variables to print
        nfail
        regterm
        objreg
        objtrue

        % They have to be defined in subclasses
        manage_op % Method deals with opFunction ?
        diagHess  % is hess a nx1 vector ?
        solver    % Solver's name
        head3     % Solver's name to print
        need_precon

        % Misc
        neg_lbnds
        eigmin_neg
        eigmax_neg
        neg_ubnds
        pos_lbnds
        eigmin_pos
        eigmax_pos
        pos_ubnds
        conds

    end

    methods (Abstract) % Should be defined in subclasses

        Print_param(o)      % Print parameters which are specific to the subclass
        Init_param(o)       % Initialize some parameters for the subclass
        Solve_Newton(o)     % Solve Newton system
        Print_results(o)    % Print some results specific to the subclass
        Reset_param(o)      % Reset some parameters specific to the subclass at the end of the loop

    end

    methods

        % Builder
        function o = pdcoO(slack, options)
            % PDCO contructor :
            % Inputs :
            %   slack         : a slackmodel or slackmodel_spot object
            %
            %   options is a structure that may contains following
            %   fields and other regarding to it subclass
            %       d1        : vector or scalar for building D1
            %                   ifndef d1 = 1.0e-4
            %       d2        : vector or scalar for building D2
            %                   ifndef d2 = 1.0e-4
            %       x0        : initial point for x
            %                   ifndef x0 = 1
            %       y0        : initial point for y
            %                   ifndef y0 = 1
            %       z0        : initial point for z
            %                   ifndef z0 = 1
            %       xsize     : scalar to scale x
            %                   ifndef xsize = 100
            %       zsize     : scalar to scale z
            %                   ifndef zsize = 100
            %       Maxiter   : scalar of max iterations in pdco
            %                   ifndef Maxiter = 30
            %       FeaTol    : scalar of accuracy in equality constraints
            %                   ifndef FeaTol = 1.0e-6
            %       OptTol    : scalar of accuracy in inequality constraints
            %                   ifndef OptTol = 1.0e-6
            %       Print     : 1 for give output, else 0
            %                   ifndef Print = 1
            %       StepTol   : scalar which control how x and z reach the
            %                   bounds
            %                   ifndef StepTol = 0.99
            %       StepSame  : 1 if stepx and stepz should be the same
            %                   ifndef StepSame = 1
            %       x0min     : Min distance between x0 and bl or bu  AFTER SCALING
            %                   ifndef x0min = 1
            %       z0min     : Min distance between abs(z0) and zero AFTER SCALING
            %                   ifndef z0min = 1
            %       mu0       : scalar of Initial mu
            %                   ifndef mu0 = 1.0e-1
            %       backtrack : specifies whether to do a backtracking linesearch
            %                   ifndef backtrack = 0

            o.slack = slack;

            o.A = slack.gcon(slack.x0);
            o.b = slack.cL;
            o.bl = slack.bL;
            o.bu = slack.bU;
            o.n = slack.n;
            o.m = slack.m;
            if isfield(options, 'd1')
                o.d1 = options.d1;
            else
                o.d1 = 1.0e-4;
            end

            if isfield(options, 'd2')
                o.d2 = options.d2;
            else
                o.d2 = 1.0e-4;
            end

            if isfield(options, 'x0')
                o.x0 = options.x0;
            else
                o.x0 = slack.x0;
            end

            if isfield(options, 'y0')
                o.y0 = options.y0;
            else
                o.y0 = ones(o.m, 1);
            end

            if isfield(options, 'z0')
                o.z0 = options.z0;
            else
                o.z0 = ones(o.n, 1);
            end

            if isfield(options, 'xsize')
                o.xsize = options.xsize;
            else
                o.xsize = 100;
            end

            if isfield(options, 'zsize')
                o.zsize = options.zsize;
            else
                o.zsize = 100;
            end

            if isfield(options, 'Maxiter')
                o.maxitn = options.Maxiter;
            else
                o.maxitn = 30;
            end

            if isfield(options, 'FeaTol')
                o.featol = options.FeaTol;
            else
                o.featol = 1.0e-6;
            end

            if isfield(options, 'OptTol')
                o.opttol = options.OptTol;
            else
                o.opttol = 1.0e-6;
            end

            if isfield(options, 'Print')
                o.Prilev = options.Print;
            else
                o.Prilev = 1;
            end

            if isfield(options, 'StepTol')
                o.steptol = options.StepTol;
            else
                o.steptol = 0.99;
            end

            if isfield(options, 'StepSame')
                o.stepSame = options.StepSame;
            else
                o.stepSame = 1;
            end

            if isfield(options, 'x0min')
                o.x0min = options.x0min;
            else
                o.x0min = 1;
            end

            if isfield(options, 'z0min')
                o.z0min = options.z0min;
            else
                o.z0min = 1;
            end

            if isfield(options, 'mu0')
                o.mu0 = options.mu0;
            else
                o.mu0 = 1.0e-1;
            end

            if isfield(options, 'backtrack')
                o.backtrack = options.backtrack;
            else
                o.backtrack = 0;
            end

            if isfield(options, 'file_id')
                o.file_id = options.file_id;
            else
                o.file_id = 1;
            end

            o.neg_lbnds = [];
            o.eigmin_neg = [];
            o.eigmax_neg = [];
            o.neg_ubnds = [];
            o.pos_lbnds = [];
            o.eigmin_pos = [];
            o.eigmax_pos = [];
            o.pos_ubnds = [];
            o.conds = [];

        end

        % TODO: slackmodel already does this!
        function categorize_bounds(o)
            % Categorize various types of bounds.
            % pos overlaps with low.
            % neg overlaps with upp.
            % two overlaps with low and upp.
            % fix and free are disjoint from all other sets.
            %
            % 26 Feb 2010: zlo and zup now point to variables with bl = 0 or bu = 0
            %              respectively, with bl < bu.  Needed to keep x
            %              strictly positive or negative.

            bigL = -9.9e+19;
            bigU = 9.9e+19;
            pos = find(o.bl == 0 & o.bu >= bigU);
            neg = find(o.bl <= bigL & o.bu == 0);
            o.low = find(o.bl > bigL & o.bl < o.bu);
            o.upp = find(o.bu < bigU & o.bl < o.bu);
            o.one = find((o.bl > bigL | o.bu < bigU) & o.bl < o.bu);  % low ou upp
            o.two = find(o.bl > bigL & o.bu < bigU & o.bl < o.bu);
            o.fix = find(o.bl == o.bu);
            o.nfix = length(o.fix);
            free = find(o.bl <= bigL & o.bu >= bigU);
            o.zlo = find(o.bl == 0 & o.bl < o.bu);
            o.zup = find(o.bu == 0 & o.bl < o.bu);

            if o.Prilev > 0
                fprintf(o.file_id, '\n\nBounds:\n  [0, inf]  [-inf, 0]');
                fprintf(o.file_id, '  Finite bl  Finite bu  Two bnds   Fixed    Free');
                fprintf(o.file_id, '\n %8g %9g %10g %10g %9g %7g %7g', ...
                    length(pos), length(neg), length(o.low),  ...
                    length(o.upp), length(o.two), length(o.fix), length(free));
                fprintf(o.file_id, '\n  [0, bu]  [bl,  0]  excluding fixed variables');
                fprintf(o.file_id, '\n %8g %9g', length(o.zlo), length(o.zup));
            end
        end

        function distrib(o, x, z)
            % distrib(x) or distrib(x, z) prints the
            % distribution of 1 or 2 vectors.
            %
            % 18 Dec 2000.  First version with 2 vectors.

            two2 = nargin > 1;
            fprintf(o.file_id, '\n\nDistribution of vector     x');
            if two2, fprintf(o.file_id, '         z'); end

            x1_d = 10^(floor(log10(max(x)+eps)) + 1);
            z1_d = 10^(floor(log10(max(z)+eps)) + 1);
            x1_d = max(x1_d, z1_d);
            kmax = 10;

            for k = 1:kmax
                x2_d = x1_d;    x1_d = x1_d/10;
                if k == kmax, x1_d = 0; end
                nx = length(find(x>= x1_d & x<x2_d));
                fprintf(o.file_id, '\n[%7.3g, %7.3g)%10g', x1_d, x2_d, nx);
                if two2
                    nz = length(find(z>= x1_d & z<x2_d));
                    fprintf(o.file_id, '%10g', nz);
                end
            end
        end

        function merit(o)
            % Evaluate the merit function for Newton's method.
            % It is the 2-norm of the three sets of residuals.

            f = [norm(o.r1)
                 norm(o.r2)
                 norm(o.rL(o.low))
                 norm(o.rU(o.upp))
                 norm(o.cL(o.low))
                 norm(o.cU(o.upp))];
            o.fmerit = norm(f);
        end

        function feasibility_resids(o)
            % Form residuals for the primal and dual equations.
            % rL, rU are output, but we input them as full vectors
            % initialized (permanently) with any relevant zeros.
            % 13 Aug 2003: z2-z1 coded more carefully
            %              (although MATLAB was doing the right thing).
            % 19 Nov 2003: r2(fix) = 0 has to be done after r2 = grad - r2;

            x_res1 = o.x;
            x_res1(o.fix) = 0;
            o.r1 = o.A * x_res1;
            o.r2 = o.A' * o.y ;

            o.r1 = o.b - o.r1 - (o.d2.^2) .* o.y;
            o.r2 = o.grad - o.r2;  % + (z2-z1);        % grad includes (d1.^2)*x
            o.r2(o.fix) = 0;
            o.r2(o.upp) = o.r2(o.upp) + o.z2(o.upp);
            o.r2(o.low) = o.r2(o.low) - o.z1(o.low);
            o.rL(o.low) = (o.bl(o.low) - o.x(o.low)) + o.x1(o.low);
            o.rU(o.upp) = (- o.bu(o.upp) + o.x(o.upp)) + o.x2(o.upp);

            o.Pinf = max([norm(o.r1, inf) norm(o.rL(o.low), inf) norm(o.rU(o.upp), inf)]);
            o.Dinf = norm(o.r2, inf);
            o.Pinf = max(o.Pinf, 1.0e-99);
            o.Dinf = max(o.Dinf, 1.0e-99);
        end

        function complementarity_resids(o)
            % Form residuals for the complementarity equations.
            % cL, cU are output, but we input them as full vectors
            % initialized (permanently) with any relevant zeros.
            % Cinf  is the complementarity residual for X1 z1 = mu e, etc.
            % Cinf0 is the same for mu = 0 (i.e., for the original problem).
            % 12 Feb 2009: Beware: If all variables are free, Cinf0 = empty.
            %                      Now changed to Cinf0 = 0;
            % 03 Apr 2010: If all variables are free, we need to define center.
            %              Set center = 1.0 arbitrarily.

            x1z1 = o.x1(o.low) .* o.z1(o.low);
            x2z2 = o.x2(o.upp) .* o.z2(o.upp);
            o.cL(o.low) = o.mu - x1z1;
            o.cU(o.upp) = o.mu - x2z2;

            maxXz = max([max(x1z1) max(x2z2)]);
            minXz = min([min(x1z1) min(x2z2)]);
            maxXz = max(maxXz, 1.0e-99);
            minXz = max(minXz, 1.0e-99);
            o.center = maxXz / minXz;
            o.Cinf = max([norm(o.cL(o.low), inf) norm(o.cU(o.upp), inf)]);
            o.Cinf0 = maxXz;
            if isempty(o.Cinf0)
                o.Cinf0 = 0;
                o.center = 1;
            end
        end

        function stepf = step_to_boundary(~, x, dx)
            % Assumes x > 0.
            % Finds the maximum step such that x + step*dx >= 0.

            stepf = 1.0e+20;
            blocking = find(dx < 0);
            if ~isempty(blocking)
                steps = x(blocking) ./ (- dx(blocking));
                stepf = min(steps);
            end
        end

        % Main function to solve the optimization problem
        function solve(o)

            o.neg_lbnds = [];
            o.eigmin_neg = [];
            o.eigmax_neg = [];
            o.neg_ubnds = [];
            o.pos_lbnds = [];
            o.eigmin_pos = [];
            o.eigmax_pos = [];
            o.pos_ubnds = [];
            o.conds = [];

            if o.Prilev > 0
               fprintf(o.file_id, '\n   --------------------------------------------------------');
               fprintf(o.file_id, '\n   pdco.m                            Version of 23 Nov 2013');
               fprintf(o.file_id, '\n   Primal-dual barrier method to minimize a convex function');
               fprintf(o.file_id, '\n   subject to linear constraints Ax + r = b,  bl <= x <= bu');
               fprintf(o.file_id, '\n                                                           ');
               fprintf(o.file_id, '\n   Michael Saunders       SOL and ICME, Stanford University');
               fprintf(o.file_id, '\n   Contributors:     Byunggyoo Kim (SOL), Chris Maes (ICME)');
               fprintf(o.file_id, '\n                     Santiago Akle (ICME), Matt Zahr (ICME)');
               fprintf(o.file_id, '\n   --------------------------------------------------------\n');
            end

            %---------------------------------------------------------------------
            % Decode A.
            %---------------------------------------------------------------------
            o.explicitA = isa(o.A, 'numeric');       % Not an opSpot.  Must be a matrix.
            if o.Prilev > 0
                if o.explicitA
                    if issparse(o.A)
                        fprintf(o.file_id, '\nA is a sparse matrix');
                        nnzA = nnz(o.A);
                        fprintf(o.file_id, '\n\nm = %8d     n = %8d      nnz(A) = %9d', o.m, o.n, nnzA);
                    else
                        fprintf(o.file_id, '\nA is a dense matrix');
                        fprintf(o.file_id, '\n\nm = %8d     n = %8d', o.m, o.n);
                    end
                else
                    fprintf(o.file_id, '\nA is a linear operator');
                    fprintf(o.file_id, '\n\nm = %8d     n = %8d', o.m, o.n);
                end
            end

            normb = norm(o.b , inf);   normx0 = norm(o.x0, inf);
            normy0 = norm(o.y0, inf);   normz0 = norm(o.z0, inf);

            if o.Prilev > 0
                fprintf(o.file_id, '\nmax |b | = %8d     max |x0| = %8.1.0e', normb , normx0);
                fprintf(o.file_id,                '      xsize = %8.1.0e', o.xsize);
                fprintf(o.file_id, '\nmax |y0| = %8d     max |z0| = %8.1.0e', normy0, normz0);
                fprintf(o.file_id,                '      zsize = %8.1.0e', o.zsize);
            end

            % Initialize.
            o.zn = zeros(o.n, 1);
            o.inner_itns = 0;
            o.inner_total = 0;
            o.inform = 0;

            % Set other parameters.
            o.eta = 1.0e-4;         % Linesearch tolerance for "sufficient descent"
            o.maxf = 10;            % Linesearch backtrack limit (function evaluations)
            o.maxfail = 1;          % Linesearch failure limit (consecutive iterations)
            o.bigcenter = 1.0e+3;   % mu is reduced if center < bigcenter
            o.gamma = max(o.d1);
            o.delta = max(o.d2);

            if o.Prilev > 0
                fprintf(o.file_id, '\n\nx0min = %8g     featol = %8.1.0e', o.x0min, o.featol);
                fprintf(o.file_id,                  '      d1max = %8.1.0e', o.gamma);
                fprintf(o.file_id,  '\nz0min = %8g     opttol = %8.1.0e', o.z0min, o.opttol);
                fprintf(o.file_id,                  '      d2max = %8.1.0e', o.delta);
                fprintf(o.file_id,  '\nmu0 = %8.1.0e     steptol = %8g', o.mu0  , o.steptol);
                fprintf(o.file_id,                  '     bigcenter = %8g'  , o.bigcenter);

                Print_param(o);
            end

            % Check for valid Method.
            o.time = cputime;
            if o.Prilev > 0
                tic
            end

            % Categorize bounds and allow for fixed variables by modifying b.
            categorize_bounds(o);

            if o.nfix > 0
                o.x1 = o.zn;
                o.x1(o.fix) = o.bl(o.fix);
                o.r1 = o.A * o.x1;
                o.b = o.b - o.r1;
                % At some stage, might want to look at normfix = norm(r1, inf);
            end

            % Scale the input data.
            % The scaled variables are
            %    xbar = x / beta,
            %    ybar = y / zeta,
            %    zbar = z / zeta.
            % Define
            %    theta = beta * zeta;
            % The scaled function is
            %    phibar = (1 / theta) fbar(beta * xbar),
            %    gradient = (beta / theta) grad,
            %    Hessian = (beta2 / theta) hess.
            o.beta = o.xsize;   if o.beta == 0, o.beta = 1; end    % beta scales b, x.
            o.zeta = o.zsize;   if o.zeta == 0, o.zeta = 1; end    % zeta scales y, z.
            o.theta = o.beta * o.zeta;                             % theta scales obj.
            % (theta could be anything, but theta = beta * zeta makes
            % scaled grad = grad / zeta = 1 approximately if zeta is chosen right.)

            o.bl(o.fix) = o.bl(o.fix) / o.beta;
            o.bu(o.fix) = o.bu(o.fix) / o.beta;
            o.bl(o.low) = o.bl(o.low) / o.beta;
            o.bu(o.upp) = o.bu(o.upp) / o.beta;
            o.d1 = o.d1 * (o.beta / sqrt(o.theta));
            o.d2 = o.d2 * (sqrt(o.theta) / o.beta);

            o.b = o.b / o.beta;   o.y0 = o.y0 / o.zeta;
            o.x0 = o.x0 / o.beta; o.z0 = o.z0 / o.zeta;

            % Initialize vectors that are not fully used if bounds are missing.
            o.rL = o.zn;   o.rU = o.zn;
            o.cL = o.zn;   o.cU = o.zn;
            o.x1 = o.zn;   o.x2 = o.zn;
            o.z1 = o.zn;   o.z2 = o.zn;
            o.dx1 = o.zn;   o.dx2 = o.zn;
            o.dz1 = o.zn;   o.dz2 = o.zn;

            % Initialize x, y, z1, z2, objective, etc.
            % 10 Aug 2003: z isn't needed here -- just at end for output.
            % 03 Jul 2008: Use pos to ensure that x = x1 for vanilla bounds.
            %              The linear constraints x(pos) - x1(pos) = bl(pos) = 0
            %              should then remain satisfied, so we'll have x(pos) > 0
            %              throughout.  At last, no danger of evaluating log(x)
            %              at negative x values.
            % 26 Feb 2010: zlo and zup now used in place of pos.
            %              See 26 Feb 2010 note above.
            o.x = o.x0;
            o.y = o.y0;
            o.x(o.fix) = o.bl(o.fix);
            o.x(o.low) = max(o.x(o.low), o.bl(o.low));
            o.x(o.upp) = min(o.x(o.upp), o.bu(o.upp));
            o.x1(o.low) = max(o.x(o.low) - o.bl(o.low), o.x0min);
            o.x2(o.upp) = max(o.bu(o.upp) - o.x(o.upp), o.x0min);
            o.z1(o.low) = max(o.z0(o.low), o.z0min);
            o.z2(o.upp) = max(-o.z0(o.upp), o.z0min);
            o.x(o.zlo) = o.x1(o.zlo);
            o.x(o.zup) = -o.x2(o.zup);

            [o.obj, o.grad, o.hess] = o.slack.obj(o.x*o.beta);
            [mH, nH] = size(o.hess);
            if o.diagHess
                if mH > 1 && nH > 1
                    fprintf('\n Warning: Using only diagonal part of hess from pdObj\n');
                    o.hess = diag(o.hess);
                end
            else
                if mH ~= o.n && nH ~= o.n
                    error('Hessian size mismatch');
                end
            end

            o.obj = o.obj / o.theta; % Scaled obj.
            o.grad = o.grad*(o.beta / o.theta) + (o.d1.^2) .* o.x; % grad includes x regularization.
            o.H = o.hess * (o.beta * o.beta / o.theta);
            if o.diagHess
                o.H = o.H + (o.d1.^2); % H includes x regularization.
            else
                o.H = o.H + sparse(1:o.n, 1:o.n, (o.d1.^2), o.n, o.n); % H includes x regularization.
            end

            % Compute primal and dual feasibility residuals:
            %    r1 = b - A*x - d2.^2*y
            %    r2 = grad - A'*y + (z2-z1)
            %    rL = bl - x + x1
            feasibility_resids(o);

            % Initialize mu and complementarity residuals:
            %    cL = mu*e - X1*z1.
            %    cU = mu*e - X2*z2.
            %
            % 25 Jan 2001: Now that b and obj are scaled (and hence x, y, z),
            %              we should be able to use mufirst = mu0 (absolute value).
            %              0.1 worked poorly on StarTest1 with x0min = z0min = 0.1.
            % 29 Jan 2001: We might as well use mu0 = x0min * z0min;
            %              so that most variables are centered after a warm start.
            % 29 Sep 2002: Use mufirst = mu0*(x0min * z0min),
            %              regarding mu0 as a scaling of the initial center.
            % 07 Aug 2003: mulast is controlled by opttol.
            %              mufirst should not be smaller.
            % 10 Aug 2003: Revert to mufirst = mu0 (absolute value).
            % 12 Jan 2006: If mu0 <= 0, reset it to the average complementarity.
            if o.mu0 <= 0
                clow = o.x1(o.low) .* o.z1( o.low);
                cupp = o.x2(o.upp) .* o.z2( o.upp);
                o.mu0 = (sum(clow) + sum(cupp)) / (length(o.low) + length(o.upp));
                o.mu0 = o.mu0*0.1;
                clear clow cupp
            end

            % mufirst = mu0 * (x0min * z0min);
            o.mufirst = o.mu0;
            o.mulast = 0.1 * o.opttol;
            o.mufirst = max(o.mufirst, o.mulast);
            o.mu = o.mufirst;
            complementarity_resids(o);

            % Compute initial merit function value.
            merit(o);

            % Initialize other things.
            o.PDitns = 0;
            converged = 0;
            Init_param(o);

            %  Iteration log.
            o.nfail = 0;
            o.regterm = norm(o.d1 .* o.x)^2 + norm(o.d2 .* o.y)^2;
            o.objreg = o.obj + 0.5 * o.regterm;
            o.objtrue = o.objreg * o.theta;

            if o.Prilev > 0
                head1 = '\n\nItn   mu stepx stepz  Pinf  Dinf';
                head2 = '  Cinf   Objective    nf  center';
                fprintf(o.file_id, [ head1 head2 o.head3 ]);
                fprintf(o.file_id, '\n%3g                 ', o.PDitns       );
                fprintf(o.file_id, '%6.1f%6.1f', log10(o.Pinf), log10(o.Dinf));
                fprintf(o.file_id, '%6.1f%15.7e', log10(o.Cinf0), o.objtrue);
                fprintf(o.file_id, '   %8.1f', o.center);
            end

            % Main loop.
            while ~converged
                o.PDitns = o.PDitns + 1;

                % record eigenvalues of o.H before Solve_Newton modifies it.
                eigsH = eig(full(o.H));
                sigmas = svd(full(o.A));

                % Compute step. This is performed in a subclass.
                Solve_Newton(o);

                % Record condition number and eigenvalues of Newton system.
                fullM = full(o.M);
                eigsM = eig(fullM);
                o.eigmin_neg = [o.eigmin_neg, min(eigsM)];
                o.eigmax_neg = [o.eigmax_neg, max(eigsM(eigsM < 0))];
                o.eigmin_pos = [o.eigmin_pos, min(eigsM(eigsM > 0))];
                o.eigmax_pos = [o.eigmax_pos, max(eigsM)];
                o.conds = [o.conds, cond(fullM)];

                eigHmin = min(eigsH);
                eigHmax = max(eigsH);
                sigmin = min(sigmas);
                sigmax = max(sigmas);

                % x1x2 = x1 for i in o.low and = x2 for i in o.upp
                xmin = min(min(o.x1(o.low)), min(o.x2(o.upp)));
                xmax = max(max(o.x1(o.low)), max(o.x2(o.upp)));

                x1z2 = zeros(o.n,1);
                x1z2(o.upp) = o.z2(o.upp);
                x1z2(o.low) = x1z2(o.low) .* o.x1(o.low);

                x2z1 = zeros(o.n,1);
                x2z1(o.low) = o.z1(o.low);
                x2z1(o.upp) = x2z1(o.upp) .* o.x2(o.upp);

                zmin = min(min(x1z2(o.low)), min(x2z1(o.low)));
                zmax = max(max(x1z2(o.upp)), max(x2z1(o.upp)));

                epsmin = (eigHmin * xmin + zmin);  % d1^2 is already in o.H.
                epsmax = (eigHmax * xmax + zmax);

                neg_lbnd = 0.5 * (o.d2^2 - epsmax - sqrt((o.d2^2 + epsmax)^2 + 4 * sigmax^2 * xmax));
                neg_ubnd = -epsmin;
                pos_lbnd = 0.5 * (o.d2^2 - epsmax + sqrt((o.d2^2 + epsmax)^2 + 4 * sigmin^2 * xmin));
                pos_ubnd = 0.5 * (o.d2^2 - epsmin + sqrt((o.d2^2 + epsmin)^2 + 4 * sigmax^2 * xmax));

                o.neg_lbnds = [o.neg_lbnds, neg_lbnd];
                o.neg_ubnds = [o.neg_ubnds, neg_ubnd];
                o.pos_lbnds = [o.pos_lbnds, pos_lbnd];
                o.pos_ubnds = [o.pos_ubnds, pos_ubnd];

                if o.inform == 4
                    break
                end

                % Find the maximum step size.
                % 13 Aug 2003: We need stepxL, stepxU also to keep x feasible
                %              so that nonlinear functions are defined.
                % 18 Nov 2003: But this gives stepx = 0 for lptest.  (??)
                stepx1 = step_to_boundary(o, o.x1(o.low), o.dx1(o.low));
                stepx2 = step_to_boundary(o, o.x2(o.upp), o.dx2(o.upp));
                stepz1 = step_to_boundary(o, o.z1(o.low), o.dz1(o.low));
                stepz2 = step_to_boundary(o, o.z2(o.upp), o.dz2(o.upp));
                %  stepxL = step_to_boundary( x(low),  dx(low));
                %  stepxU = step_to_boundary( x(upp),  dx(upp));
                %  stepx = min([stepx1,   stepx2,   stepxL,   stepxU]);
                stepx = min([stepx1, stepx2]);
                stepz = min([stepz1, stepz2]);
                stepx = min([o.steptol * stepx, 1]);
                stepz = min([o.steptol * stepz, 1]);
                if o.stepSame                      % For NLPs, force same step
                    stepx = min(stepx, stepz);   % (true Newton method)
                    stepz = stepx;
                end

                % Backtracking linesearch.
                fail = true;
                o.nf = 0;

                while o.nf < o.maxf
                    o.nf = o.nf + 1;
                    o.x1(o.low) = o.x1(o.low) + stepx * o.dx1(o.low);
                    o.x2(o.upp) = o.x2(o.upp) + stepx * o.dx2(o.upp);
                    o.z1(o.low) = o.z1(o.low) + stepz * o.dz1(o.low);
                    o.z2(o.upp) = o.z2(o.upp) + stepz * o.dz2(o.upp);
                    o.x = o.x + stepx * o.dx;
                    o.y = o.y + stepz * o.dy;
                    o.x(o.zlo) = o.x1(o.zlo);
                    o.x(o.zup) = -o.x2(o.zup);

                    [o.obj, o.grad, o.hess] = o.slack.obj(o.x * o.beta);
                    if o.diagHess
                        if mH>1 && nH>1
                            o.hess = diag(o.hess);
                        end
                    end
                    o.obj = o.obj / o.theta;
                    o.grad = o.grad * (o.beta / o.theta) + (o.d1.^2) .* o.x;
                    o.H = o.hess * (o.beta * o.beta / o.theta);
                    if o.diagHess
                        o.H = o.H + (o.d1.^2); % H includes x regularization.
                    else
                        o.H = o.H + sparse(1:o.n, 1:o.n, (o.d1.^2), o.n, o.n); % H includes x regularization.
                    end

                    feasibility_resids(o);
                    complementarity_resids(o);
                    merit(o);
                    o.step = min(stepx, stepz);

                    if ~o.backtrack || fmeritnew <= (1 - o.eta*o.step)*o.fmerit
                        fail = false;
                        break;
                    end

                    % Merit function didn't decrease.
                    % Restore variables to previous values.
                    % (This introduces a little error, but save lots of space.)
                    o.x = o.x - stepx * o.dx;
                    o.y = o.y - stepz * o.dy;
                    o.x1(o.low) = o.x1(o.low) - stepx * o.dx1(o.low);
                    o.x2(o.upp) = o.x2(o.upp) - stepx * o.dx2(o.upp);
                    o.z1(o.low) = o.z1(o.low) - stepz * o.dz1(o.low);
                    o.z2(o.upp) = o.z2(o.upp) - stepz * o.dz2(o.upp);
                    o.x(o.zlo) = o.x1(o.zlo);
                    o.x(o.zup) = -o.x2(o.zup);

                    % Back-track.
                    % If it's the first time,
                    % make stepx and stepz the same.
                    if o.nf == 1 && stepx ~= stepz
                        stepx = o.step;
                    elseif o.nf < o.maxf
                        stepx = stepx / 2;
                    end;
                        stepz = stepx;
                end

                if fail
                    fprintf(o.file_id, '\n     Linesearch failed (nf too big)');
                    o.nfail = o.nfail + 1;
                else
                    o.nfail = 0;
                end

                % Set convergence measures.
                o.regterm = norm(o.d1.*o.x)^2 + norm(o.d2.*o.y)^2;
                o.objreg = o.obj + 0.5 * o.regterm;
                o.objtrue = o.objreg * o.theta;

                primalfeas = o.Pinf <= o.featol;
                dualfeas = o.Dinf <= o.featol;
                complementary = o.Cinf0 <= o.opttol;
                enough = o.PDitns >= 4; % Prevent premature termination.
                converged = primalfeas & dualfeas & complementary & enough;

                % Iteration log.
                if o.Prilev > 0
                    str1 = sprintf('\n%3g%5.1f', o.PDitns, log10(o.mu));
                    str2 = sprintf('%6.3f%6.3f', stepx, stepz);
                    if stepx < 0.001 || stepz < 0.001
                        str2 = sprintf('%6.1.0e%6.1.0e', stepx, stepz);
                    end
                    str3 = sprintf('%6.1f%6.1f', log10(o.Pinf), log10(o.Dinf));
                    str4 = sprintf('%6.1f%15.7e', log10(o.Cinf0), o.objtrue);
                    str5 = sprintf('%3g%8.1f', o.nf, o.center);
                    if o.center > 99999
                        str5 = sprintf('%3g%8.1.0e', o.nf, o.center);
                    end
                    fprintf(o.file_id, [str1 str2 str3 str4 str5]);
                    Print_results(o);
                end

                % Test for termination.
                if converged
                    if o.Prilev > 0
                        fprintf(o.file_id, '\n   Converged');
                    end
                elseif o.PDitns >= o.maxitn
                    if o.Prilev > 0
                        fprintf(o.file_id, '\n   Too many iterations');
                    end
                    o.inform = 1;
                    break
                elseif o.nfail >= o.maxfail
                    if o.Prilev > 0
                        fprintf(o.file_id, '\n   Too many linesearch failures');
                    end
                    o.inform = 2;
                    break
                elseif o.step <= 1.0e-10
                    if o.Prilev > 0
                        fprintf(o.file_id, '\nStep lengths too small');
                    end
                    o.inform = 3;
                    break
                else

                    % Reduce mu, and reset certain residuals.
                    stepmu = min(stepx, stepz);
                    stepmu = min(stepmu, o.steptol);
                    muold = o.mu;
                    mumin = 0.1 * max([o.Pinf o.Dinf o.Cinf]); % Target mu shouldn't be too small.
                    mumin = min(o.mu, mumin); % mu should be monotonic.

                    o.mu = o.mu - stepmu * o.mu; % mu being treated as a variable.
                    if o.center >= o.bigcenter, o.mu = muold; end  % Keep old mu if far from center.
                    o.mu = max(o.mu, mumin); % 04 May 2008: No smaller than target.
                    o.mu = max(o.mu, o.mulast); % 13 Jun 1998: No need for smaller mu.

                    % mutrad = mu0 * (sum(Xz) / n); % 24 May 1998: Traditional value, but
                    % mu = min(mu, mutrad); % it seemed to decrease mu too much.

                    o.mu = max(o.mu, o.mulast);  % 13 Jun 1998: No need for smaller mu.
                    complementarity_resids(o);
                    merit(o);
                    Reset_param(o);
                end
            end
            %---------------------------------------------------------------------
            % End of main loop.
            %---------------------------------------------------------------------

            % Reconstruct z.
            o.x(o.fix) = 0; % Exclude x(fix) temporarily from |x|.
            o.z = zeros(o.n, 1); % Exclude z(fix) also.
            o.z(o.low) = o.z1(o.low);
            o.z(o.upp) = o.z(o.upp) - o.z2(o.upp);

            % Print statistics.
            if o.Prilev > 0
                fprintf(o.file_id, '\n\nmax |x| = %10.3f', norm(o.x, inf));
                fprintf(o.file_id, '    max |y| = %10.3f', norm(o.y, inf));
                fprintf(o.file_id, '    max |z| = %10.3f', norm(o.z, inf));  % excludes z(fix)
                fprintf(o.file_id, '   scaled');
            end

            o.bl(o.fix) = o.bl(o.fix) * o.beta; % Unscale bl, bu, x, y, z.
            o.bu(o.fix) = o.bu(o.fix) * o.beta;
            o.bl(o.low) = o.bl(o.low) * o.beta;
            o.bu(o.upp) = o.bu(o.upp) * o.beta;

            o.x = o.x * o.beta; o.y = o.y * o.zeta; o.z = o.z * o.zeta;

            if o.Prilev > 0
                fprintf(o.file_id,  '\nmax |x| = %10.3f', norm(o.x, inf));
                fprintf(o.file_id, '    max |y| = %10.3f', norm(o.y, inf));
                fprintf(o.file_id, '    max |z| = %10.3f', norm(o.z, inf));  % excludes z(fix)
                fprintf(o.file_id, ' unscaled');
            end

            % Reconstruct x(fix) and z(fix).
            o.r2 = o.A' * o.y ;
            if o.nfix > 0
                o.x(o.fix) = o.bl(o.fix);
                o.z(o.fix) = o.grad(o.fix) - o.r2(o.fix); % z = grad - A'y
                o.z1(o.fix) = max(o.z(o.fix), 0);
                o.z2(o.fix) = max(-o.z(o.fix), 0);
            end

            % Reconstruct b.
            o.b = o.b * o.beta;
            if o.nfix > 0
                o.x1 = zeros(o.n, 1);
                o.x1(o.fix) = o.bl(o.fix);
                o.r1 = o.A *o.x1 ;
                o.b = o.b + o.r1;
                if o.Prilev > 0
                    fprintf(o.file_id, '\nmax |x| and max |z| exclude fixed variables');
                end
            end

            % Evaluate function at final point.
            % Recompute all of z.
            % 03 Apr 2010: z now includes (d1.^2).*x from the primal regularization.
            [o.obj, o.grad, o.hess] = o.slack.obj(o.x);
            o.z = o.grad + (o.d1.^2) .* o.x - o.r2; % z = grad + (d1.^2).*x - A'y

            o.time = cputime - o.time;
            if o.Prilev > 0
                fprintf(o.file_id, ...
                        '\nPDitns = %10g %sitns = %10g    cputime = %10.1f', ...
                        o.PDitns, o.solver, o.inner_total, o.time);
                distrib(o, abs(o.x), abs(o.z));
                toc
            end
        end
    end
end
