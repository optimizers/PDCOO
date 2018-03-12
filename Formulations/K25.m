classdef K25 < handle

    properties
        M
        rhs
        sol
    end

    methods (Abstract)
        Solver(o)
    end

    methods
        function o = K25(options)
            o.diagHess = false;
        end

        function y = opK25(x, ~)
            error('Implement me!')
            % dx = x(1:o.n);
            % dy = x(o.n+1:o.n+o.m);
            % u = -o.H*dx + o.A'*dy;
            % v = o.A*dx + o.d2.^2 .* dy;
            % y = [u;v];
        end

        function Solve_Newton(o)
            %-----------------------------------------------------------------
            %  Solve (*) for [dx ; dy].
            %-----------------------------------------------------------------
            %  Define a damped Newton iteration for solving f = 0,
            %  keeping  x1, x2, z1, z2 > 0.  We eliminate dx1, dx2, dz1, dz2
            %  to obtain the system
            %
            %  [-H2  A' ] [dx] = [w ]     (*),   H2 = H + D1^2 + X1inv Z1 + X2inv Z2,
            %  [ A  D2^2] [dy] = [r1]            w  = r2 - X1inv(cL + Z1 rL)
            %                                       + X2inv(cU + Z2 rU),
            %
            %----------------------------------------------------------------

            x1s = sqrt(o.x1); X1s = spdiags(x1s, 0, o.n, o.n);
            x2s = sqrt(o.x2); X2s = spdiags(x2s, 0, o.n, o.n);

            % Suboptimal!!
            o.H(o.upp, :) = X2s(o.upp, o.upp) * o.H(o.upp, :);
            o.H(o.low, :) = X1s(o.low, o.low) * o.H(o.low, :);
            o.H(:, o.low) = o.H(:, o.low) * X1s(o.low, o.low);
            o.H(:, o.upp) = o.H(:, o.upp) * X2s(o.upp, o.upp);

            %o.H      = o.H + sparse(o.two,o.two, o.z1(o.two).*o.x2(o.two), o.n, o.n);
            %o.H      = o.H + sparse(o.two,o.two, o.z2(o.two).*o.x1(o.two), o.n, o.n);

            x1z2 = zeros(o.n,1);
            x1z2(o.upp) = o.z2(o.upp);
            x1z2(o.low) = x1z2(o.low) .* o.x1(o.low);

            x2z1 = zeros(o.n,1);
            x2z1(o.low) = o.z1(o.low);
            x2z1(o.upp) = x2z1(o.upp) .* o.x2(o.upp);

            o.H = o.H + spdiags(x1z2, 0, o.n, o.n) + spdiags(x2z1, 0, o.n, o.n);

            A = o.A;
            A(:, o.low) = A(:, o.low) * X1s(o.low, o.low);
            A(:, o.upp) = A(:, o.upp) * X2s(o.upp, o.upp);

            r3 = zeros(o.n,1);
            r3(o.low) = o.cL(o.low) + o.z1(o.low) .* o.rL(o.low);
            r3(o.low) = r3(o.low) ./ x1s(o.low);
            r3(o.upp) = r3(o.upp) .* x2s(o.upp);

            r4 = zeros(o.n,1);
            r4(o.upp) = o.cU(o.upp) + o.z2(o.upp) .* o.rU(o.upp);
            r4(o.low) = r4(o.low) .* x1s(o.low);
            r4(o.upp) = r4(o.upp) ./ x2s(o.upp);

            % Assemble rhs.
            w = o.r2;
            w(o.low) = w(o.low) .* x1s(o.low);
            w(o.upp) = w(o.upp) .* x2s(o.upp);

            w = w - r3 + r4;

            if o.nfix > 0 && o.explicitA
                [ih, jh, vh] = find(o.H);
                for k = o.fix'
                    vh(ih == k & ih ~= jh) = 0;
                    vh(jh == k & ih ~= jh) = 0;
                end
                o.H = sparse(ih, jh, vh);
            end

            if ~o.explicitA
                o.M = opFunction(o.m + o.n, o.m + o.n, @opK25);
            else
                if o.nfix == 0
                    o.M = [ -o.H  A'
                             A    sparse(1:o.m, 1:o.m, o.d2.^2, o.m,o.m)];
                else
                    A(:, o.fix) = 0;
                    o.M = [ -o.H  A'
                             A    sparse(1:o.m, 1:o.m, o.d2.^2, o.m,o.m)];
                end
            end

            if o.need_precon
                o.precon = speye(size(o.M));
            end

            o.rhs = [w; o.r1];
            o.rhs(o.fix) = 0;

            Solver(o);

            o.dx  = o.sol(1:o.n);
            o.dx(o.low) = o.dx(o.low) .* x1s(o.low);
            o.dx(o.upp) = o.dx(o.upp) .* x2s(o.upp);
            o.dy  = o.sol(o.n+1:o.n+o.m);

            o.dx1(o.low) = -o.rL(o.low) + o.dx(o.low);
            o.dx2(o.upp) = -o.rU(o.upp) - o.dx(o.upp);

            o.dz1(o.low) = (o.cL(o.low) - o.z1(o.low) .* o.dx1(o.low)) ./ o.x1(o.low);
            o.dz2(o.upp) = (o.cU(o.upp) - o.z2(o.upp) .* o.dx2(o.upp)) ./ o.x2(o.upp);
        end
    end
end
