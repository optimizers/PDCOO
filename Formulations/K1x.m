classdef K1x < handle

    properties
        M
        rhs
        sol
        norm_A
    end

    methods (Abstract)
        Solver(o)
    end

    methods
        function o = K1x(options)
            o.diagHess = true;
        end

        function y = opK1y(x, ~)
            t = o.A * x;
            t = t ./ (o.d2.^2);
            t = o.A' * t;
            y = t + o.H * x;
        end

        function Solve_Newton(o)
            %-----------------------------------------------------------------
            %  Solve (*) for dx.
            %-----------------------------------------------------------------
            %  Define a damped Newton iteration for solving f = 0,
            %  keeping  x1, x2, z1, z2 > 0.  We eliminate dx1, dx2, dz1, dz2
            %  to obtain the system
            %
            %  [-H2  A' ] [dx] = [w ],   H2 = H + D1^2 + X1inv Z1 + X2inv Z2,
            %  [ A  D2^2] [dy] = [r1]    w  = r2 - X1inv(cL + Z1 rL)
            %                                    + X2inv(cU + Z2 rU),
            %
            %  Then we eliminate dy to obtain :
            %
            %     (A'*D2^-2*A + H2) dx = A'*D2^-2*r1 - w   (*)
            %-----------------------------------------------------------------
            % For this method to work, H must be diagonal
            o.H(o.low) = o.H(o.low) + o.z1(o.low) ./ o.x1(o.low);
            o.H(o.upp) = o.H(o.upp) + o.z2(o.upp) ./ o.x2(o.upp);
            o.H(o.fix) = Inf;
            w = o.r2;
            w(o.low) = w(o.low) - (o.cL(o.low) + o.z1(o.low) .* o.rL(o.low)) ./ o.x1(o.low);
            w(o.upp) = w(o.upp) + (o.cU(o.upp) + o.z2(o.upp) .* o.rU(o.upp)) ./ o.x2(o.upp);

            if o.PDitns == 1
                o.norm_A = normest(o.A, 1.0e-2);
            end

            if o.need_precon
                o.precon = 1 ./ (o.norm_A^2 / o.d2^2 + 1./o.H);
                o.precon = diag(sparse(o.precon));
            end

            if o.explicitA
                o.M = o.A' * sparse(1:o.m, 1:o.m, 1 ./ (o.d2.^2), o.m, o.m) * o.A;
                o.M = o.M + sparse(1:o.n, 1:o.n, o.H, o.n, o.n);
            else
                o.M = opFunction(o.m, o.m, @opK1y);
            end

            o.rhs = o.A' * (o.r1 ./ (o.d2.^2)) - w;

            Solver(o);

            o.dx = o.sol;

            % dx is now known.  Get dy, dx1, dx2, dz1, dz2.
            o.dy = (o.r1 - o.A * o.dx) ./ (o.d2.^2);
            o.dx1(o.low) = -o.rL(o.low) + o.dx(o.low);
            o.dx2(o.upp) = -o.rU(o.upp) - o.dx(o.upp);
            o.dz1(o.low) = (o.cL(o.low) - o.z1(o.low) .* o.dx1(o.low)) ./ o.x1(o.low);
            o.dz2(o.upp) = (o.cU(o.upp) - o.z2(o.upp) .* o.dx2(o.upp)) ./ o.x2(o.upp);
        end
    end
end
