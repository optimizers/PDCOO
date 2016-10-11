classdef K2 < pdcoO

    properties
        M
        rhs
        sol
        Afree
    end

    methods (Abstract)
        Solver(o)
    end

    methods
        function o = K2(slack, options)
            o = o@pdcoO(slack, options);
            o.diagHess = false;
        end

        function y = opK2(x, ~)
            dx = x(1:o.n);
            dy = x(o.n+1:o.n+o.m);
            u = -o.H * dx + o.A' * dy;
            v = o.A * dx + o.d2.^2 .* dy;
            y = [u; v];
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
            o.H      = o.H + sparse(o.low,o.low, o.z1(o.low)./o.x1(o.low), o.n, o.n);
            o.H      = o.H + sparse(o.upp,o.upp, o.z2(o.upp)./o.x2(o.upp), o.n,o.n);
            w      = o.r2;
            w(o.low) = w(o.low) - (o.cL(o.low) + o.z1(o.low).*o.rL(o.low))./o.x1(o.low);
            w(o.upp) = w(o.upp) + (o.cU(o.upp) + o.z2(o.upp).*o.rU(o.upp))./o.x2(o.upp);

            if o.nfix > 0 && o.explicitA
                [ih, jh, vh] = find(o.H);
                for k = o.fix'
                    vh(ih == k & ih ~= jh) = 0;
                    vh(jh == k & ih ~= jh) = 0;
                end
                o.H = sparse(ih, jh, vh);
            end

            if ~o.explicitA
                o.M = opFunction(o.m + o.n, o.m + o.n, @opK2);
            else
                if o.nfix == 0
                    o.M = [ -o.H  o.A'
                             o.A  sparse(1:o.m, 1:o.m, o.d2.^2, o.m,o.m)];
                else
                    if o.PDitns == 1
                        o.Afree = o.A;
                        o.Afree(:, o.fix) = 0;
                    end
                    o.M = [ -o.H      o.Afree'
                             o.Afree  sparse(1:o.m, 1:o.m, o.d2.^2, o.m,o.m)];
                end
            end

            if o.need_precon
                o.precon = speye(size(o.M));
            end

            o.rhs = [w; o.r1];
            o.rhs(o.fix) = 0;

            Solver(o);

            o.dx = o.sol(1:o.n);
            o.dy = o.sol(o.n+1:o.n+o.m);

            o.dx1(o.low) = -o.rL(o.low) + o.dx(o.low);
            o.dx2(o.upp) = -o.rU(o.upp) - o.dx(o.upp);
            o.dz1(o.low) = (o.cL(o.low) - o.z1(o.low).*o.dx1(o.low)) ./ o.x1(o.low);
            o.dz2(o.upp) = (o.cU(o.upp) - o.z2(o.upp).*o.dx2(o.upp)) ./ o.x2(o.upp);
        end
    end
end
