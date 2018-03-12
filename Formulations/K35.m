classdef K35 < handle

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
        function o = K35(options)
            o.diagHess = false;
        end

        function Solve_Newton(o)
            %-----------------------------------------------------------------
            %  Solve (*) for [dx ; dy ; dz1* ; dz2*].
            %-----------------------------------------------------------------
            %  Define a damped Newton iteration for solving f = 0,
            %  keeping  x1, x2 > 0.  We eliminate dx1, dx2
            %  to obtain the system
            %
            %  [-H1     A'    Z1^1/2   -Z2^1/2] [dx]   =                     [r2]     (*),   H1 = H + D1^2
            %  [ A     D2^2     0          0  ] [dy]   =                     [r1]          dz1* = Z1^-1/2 * dz1
            %  [Z1^1/2   0      X1         0  ] [dz1*] = [Z1^-1/2*cl + Z1^1/2*rl]          dz2* = Z2^-1/2 * dz2
            %  [-Z2^1/2  0      0          X2 ] [dz2*] = [Z2^-1/2*cu + Z2^1/2*ru]
            %----------------------------------------------------------------

            nlow = length(o.low) ; nupp = length(o.upp);

            % Z1s = √Z1, Z2s = √Z2 but are rectangular.
            z1s = sqrt(o.z1(o.low)); Z1s = sparse(1:nlow, o.low, z1s, nlow, o.n);
            z2s = sqrt(o.z2(o.upp)); Z2s = sparse(1:nupp, o.upp, z2s, nupp, o.n);

            if o.nfix > 0
                [ih, jh, vh] = find(o.H);
                for k = o.fix'
                    vh(ih == k & ih ~= jh) = 0;
                    vh(jh == k & ih ~= jh) = 0;
                end
                o.H = sparse(ih, jh, vh, o.n, o.n);
            end

            % X1 and X2 are square.
            X1 = sparse(1:nlow, 1:nlow, o.x1(o.low), nlow, nlow);
            X2 = sparse(1:nupp, 1:nupp, o.x2(o.upp), nupp, nupp);

            D2reg = sparse(1:o.m, 1:o.m, o.d2.^2, o.m,o.m);

            if o.explicitA

                if o.nfix == 0
                    o.M = [ -o.H    o.A'               Z1s'                -Z2s'
                             o.A    D2reg              sparse(o.m,nlow)     sparse(o.m,nupp)
                             Z1s    sparse(nlow,o.m)   X1                   sparse(nlow,nupp)
                            -Z2s    sparse(nupp,o.m)   sparse(nupp,nlow)    X2                ];
                else
                    if o.PDitns == 1
                        o.Afree = o.A;
                        o.Afree(:, o.fix) = 0;
                    end
                    o.M = [ -o.H      o.Afree'           Z1s'                -Z2s'
                             o.Afree  D2reg              sparse(o.m,nlow)     sparse(o.m,nupp)
                             Z1s      sparse(nlow,o.m)   X1                   sparse(nlow,nupp)
                            -Z2s      sparse(nupp,o.m)   sparse(nupp,nlow)    X2                ];
                end

            else

                o.M = [ -o.H    o.A'                Z1s'                -Z2s'
                         o.A    D2reg               sparse(o.m,nlow)     sparse(o.m,nupp)
                         Z1s    sparse(nlow,o.m)    X1                   sparse(nlow,nupp)
                        -Z2s    sparse(nupp,o.m)    sparse(nupp,nlow)    X2                ];

            end

            rmu1 = (o.cL(o.low) + o.z1(o.low) .* o.rL(o.low)) ./ z1s;
            rmu2 = (o.cU(o.upp) + o.z2(o.upp) .* o.rU(o.upp)) ./ z2s;

            o.rhs = [ o.r2 ; o.r1 ; rmu1 ; rmu2 ];

            if o.explicitA
                o.rhs(o.fix) = 0;
            end

            if o.need_precon
                o.precon = speye(size(o.M));
            end

            Solver(o);

            o.dx  = o.sol(1 : o.n);
            o.dy  = o.sol(o.n+1 : o.n+o.m);
            o.dz1(o.low) = o.sol(o.n+o.m+1 : o.n+o.m+nlow) .* z1s;
            o.dz2(o.upp) = o.sol(o.n+o.m+nlow+1 : o.n+o.m+nlow+nupp) .* z2s;

            o.dx1(o.low)  = -o.rL(o.low) + o.dx(o.low);
            o.dx2(o.upp)  = -o.rU(o.upp) - o.dx(o.upp);
        end
    end
end
