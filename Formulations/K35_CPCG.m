classdef K35_CPCG < pdcoO

    properties
        C
        P
        rhs
        sol
        Afree
    end

    methods (Abstract)
        Solver(o)
    end

    methods
        function o = K35_CPCG(slack, options)
            o = o@pdcoO(slack, options);
            o.diagHess = false;
        end

        function Solve_Newton(o)

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

            rmu1 = (o.cL(o.low) + o.z1(o.low) .* o.rL(o.low)) ./ z1s;
            rmu2 = (o.cU(o.upp) + o.z2(o.upp) .* o.rU(o.upp)) ./ z2s;

            if o.nfix==0
                J = [ o.A
                      Z1s
                     -Z2s   ];
            else
                if o.PDitns == 1
                    o.Afree = o.A;
                    o.Afree(:, o.fix) = 0;
                end
                J = [ o.Afree
                      Z1s
                     -Z2s     ];
            end

            dH = spdiags(diag(o.H), 0, o.n, o.n);
            o.C  = spdiags([o.d2.^2 * ones(o.m,1) ; o.x1(o.low) ; o.x2(o.upp)], 0, o.m+nlow+nupp, o.m+nlow+nupp);
            G  = [ dH  J'
                    J  -o.C ];

            % Shift system.
            % TODO: use alternative starting point.
            yz0 = o.C \ [o.r1 ; rmu1 ; rmu2];
            o.rhs = o.r2 - J' * yz0;

            % Full LDL of constraint preconditioner.
            o.P = opLDL(G, 0);

            Solver(o);

            o.dx = o.sol;

            % Recover dyz.
            dyz = o.C \ (J * o.dx);
            dyz = yz0 - dyz;
            o.dy  = dyz(1:o.m);
            o.dz1(o.low) = dyz(o.m+1:o.m+nlow) .* z1s;
            o.dz2(o.upp) = dyz(o.m+nlow+1:o.m+nlow+nupp) .* z2s;

            o.dx1(o.low)  = -o.rL(o.low) + o.dx(o.low);
            o.dx2(o.upp)  = -o.rU(o.upp) - o.dx(o.upp);
        end
    end
end
