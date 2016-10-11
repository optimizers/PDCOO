classdef K2_CPCG < pdcoO
    
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
        function o = K2_CPCG(slack,options)
            o = o@pdcoO(slack,options);
            o.diagHess = false;
        end
        
        function Solve_Newton(o)
            o.H = o.H + sparse(o.low, o.low, o.z1(o.low) ./ o.x1(o.low), o.n, o.n); 
            o.H = o.H + sparse(o.upp, o.upp, o.z2(o.upp) ./ o.x2(o.upp), o.n, o.n); 
            w = o.r2; 
            w(o.low) = w(o.low) - (o.cL(o.low) + o.z1(o.low) .* o.rL(o.low)) ./ o.x1(o.low);
            w(o.upp) = w(o.upp) + (o.cU(o.upp) + o.z2(o.upp) .* o.rU(o.upp)) ./ o.x2(o.upp);
            
            if o.nfix > 0 && o.explicitA
                [ih, jh, vh] = find(o.H);
                for k = o.fix'
                    vh(ih == k & ih ~= jh) = 0;
                    vh(jh == k & ih ~= jh) = 0;
                end 
                o.H = sparse(ih, jh, vh); 
            end
            
            if o.PDitns == 1
                o.Afree = o.A;
                o.Afree(:, o.fix) = 0;
            end
            
            dH = spdiags(diag(o.H), 0, o.n, o.n);
            o.C = sparse(1:o.m, 1:o.m, o.d2.^2, o.m, o.m);

            if o.nfix == 0
                J = o.A;
            else
                J = o.Afree;
            end

            G  = [ dH  J'
                    J  -o.C ];
            % Full LDL of constraint preconditioner.
            o.P = opLDL(G, 0);
            
            rhs1 = [w; o.r1]; 
            rhs1(o.fix) = 0;
            
            % Shift rhs.
            % TODO: use alternative starting point.
            y0 = o.C \ rhs1(o.n+1:o.n+o.m);
            o.rhs = rhs1(1:o.n) - J' * y0;
            
            Solver(o);
            
            o.dx = o.sol;
            
            % Recover y.
            o.dy = o.C \ (J * o.dx);
            o.dy = y0 - o.dy;

            o.dx1(o.low) = -o.rL(o.low) + o.dx(o.low);
            o.dx2(o.upp) = -o.rU(o.upp) - o.dx(o.upp);
            o.dz1(o.low) = (o.cL(o.low) - o.z1(o.low) .* o.dx1(o.low)) ./ o.x1(o.low);
            o.dz2(o.upp) = (o.cU(o.upp) - o.z2(o.upp) .* o.dx2(o.upp)) ./ o.x2(o.upp);
        end
    end
end
