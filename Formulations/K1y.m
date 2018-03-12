classdef K1y < handle
    
    properties
        M
        rhs
        sol
    end
    
    methods (Abstract)
        Solver(o)
    end
    
    methods
        function o = K1y(options)
            o.diagHess = true;
        end
        
        function y = opK1y(x, ~)
            t = o.A'*x;
            t = t ./ o.H;
            t = o.A * t;
            y = t + o.d2.^2 .* x;
        end
        
        function Solve_Newton(o)
            %-----------------------------------------------------------------
            %  Solve (*) for dy.
            %-----------------------------------------------------------------
            %  Define a damped Newton iteration for solving f = 0,
            %  keeping  x1, x2, z1, z2 > 0.  We eliminate dx1, dx2, dz1, dz2
            %  to obtain the system
            %
            %  [-H2  A' ] [dx] = [w ],   H2 = H + D1^2 + X1inv Z1 + X2inv Z2,
            %  [ A  D2^2] [dy] = [r1]    w  = r2 - X1inv(cL + Z1 rL)
            %                                    + X2inv(cU + Z2 rU),
            %
            %  Then we eliminate dx to obtain :
            %
            %     (A*H2^-1*A' + D2^2) dy = A*H2^-1*w + r1   (*)
            %-----------------------------------------------------------------
            % For this method to work, H must be diagonal
            o.H(o.low) = o.H(o.low) + o.z1(o.low) ./ o.x1(o.low);
            o.H(o.upp) = o.H(o.upp) + o.z2(o.upp) ./ o.x2(o.upp);
            o.H(o.fix) = Inf;
            w = o.r2;
            w(o.low) = w(o.low) - (o.cL(o.low) + o.z1(o.low) .* o.rL(o.low)) ./ o.x1(o.low);
            w(o.upp) = w(o.upp) + (o.cU(o.upp) + o.z2(o.upp) .* o.rU(o.upp)) ./ o.x2(o.upp);

            if o.need_precon
                if o.explicitA
                    D = sqrt(o.H);
                    AD = o.A*diag(sparse(D));
                    AD = AD.^2;
                    wD = sum(AD,2); % Sum of sqrs of each row of AD.  %(Sparse)
                    wD = sqrt(full(wD) + (o.d2.^2)); %(Dense)
                    o.precon = 1 ./ wD;
                    clear AD wD
                else
                    o.precon = ones(o.m, 1);
                end
                o.precon = diag(sparse(o.precon));
            end
            
            if o.explicitA
                o.M = o.A * sparse(1:o.n, 1:o.n, 1 ./ o.H, o.n, o.n) * o.A';
                o.M = o.M + sparse(1:o.m, 1:o.m, o.d2.^2, o.m, o.m);
            else
                o.M = opFunction(o.m, o.m, @opK1y);
            end
            
            o.rhs = o.A * (w ./ o.H) + o.r1;
            
            Solver(o);

            o.dy = o.sol;

            % dy is now known.  Get dx, dx1, dx2, dz1, dz2.
            o.grad = o.A' * o.dy;
            o.grad(o.fix) = 0; % Is this needed?      % grad is a work vector
            o.dx = (o.grad - w) ./ o.H;
            o.dx1(o.low) = -o.rL(o.low) + o.dx(o.low);
            o.dx2(o.upp) = -o.rU(o.upp) - o.dx(o.upp);
            o.dz1(o.low) = (o.cL(o.low) - o.z1(o.low) .* o.dx1(o.low)) ./ o.x1(o.low);
            o.dz2(o.upp) = (o.cU(o.upp) - o.z2(o.upp) .* o.dx2(o.upp)) ./ o.x2(o.upp);
        end
    end
end
