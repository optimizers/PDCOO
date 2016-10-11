classdef K2_ls < pdcoO
    
    properties
        M
        N
        rhs
        sol
    end
    
    methods (Abstract)
        Solver(o)
    end
    
    methods
        function o = K2_ls(slack, options)
            o = o@pdcoO(slack, options);
            o.diagHess = true;
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
            %  which is equivalent to the least-squares formultion
            %  [H2    A  ][-dx] = [rhs]  rhs = w - A'*D2^-2 * r1
            %  [A'  -D2^2][ dy] = [ 0 ]
            %
            %-----------------------------------------------------------------
            % For this method to work, H must be diagonal
            
            o.H(o.low) = o.H(o.low) + o.z1(o.low) ./o.x1 (o.low);
            o.H(o.upp) = o.H(o.upp) + o.z2(o.upp) ./o.x2 (o.upp);
            o.H(o.fix) = Inf;
            w = o.r2;
            w(o.low) = w(o.low) - (o.cL(o.low) + o.z1(o.low) .* o.rL(o.low)) ./ o.x1(o.low);
            w(o.upp) = w(o.upp) + (o.cU(o.upp) + o.z2(o.upp) .* o.rU(o.upp)) ./ o.x2(o.upp);
            
            o.M = opFunction(o.n, o.n, @(x,mode) x ./ o.H);
            o.N = opFunction(o.m, o.m, @(x,mode) x ./ (o.d2.^2));
            
            dy0 = o.N * o.r1;
            o.rhs = w - o.A' * dy0;
            
            Solver(o);

            o.dy = o.sol + dy0;
            % dy is now known.  Get dx, dx1, dx2, dz1, dz2.

            o.grad = o.A' * o.dy;
            o.grad(o.fix) = 0; % Is this needed?
            o.dx = (o.grad - w) ./ o.H;
            o.dx1(o.low) = -o.rL(o.low) + o.dx(o.low);
            o.dx2(o.upp) = -o.rU(o.upp) - o.dx(o.upp);
            o.dz1(o.low) = (o.cL(o.low) - o.z1(o.low) .* o.dx1(o.low)) ./ o.x1(o.low);
            o.dz2(o.upp) = (o.cU(o.upp) - o.z2(o.upp) .* o.dx2(o.upp)) ./ o.x2(o.upp);
        end
    end
end
