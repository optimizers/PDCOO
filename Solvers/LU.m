classdef LU < handle
    
    properties
        L
        pSQD
    end
    
    methods
        function o = LU(options)
            o.manage_op = false;
            o.need_precon = false;
            o.solver  = '    LU';
            o.head3 = '       LU';
        end
        
        function Print_param(o)
            fprintf(o.file_id, '\n\nMethod   = LU\n');
        end
        
        function Init_param(~)
        end
        
        function Solver(o)
            
            if o.PDitns == 1
                o.pSQD = symamd(o.M); 
            end % Do ordering only once.
            
            % Use Matlab's old sparse LU (Gilbert and Peierls) on K. 
            % Note that K is symmetric quasi-definite (SQD).
            % If delta isn't too small, we can suppress row permutations
            % and still have a sufficiently stable method.

            thresh = eps;        % eps ~= 2e-16 suppresses partial pivoting
            [o.L, U, perm] = lu(o.M(o.pSQD, o.pSQD), thresh, 'vector');
            if any(perm ~= 1:size(o.M, 1))  % We trust this gives perm=1:n
                error('[L, U, perm]=lu(K(pSQD, pSQD), ...) gave non-identity perm')
            end
            o.sol = U \ (o.L \ o.rhs(o.pSQD));
            o.sol(o.pSQD) = o.sol;
        end  
        
        function Print_results(o)
            if o.PDitns == 1, fprintf(o.file_id, ' %8g x2', nnz(o.L)); end
        end
        
        function Reset_param(~)
        end
    end
end
