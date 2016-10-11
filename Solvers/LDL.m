classdef LDL < pdcoO

    properties
        L
    end

    methods
        function o = LDL(slack, options)
            o = o@pdcoO(slack, options);
            o.manage_op = false;
            o.need_precon = false;
            o.solver  = '    LDL';
            o.head3 = '       LDL';
        end

        function Print_param(o)
            fprintf(o.file_id,'\n\nMethod   = LDL\n');
        end

        function Init_param(~)
        end

        function Solver(o)
            % Use MA57 via Matlab's sparse ldl.
            % MA57 uses multiple cores if available.
            thresh = 0;      % tells MA57 to keep its sparsity-preserving order
            [o.L, D, P, S] = ldl(o.M, thresh);
            if nnz(D) ~= size(o.M, 1);
                error('[L, D, P, S] = ldl(K, 0) gave non-diagonal D')
            end
            o.sol = S * (P * (o.L' \ (D \ (o.L \ (P' * (S * o.rhs))))));
        end

        function Print_results(o)
            if o.PDitns == 1, fprintf(o.file_id, ' %8g x2', nnz(o.L)); end
        end

        function Reset_param(~)
        end
    end
end
