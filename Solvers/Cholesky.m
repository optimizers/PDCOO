classdef Cholesky < handle

    properties
        R
        P
    end

    methods
        function o = Cholesky(slack, options)
            o.manage_op = false;
            o.need_precon = false;
            o.solver = '  Chol';  o.head3 = '     Chol';
        end

        function Print_param(o)
            fprintf(o.file_id, '\n\nMethod   = Cholesky\n');
        end

        function Init_param(~)
        end

        function Solver(o)

            if o.PDitns==1, o.P = symamd(o.M); end % Do ordering only once.

            [o.R, indef] = chol(o.M(o.P, o.P));
            if indef
                fprintf(o.file_id, '\n\n   chol says M is not pos def');
                fprintf(o.file_id, '\n   Use bigger d2, or choose another method');
                o.inform = 4;
            end

            % sol = M \ rhs;
            o.sol = o.R \ (o.R' \ o.rhs(o.P));
            o.sol(o.P) = o.sol;
        end

        function Print_results(o)
            if o.PDitns == 1, fprintf(o.file_id, ' %8g', nnz(o.R)); end
        end

        function Reset_param(~)
        end
    end
end
