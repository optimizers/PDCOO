classdef QR < handle
    
    properties
        R
        P
    end
    
    methods
        function o = QR(options)
            
            o.manage_op = false;
            o.need_precon = false;
            o.solver  = '    QR';
            o.head3 = '       QR';
        end
        
        function Print_param(o)
            fprintf(o.file_id, '\n\nMethod   = QR\n');
        end
        
        function Init_param(~)
        end
        
        function Solver(o)
            if o.PDitns == 1, o.P = colamd(o.M); end % Do ordering only once.
            
            % sol = M \ rhs;
            [rhs, o.R] = qr(o.M(:,o.P), o.rhs, 0);
            o.sol = o.R \ rhs;
            o.sol(o.P) = o.sol;
        end
        
        function Print_results(o)
            if o.PDitns == 1, fprintf(o.file_id, ' %8g', nnz(o.R)); end
        end
        
        function Reset_param(~)
        end
    end
end
