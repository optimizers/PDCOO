classdef PCG < handle
    properties
        itnlim
        atol1
        atol2
        atolmin
        atol
        r3ratio
        atolold
        precon
    end

    methods
        function o = PCG(options)

            if isfield(options, 'atol1')
                o.atol1 = options.atol1;
            else
                o.atol1 = 1e-5;
            end

            if isfield(options, 'atol2')
                o.atol2 = options.atol2;
            else
                o.atol2 = 1e-10;
            end

            if isfield(options, 'itnlim')
                o.itnlim = options.itnlim * min(o.m, o.n);
            else
                o.itnlim = 10 * min(o.m, o.n);
            end

            o.manage_op = true;
            o.need_precon = true;
            o.solver  = '  PCG';
            o.head3 = '  atol   PCG Inexact';

            % Parameters for LSMR and MINRES.
            o.atolmin = eps;    % Smallest atol if linesearch back-tracks
        end

        function Print_param(o)
            fprintf(o.file_id, '\n\nPCG:');
            fprintf(o.file_id, '\natol1    = %8.1e     atol2    = %8.1e', o.atol1, o.atol2);
            fprintf(o.file_id, '\n                     itnlim   = %8g', o.itnlim);
            fprintf(o.file_id, '\n\nMethod   = PCG\n');
        end

        function Init_param(o)
            o.atol = o.atol1;
            o.atol2 = max(o.atol2, o.atolmin);
            o.atolmin = o.atol2;
        end

        function Solver(o)
            % 31 Jan 2001: Set atol according to progress, a la Inexact Newton.
            % 07 Feb 2001: 0.1 not small enough for Satellite problem.  Try 0.01.
            % 25 Apr 2001: 0.01 seems wasteful for Star problem.
            %              Now that starting conditions are better, go back to 0.1.

            r3norm = max([o.Pinf, o.Dinf, o.Cinf]);
            o.atol   = min([o.atol, r3norm * 0.1]);
            o.atol   = max([o.atol, o.atolmin]);

            [o.sol, cg_flag, cg_relres, o.inner_itns] = ...
                pcg(o.M, o.rhs, o.atol, o.itnlim, o.precon);

            if cg_flag == 1 || cg_flag == 3 || cg_flag == 4   % conlim or itnlim
                fprintf('\n    CG failed to converge:  cg_flag = %3d', cg_flag)
            end

            o.atolold = o.atol;
            o.r3ratio = cg_relres * norm(o.rhs) / o.fmerit;
            o.inner_total = o.inner_total + o.inner_itns;

            % See if atol need to be reduced.
            if o.atol > o.atolmin
                if o.r3ratio >= 0.001     % Accept dy but make next one more accurate.
                    o.atol = o.atol * 0.1;
                end
            end
        end

        function Print_results(o)
            fprintf(o.file_id, ' %5.1f%7g%7.3f', log10(o.atolold), o.inner_itns, o.r3ratio);
        end

        function Reset_param(o)
            % Reduce atol for LSMR (and MINRES).
            % NOW DONE AT TOP OF LOOP.

            o.atolold = o.atol;
            % if atol > atol2
            %   atolfac = (mu/mufirst)^0.25;
            %   atol    = max( atol*atolfac, atol2 );
            % end

            % atol = min( atol, mu );     % 22 Jan 2001: a la Inexact Newton.
            % atol = min( atol, 0.5*mu ); % 30 Jan 2001: A bit tighter

            % If the linesearch took more than one function (nf > 1),
            % we assume the search direction needed more accuracy
            % (though this may be true only for LPs).
            % 12 Jun 1998: Ask for more accuracy if nf > 2.
            % 24 Nov 2000: Also if the steps are small.
            % 30 Jan 2001: Small steps might be ok with warm start.
            % 06 Feb 2001: Not necessarily.  Reinstated tests in next line.
            if o.nf > 2 || o.step <= 0.01
                o.atol = o.atolold * 0.1;
            end
        end
    end
end
