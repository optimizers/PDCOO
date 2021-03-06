classdef lpmodel < model.qpmodel
    % nlpmodel where 
    %   objective is linear: f(x) = c' * x
    %   constraints are linear: c(x) = A  * x

    methods (Sealed = true)

        function o = lpmodel(name, x0, cL, cU, bL, bU, A, c)
            n = size(x0, 1);
            o = o@model.qpmodel(name, x0, cL, cU, bL, bU, A, c, sparse(n, n));
        end
    end
end
