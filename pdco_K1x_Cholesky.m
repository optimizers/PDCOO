classdef pdco_K1x_Cholesky < pdcoO & K1x & Cholesky
  properties
  end

  methods
  function o = pdco_K1x_Cholesky(slack, options_pdco, options_form, options_solv)
    o = o@pdcoO(slack, options_pdco);
    o = o@K1x(options_form);
    o = o@Cholesky(options_solv);
  end
  end
end
