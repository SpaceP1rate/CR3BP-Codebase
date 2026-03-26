function [alpha, beta] = cr3bp_sys_stabilparams(Phi)

    alpha = 2 - trace(Phi);
    beta = 0.5*(alpha^2 + 2 - trace(Phi^2));

end