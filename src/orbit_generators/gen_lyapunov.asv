function [Corrected_IC, T_half, varargout] = gen_lyapunov(C_target, LP, mu)
% GEN_LYAPUNOV: Generates a Planar Lyapunov orbit around L1/L2 for a 
% target Jacobi constant.
%
%   [Corrected_IC, T_half]              = gen_lyapunov(C_target, LP, mu)   
%   [Corrected_IC, T_half, Jacobi]      = gen_lyapunov(C_target, LP, mu)
%   [Corrected_IC, T_half, Jacobi, nu]  = gen_lyapunov(C_target, LP, mu)
%
% -------------------------------------------------------------------------
% Inputs:
%   C_target     - Target Jacobi constant (scalar double)
%   LP           - String identifying the Lagrange point ('l1', 'l2')
%   mu           - Mass ratio parameter (scalar double, 0 < mu < 0.5)
% -------------------------------------------------------------------------
% Outputs:
%   Corrected_IC - Corrected initial state [x, 0, 0, 0, vy, 0]
%   T_half       - Half of the orbital period (double)
%   Jacobi       - Final Jacobi constant of the converged orbit (double)
%   nu           - Vertical stability index of the converged orbit (double)
% -------------------------------------------------------------------------
% Notes:
%   - Uses a single-shooting differential correction scheme.
%   - Shoots to match Jacobi constant with desired one.
%   - Assumes normalized CR3BP units in the rotating frame.
% -------------------------------------------------------------------------
% Version: 1.2 | Created: January 24, 2026 | Modified: February 3, 2026 
%
% Version History:
%   1.0 - Initial release.
%
%   1.1 - Modified method to increase convergence, added dynamic scaling 
%         multiplier to solver update step tuned for each lagrange point.
%
%   1.2 - Added variable argument output, defaulting only to IC and Half
%         period. If prompted will also return the Jacobi constant and the
%         vertical stability index nu of the converged orbit.
% -------------------------------------------------------------------------
% Copyright (c) 2026 Aleksandr Hakobyan
% Laboratory of Spaceflight and Planetary Exploration
% Worcester Polytechnic Institute
%
% ACADEMIC USE LICENSE
%
% Permission is granted, free of charge, to use, copy, modify, and
% distribute this software for research and educational purposes only,
% subject to the following conditions:
%
%   1. This notice must appear in all copies or substantial portions of
%      the software.
%   2. The software may not be used for commercial purposes without prior
%      written permission from the copyright holder.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED.

L = cr3bp_sys_lagpoints(mu, LP);
xL = L(1);
nStates = 6;
tol = 1e-8;
maxIter = 15000; 
opts = odeset('Events', @ode_event_xcross, 'RelTol', 1e-5, 'AbsTol', 1e-5);

% Linearized analytical estimate for initial Vy guess
% Based on the characteristic equation of CR3BP linearized at L1/L2
gamma = abs(xL - (1-mu)); % Distance from secondary to LP
c2 = (mu/gamma^3) + ((1-mu)/(1-gamma)^3); % Approximate for L1/L2
lambda = sqrt(0.5*(c2 - 2 + sqrt(9*c2^2 - 8*c2)));
k = (lambda^2 + 1 + 2*c2) / (2*lambda);

switch LP
    case 'l1'
        x_start = xL + 0.02;   % Small initial perturbation
        v_direction = -1;      % Prograde
        tf = 4;
        vy0 = v_direction * abs(k * lambda * (x_start - xL));
        alpha = [0.37, 0.45];  % Dynamic scalar tuned for L1
    case 'l2'
        x_start = xL + 0.004;
        v_direction = -1;
        tf = 4;
        vy0 = v_direction * abs(k * lambda * (x_start - xL));
        alpha = [0.1245, 0.145]; % Dynamic scalar tuned for L2
end

fprintf('Shooting for Lyapunov orbit around %s.\n%s\n',LP,repmat('-', 1, 44))
for iter = 1:maxIter

    X0 = [x_start; 0; 0; 0; vy0; 0];
    Jacobi = cr3bp_sys_jacconst(X0', mu);
    
    Y0 = [X0; reshape(eye(nStates), [], 1)];
    [~, ~, TE, YE, ~] = ode45(@(t,Y) cr3bp_eom_stm(t,Y,mu), [0 tf], Y0, opts);
    
    X_end = YE(end, 1:6);
    Phi   = reshape(YE(end, nStates+1:end), nStates, nStates);
    
    Fx = X_end(4);               % vx(tf)
    Fc = Jacobi - C_target;      % Jacobi error
    
    Cx = 2 * cr3bp_sys_dUdx(x_start,0,0,mu);
    Cy = -2 * vy0;
    
    J = [ Phi(4,1), Phi(4,5);
          Cx,        Cy ];
    
    delta = -J \ [Fx; Fc];
    
    x_start = x_start + alpha(1)*delta(1);
    vy0     = vy0     + alpha(2)*delta(2);
    
    if mod(iter, 20) == 0 || iter == 1
        fprintf('Iter %d: x = %f, error vx = %e\n', iter, x_start, abs(Fx));
    end
    
    Corrected_IC = [x_start, 0, 0, 0, vy0, 0];
    Jacobi = cr3bp_sys_jacconst(Corrected_IC, mu);

    if abs(Fx) < tol
        fprintf('%s\nConverged in %d iterations.\n',repmat('-', 1, 44), iter);
        fprintf('Jacobi Constant Error: %.3e\n\n',abs(Jacobi-C_target));
        break;
    end

    if iter == maxIter
        fprintf('%s\nNot Converged after max iterations.\n',repmat('-', 1, 44));
        fprintf('Jacobi Constant Error: %.3e\n\n',abs(Jacobi-C_target));
    end
end

Corrected_IC = [x_start, 0, 0, 0, vy0, 0];
T_half = TE(end);

switch nargout
    case 3
        varargout{1} = cr3bp_sys_jacconst(Corrected_IC, mu);
    case 4
        varargout{1} = cr3bp_sys_jacconst(Corrected_IC, mu); % Jacobi
        varargout{2} = 0.5 * trace(Phi([3 6],[3 6])); % 0.5 trace(Phi_z) or nu
end

end