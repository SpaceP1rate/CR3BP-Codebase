function dYdt = cr3bp_eom_stm(~, Y, mu)
% CR3BP_EOM_STM: Computes time derivatives for planar or spatial CR3BP and propagates STM.
%
%   dYdt = cr3bp_eom_stm(t, Y, mu)
%
% -------------------------------------------------------------------------
% Inputs:
%     t  - time (unused, for ODE compatibility)
%     Y  - augmented state vector [X; Phi(:)]
%     X = state vector (4 or 6 elements)
%     Phi = STM matrix reshaped as vector (16 or 36 elements)
%     mu - mass parameter
% -------------------------------------------------------------------------
% Outputs:
%     dYdt - time derivative of augmented state [dX/dt; dPhi/dt(:)]
% -------------------------------------------------------------------------
% Notes:
%     - Units are normalized CR3BP units:
%     - Distance unit  = primary-secondary distance
%     - Time unit      = one rotation of secondary and primary about
%                        their baricenter
%     - Mass parameter = mu = m2 / (m1 + m2)
% -------------------------------------------------------------------------
% Version: 1.0
% Date:    August 12, 2025
%
% Version History:
%   1.0 - Initial release.
% -------------------------------------------------------------------------
% Copyright (c) 2025 Aleksandr Hakobyan
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

    n = length(Y);
    if n == 4 + 16
        nStates = 4;
    elseif n == 6 + 36
        nStates = 6;
    else
        error('Input length must be 20 (4+16) or 42 (6+36).');
    end

    % Extract state and STM
    X = Y(1:nStates);
    Phi = reshape(Y(nStates+1:end), nStates, nStates);

    % Unpack state
    x = X(1);
    y = X(2);
    if nStates == 6
        z = X(3);
        vx = X(4);
        vy = X(5);
        vz = X(6);
    else
        z = 0;
        vx = X(3);
        vy = X(4);
        vz = 0;
    end

    % Compute distances to primaries
    r1 = sqrt((x + mu)^2 + y^2 + z^2);
    r2 = sqrt((x - 1 + mu)^2 + y^2 + z^2);

    % Accelerations
    ax = 2*vy + x - (1 - mu)*(x + mu)/r1^3 - mu*(x - 1 + mu)/r2^3;
    ay = -2*vx + y - (1 - mu)*y/r1^3 - mu*y/r2^3;
    az = -(1 - mu)*z/r1^3 - mu*z/r2^3;

    % State derivatives
    if nStates == 6
        dXdt = [vx; vy; vz; ax; ay; az];
    else
        dXdt = [vx; vy; ax; ay];
    end

    % Construct Jacobian matrix A
    A = zeros(nStates);

    % Top-right block identity (position derivatives wrt velocity)
    A(1:nStates/2, nStates/2+1:end) = eye(nStates/2);

    % Precompute powers
    r1_3 = r1^3; r1_5 = r1^5;
    r2_3 = r2^3; r2_5 = r2^5;

    dx1 = x + mu;
    dx2 = x - 1 + mu;

    if nStates == 4
        Uxx = 1 - (1 - mu)*(1/r1_3 - 3*dx1^2/r1_5) - mu*(1/r2_3 - 3*dx2^2/r2_5);
        Uyy = 1 - (1 - mu)*(1/r1_3 - 3*y^2/r1_5) - mu*(1/r2_3 - 3*y^2/r2_5);
        Uxy = 3*((1 - mu)*dx1*y/r1_5 + mu*dx2*y/r2_5);

        A(3,1) = Uxx; A(3,2) = Uxy;
        A(4,1) = Uxy; A(4,2) = Uyy;

        A(3,4) = 2; A(4,3) = -2;

    else
        Uxx = 1 - (1 - mu)*(1/r1_3 - 3*dx1^2/r1_5) - mu*(1/r2_3 - 3*dx2^2/r2_5);
        Uyy = 1 - (1 - mu)*(1/r1_3 - 3*y^2/r1_5) - mu*(1/r2_3 - 3*y^2/r2_5);
        Uzz = - (1 - mu)*(1/r1_3 - 3*z^2/r1_5) - mu*(1/r2_3 - 3*z^2/r2_5);

        Uxy = 3*((1 - mu)*dx1*y/r1_5 + mu*dx2*y/r2_5);
        Uxz = 3*((1 - mu)*dx1*z/r1_5 + mu*dx2*z/r2_5);
        Uyz = 3*((1 - mu)*y*z/r1_5 + mu*y*z/r2_5);

        A(4,1) = Uxx; A(4,2) = Uxy; A(4,3) = Uxz;
        A(5,1) = Uxy; A(5,2) = Uyy; A(5,3) = Uyz;
        A(6,1) = Uxz; A(6,2) = Uyz; A(6,3) = Uzz;

        A(4,5) = 2; A(5,4) = -2;
    end

    % STM derivative
    dPhidt = A * Phi;
    dPhi_vec = reshape(dPhidt, nStates^2, 1);

    % Combine derivatives
    dYdt = [dXdt; dPhi_vec];

end
