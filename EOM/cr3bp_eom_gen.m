function dXdt = cr3bp_eom_gen(~, X, mu)
%CR3BP_EOM_GEN: Computes time derivatives for planar or spatial CR3BP
%
%   dXdt = cr3bp_eom_gen(t, X, mu)
%
%--------------------------------------------------------------------------
% Inputs:
%     t  - time (unused but kept for ODE compatibility)
%     X  - state vector
%          Planar case: [x, y, vx, vy]         (4x1)
%          Spatial case: [x, y, z, vx, vy, vz] (6x1)
%     mu - mass parameter (m2 / (m1 + m2))
%--------------------------------------------------------------------------
% Outputs:
%     dXdt - time derivative of state (same dimension as X)
%--------------------------------------------------------------------------
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

    n = length(X); % detect state dimension
    
    if n ~= 4 && n ~= 6
        error('State vector must have 4 elements (planar) or 6 elements (spatial).');
    end
    
    % Unpack state
    x = X(1);
    y = X(2);
    if n == 6
        z = X(3);
        vx = X(4);
        vy = X(5);
        vz = X(6);
    else
        vx = X(3);
        vy = X(4);
        z = 0;  % planar case
        vz = 0;
    end
    
    % Distances to primaries
    r1 = sqrt( (x + mu)^2 + y^2 + z^2 );
    r2 = sqrt( (x - 1 + mu)^2 + y^2 + z^2 );
    
    % Accelerations
    ax = 2*vy + x - (1 - mu)*(x + mu)/(r1^3) - mu*(x - 1 + mu)/(r2^3);
    ay = -2*vx + y - (1 - mu)*y/(r1^3) - mu*y/(r2^3);
    az = -(1 - mu)*z/(r1^3) - mu*z/(r2^3);
    
    % Assemble derivative vector
    if n == 6
        dXdt = [vx; vy; vz; ax; ay; az];
    else
        dXdt = [vx; vy; ax; ay];
    end
end
