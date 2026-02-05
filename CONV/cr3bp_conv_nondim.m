function [X_nd, T_nd] = cr3bp_conv_nondim(X_dim, t_dim,a, method, value)
% CR3BP_CONV_DIM: Converts dimensional CR3BP state and time into
%                 nondimensional  vectors.
%
%   [X_dim, T_nd] = cr3bp_conv_nondim(X_dim, t_dim, a, method, value)
%           converts dimensional CR3BP state and time vectors into 
%           nondimensional  units based on the specified 
%           nondimensionalization method.
%
% -------------------------------------------------------------------------
% Inputs:
%   X_dim   - dimensional state vector(s) in rotating frame [m]
%             Planar case: [x, y, vx, vy]         (n x 4)
%             Spatial case: [x, y, z, vx, vy, vz] (n x 6)
%             Rows correspond to different times
%   t_dim   - dimensional time vector (n x 1)               [s]
%   a       - Distance between primary and secondary bodies [m]
%   method  - String specifying dimensionalization method:
%               'Synodic Period' - uses synodic period value
%               'Primary Mass'   - uses primary mass and mass ratio mu
%               'Secondary Mass' - uses secondary mass and mass ratio mu
%   value   - Numeric value corresponding to chosen method:
%               - 'Synodic Period': Psyn     [s]
%               - 'Primary Mass'  : [m1, mu] [kg, nondim (0 < mu < 0.5)]
%               - 'Secondary Mass': [m2, mu] [kg, nondim (0 < mu < 0.5)]
% -------------------------------------------------------------------------
% Outputs:
%   X_dim    - nondimensional state vector(s) 
%             Planar case: [x, y, vx, vy]        
%             Spatial case: [x, y, z, vx, vy, vz] 
%   T_nd    - dimensional time vector [s]
% -------------------------------------------------------------------------
% Notes:
%   - Supports both planar (4-element) and spatial (6-element) states.
%   - Uses gravitational constant G = 6.67430e-11 m^3/kg/s^2.
%   - Positions are scaled by 'a', velocities by L*n, and time by 1/n,
%     where n is computed from the chosen method.
% -------------------------------------------------------------------------
% Version: 1.0 | Date: August 15, 2025
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

    G = 6.67430e-11;     % Gravitational Constant

    switch method
        case 'Synodic Period'
            n = (2 * pi) / value;

        case 'Primary Mass'
            n = sqrt(G * value(1) / ((1-value(2)) * a^3)); %n = sqrt((G*m1/((1-mu)*a^3))) value = [m1, mu]

        case 'Secondary Mass'
            n = sqrt(G * value(1) / (value(2) * a^3));     %n = sqrt((G*m2/(mu*a^3))) value = [m2, mu]
    end

    [~, nCols] = size(X_dim);
    
    % Unpack State vector
    x = X_dim(:,1);
    y = X_dim(:,2);
    
    % Scaling Coefficients

    L = a;       % position scaling
    V = L * n;   % velocity scaling
    T = 1/n;     % temporal scaling

    if nCols == 6
        z  = X_dim(:,3);
        vx = X_dim(:,4);
        vy = X_dim(:,5);
        vz = X_dim(:,6);

        % Dimensionalization Spatial Case
        X_nd(:,1) = x / L;
        X_nd(:,2) = y / L;
        X_nd(:,3) = z / L;
        X_nd(:,4) = vx / V;
        X_nd(:,5) = vy / V;
        X_nd(:,6) = vz / V;
        T_nd = t_dim / T;

    elseif nCols == 4
       vx = X_dim(:,3);
       vy = X_dim(:,4);

        % Dimensionalization Planar Case
        X_nd(:,1) = x / L;
        X_nd(:,2) = y / L;
        X_nd(:,3) = vx / V;
        X_nd(:,4) = vy / V;
        T_nd = t_dim / T;
    else 
        error('State vector must have 4 columns (planar) or 6 columns (spatial).');
    end

end