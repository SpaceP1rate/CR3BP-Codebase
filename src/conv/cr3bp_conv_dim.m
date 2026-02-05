function [X_dim, T_dim] = cr3bp_conv_dim(X_nd, t_nd,a, method, value)
% CR3BP_CONV_DIM: Converts nondimensional CR3BP state and time vectors to 
%                 dimensional units.
%
%   [X_dim, T_dim] = cr3bp_conv_dim(X_nd, t_nd, a, method, value)
%           converts nondimensional CR3BP state and time vectors into 
%           dimensional  units based on the specified dimensionalization 
%           method.
%
% -------------------------------------------------------------------------
% Inputs:
%   X_nd    - state vector(s) in rotating frame
%             Planar case: [x, y, vx, vy]         (n x 4)
%             Spatial case: [x, y, z, vx, vy, vz] (n x 6)
%             Rows correspond to different times
%   t_nd    - nondimensional time vector (n x 1)
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
%   X_dim   - dimensional state vector(s) in SI units
%             Planar case: [x, y, vx, vy]         (m, m/s)
%             Spatial case: [x, y, z, vx, vy, vz] (m, m/s)
%   T_dim   - dimensional time vector [s]
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

    [~, nCols] = size(X_nd);
    
    % Unpack State vector
    x = X_nd(:,1);
    y = X_nd(:,2);
    
    % Scaling Coefficients

    L = a;       % position scaling
    V = L * n;   % velocity scaling
    T = 1/n;     % temporal scaling

    if nCols == 6
        z  = X_nd(:,3);
        vx = X_nd(:,4);
        vy = X_nd(:,5);
        vz = X_nd(:,6);

        % Dimensionalization Spatial Case
        X_dim(:,1) = L * x;
        X_dim(:,2) = L * y;
        X_dim(:,3) = L * z;
        X_dim(:,4) = V * vx;
        X_dim(:,5) = V * vy;
        X_dim(:,6) = V * vz;
        T_dim = T * t_nd;

    elseif nCols == 4
       vx = X_nd(:,3);
       vy = X_nd(:,4);

        % Dimensionalization Planar Case
        X_dim(:,1) = L * x;
        X_dim(:,2) = L * y;
        X_dim(:,3) = V * vx;
        X_dim(:,4) = V * vy;
        T_dim = T * t_nd;
    else 
        error('State vector must have 4 columns (planar) or 6 columns (spatial).');
    end

end