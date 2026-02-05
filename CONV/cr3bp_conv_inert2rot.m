function X_rot = cr3bp_conv_inert2rot(X_inertial, t)
% CR3BP_CONV_INERT2ROT: Converts CR3BP states from inertial to 
% rotating frame
%
%   X_rot = cr3bp_conv_inert2rot(X_inertial, t)
%
%--------------------------------------------------------------------------
% Inputs:
%   X_inertial - state vector(s) in inertial frame
%                Planar case: [x, y, vx, vy]         (n x 4)
%                Spatial case: [x, y, z, vx, vy, vz] (n x 6)
%                Rows correspond to different times
%   t          - time(s), same length as rows in X_rot
%--------------------------------------------------------------------------
% Outputs:
%   X_rot - state vector(s) in rotating frame
%           same dimension as X_inertial
%--------------------------------------------------------------------------
% Notes:
%   - Assumes CR3BP normalized units
%   - Angular velocity of rotating frame = 1 rad/time unit
%--------------------------------------------------------------------------
% Version: 1.0
% Date:    August 13, 2025
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

    [nRows, nCols] = size(X_inertial);

    if nCols ~= 4 && nCols ~= 6
        error('State vector must have 4 columns (planar) or 6 columns (spatial).');
    end

    if numel(t) ~= nRows
        error('Length of time vector must match number of state rows.');
    end

    % Unpack rotating frame states
    x  = X_inertial(:, 1);
    y  = X_inertial(:, 2);

    if nCols == 6
        z  = X_inertial(:, 3);
        vx = X_inertial(:, 4);
        vy = X_inertial(:, 5);
        vz = X_inertial(:, 6);
    else
        z  = zeros(nRows, 1);
        vx = X_inertial(:, 3);
        vy = X_inertial(:, 4);
        vz = zeros(nRows, 1);
    end

    % Rotation angle (omega = 1 in normalized units)
    theta = t; % since omega = 1, theta = t

    % Precompute sin/cos
    cth = cos(theta);
    sth = sin(theta);

    % Transform positions
    X_rot =  cth .* x + sth .* y;
    Y_rot =  -sth .* x + cth .* y;
    Z_rot =  z; % unchanged in rotation

    % Transform velocities

    VX_rot =  (cth .* vx + sth .* vy) + Y_rot;
    VY_rot =  (-sth .* vx + cth .* vy) - X_rot;
    VZ_rot =  vz; % unchanged in rotation

    % Assemble output
    if nCols == 6
        X_rot = [X_rot, Y_rot, Z_rot, VX_rot, VY_rot, VZ_rot];
    else
        X_rot = [X_rot, Y_rot, VX_rot, VY_rot];
    end
end
