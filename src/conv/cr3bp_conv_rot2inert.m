function X_inertial = cr3bp_conv_rot2inert(X_rot, t)
% CR3BP_CONV_ROT2INERT: Converts CR3BP states from rotating to 
% inertial frame
%
%   X_inertial = cr3bp_conv_rot2inert(X_rot, t)
%
%--------------------------------------------------------------------------
% Inputs:
%   X_rot - state vector(s) in rotating frame
%           Planar case: [x, y, vx, vy]         (n x 4)
%           Spatial case: [x, y, z, vx, vy, vz] (n x 6)
%           Rows correspond to different times
%   t     - time(s), same length as rows in X_rot
%--------------------------------------------------------------------------
% Outputs:
%   X_inertial - state vector(s) in inertial frame
%                same dimension as X_rot
%--------------------------------------------------------------------------
% Notes:
%   - Assumes CR3BP normalized units
%   - Angular velocity of rotating frame = 1 rad/time unit
%--------------------------------------------------------------------------
% Version: 1.0 | Date:August 12, 2025
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

    [nRows, nCols] = size(X_rot);

    if nCols ~= 4 && nCols ~= 6
        error('State vector must have 4 columns (planar) or 6 columns (spatial).');
    end

    if numel(t) ~= nRows
        error('Length of time vector must match number of state rows.');
    end

    % Unpack rotating frame states
    x  = X_rot(:, 1);
    y  = X_rot(:, 2);

    if nCols == 6
        z  = X_rot(:, 3);
        vx = X_rot(:, 4);
        vy = X_rot(:, 5);
        vz = X_rot(:, 6);
    else
        z  = zeros(nRows, 1);
        vx = X_rot(:, 3);
        vy = X_rot(:, 4);
        vz = zeros(nRows, 1);
    end

    % Rotation angle (omega = 1 in normalized units)
    theta = t; % since omega = 1, theta = t

    % Precompute sin/cos
    cth = cos(theta);
    sth = sin(theta);

    % Transform positions
    X_pos =  cth .* x - sth .* y;
    Y_pos =  sth .* x + cth .* y;
    Z_pos =  z; % unchanged in rotation

    % Transform velocities
    % Inertial velocity = rotation_matrix * (v_rot + omega × r_rot)
    % Here omega = [0, 0, 1] in normalized units
    vxi_rot = vx - y; % add omega x r in rotating frame
    vyi_rot = vy + x;

    VX_inertial =  cth .* vxi_rot - sth .* vyi_rot;
    VY_inertial =  sth .* vxi_rot + cth .* vyi_rot;
    VZ_inertial =  vz; % unchanged in rotation

    % Assemble output
    if nCols == 6
        X_inertial = [X_pos, Y_pos, Z_pos, VX_inertial, VY_inertial, VZ_inertial];
    else
        X_inertial = [X_pos, Y_pos, VX_inertial, VY_inertial];
    end
end
