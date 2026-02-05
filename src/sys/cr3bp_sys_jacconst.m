function J = cr3bp_sys_jacconst(X,mu)
% CR3BP_SYS_JACCONST: Calculates the Jacobi Constant in the CR3BP system.
%
%   J = cr3bp_sys_jacconst(X,mu) returns the Jacobi constant at every
%                                point in time.
%
% -------------------------------------------------------------------------
% Inputs:
%   X       - state vector(s) in rotating frame
%             Planar case: [x, y, vx, vy]         (n x 4)
%             Spatial case: [x, y, z, vx, vy, vz] (n x 6)
%             Rows correspond to different times
%   mu      - Mass ratio parameter (scalar double, 0 < mu < 0.5)
% -------------------------------------------------------------------------
% Outputs:
%   J       - (n x 1) vector of Jacobi Consant values
% -------------------------------------------------------------------------
% Notes:
%   - Assumes normalized CR3BP units in the rotating frame.
% -------------------------------------------------------------------------
% Version: 1.0 | Date: August 14, 2025
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

    [nRows, nCols] = size(X);

    % Unpacking the state vector
    x = X(:,1);
    y = X(:,2);
    
    if nCols == 6
       z  = X(:,3);
       vx = X(:,4);
       vy = X(:,5);
       vz = X(:,6);
    elseif nCols == 4
       z  = zeros(nRows, 1);
       vx = X(:,3);
       vy = X(:,4);
       vz = zeros(nRows, 1);
    else 
        error('State vector must have 4 columns (planar) or 6 columns (spatial).');
    end

    % Distances to the two primaries
    r1 = sqrt((x + mu).^2 + y.^2 + z.^2);
    r2 = sqrt((x - 1 + mu).^2 + y.^2 + z.^2);
    
    % Effective potential
    U = ((1 - mu) ./ r1) + (mu ./ r2) + (0.5 .* (x.^2 + y.^2));

    % Jacobi constant calculaltion
    J = 2 .* U - (vx.^2 + vy.^2 + vz.^2);

end