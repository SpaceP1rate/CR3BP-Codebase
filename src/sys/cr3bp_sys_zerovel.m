function varargout = cr3bp_sys_zerovel(mu, J, bounds, res, option)
% CR3BP_SYS_ZEROVEL: Compute Zero-Velocity Curves/Surfaces and Forbidden 
%                    Regions in the CR3BP system.
%
%   [X, Y, Phi, Forbidden_mask] = cr3bp_sys_zerovel(mu, J, bounds, res)
%       computes the planar (z = 0) zero-velocity function Phi and the
%       forbidden region mask on a 2D grid.
%
%   [X, Y, Z, Phi, Forbidden_mask] = cr3bp_sys_zerovel(mu, J, bounds, res, 
%       'spatial') computes the spatial (3D) zero-velocity surface on a 
%       3D grid.
%
% -------------------------------------------------------------------------
% Inputs:
%   mu      - Mass ratio parameter (scalar double, 0 < mu < 0.5)
%   J       - Jacobi constant (scalar double)
%   bounds  - Vector of axis half-ranges:
%             planar: [xlim, ylim]
%             spatial: [xlim, ylim, zlim], the grid spans [-lim, +lim].
%   res     - Grid resolution:
%             scalar: n  -> creates [n n] (planar) or [N N N] (spatial)
%             vector: [nx ny] (planar) or [nx ny nz] (spatial)
%   option  - (optional) 'planar' (default) or 'spatial'
% -------------------------------------------------------------------------
% Outputs:
%   Planar call:
%       X, Y            - 2D meshgrid coordinates
%       Phi             - Zero-velocity function, Phi = 2*U - J
%       Forbidden_mask  - Logical mask of Phi < 0
%
%   Spatial call:
%       X, Y, Z         - 3D ndgrid coordinates
%       Phi             - Zero-velocity function, Phi = 2*U - J
%       Forbidden_mask  - Logical mask of Phi < 0
% -------------------------------------------------------------------------
% Notes:
%   - Assumes normalized CR3BP units in the rotating frame.
%   - Effective potential: U = (1-mu)/r1 + mu/r2 + 0.5*(x^2 + y^2).
%   - Zero-velocity curves/surfaces satisfy Phi = 0.
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
% -------------------------------------------------------------------------
    
    if nargin < 5
        option = 'planar';
    end
   
    option = lower(string(option));
    if ~(option == "planar" || option == "spatial")
        error('Invalid option. Choose ''planar'' or ''spatial''.');
    end

    % Handle scalar vs vector resolution
    if isscalar(res)
        switch option
            case 'spatial'
                res = [res res res];
            case 'planar'
            res = [res res];
        end
    end
    
    % Define grids
    x = linspace(-bounds(1), bounds(1), res(1));
    y = linspace(-bounds(2), bounds(2), res(2));

    switch option
        case 'planar'
        Z = zeros(res(2), res(1)); 
        [X, Y] = meshgrid(x, y);

        % Distances to primaries planar case
        r1 = sqrt((X + mu).^2 + Y.^2 + Z.^2);
        r2 = sqrt((X - 1 + mu).^2 + Y.^2 + Z.^2);

        case 'spatial'
        z = linspace(-bounds(3), bounds(3), res(3));
        [X, Y, Z] = ndgrid(x, y, z);
        
        % Distances to primaries spatial case
        r1 = sqrt((X + mu).^2 + Y.^2 + Z.^2);
        r2 = sqrt((X - 1 + mu).^2 + Y.^2 + Z.^2);
    end

    % Calculate Effective Potential
    U = (1 - mu) ./ r1 + mu ./ r2 + 0.5*(X.^2 + Y.^2);

    % Zero-Velocity function
    Phi = 2*U - J;

    % Masks
    Forbidden_mask = Phi < 0;          % true where forbidden

    switch option
        case 'planar'    
        varargout = {X, Y, Phi, Forbidden_mask};

        case 'spatial'
        varargout = {X, Y, Z, Phi, Forbidden_mask};
    end

end
