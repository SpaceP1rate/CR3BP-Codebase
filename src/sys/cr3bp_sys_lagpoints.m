function L = cr3bp_sys_lagpoints(mu, option)
% CR3BP_SYS_LAGPOINTS: Finds Lagrange Points in the CR3BP system.
%
%   L = cr3bp_sys_lagpoints(mu) returns all 5 Lagrange points [x,y].
%
%   L = cr3bp_sys_lagpoints(mu, option) returns points based on option:
%       'all'        - all five points (default)
%       'collinear'  - L1, L2, L3
%       'triangular' - L4, L5
%       'l1','l2','l3','l4','l5' - return single specified point
%
% -------------------------------------------------------------------------
% Inputs:
%   mu      - Mass ratio parameter (scalar double, 0 < mu < 0.5)
%   option  - (optional) String specifying which points to return,
%             default is 'all'
% -------------------------------------------------------------------------
% Outputs:
%   L       - (n x 2) matrix of [x, y] coordinates
%   version - (optional) version string of this function
% -------------------------------------------------------------------------
% Notes:
%   - Assumes normalized CR3BP units in the rotating frame.
%   - Collinear points found numerically, triangular points analytically.
% -------------------------------------------------------------------------
% Version: 1.0 | Date: August 12, 2025
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

    if nargin < 2
        option = 'all';
    end

    % Function for root-finding of collinear points along x-axis
    function val = collinearEq(x, mu)
        r1 = abs(x + mu);
        r2 = abs(x - (1 - mu));
        val = x - (1 - mu)*(x + mu)/r1^3 - mu*(x - (1 - mu))/r2^3;
    end

    % Compute all points so we can return individual easily
    % (only compute collinear once when needed)
    computeCollinear = @(mu) [ ...
        fzero(@(x) collinearEq(x, mu), [0.5 - mu, 1 - mu - 1e-6]);  % L1
        fzero(@(x) collinearEq(x, mu), [1 - mu + 1e-6, 2]);          % L2
        fzero(@(x) collinearEq(x, mu), [-2, -mu - 1e-6])             % L3
    ];
    % Triangular points (analytical)
    L4 = [0.5 - mu,  sqrt(3)/2];
    L5 = [0.5 - mu, -sqrt(3)/2];

    option = lower(option);

    switch option
        case 'collinear'
            colPts = computeCollinear(mu);
            L = [colPts(1), 0;
                 colPts(2), 0;
                 colPts(3), 0];

        case 'triangular'
            L = [L4; L5];

        case 'all'
            colPts = computeCollinear(mu);
            L = [colPts(1), 0;
                 colPts(2), 0;
                 colPts(3), 0;
                 L4;
                 L5];

        case 'l1'
            colPts = computeCollinear(mu);
            L = [colPts(1), 0];

        case 'l2'
            colPts = computeCollinear(mu);
            L = [colPts(2), 0];

        case 'l3'
            colPts = computeCollinear(mu);
            L = [colPts(3), 0];

        case 'l4'
            L = L4;

        case 'l5'
            L = L5;

        otherwise
            error('Invalid option. Choose ''all'', ''collinear'', ''triangular'', or ''l1''-''l5''.');
    end
end
