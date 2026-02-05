function dUdx = cr3bp_dUdx(x, y, z, mu)
% CR3BP_DUDX  Partial derivative of the CR3BP pseudo-potential w.r.t. x
%
%   dUdx = cr3bp_dUdx(x, y, z, mu)
%
% Inputs:
%   x  - x-position in rotating frame
%   y  - y-position 
%   z  - z-position 
%   mu - mass parameter (0 < mu < 0.5)
%
% Output:
%   dUdx - partial derivative dU/dx
%
% Normalized CR3BP units

    r1 = sqrt((x + mu)^2       + y^2 + z^2);
    r2 = sqrt((x - 1 + mu)^2   + y^2 + z^2);

    dUdx = x - (1 - mu) * (x + mu)   / r1^3 ...
         - mu       * (x - 1 + mu)   / r2^3;

end