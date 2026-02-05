function dUdz = cr3bp_dUdz(x, y, z, mu)
% CR3BP_DUDZ  Partial derivative of the CR3BP pseudo-potential w.r.t. z
%
%   dUdz = cr3bp_dUdz(x, y, z, mu)
%
% Inputs:
%   x  - x-position in rotating frame
%   y  - y-position 
%   z  - z-position 
%   mu - mass parameter (0 < mu < 0.5)
%
% Output:
%   dUdz - partial derivative dU/dz
%
% Normalized CR3BP units

r1 = sqrt((x + mu)^2 + y^2 + z^2);
r2 = sqrt((x - 1 + mu)^2 + y^2 + z^2);

dUdz = - (1 - mu) * z / r1^3 ...
       - mu * z / r2^3;
end