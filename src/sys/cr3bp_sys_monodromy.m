function [M, vU, vS] = cr3bp_sys_monodromy(Y)
% CR3BP_MONODROMY_SU
% Compute monodromy and stable/unstable eigenvectors from augmented history.
%
% Input:
%   Y : N x (n + n^2) array [X(1:n), Phi(:)] over one full period
%       (n=4 for planar, n=6 for spatial CR3BP)
%
% Output:
%   M  : n x n monodromy matrix (STM over one period)
%   vU : unstable eigenvector (unit, real part)
%   vS : stable   eigenvector (unit, real part)

    % infer state dimension n
    m = size(Y,2);
    n = round((-1 + sqrt(1 + 4*m)) / 2);

    % monodromy
    M = reshape(Y(end, n+1:end), n, n);

    % eigendecomposition
    [V, D] = eig(M);
    lam = diag(D);
    
    % unstable = largest |λ|, stable = smallest |λ|
    [~, idxU] = max(abs(lam));
    [~, idxS] = min(abs(lam));

    vU = real(V(:, idxU)); vU = vU / norm(vU);
    vS = real(V(:, idxS)); vS = vS / norm(vS);
end
