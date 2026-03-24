function acceleration = cr3bp_sys_accel(X,mu)
%CR3BP_SYS_ACCEL Summary of this function goes here
%   Detailed explanation goes here
    n = length(X); % detect state dimension
    
    if n ~= 4 && n ~= 6
        error('State vector must have 4 elements (planar) or 6 elements (spatial).');
    end
    
    % Unpack state
    x = X(1);
    y = X(2);
    if n == 6
        z = X(3);
        vx = X(4);
        vy = X(5);
    else
        vx = X(3);
        vy = X(4);
        z = 0;  % planar case
    end
    
    % Distances to primaries
    r1 = sqrt( (x + mu)^2 + y^2 + z^2 );
    r2 = sqrt( (x - 1 + mu)^2 + y^2 + z^2 );
    
    % Accelerations
    acceleration(1) = 2*vy + x - (1 - mu)*(x + mu)/(r1^3) - mu*(x - 1 + mu)/(r2^3);
    acceleration(2) = -2*vx + y - (1 - mu)*y/(r1^3) - mu*y/(r2^3);
    acceleration(3) = -(1 - mu)*z/(r1^3) - mu*z/(r2^3);
end

