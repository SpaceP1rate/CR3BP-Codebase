function [Corrected_IC, T_half, varargout] = gen_lyapunov(deltaX, LP, mu)

    L = cr3bp_sys_lagpoints(mu, LP);
    xL = L(1);
    nStates = 6;
    tol = 1e-11;
    maxIter = 500; 
    opts = odeset('Events', @ode_event_xcross, 'RelTol', 1e-5, 'AbsTol', 1e-5);
    
    % Linearized analytical estimate for initial Vy guess
    % Based on the characteristic equation of CR3BP linearized at L1/L2
    gamma = abs(xL - (1-mu)); % Distance from secondary to LP
    c2 = (mu/gamma^3) + ((1-mu)/(1-gamma)^3); % Approximate for L1/L2
    lambda = sqrt(0.5*(c2 - 2 + sqrt(9*c2^2 - 8*c2)));
    k = (lambda^2 + 1 + 2*c2) / (2*lambda);
    
    switch LP
        case 'l1'
            x_start = xL + deltaX;
            v_direction = -1; 
            tf = 3;
            vy0 = v_direction * abs(k * lambda * (deltaX));
        case 'l2'
            x_start = xL - deltaX;
           v_direction = 1; 
            tf = 3;
            vy0 = v_direction * abs(k * lambda * (deltaX));
    end
    
    fprintf('Shooting for Lyapunov orbit around %s.\n%s\n',LP,repmat('-', 1, 46))
    for iter = 1:maxIter
    
        X0 = [x_start; 0; 0; 0; vy0; 0];
        
        Y0 = [X0; reshape(eye(nStates), [], 1)];
        [T, Xtemp, TE, YE, ~] = ode45(@(t,Y) cr3bp_eom_stm(t,Y,mu), [0 tf], Y0, opts);
        
        X_end = YE(end, 1:6);
        Phi   = reshape(YE(end, nStates+1:end), nStates, nStates);
        
        dt= T(end) - T(end-1);
        yddot = gradient(Xtemp(:,5),dt);
    
        idx = find(T == TE(end));

        F = [0;X_end(4)];   % [0, delta_vx]
     
        J = [Phi(2,5), Xtemp(idx,5);
             Phi(4,5), yddot(idx)];
        
        delta = -J \ F;
        
        vy0 = vy0 + 0.6*delta(1);
        tf = tf + 0.6*delta(2);

        if mod(iter, 50) == 0 || iter == 1
            fprintf('Iter %d: vy = %f, error vx = %e\n', iter, vy0, abs(F(2)));
        end
    
        if abs(F(2)) < tol
            fprintf('%s\nConverged in %d iterations.\n',repmat('-', 1, 44), iter);
            break;
        end

        if iter == maxIter
            fprintf('%s\nNot Converged after max iterations.\n',repmat('-', 1, 44));
        end
    end
    
    Corrected_IC = [x_start, 0, 0, 0, vy0, 0];
    T_half = TE(end);
    
    switch nargout
        case 3
            varargout{1} = cr3bp_sys_jacconst(Corrected_IC, mu);
        case 4
            varargout{1} = cr3bp_sys_jacconst(Corrected_IC, mu); % Jacobi
            varargout{2} = 0.5 * trace(Phi([3 6],[3 6])); % 0.5 trace(Phi_z) or nu
    end
end