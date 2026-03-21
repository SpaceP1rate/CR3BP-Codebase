function [Corrected_IC, T_half, Z_new, Flag] = gen_halo_pac(IC_guess, IC_prev, Z_prev, ds,type,LP,mu,alpha)
    % IC_guess: Predictor [x, 0, z, 0, vy, 0]
    % T_guess:  Predictor half-period
    % Z_prev:   The 4x1 tangent vector from the previous orbit
    
    nStates = 6;
    tol = 1e-10;
    maxIter = 1000;
    opts = odeset('Events', @ode_event_xcross, 'RelTol', 1e-8, 'AbsTol', 1e-8);
   

    x0  = IC_guess(1);
    z0  = IC_guess(3);
    vy0 = IC_guess(5);
    
    X_prev_vec = [IC_prev(1); IC_prev(3); IC_prev(5)];
    tf = 5;
    fprintf('Shooting for %s Halo orbit around %s | Pseudo Arc-Length.\n%s\n',type,LP,repmat('-', 1, 105))
    for iter = 1:maxIter
     
        X0 = [x0; 0; z0; 0; vy0; 0];
        Y0 = [X0; reshape(eye(nStates), [], 1)];

        [~, ~, TE, YE, ~] = ode89(@(t,Y) cr3bp_eom_stm(t,Y,mu), [0 tf], Y0, opts);
        
        X_f = YE(end, 1:6);
        Phi = reshape(YE(end, nStates+1:end), nStates, nStates);
       

        J = [Phi(4,1),Phi(4,3), Phi(4,5);
             Phi(6,1),Phi(6,3), Phi(6,5)];
        
        X_curr_vec = [x0; z0; vy0];

        G = [X_f(4);X_f(6);(X_curr_vec - X_prev_vec)'*Z_prev - ds];
         
        DG = [J; Z_prev'];

        delta = -DG \ G;
        
        x0  = x0  + alpha*delta(1);
        z0  = z0  + alpha*delta(2);
        vy0 = vy0 + alpha*delta(3);

        if max(abs(G)) < tol
            fprintf('%s\nConverged in %d iterations.\n\n',repmat('-', 1, 105), iter);
            Flag = true;
            break;
        end

        if iter == maxIter
            fprintf('%s\nNot Converged after max iterations.\n\n',repmat('-', 1, 105));
            Flag = false;
        end

        if mod(iter, 5) == 0 || iter == 1
            fprintf('Iter %d: x0 = %f, z0 = %f, vy0 = %f, error vx = %e, error vz = %e\n',...
                    iter, x0,z0,vy0, abs(G(1)),abs(G(2)));
        end
    end
    
    Corrected_IC = [x0, 0, z0, 0, vy0, 0];
    T_half = TE(end);

    Z_new = null(J); 
    if dot(Z_new, Z_prev) < 0, Z_new = -Z_new; end
end