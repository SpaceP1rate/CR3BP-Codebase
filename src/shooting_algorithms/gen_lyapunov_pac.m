function [Corrected_IC, T_half, nu, Z_new, Flag] = gen_lyapunov_pac(IC_guess, IC_prev, Z_prev,T_cur,T_prev, ds,LP,mu)

    nStates = 6;
    tol = 1e-8;
    maxIter = 500; 
    opts = odeset('Events', @ode_event_xcross, 'RelTol', 1e-8, 'AbsTol', 1e-8);
    
    x0  = IC_guess(1);
    vy0 = IC_guess(5);
    
    X_prev_vec = [IC_prev(1); IC_prev(5); T_prev];
    tf = T_cur;

    fprintf('Shooting for Lyapunov orbit around %s | Pseudo Arc-Length.\n%s\n',LP,repmat('-', 1, 46))
    for iter = 1:maxIter
    
        X0 = [x0; 0; 0; 0; vy0; 0];
        
        Y0 = [X0; reshape(eye(nStates), [], 1)];
        [~,~, TE, YE, ~] = ode89(@(t,Y) cr3bp_eom_stm(t,Y,mu), [0 tf], Y0, opts);
        
        X_f = YE(end, 1:6);
        Phi   = reshape(YE(end, nStates+1:end), nStates, nStates);
     
        Xddot = cr3bp_sys_accel(YE(end,1:6),mu);

        J = [Phi(2,1),Phi(2,5),YE(end,5)
             Phi(4,1),Phi(4,5),Xddot(1)];
        
        X_curr_vec = [x0; vy0; tf];

        G = [0;X_f(4);(X_curr_vec - X_prev_vec)'*Z_prev - ds];

        DG = [J; Z_prev'];

        delta = -pinv(DG)*G;
        
        x0 = x0 + 0.8*delta(1);
        vy0 = vy0 + 0.8*delta(2);
        tf = tf + 0.8*delta(3);

        if mod(iter, 5) == 0 || iter == 1
            fprintf('Iter %d: x0 = %f, vy0 = %f, error vx = %e\n', iter,x0, vy0, abs(G(2)));
        end
    
        if max(abs(G)) < tol
            fprintf('%s\nConverged in %d iterations.\n\n',repmat('-', 1, 105), iter);
            Flag = true;
            break;
        end

        if iter == maxIter
            fprintf('%s\nNot Converged after max iterations.\n\n',repmat('-', 1, 105));
            Flag = false;
        end
    end
    
    Corrected_IC = [x0, 0, 0, 0, vy0, 0];
    T_half = TE(end);

    Y0 = [Corrected_IC'; reshape(eye(nStates), [], 1)];

    [~, Yf] = ode78(@(t,X) cr3bp_eom_stm(t,X,mu), [0,2*T_half], Y0');
    
    M = reshape(Yf(end, nStates+1:end), nStates, nStates);

    [~,D] = eig(M);

    lambda_max = max(diag(D));

    nu = 0.5*(lambda_max + 1/lambda_max);
    
    Z_new = null(J); 
    if dot(Z_new, Z_prev) < 0, Z_new = -Z_new; end   
end