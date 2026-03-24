function [Corrected_IC, T_half, nu] = gen_halo_z(IC,type,LP, mu)

    nStates = 6;
    tol = 1e-8;
    maxIter = 500; 
    opts = odeset('Events', @ode_event_xcross, 'RelTol', 1e-8, 'AbsTol', 1e-8);
    
    x_start = IC(1);
    vy0 = IC(5);
    
    switch type
        case 'south'
            z0 = -IC(3); 
        case 'north'
            z0 = IC(3); 
    end
       
    tf = 5;

    fprintf('Shooting for %s Halo orbit around %s | Fixed Variable z0.\n%s\n',type,LP,repmat('-', 1, 89))
    for iter = 1:maxIter
    
        X0 = [x_start; 0; z0; 0; vy0; 0];
        
        Y0 = [X0; reshape(eye(nStates), [], 1)];
        [~ ,~ , TE, YE, ~] = ode78(@(t,Y) cr3bp_eom_stm(t,Y,mu), [0 tf], Y0, opts);
        
        X_end = YE(end, 1:6);
        Phi   = reshape(YE(end, nStates+1:end), nStates, nStates);
        

        F = [X_end(4);X_end(6)]; 


        J = [Phi(4,1), Phi(4,5);
             Phi(6,1), Phi(6,5)];
        
        delta = -pinv(J)*F;

        x_start = x_start + 0.1*delta(1);
        vy0 = vy0 + 0.1*delta(2);

        if mod(iter, 5) == 0 || iter == 1
            fprintf('Iter %d: x0 = %f, vy0 = %f, error vx = %e, error vz = %e\n',...
                iter, x_start,vy0, abs(F(1)),abs(F(2)));
        end
    
        if abs(norm(F)) < tol
            fprintf('%s\nConverged in %d iterations.\n\n',repmat('-', 1, 89), iter);
            break;
        end

        if iter == maxIter
            fprintf('%s\nNot Converged after max iterations.\n\n',repmat('-', 1, 89));
        end
    end
    
    Corrected_IC = [x_start, 0, z0, 0, vy0, 0];
    T_half = TE(end);
    
    Y0 = [Corrected_IC'; reshape(eye(nStates), [], 1)];

    [~, Yf] = ode89(@(t,X) cr3bp_eom_stm(t,X,mu), [0,2*T_half], Y0');
    
    M = reshape(Yf(end, nStates+1:end), nStates, nStates);

    [~,D] = eig(M);

    lambda_max = max(diag(D));

    nu = 0.5*(lambda_max + 1/lambda_max);

end