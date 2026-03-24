function [Corrected_IC, T_half, nu] = gen_halo_x(IC,type,LP, mu)

    nStates = 6;
    tol = 1e-10;
    maxIter = 500; 
    opts = odeset('Events', @ode_event_xcross, 'RelTol', 1e-10, 'AbsTol', 1e-12);
    
    x0 = IC(1);
    vy0 = IC(5);
    
    switch type
        case 'north'
            switch LP
                case 'l1'
                z_direction = 1;
                case 'l2'
                z_direction = -1;
            end
        case 'south'
            switch LP
                case 'l1'
                z_direction = 1;
                case 'l2'
                z_direction = 1;
            end
    end

    z_start = IC(3)*z_direction;    
    tf = 5;

    fprintf('Shooting for %s Halo orbit around %s.\n%s\n',type,LP,repmat('-', 1, 46))
    for iter = 1:maxIter
    
        X0 = [x0; 0; z_start; 0; vy0; 0];
        
        Y0 = [X0; reshape(eye(nStates), [], 1)];
        [~ ,~ , TE, YE, ~] = ode78(@(t,Y) cr3bp_eom_stm(t,Y,mu), [0 tf], Y0, opts);
        
        X_end = YE(end, 1:6);
        Phi   = reshape(YE(end, nStates+1:end), nStates, nStates);
        

        F = [X_end(4);X_end(6)]; 


        J = [Phi(4,3), Phi(4,5);
             Phi(6,3), Phi(6,5)];
        
        delta = -J \ F;

        z_start = z_start + 0.3*delta(1);
        vy0 = vy0 + 0.3*delta(2);

        if mod(iter, 5) == 0 || iter == 1
            fprintf('Iter %d: z0 = %f, vy0 = %f, error vx = %e, error vz = %e\n',...
                iter, z_start,vy0, abs(F(1)),abs(F(2)));
        end
    
        if abs(norm(F)) < tol
            fprintf('%s\nConverged in %d iterations.\n\n',repmat('-', 1, 44), iter);
            break;
        end

        if iter == maxIter
            fprintf('%s\nNot Converged after max iterations.\n\n',repmat('-', 1, 44));
        end
    end
    
    switch type
        case 'north'
            if z_start > 0
                Corrected_IC = [x0, 0, z_start, 0, vy0, 0];
            else
                Corrected_IC = [x0, 0, -z_start, 0, vy0, 0];
            end
        case 'south'
            Corrected_IC = [x0, 0, z_start, 0, vy0, 0];
    end

    T_half = TE(end);
    
    Y0 = [Corrected_IC'; reshape(eye(nStates), [], 1)];

    [~, Yf] = ode89(@(t,X) cr3bp_eom_stm(t,X,mu), [0,2*T_half], Y0');
    
    M = reshape(Yf(end, nStates+1:end), nStates, nStates);

    [~,D] = eig(M);

    lambda_max = max(diag(D));

    nu = 0.5*(lambda_max + 1/lambda_max);
end