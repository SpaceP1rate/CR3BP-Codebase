function [Corrected_IC, T_half, varargout] = gen_halo(IC,type,LP, mu)

    nStates = 6;
    tol = 1e-12;
    maxIter = 500; 
    opts = odeset('Events', @ode_event_xcross, 'RelTol', 1e-5, 'AbsTol', 1e-5);
    
    x_start = IC(1);
    vy0 = IC(5);
    
    switch type
        case 'north'
            switch LP
                case 'l1'
                z_direction = -1;
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

    z0 = z_direction * IC(3);    
    tf = 5;

    fprintf('Shooting for %s Halo orbit around %s.\n%s\n',type,LP,repmat('-', 1, 46))
    for iter = 1:maxIter
    
        X0 = [x_start; 0; z0; 0; vy0; 0];
        
        Y0 = [X0; reshape(eye(nStates), [], 1)];
        [~ ,~ , TE, YE, ~] = ode45(@(t,Y) cr3bp_eom_stm(t,Y,mu), [0 tf], Y0, opts);
        
        X_end = YE(end, 1:6);
        Phi   = reshape(YE(end, nStates+1:end), nStates, nStates);
        

        F = [X_end(4);X_end(6)]; 


        J = [Phi(4,1), Phi(4,5);
             Phi(6,1), Phi(6,5)];
        
        delta = -J \ F;
        
        x_start = x_start + 0.3*delta(1);
        vy0 = vy0 + 0.3*delta(2);

        if mod(iter, 5) == 0 || iter == 1
            fprintf('Iter %d: x0 = %f, vy0 = %f, error vx = %e, error vz = %e\n',...
                iter, x_start,vy0, abs(F(1)),abs(F(2)));
        end
    
        if abs(F(2)) < tol
            fprintf('%s\nConverged in %d iterations.\n\n',repmat('-', 1, 44), iter);
            break;
        end

        if iter == maxIter
            fprintf('%s\nNot Converged after max iterations.\n\n',repmat('-', 1, 44));
        end
    end
    
    Corrected_IC = [x_start, 0, z0, 0, vy0, 0];
    T_half = TE(end);
    
    switch nargout
        case 3
            varargout{1} = cr3bp_sys_jacconst(Corrected_IC, mu);
        case 4
            varargout{1} = cr3bp_sys_jacconst(Corrected_IC, mu); % Jacobi
            varargout{2} = 0.5 * trace(Phi([3 6],[3 6])); % 0.5 trace(Phi_z) or nu
    end
end