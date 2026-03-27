function lyapunov_orbit_continuation(mu,LP,num_orbits)

    ds = 1e-3;
    deltaX = 0.0001; 
    
    % Lyapunov Orbit Generator
     
    [Xl0,T_half,stability] = gen_lyapunov(deltaX,LP,mu);
    X0L(1,:)= Xl0;
    T_Lfam(1) = 2*T_half;
    stab_fam(1,:) = stability;
    
    Y0 = [X0L'; reshape(eye(6), [], 1)];
    [~, YE] = ode78(@(t,Y) cr3bp_eom_stm(t,Y,mu), [0 T_Lfam], Y0);
    Phi   = reshape(YE(end, 7:end),6,6);
    
    Xddot = cr3bp_sys_accel(YE(end,1:6),mu);
    
    J = [Phi(2,1),Phi(2,5),YE(end,5)
         Phi(4,1),Phi(4,5),Xddot(1)];
    
    Xstar = null(J);
    
    Xl0_pac = Xl0;
    Xl0_pac(1) = Xl0(1) + Xstar(1)*ds;
    Xl0_pac(5) = Xl0(5) + Xstar(2)*ds;
    T          = T_Lfam + Xstar(3)*ds;
    
    tic
    for i=2:num_orbits
        [Xl0, T_halfl, stability, Xstar_new, flag] = gen_lyapunov_pac(Xl0_pac,Xl0,Xstar,T,T_Lfam(i-1),ds,LP,mu);
        X0L(i,:) = Xl0;
        T_Lfam(i) = 2 * T_halfl;
        stab_fam(i,:) = stability;
    
        Xstar = Xstar_new;
        Xl0_pac = Xl0;
        Xl0_pac(1) = Xl0(1) + Xstar(1)*ds;
        Xl0_pac(5) = Xl0(5) + Xstar(2)*ds;
        T          = T_Lfam(i) + Xstar(3)*ds;
        if ~flag
            fprintf("Stopping Continuation due to Failed Convergence\n\n")
            break;   
        end
    end

    family.(LP).IC = X0L;
    family.(LP).T = T_Lfam;
    family.(LP).params = stab_fam;

    fprintf("Continuation took %.2f minutes!\n\n",toc/60)
    save([LP,'_lyapunov_ic.mat'],"family")
end