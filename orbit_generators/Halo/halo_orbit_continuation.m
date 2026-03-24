function halo_orbit_continuation(mu,LP,type,num_orbits)

    ds = 1e-4;
    
    % Lyapunov Orbit Generator
    deltaX = 0.02;  
    [X0,~,~] = gen_lyapunov(deltaX,LP,mu);
    X0L= X0;
    
    % Halo Seed
    X0h = X0L;
    X0h(3) = 1e-3;
    
    % Halo Orbit Generator
    [Xh0, T_halfh, nu] = gen_halo_z(X0h, type, LP, mu);
    
    X0H(1,:) = Xh0;
    T_Hfam(1) = 2 * T_halfh;
    nu_fam(1) = nu;
    
    Y0 = [X0H'; reshape(eye(6), [], 1)];
    [~, YE] = ode78(@(t,Y) cr3bp_eom_stm(t,Y,mu), [0 T_Hfam], Y0);
    Phi   = reshape(YE(end, 7:end),6,6);
    
    Xddot = cr3bp_sys_accel(YE(end,1:6),mu);
    
    J = [Phi(2,1),Phi(2,3), Phi(2,5), YE(end,5)
         Phi(4,1),Phi(4,3), Phi(4,5), Xddot(1);
         Phi(6,1),Phi(6,3), Phi(6,5), Xddot(3)];
    
    Xstar = null(J);
    
    switch type
        case 'south'
            if Xstar(2)>0 
                Xstar(2) = -Xstar(2);
            end
        case 'north'
            if Xstar(2)<0 
                Xstar(2) = -Xstar(2);
            end
    end
    
    Xh0_pac = Xh0;
    Xh0_pac(1) = Xh0(1) + Xstar(1)*ds;
    Xh0_pac(3) = Xh0(3) + Xstar(2)*ds;
    Xh0_pac(5) = Xh0(5) + Xstar(3)*ds;
    T = T_Hfam + Xstar(4)*ds;
    alpha = 1;
    tic
    
    for i=2:num_orbits
        [Xh0, T_halfh, nu, Xstar_new, flag] = gen_halo_pac(Xh0_pac,Xh0,Xstar,T,T_Hfam(i-1),ds,type,LP,mu,alpha);
        X0H(i,:) = Xh0;
        T_Hfam(i) = 2 * T_halfh;
        nu_fam(i) = nu;
    
        if abs(Xh0(1)/(1-mu) - 1)*100 < 15
            alpha = 1;  
        end
        fprintf("Alpha: %f\n\n",alpha)
    
        Xstar = Xstar_new;
        Xh0_pac = Xh0;
        Xh0_pac(1) = Xh0(1) + Xstar(1)*ds;
        Xh0_pac(3) = Xh0(3) + Xstar(2)*ds;
        Xh0_pac(5) = Xh0(5) + Xstar(3)*ds;
        T       = T_Hfam(i) + Xstar(4)*ds;
    
        if ~flag
            fprintf("Stopping Continuation due to Failed Convergence\n\n")
            break;   
        end
    
    end
    
    family.([LP,'_',type]).IC = X0H;
    family.([LP,'_',type]).T = T_Hfam;

    fprintf("Continuation took %.2f minutes!\n\n",toc/60)
    save([type,'_halo_',LP,'_ic.mat'],"family")
end
