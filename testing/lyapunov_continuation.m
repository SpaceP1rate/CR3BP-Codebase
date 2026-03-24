clear; close; clc;
%% Initialization
addpath(genpath('..\src'))

mu = 0.0121505;

LP = 'l2';
num_orbits = 2000;
ds = 1e-3;
deltaX = 0.0001; 
opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

L = cr3bp_sys_lagpoints(mu, LP);
LP_Jac = cr3bp_sys_jacconst([L(1), 0, 0, 0, 0, 0], mu);

%% Lyapunov Orbit Generator
 
[Xl0,T_half,nu] = gen_lyapunov(deltaX,LP,mu);
X0L(1,:)= Xl0;
T_Lfam(1) = 2*T_half;
nu_fam(1) = nu;

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
    [Xl0, T_halfl, nu, Xstar_new, flag] = gen_lyapunov_pac(Xl0_pac,Xl0,Xstar,T,T_Lfam(i-1),ds,LP,mu);
    X0L(i,:) = Xl0;
    T_Lfam(i) = 2 * T_halfl;
    nu_fam(i,:) = nu;

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
gen_num = i;
fprintf("Continuation took %.2f minutes!\n\n",toc/60)
%% Plotting
clf;
figure(1);hold on
scatter3(1-mu,0,0,'ko','filled','DisplayName','Moon')
scatter3(L(1),0,0,'gx','DisplayName','L_1')
for i = 1:50:gen_num
    [tl, Xl] = ode78(@(t,X) cr3bp_eom_gen(t,X,mu), [0,T_Lfam(i)], X0L(i,:)',opts);
    p2=plot3(Xl(:,1), Xl(:,2),Xl(:,3),'LineWidth',1.2,'DisplayName','Lyapunov','Color','r');
end

axis equal
xlim([1-mu-0.3,1-mu+0.3])
ylim([-0.3,0.3])
hold off; grid on;
xlabel('x'); ylabel('y');
title('Periodic Orbits in CR3BP');

figure(2)
plot(X0L(:,1),nu_fam)