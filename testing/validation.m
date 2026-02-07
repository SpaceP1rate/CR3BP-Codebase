clear; close; clc;
%% Initialization
addpath(genpath('..\src'))

mu = 0.0121505;
LP = 'l2';
L = cr3bp_sys_lagpoints(mu, LP);

LP_Jac = cr3bp_sys_jacconst([L(1), 0, 0, 0, 0, 0], mu);

C_target = linspace(LP_Jac - 0.01,LP_Jac - 0.24,4);

%% Lyapunov Orbit Generator
for i=1:length(C_target)
    [X0,T_half,JC,nu] = gen_lyapunov(C_target(i),LP,mu);
    X0f(i,:) = X0;
    T_fam(i,1) = 2*T_half;
    nu_fam(i,1) = nu;
end

%% Plotting
a = 3.844 * 1e8;
m1 = 5.97219*1e24;

figure(1);hold on
l1p_str = sprintf('L_2 Lyapunov (JC=%.4f)',JC);
scatter((1-mu)*a,0,40,'o','filled')
opts = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
for i = 1:length(C_target)
    try
        [t, X] = ode89(@(t,X) cr3bp_eom_gen(t,X,mu), 0:0.0001:T_fam(i), X0f(i,:)');
        X_In = cr3bp_conv_rot2inert(X,t);
        [Xdim, tdim] = cr3bp_conv_dim(X,t,a, 'Primary Mass', [m1, mu]);
    catch
        fprintf("Not possible to resolve\n");
    end
    p1=plot3(Xdim(:,1), Xdim(:,2),Xdim(:,3),'LineWidth',1.5,'DisplayName',l1p_str);
    axis("equal")
end
hold off
xlabel('x'); ylabel('y');
title('Planar Lyapunov Orbit [Rotating]');

figure(2); hold on
plot(C_target,nu_fam,'DisplayName','Vertical Stability Index')
yline(1,'DisplayName','Bifurcation Point')
xlabel("Jacobi Constant")
ylabel("Vertical Stability Index \nu")
legend('Location','northeast')
title("Jacobi Constant vs Vertical Stability Index")


