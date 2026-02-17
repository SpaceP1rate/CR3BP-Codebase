clear; close; clc;
%% Initialization
addpath(genpath('..\src'))

mu = 0.0121505;
LP = 'l1';
L = cr3bp_sys_lagpoints(mu, LP);
LP_Jac = cr3bp_sys_jacconst([L(1), 0, 0, 0, 0, 0], mu);

deltaX = linspace(0.01,0.11,20); % max 0.13 for L1, 0.15 for L2 
deltaZ = linspace(1e-3,1e-2,20);
%% Lyapunov Orbit Generator
for i=1:length(deltaX)
    [X0,T_half,JC,nu] = gen_lyapunov(deltaX(i),LP,mu);
    X0L(i,:) = X0;
    T_Lfam(i,1) = 2*T_half;
    nu_Lfam(i,1) = nu;
end

X0h = X0L(2,:);
X0h(3) = 1e-4;
type = 'south';
for i =1:length(deltaZ)
    [Xh0,T_halfh,JCh,nuh] = gen_halo(X0h,type,LP,mu);
    X0H(i,:) = Xh0;
    T_Hfam(i,1) = 2*T_halfh;
    nu_Hfam(i,1) = nuh;
    X0h = Xh0;
    X0h(3) = deltaZ(i);
end

%% Plotting
a = 3.844 * 1e8;
m1 = 5.97219*1e24;

figure(1);hold on
l1p_str = sprintf('L_2 Lyapunov (JC=%.4f)',JC);
scatter((1-mu),0,40,'o','filled')
opts = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
for i = 1:length(deltaX)
    try
        [t, X] = ode89(@(t,X) cr3bp_eom_gen(t,X,mu), 0:0.0001:T_Lfam(i), X0L(i,:)');
        [Xdim, tdim] = cr3bp_conv_dim(X,t,a, 'Primary Mass', [m1, mu]);
        [th, Xh] = ode89(@(t,X) cr3bp_eom_gen(t,X,mu), 0:0.0001:T_Hfam(i), X0H(i,:)');
        [Xdimh, tdimh] = cr3bp_conv_dim(Xh,th,a, 'Primary Mass', [m1, mu]);

    catch
        fprintf("Not possible to resolve\n");
    end
    p1=plot3(X(:,1), X(:,2),X(:,3),'LineWidth',1.5,'DisplayName',l1p_str,'Color','b');
    p2=plot3(Xh(:,1), Xh(:,2),Xh(:,3),'LineWidth',1.5,'DisplayName','North Halo','Color','r');
    % axis('auto')
end

zlim([-1e-1,1e-1])
xlim([1-mu-0.3,1-mu+0.3])
ylim([-0.3,0.3])
hold off
xlabel('x'); ylabel('y');
title('Planar Lyapunov Orbit [Rotating]');
view(3)
figure(2); hold on
plot(deltaX,nu_Lfam,'DisplayName','Vertical Stability Index Lyapunov')
plot(deltaZ,nu_Hfam,'DisplayName','Vertical Stability Index Halo')
yline(1,'DisplayName','Bifurcation Point')
xlabel("X Displacement")
ylabel("Vertical Stability Index \nu")
legend('Location','northeast')
title("Jacobi Constant vs Vertical Stability Index")


