clear; close; clc;
%% Initialization
addpath(genpath('..\src'))

mu = 0.0121505;
LP = 'l1';
L = cr3bp_sys_lagpoints(mu, LP);
LP_Jac = cr3bp_sys_jacconst([L(1), 0, 0, 0, 0, 0], mu);

deltaX = 0.03; % max 0.13 for L1, 0.15 for L2 
%% Lyapunov Orbit Generator

[X0,T_half,JC,Phi] = gen_lyapunov(deltaX,LP,mu);
X0L= X0;
T_Lfam = 2*T_half;

X0h = X0L;
X0h(3) = 3e-2;
type = 'south';

[Xh0,T_halfh,JCh,nuh] = gen_halo(X0h,type,LP,mu);
X0H = Xh0;
T_Hfam = 2*T_halfh;
nu_Hfam = nuh;

Y0 = [X0H'; reshape(eye(6), [], 1)];

[~,Yp] = ode45(@(t,Y) cr3bp_eom_stm(t,Y,mu), [0 T_Lfam], Y0);
M = reshape(Yp(end, 7:end), 6, 6);

X_halo1 = Yp(45,1:6);
X_halo2 = Yp(75,1:6);
[V, D] = eig(M);

Xp = X0H' - 1e-3* V(:,2);
Xp1 = X_halo1' - 1e-3 * V(:,2);
Xp2 = X_halo2' - 1e-3 * V(:,2);
%% Plotting
clf;
figure(1);hold on
scatter3(1-mu,0,0,'ko','filled','DisplayName','Moon')
scatter3(L(1),0,0,'gx','DisplayName','L_1')
[t, X] = ode89(@(t,X) cr3bp_eom_gen(t,X,mu), 0:0.0001:T_Lfam, X0L');
[th, Xh] = ode89(@(t,X) cr3bp_eom_gen(t,X,mu), 0:0.0001:T_Hfam, X0H');
[thu, Xhu] = ode89(@(t,X) cr3bp_eom_gen(t,X,mu), 0:-0.0001:-1.5*T_Hfam, Xp);
[thu1, Xhu1] = ode89(@(t,X) cr3bp_eom_gen(t,X,mu), 0:-0.0001:-1.5*T_Hfam, Xp1);
[thu2, Xhu2] = ode89(@(t,X) cr3bp_eom_gen(t,X,mu), 0:-0.0001:-1.5*T_Hfam, Xp2);
p1=plot3(X(:,1), X(:,2),X(:,3),'LineWidth',1.2,'DisplayName','L_1 Lyapunov','Color','b');
p2=plot3(Xh(:,1), Xh(:,2),Xh(:,3),'LineWidth',1.2,'DisplayName','South Halo','Color','k');
p3=plot3(Xhu(:,1), Xhu(:,2),Xhu(:,3),'LineWidth',0.7,'DisplayName','South Halo Unstable Manifold','Color','r');
p4=plot3(Xhu1(:,1), Xhu1(:,2),Xhu1(:,3),'LineWidth',0.7,'DisplayName','South Halo Unstable Manifold','Color','r');
p5=plot3(Xhu2(:,1), Xhu2(:,2),Xhu2(:,3),'LineWidth',0.7,'DisplayName','South Halo Unstable Manifold','Color','r');

axis equal
zlim([-1e-1,1e-1])
xlim([1-mu-0.6,1-mu+0.1])
ylim([-0.3,0.3])
hold off; grid on;
xlabel('x'); ylabel('y');
title('Periodic Orbits in CR3BP');
legend()
view(3)



