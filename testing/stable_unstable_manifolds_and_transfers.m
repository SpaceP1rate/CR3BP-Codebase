clear; close; clc;
%% Initialization
addpath(genpath('..\src'))
mu = 0.0121505;
n = 8; 
opts1 = odeset('Events', @(t,X) poincare_neg_90L1(t,X),'RelTol',1e-12,'AbsTol',1e-12);
opts2 = odeset('Events', @(t,X) poincare_neg_90L2(t,X),'RelTol',1e-12,'AbsTol',1e-12);
% max 0.13 for L1, 0.15 for L2 
%%  Lyapunov Orbit around L2
LP = 'l2';
deltaX = 0.007; 
[X0_L2,T_half_L2] = gen_lyapunov(deltaX,LP,mu);

%% Lyapunov Orbit around L1
LP = 'l1';  
deltaX = 0.007;
[X0_L1,T_half_L1] = gen_lyapunov(deltaX,LP,mu);

%% Propagating Orbits with STM

Y0L2 = [X0_L2'; reshape(eye(6), [], 1)];
Y0L1 = [X0_L1'; reshape(eye(6), [], 1)];

[~,Y_L2] = ode45(@(t,Y) cr3bp_eom_stm(t,Y,mu), linspace(0,2*T_half_L2,n), Y0L2);
[~,Y_L1] = ode45(@(t,Y) cr3bp_eom_stm(t,Y,mu), linspace(0,2*T_half_L1,n), Y0L1);

%% Calculating Unstable and Stable Manifold using Monodromy Matrix

ML2 = reshape(Y_L2(end, 7:end), 6, 6);
ML1 = reshape(Y_L1(end, 7:end), 6, 6);

[VL1, DL1] = eig(ML1);
[VL2, DL2] = eig(ML2);

X_L2 = Y_L2(:,1:6);
X_L1 = Y_L1(:,1:6);

X_L2_unstable = X_L2 + 1e-3*VL2(:,1)';
X_L1_stable   = X_L1 - 1e-3*VL1(:,2)';

%% Plotting
clf;
figure(1); hold on
for i = 1:n
    [tL2u, XL2u,TPL2,XPL2d,~] = ode89(@(t,X) cr3bp_eom_gen(t,X,mu), 0:0.0001:4*T_half_L2, X_L2_unstable(i,:)',opts2);
    p1=plot3(XL2u(:,1), XL2u(:,2), XL2u(:,3),'LineWidth',0.7,'DisplayName','L_2 Unstable','Color','r');
    [tL1u, XL1u,TPL1,XPL1d,~] = ode89(@(t,X) cr3bp_eom_gen(t,X,mu), 0:-0.0001:-4*T_half_L1, X_L1_stable(i,:)',opts1);
    p2=plot3(XL1u(:,1), XL1u(:,2), XL1u(:,3),'LineWidth',0.7,'DisplayName','L_1 Stable','Color','g');

    XPL2(i,:) = XPL2d;
    XPL1(i,:) = XPL1d;
end
scatter(1-mu,0,'ko','filled','DisplayName','Moon')
axis equal

figure(2)
subplot(1,2,1); hold on
scatter(XPL1(:,2),XPL1(:,4),'go','filled')
scatter(XPL2(:,2),XPL2(:,4),'ro','filled')
xlabel("Y [LU]")
ylabel("V_x [LU/TU]")
axis equal
subplot(1,2,2); hold on
scatter(XPL1(:,2),XPL1(:,5),'go','filled')
scatter(XPL2(:,2),XPL2(:,5),'ro','filled')
xlabel("Y [LU]")
ylabel("V_y [LU/TU]")
axis equal
