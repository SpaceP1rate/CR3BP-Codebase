clear; close; clc;
%% Initialization
addpath(genpath('..\src'))

mu = 0.0121505;
LP = 'l1';
L = cr3bp_sys_lagpoints(mu, LP);
LP_Jac = cr3bp_sys_jacconst([L(1), 0, 0, 0, 0, 0], mu);

deltaX = 0.018; % max 0.13 for L1, 0.15 for L2 
num_orbits = 1200;
%% Lyapunov Orbit Generator

[X0,T_half,Phi] = gen_lyapunov(deltaX,LP,mu);
X0L= X0;
T_Lfam = 2*T_half;

X0h = X0L;

type = 'north';

switch type
    case 'north'
        dir = -1;
    case 'south'
        dir = 1;
end
X0h(3) = 1e-4*dir;
%% Halo Generator
[Xh0, T_halfh, JCh, nuh] = gen_halo_z(X0h, type, LP, mu);
    
X0H(1,:) = Xh0;
T_Hfam(1) = 2 * T_halfh;
Y0 = [X0H'; reshape(eye(6), [], 1)];
[~, YE] = ode78(@(t,Y) cr3bp_eom_stm(t,Y,mu), [0 T_Hfam], Y0);
Phi   = reshape(YE(end, 7:end),6,6);

J = [Phi(4,1),Phi(4,3), Phi(4,5);
     Phi(6,1),Phi(6,3), Phi(6,5)];

Xstar = null(J);

switch type
    case 'south'

        if Xstar(2)<0 && Xstar(1) <0
            Xstar = -Xstar;
        end
    case 'north'
        if Xstar(2)<0 && Xstar(1) <0
            Xstar = null(J);
        end
end

ds = 1e-3;

Xh0_pac = Xh0;
Xh0_pac(1) = Xh0(1) + Xstar(1)*ds;
Xh0_pac(3) = Xh0(3) + Xstar(2)*ds;
Xh0_pac(5) = Xh0(5) + Xstar(3)*ds;
alpha = 0.3;
for i=2:num_orbits
    [Xh0, T_halfh, Xstar_new] = gen_halo_pac(Xh0_pac,Xh0,Xstar,ds,type,LP,mu,alpha);
    X0H(i,:) = Xh0;
    T_Hfam(i) = 2 * T_halfh;
    
    if Xh0(1) > 0.97
        alpha = 0.05;
        fprintf("Refined Alpha ds = %f\n\n",alpha)
    elseif Xh0(1) > 0.95
        % ds = 1e-4;
        alpha = 0.2;
        fprintf("Moderately Alpha to ds = %f\n\n",alpha)
    else
        % ds = 1e-3;
        alpha = 0.3;
        fprintf("Coarse Alpha to ds = %f\n\n",alpha)
    end
    Xstar = Xstar_new;
    Xh0_pac = Xh0;
    Xh0_pac(1) = Xh0(1) + Xstar(1)*ds;
    Xh0_pac(3) = Xh0(3) + Xstar(2)*ds;
    Xh0_pac(5) = Xh0(5) + Xstar(3)*ds;
    if Xh0(1) > 0.96 && abs(Xh0(3)) < 1e-3
        break;
        
    end
end
gen_num = i;
%% Plotting
clf;
figure(1);hold on
scatter3(1-mu,0,0,'ko','filled','DisplayName','Moon')
scatter3(L(1),0,0,'gx','DisplayName','L_1')
[t, X] = ode78(@(t,X) cr3bp_eom_gen(t,X,mu), 0:0.0001:T_Lfam, X0L');
for i = 1:5:gen_num
    [th, Xh] = ode78(@(t,X) cr3bp_eom_gen(t,X,mu), 0:0.0001:T_Hfam(i), X0H(i,:)');
    p2=plot3(Xh(:,1), Xh(:,2),Xh(:,3),'LineWidth',1.2,'DisplayName','South Halo','Color','r');
end
p1=plot3(X(:,1), X(:,2),X(:,3),'LineWidth',1.2,'DisplayName','L_1 Lyapunov','Color','b');

axis equal
zlim([-1e-1,1e-1])
xlim([1-mu-0.6,1-mu+0.1])
ylim([-0.3,0.3])
hold off; grid on;
xlabel('x'); ylabel('y');
title('Periodic Orbits in CR3BP');
view(3)



