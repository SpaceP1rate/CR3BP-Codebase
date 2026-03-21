clear; close; clc;
%% Initialization
addpath(genpath('..\src'))

mu = 0.0121505;

LP = 'l2';
type = 'south';
num_orbits = 28000;
ds = 5e-5;
opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

L = cr3bp_sys_lagpoints(mu, LP);
LP_Jac = cr3bp_sys_jacconst([L(1), 0, 0, 0, 0, 0], mu);
%% Lyapunov Orbit Generator
deltaX = 0.028;  
[X0,T_half,JC,~] = gen_lyapunov(deltaX,LP,mu);
X0L= X0;
T_Lfam = 2*T_half;

X0h = X0L;
X0h(3) = 1e-4;
%% Halo Orbit Generator
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

        if Xstar(2)<0 
            Xstar(2) = -Xstar(2);
        end
    case 'north'
        if Xstar(2)>0 
            Xstar(2) = -Xstar(2);
        end
end

Xh0_pac = Xh0;
Xh0_pac(1) = Xh0(1) + Xstar(1)*ds;
Xh0_pac(3) = Xh0(3) + Xstar(2)*ds;
Xh0_pac(5) = Xh0(5) + Xstar(3)*ds;
alpha = 0.3;
tic
for i=2:num_orbits
    [Xh0, T_halfh, Xstar_new, flag] = gen_halo_pac(Xh0_pac,Xh0,Xstar,ds,type,LP,mu,alpha);
    X0H(i,:) = Xh0;
    T_Hfam(i) = 2 * T_halfh;
    
    if abs(Xh0(1)/(1-mu) - 1)*100 < 8e-2
        alpha = 0.08;
        fprintf("Refined Alpha: %f\n\n",alpha)
    elseif abs(Xh0(1)/(1-mu) - 1)*100 < 8e-3
        alpha = 0.2;
        fprintf("Moderate Alpha: %f\n\n",alpha)
    else
        fprintf("Coarse Alpha: %f\n\n",alpha)
    end
    Xstar = Xstar_new;
    Xh0_pac = Xh0;
    Xh0_pac(1) = Xh0(1) + Xstar(1)*ds;
    Xh0_pac(3) = Xh0(3) + Xstar(2)*ds;
    Xh0_pac(5) = Xh0(5) + Xstar(3)*ds;
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
[t, X] = ode78(@(t,X) cr3bp_eom_gen(t,X,mu), [0,T_Lfam], X0L');
for i = 1:500:gen_num
    [th, Xh] = ode89(@(t,X) cr3bp_eom_gen(t,X,mu), [0,T_Hfam(i)], X0H(i,:)',opts);
    p2=plot3(Xh(:,1), Xh(:,2),Xh(:,3),'LineWidth',1.2,'DisplayName','South Halo','Color','r');
end
p1=plot3(X(:,1), X(:,2),X(:,3),'LineWidth',1.2,'DisplayName','L_1 Lyapunov','Color','b');

axis equal
zlim([-2.5e-1,2.5e-1])
xlim([1-mu-0.1,1-mu+0.3])
ylim([-0.3,0.3])
hold off; grid on;
xlabel('x'); ylabel('y');
title('Periodic Orbits in CR3BP');
view(3)



