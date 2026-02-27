clear; close; clc;
%% Initialization
addpath(genpath('..\src'))

mu = 0.0121505;
LP = 'l1';
L = cr3bp_sys_lagpoints(mu, LP);
LP_Jac = cr3bp_sys_jacconst([L(1), 0, 0, 0, 0, 0], mu);

deltaX = 0.03; % max 0.13 for L1, 0.15 for L2 
num_orbits = 150;
%% Lyapunov Orbit Generator

[X0,T_half,JC,Phi] = gen_lyapunov(deltaX,LP,mu);
X0L= X0;
T_Lfam = 2*T_half;

X0h = X0L;
X0h(3) = 1e-3;
type = 'south';
cond = true;
step_mode = 'z-step';
ds_z = 0.001; % Step size for Z
ds_x = 0.0015; % Step size for X

%% Halo Generator
for i = 1:num_orbits
    % --- 1. Correct (Correction Step) ---
    if strcmp(step_mode, 'z-step')
        [Xh0, T_halfh, JCh, nuh] = gen_halo_z(X0h, type, LP, mu);
    else
        [Xh0, T_halfh, JCh, nuh] = gen_halo_x(X0h, type, LP, mu);
    end
    
    % --- 2. Store Data ---
    X0H(i,:) = Xh0;
    T_Hfam(i) = 2 * T_halfh;
    
    % --- 3. Update Mode (Selection Step) ---
    if i > 1
        dx = abs(X0H(i,1) - X0H(i-1,1));
        dz = abs(X0H(i,3) - X0H(i-1,3));
        
        if dz < 0.8* dx
            step_mode = 'x-step';
        elseif dz >= 1.2 * dx
            step_mode = 'z-step';
        end
    end
    
    % --- 4. Predict Next Guess (Prediction Step) ---
    X0h = Xh0; % Always seed next guess with current solution
    if strcmp(step_mode, 'z-step')
        % Step away from the plane: for South start @ North crossing, z increases
        X0h(3) = Xh0(3) + sign(Xh0(3)) * ds_z; 
        
    else
        % Step toward the Moon: for L1, x increases
        x_dir = 1; if strcmp(LP, 'l2'), x_dir = -1; end
        X0h(1) = Xh0(1) + x_dir * ds_x;
    end
end
%% Plotting
clf;
figure(1);hold on
scatter3(1-mu,0,0,'ko','filled','DisplayName','Moon')
scatter3(L(1),0,0,'gx','DisplayName','L_1')
[t, X] = ode89(@(t,X) cr3bp_eom_gen(t,X,mu), 0:0.0001:T_Lfam, X0L');
for i = 1:num_orbits
    [th, Xh] = ode89(@(t,X) cr3bp_eom_gen(t,X,mu), 0:0.0001:T_Hfam(i), X0H(i,:)');
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



