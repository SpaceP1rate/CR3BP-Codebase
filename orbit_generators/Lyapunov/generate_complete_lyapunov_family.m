clear; close; clc;clf;
addpath(genpath('..\..\src'))

mu = 0.0121505;
num_orbits = 2000;
 
opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-13);

% Generate Complete Family of Halo Orbits
if ~exist('l1_lyapunov_ic.mat')
    LP = 'l1';
    lyapunov_orbit_continuation(mu,LP,num_orbits)
    sendTelegram(sprintf('Simulation Complete for %s lyapunov family!',LP))
end
if ~exist('l2_lyapunov_ic.mat')
    LP = 'l2';
    lyapunov_orbit_continuation(mu,LP,num_orbits)
    sendTelegram(sprintf('Simulation Complete for %s lyapunov family!',LP))
end

load('l1_lyapunov_ic.mat')
L = load('l2_lyapunov_ic.mat');
family.l2 = L.family.l2;
clear("L")

beta = linspace(-5500,5500,num_orbits);

[~,idx_t_l1]= find(abs(family.l1.params(:,1)-(beta + 2)/(-2))<1& abs(family.l1.params(:,2)-beta)<1,2,'last');
[~,idx_t_l2]= find(abs(family.l2.params(:,1)-(beta + 2)/(-2))<7e-1 & abs(family.l2.params(:,2)-beta)<7e-1,2,'first');

figure(1);clf;hold on

for i=1:50:num_orbits
    [~,XL1] = ode89(@(t,X) cr3bp_eom_gen(t,X,mu), [0,family.l1.T(i)], family.l1.IC(i,:)',opts);
    [~,XL2] = ode89(@(t,X) cr3bp_eom_gen(t,X,mu), [0,family.l2.T(i)], family.l2.IC(i,:)',opts);

    if i ==1
        plot3(XL1(:,1), XL1(:,2),XL1(:,3),'LineWidth',0.1,'DisplayName','L_1 Lyapunov','Color','b');
        plot3(XL2(:,1), XL2(:,2),XL2(:,3),'LineWidth',0.1,'DisplayName','L_2 Lyapunov','Color','r');

    else
        plot3(XL1(:,1), XL1(:,2),XL1(:,3),'LineWidth',0.1,'HandleVisibility','off','Color','b');
        plot3(XL2(:,1), XL2(:,2),XL2(:,3),'LineWidth',0.1,'HandleVisibility','off','Color','r');

    end
end

% Tangent Bifurcations
for i = 1:2
    [~,XL1] = ode89(@(t,X) cr3bp_eom_gen(t,X,mu), [0,family.l1.T(idx_t_l1(i))], family.l1.IC(idx_t_l1(i),:)',opts);
    [~,XL2] = ode89(@(t,X) cr3bp_eom_gen(t,X,mu), [0,family.l2.T(idx_t_l2(i))], family.l2.IC(idx_t_l2(i),:)',opts);
    
    if i ==1
        plot3(XL1(:,1), XL1(:,2),XL1(:,3),'--k','LineWidth',0.5,'DisplayName','L_1 Lyapunov Tangent Bifurcation');
        plot3(XL2(:,1), XL2(:,2),XL2(:,3),'-.k','LineWidth',0.5,'DisplayName','L_2 Lyapunov Tangent Bifurcation'); 
    else
        plot3(XL1(:,1), XL1(:,2),XL1(:,3),'--k','LineWidth',0.5,'HandleVisibility','off');
        plot3(XL2(:,1), XL2(:,2),XL2(:,3),'-.k','LineWidth',0.5,'HandleVisibility','off'); 
    end
end

scatter3(1-mu,0,0,10,'ko','filled','DisplayName','Moon')

axis('equal')

grid('on')
legend('Location','northeast')

xlabel('X [LU]'); ylabel('Y [LU]')
title('L_1/L_2 Lyapunov Family')

figure(2); hold on
plot(family.l1.IC(1:end,1),family.l1.T(1:end),'-.k','LineWidth',1.2,'DisplayName','L_1 Lyapunov')
plot(family.l2.IC(1:end,1),family.l2.T(1:end),'--k','LineWidth',1.2,'DisplayName','L_2 Lyapunov')

title("Orbital Period of Lyapunov Families vs X Displacement")
xlabel('X [LU]'); ylabel("Period [TU]")
axis auto
legend('Location','northeast')
grid on


id = 1:num_orbits;


figure(3);hold on
plot(family.l1.params(id,1),family.l1.params(id,2),'-b','DisplayName','L_1 Lyapunov')
plot(family.l2.params(id,1),family.l2.params(id,2),'-r','DisplayName','L_2 Lyapunov')
plot((beta + 2)/(-2),beta,'--k','DisplayName','Tangent Bifurcation')
plot((beta + 2)/(2),beta,'-.k','DisplayName','Period-Doubling Bifurcation')


title("Broucke Stability Diagram")
xlabel('\alpha'); ylabel("\beta")
axis tight
legend('Location','northeast')
grid on
