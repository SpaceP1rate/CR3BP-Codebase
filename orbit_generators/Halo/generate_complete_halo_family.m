clear; close; clc;
addpath(genpath('..\src'))

mu = 0.0121505;
num_orbits = 12000;

opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

% Generate Complete Family of Halo Orbits
if ~exist('south_halo_l1_ic.mat')
    LP = 'l1'; type = 'south';
    halo_orbit_continuation(mu,LP,type,num_orbits)
    sendTelegram(sprintf('Simulation Complete for %s %s halo family!',type,LP))
end
if ~exist('north_halo_l1_ic.mat')
    LP = 'l1'; type = 'north';
    halo_orbit_continuation(mu,LP,type,num_orbits)
    sendTelegram(sprintf('Simulation Complete for %s %s halo family!',type,LP))
end
if ~exist('south_halo_l2_ic.mat')
    LP = 'l2'; type = 'south';
    halo_orbit_continuation(mu,LP,type,num_orbits) 
    sendTelegram(sprintf('Simulation Complete for %s %s halo family!',type,LP))
end
if ~exist('north_halo_l2_ic.mat')
    LP = 'l2'; type = 'north';
    halo_orbit_continuation(mu,LP,type,num_orbits) 
    sendTelegram(sprintf('Simulation Complete for %s %s halo family!',type,LP))
end

load('south_halo_l1_ic.mat')
L = load('north_halo_l1_ic.mat');
family.l1_north = L.family.l1_north;
L = load('south_halo_l2_ic.mat');
family.l2_south = L.family.l2_south;
L = load('north_halo_l2_ic.mat');
family.l2_north = L.family.l2_north;
clear("L")


figure(1)
s1 = subplot(1,1,1); hold on
figure(2)
s2 = subplot(1,1,1); hold on
for i=1:200:num_orbits
    [~,XL1s] = ode78(@(t,X) cr3bp_eom_gen(t,X,mu), [0,family.l1_south.T(i)], family.l1_south.IC(i,:)',opts);
    [~,XL1n] = ode78(@(t,X) cr3bp_eom_gen(t,X,mu), [0,family.l1_north.T(i)], family.l1_north.IC(i,:)',opts);
    [~,XL2s] = ode78(@(t,X) cr3bp_eom_gen(t,X,mu), [0,family.l2_south.T(i)], family.l2_south.IC(i,:)',opts);
    [~,XL2n] = ode78(@(t,X) cr3bp_eom_gen(t,X,mu), [0,family.l2_north.T(i)], family.l2_north.IC(i,:)',opts);
    if i ==1
        plot3(s1,XL1s(:,1), XL1s(:,2),XL1s(:,3),'LineWidth',0.1,'DisplayName','L_1 South','Color','r');
        plot3(s1,XL1n(:,1), XL1n(:,2),XL1n(:,3),'LineWidth',0.1,'DisplayName','L_1 North','Color','b');
        plot3(s2,XL2s(:,1), XL2s(:,2),XL2s(:,3),'LineWidth',0.1,'DisplayName','L_2 South','Color','r');
        plot3(s2,XL2n(:,1), XL2n(:,2),XL2n(:,3),'LineWidth',0.1,'DisplayName','L_2 North','Color','b');
    else
        plot3(s1,XL1s(:,1), XL1s(:,2),XL1s(:,3),'LineWidth',0.1,'HandleVisibility','off','Color','r');
        plot3(s1,XL1n(:,1), XL1n(:,2),XL1n(:,3),'LineWidth',0.1,'HandleVisibility','off','Color','b');
        plot3(s2,XL2s(:,1), XL2s(:,2),XL2s(:,3),'LineWidth',0.1,'HandleVisibility','off','Color','r');
        plot3(s2,XL2n(:,1), XL2n(:,2),XL2n(:,3),'LineWidth',0.1,'HandleVisibility','off','Color','b');
    end
end
scatter3(s1,1-mu,0,0,10,'ko','filled','DisplayName','Moon')
scatter3(s2,1-mu,0,0,10,'ko','filled','DisplayName','Moon')
view(s1,66,20); view(s2,-66,20)
axis(s1,'equal'); axis(s2,'equal')
% zlim(s2,zlim(s1));
grid(s1,'on'); grid(s2,'on')
legend(s1,'Location','northeast')
legend(s2,'Location','northeast')
xlabel(s1,'X [LU]'); ylabel(s1,'Y [LU]'); zlabel(s1,'Z [LU]')
xlabel(s2,'X [LU]'); ylabel(s2,'Y [LU]'); zlabel(s2,'Z [LU]')
title(s1,'L_1 North/South Halo Families'); title(s2,'L_2 North/South Halo Families')

figure(3)
subplot(1,2,1);hold on
plot(family.l1_north.IC(1:end,1),family.l1_north.T(1:end),'-.k','LineWidth',1.2,'DisplayName','L_1 North Halo')
plot(family.l2_north.IC(1:end,1),family.l2_north.T(1:end),'--k','LineWidth',1.2,'DisplayName','L_2 North Halo')

title("Orbital Period of Halo Families vs X Displacement")
xlabel('X [LU]'); ylabel("Period [TU]")
axis auto
legend('Location','northeast')
grid on

subplot(1,2,2);hold on

plot(family.l1_north.IC(1:end,3),family.l1_north.T(1:end),'-.k','LineWidth',1.2,'DisplayName','L_1 North Halo')
plot(family.l2_north.IC(1:end,3),family.l2_north.T(1:end),'--k','LineWidth',1.2,'DisplayName','L_2 North Halo')

title("Orbital Period of Halo Families vs Z Displacement")
xlabel('Z [LU]'); ylabel("Period [TU]")
axis tight
legend('Location','northeast')
grid on
