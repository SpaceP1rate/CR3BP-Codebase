%% CR3BP Final 5-Plot Orbit Set (Refined View Angles)
clear; close all; clc;
%% 1. Parameters & Setup
mu = 0.0121505;
colors = struct( ...
    'L1',    [0 0.447 0.741], ...
    'L2',    [0.850 0.325 0.098], ...
    'North', [0.150 0.650 0.150], ... 
    'South', [0.850 0.650 0.100], ... 
    'Seed',  [0 0 0]);
dx_steps = linspace(0.01, 0.11, 8);
dz_steps = linspace(0.005, 0.08, 12); 
moon_x = 1 - mu;

%% ==========================================================
%% 2. FIGURE 1 — PLANAR LYAPUNOV FAMILIES (Generates Seeds)
%% ==========================================================
fig1 = figure('Name','Planar Lyapunov Families','Units','normalized','Position',[0.05 0.55 0.4 0.35]);
hold on; grid on; box on;
[X0_seed1, X_bif1] = plot_lyapunov_family('l1', mu, dx_steps, colors.L1, colors.Seed, 2);
[X0_seed2, X_bif2] = plot_lyapunov_family('l2', mu, dx_steps, colors.L2, colors.Seed, 3);
scatter(moon_x, 0, 40, 'k', 'filled', 'DisplayName', 'Moon');
title('Planar Lyapunov Families'); 
xlabel('x [LU]'); ylabel('y [LU]'); axis equal;
legend('Location', 'northeast');

%% ==========================================================
%% 3. FIGURE 2 — L1 3D OVERVIEW (Linked to Seed 1)
%% ==========================================================
fig2 = figure('Name','L1 3D Overview','Units','normalized','Position',[0.5 0.55 0.4 0.35]);
hold on; grid on;
L1_coords = cr3bp_sys_lagpoints(mu, 'l1');
plot3(X_bif1(:,1), X_bif1(:,2), X_bif1(:,3), 'Color', colors.Seed, 'LineWidth', 2, 'DisplayName', 'L1 Seed');
[~, Xn1] = get_single_halo(X0_seed1, 'north', 'l1', mu, dz_steps(5));
[~, Xs1] = get_single_halo(X0_seed1, 'south', 'l1', mu, dz_steps(5));
plot3(Xn1(:,1), Xn1(:,2), Xn1(:,3), 'Color', colors.North, 'LineWidth', 1.2, 'DisplayName', 'North Halo');
plot3(Xs1(:,1), Xs1(:,2), Xs1(:,3), 'Color', colors.South, 'LineWidth', 1.2, 'DisplayName', 'South Halo');
scatter3(L1_coords(1), 0, 0, 40, 'r', 'x', 'LineWidth', 1.5, 'DisplayName', 'L1 Point');
view_width = 0.1; 
xlim([L1_coords(1)-view_width, L1_coords(1)+view_width]);
ylim([-view_width/2, view_width/2]); zlim([-view_width/2, view_width/2]);
title('L1 Bifurcation Overview'); 
xlabel('x [LU]'); ylabel('y [LU]'); zlabel('z [LU]');
axis equal; 
view([15, 15]); % Improved perspective for slenderness
legend('Location', 'best');

%% ==========================================================
%% 4. FIGURE 3 — L2 3D OVERVIEW (Linked to Seed 2)
%% ==========================================================
fig3 = figure('Name','L2 3D Overview','Units','normalized','Position',[0.5 0.1 0.4 0.35]);
hold on; grid on;
L2_coords = cr3bp_sys_lagpoints(mu, 'l2');
plot3(X_bif2(:,1), X_bif2(:,2), X_bif2(:,3), 'Color', colors.Seed, 'LineWidth', 2, 'DisplayName', 'L2 Seed');
[~, Xn2] = get_single_halo(X0_seed2, 'north', 'l2', mu, dz_steps(5));
[~, Xs2] = get_single_halo(X0_seed2, 'south', 'l2', mu, dz_steps(5));
plot3(Xn2(:,1), Xn2(:,2), Xn2(:,3), 'Color', colors.North, 'LineWidth', 1.2, 'DisplayName', 'North Halo');
plot3(Xs2(:,1), Xs2(:,2), Xs2(:,3), 'Color', colors.South, 'LineWidth', 1.2, 'DisplayName', 'South Halo');
scatter3(L2_coords(1), 0, 0, 40, 'r', 'x', 'LineWidth', 1.5, 'DisplayName', 'L2 Point');
xlim([L2_coords(1)-view_width, L2_coords(1)+view_width]);
ylim([-view_width/2, view_width/2]); zlim([-view_width/2, view_width/2]);
title('L2 Bifurcation Overview'); 
xlabel('x [LU]'); ylabel('y [LU]'); zlabel('z [LU]');
axis equal; 
view([-15, 15]); % Flipped for L2 perspective
legend('Location', 'best');

%% ==========================================================
%% 5. FIGURE 4 — L1 HALO FAMILY (Uses Seed 1)
%% ==========================================================
fig4 = figure('Name','L1 Halo Family','Units','normalized','Position',[0.05 0.1 0.4 0.35]);
hold on; grid on;
plot_halo_family(X0_seed1, 'north', 'l1', mu, dz_steps, colors.North, 'L1 North');
plot_halo_family(X0_seed1, 'south', 'l1', mu, dz_steps, colors.South, 'L1 South');
format_halo_plot('L1 Halo Family', mu, L1_coords, moon_x, 0.12);
view([25, 10]); % Consistent family view

%% ==========================================================
%% 6. FIGURE 5 — L2 HALO FAMILY (Uses Seed 2)
%% ==========================================================
fig5 = figure('Name','L2 Halo Family','Units','normalized','Position',[0.25 0.05 0.4 0.35]);
hold on; grid on;
plot_halo_family(X0_seed2, 'north', 'l2', mu, dz_steps, colors.North, 'L2 North');
plot_halo_family(X0_seed2, 'south', 'l2', mu, dz_steps, colors.South, 'L2 South');
format_halo_plot('L2 Halo Family', mu, L2_coords, moon_x, 0.18);
view([-25, 10]); % Consistent family view

%% ---------------- PRINTING COMMANDS (600 DPI) -------------
exportgraphics(fig1, 'Planar_Lyapunov.png', 'Resolution', 600);
exportgraphics(fig2, 'L1_Overview.png', 'Resolution', 600);
exportgraphics(fig3, 'L2_Overview.png', 'Resolution', 600);
exportgraphics(fig4, 'L1_Halo_Family.png', 'Resolution', 600);
exportgraphics(fig5, 'L2_Halo_Family.png', 'Resolution', 600);

%% ---------------- HELPER FUNCTIONS ------------------------
function [X0_seed, X_bif] = plot_lyapunov_family(LP, mu, dx_steps, col, seed_col, target_idx)
    L = cr3bp_sys_lagpoints(mu, LP);
    for i = 1:length(dx_steps)
        [X0, Th] = gen_lyapunov(dx_steps(i), LP, mu);
        [~, X] = ode89(@(t,X) cr3bp_eom_gen(t,X,mu), [0 2*Th], X0);
        if i == target_idx
            plot(X(:,1), X(:,2), 'Color', seed_col, 'LineWidth', 2, 'DisplayName', sprintf('%s Seed', upper(LP)));
            X0_seed = X0; X_bif = X;
        else
            plot(X(:,1), X(:,2), 'Color', col, 'LineWidth', 0.5, 'HandleVisibility', 'off');
        end
    end
    plot(NaN, NaN, 'Color', col, 'DisplayName', sprintf('%s Lyapunovs', upper(LP)));
    scatter(L(1), 0, 30, 'r', 'x', 'HandleVisibility', 'off');
end
function plot_halo_family(seed, type, LP, mu, dz_array, col, label)
    h_seed = seed;
    for j = 1:length(dz_array)
        h_seed(3) = dz_array(j);
        try
            [Xh0, Th] = gen_halo(h_seed, type, LP, mu);
            [~, Xh] = ode89(@(t,X) cr3bp_eom_gen(t,X,mu), [0 2*Th], Xh0);
            opacity = 0.2 + (0.8 * (j/length(dz_array)));
            plot3(Xh(:,1), Xh(:,2), Xh(:,3), 'Color', [col, opacity], 'LineWidth', 0.8, ...
                'HandleVisibility', ifexpr(j == length(dz_array), 'on', 'off'), ...
                'DisplayName', label);
            h_seed = Xh0;
        catch, continue; end
    end
end
function [X0_h, X_traj] = get_single_halo(seed, type, LP, mu, dz)
    s = seed; s(3) = dz;
    [X0_h, Th] = gen_halo(s, type, LP, mu);
    [~, X_traj] = ode89(@(t,X) cr3bp_eom_gen(t,X,mu), [0 2*Th], X0_h);
end
function format_halo_plot(lbl, mu, L, moon_x, width)
    title(lbl); 
    xlabel('x [LU]'); ylabel('y [LU]'); zlabel('z [LU]');
    scatter3(moon_x, 0, 0, 45, 'k', 'filled', 'DisplayName', 'Moon');
    scatter3(L(1), 0, 0, 40, 'r', 'x', 'LineWidth', 1.5, 'DisplayName', 'LP');
    midpoint = (L(1) + moon_x) / 2;
    xlim([midpoint - width, midpoint + width]);
    axis equal; grid on; legend('Location', 'best');
end
function val = ifexpr(cond, v1, v2)
    if cond, val = v1; else, val = v2; end
end