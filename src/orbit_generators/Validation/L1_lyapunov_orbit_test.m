clear; close; clc;

addpath(genpath('Codebase'))
addpath('planet3D')


mu = 0.01216;
% Lagrange point calculation test
LP = cr3bp_sys_lagpoints(mu,'l1');
L1 = LP(1);

opts = odeset('Events', @(t,X) ode_event_xcrossneg(t,X),'RelTol',1e-12,'AbsTol',1e-12);
x0 = L1 - 0.1; y0 = 0;
vx0 = 0; vy0 = 2.2;
X0 = [x0,y0,x0,vy0]';
nStates = size(X0,1);

Phi0 = eye(nStates);

Y0 = [X0; Phi0(:)];
tol = 1e-13;          % Convergence tolerance
maxIter = 2000;         % Maximum iterations
alpha = 0.9;            % Correction coefficient of shooting slope
% Lyapunov Shooting Algorithm tests
for iter = 1:maxIter
    % Integrate until first y=0 crossing (half-period)
    [~,Y,TE,YE,~] = ode89(@(t,Y) cr3bp_eom_stm(t,Y,mu), [0 10], Y0, opts);
    
    % Extract state at crossing
    X_cross = YE(:,1:nStates)';
    dx_cross = X_cross(3,end);                  % dx at crossing
    Phi_cross = reshape(YE(end,nStates+1:end), nStates, nStates);
    
    % STM element: d(dx_cross)/d(dy0)
    dF_dvy0 = Phi_cross(3,4);
    
    % Newton-Raphson update for vy0
    vy0_new = vy0  - alpha* dx_cross / dF_dvy0;
   
    fprintf('Iter %d: vy0 = %.10f, dx_cross = %.3e\n', iter, vy0_new, dx_cross);
    
    % Check convergence
    if abs(dx_cross) < tol
        fprintf('Converged after %d iterations.\n', iter);
        vy0 = vy0_new;
        break;

    end
    vy0 = vy0_new;
    % Update for next iteration
    Y0 = [x0; y0; vx0; vy0; Phi0(:)];
end
Phi = permute(reshape(Y(:,5:end), 4, 4, []), [3 1 2]);  % Nt × 4 × 4

fprintf('Final initial conditions for Lyapunov orbit:\n');
fprintf('x0 = %.10f, y0 = %.10f, vx0 = %.10f, vy0 = %.10f\n\n', x0, y0, vx0, vy0);
[M, vU, vS] = cr3bp_sys_monodromy(Y);

% Integrate full orbit
[t, X] = ode89(@(t,X) cr3bp_eom_gen(t,X,mu), 0:0.01:2*TE(end), Y0(1:4)');

% for i = 1:length(t)
%     Uv(i,:) = (squeeze(Phi(i,:,:))*vU)';
%     X_s(i,:) = X(i,:) + 2e-5 * Uv(i,:);
%     % Integrate Stable Manifold
%     [t_tmp, X_tmp] = ode45(@(t,X) cr3bp_eom_gen(t,X,mu), 0:0.01:5, X_s(i,:)');
%     fprintf("Integrating Stable manifold at initial condition: \n",string(i));
% 
%     t_sm(i,1:length(t_tmp)) = reshape(t_tmp(:)',[1,length(t_tmp)]);
%     X_sm(i,1:size(X_tmp,1),1:size(X_tmp,2)) = reshape(X_tmp(:,:),[1,length(t_tmp),4]);
% end
% figure(6)
% hold on
% for i = 1:length(t)
%     plot(X_sm(i,:,1),X_sm(i,:,2), 'Color','b')
% end
% axis square

% Rotating to Inertial to Rotating Test
X_In = cr3bp_conv_rot2inert(X,t);
X_rrot = cr3bp_conv_inert2rot(X_In,t);

% Dimensionalization Test
m1 = 5.97219*1e24;
a = 3.844 * 1e8;
[X_rotdim, t_dim] = cr3bp_conv_dim(X_rrot,t,a, 'Primary Mass', [m1, mu]);
[X_rend, t_rend]  = cr3bp_conv_nondim(X_rotdim, t_dim, a, 'Primary Mass', [m1,mu]);

figure(8)
background('Milky Way')
hold on

% opts.Position = [-mu*a;0;0];

moon_opts.Position = [(1-mu)*a;0;0];
moon_opts.Clipping = 'on';
moon_opts.RefPlane = 'equatorial';
% planet3D('Earth Cloudy',opts)
planet3D('Moon',moon_opts)
light('Position',[1,-1,0]);
plot3(X_rotdim(:,1),X_rotdim(:,2),zeros(size(X_rotdim,1),1))

hold off



% Minor differences after redimensionalization noticed, possibly due to
% machine percision, or cutting/rounding error compounding from all the
% operations. Seems to be of acceptable magnitude. Largest Order of 1e-17 
% for X and t.  

abs_err_ndrd_X = mean(abs(X_In - X_rend),1);
abs_err_ndrd_t = mean(abs(t - t_rend));

fprintf("Nondim/Dim conversion erros.\n")
fprintf("Absolute error of %.3e [m] for X\n", abs_err_ndrd_X)
fprintf("\n")
fprintf("Absolute error of %.3e [s] for t vector\n", abs_err_ndrd_t)


% Plotting
figure(1);
subplot(1,3,1)
plot(X(:,1), X(:,2), 'b-');
axis equal;
xlabel('x'); ylabel('y');
title('Corrected Planar Lyapunov Orbit [Rotating]');

subplot(1,3,2)
plot(X_In(:,1), X_In(:,2), 'b-');
axis equal;
xlabel('x'); ylabel('y');
title('Corrected Planar Lyapunov Orbit [Inertial]');

subplot(1,3,3)
plot(X_rrot(:,1), X_rrot(:,2), 'b-');
axis equal;
xlabel('x'); ylabel('y');
title('Corrected Planar Lyapunov Orbit [Rotating Reconstructed]');

% Jacobi Constant Calculation
Cj = cr3bp_sys_jacconst(X,mu);
figure(2);
plot(t, Cj, 'b-');
axis equal;
xlabel('Time'); ylabel('Jacobi Constant');
title('Jacobi Constant of Lyapunov Orbit');

% Zero-Velocity Curve calculation
[x,y,Phi,Forb] = cr3bp_sys_zerovel(mu,Cj(1),[1.5,1.5],3000);

figure(3);hold on
% Forbidden region
contour(x, y, Forb, 1,'DisplayName','Zero-Velocity Curve'); colormap([1,1,1;0.8,0.2,0]);

% More Plotting
plot(X(:,1), X(:,2), 'b-','DisplayName','Trajectory');
plot(-mu,0,'go','MarkerFaceColor','g','DisplayName', 'Earth');
plot(1-mu,0,'ko','MarkerFaceColor','k','DisplayName', 'Moon');
scatter(LP(1),LP(2),'rd','MarkerFaceColor','r','DisplayName','L1')
legend('Location','best')
xlabel('x'); ylabel('y'); title('Planar ZVC and Forbidden Region');
axis square


