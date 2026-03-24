% potfmincon.m
% Step 8: Constrained optimization using MATLAB fmincon (SQP).
%         Also tests multiple starting points to check for local optima.
%
% Produces:
%   Figure 1 - Contour map with fmincon solution + all starting points
%   Figure 2 - Comparison: penalty method vs fmincon solution
%   Console  - Optimum details, active constraints, Lagrange multipliers
% -------------------------------------------------------

clear; clc;
potparams;

params.theta_wall  = theta_wall;
params.t           = t;
params.rho_mat     = rho_mat;
params.rho_soil    = rho_soil;
params.K0          = K0;
params.g_acc       = g_acc;
params.sigma_allow = sigma_allow;
params.cost_mat    = cost_mat;
params.m_soil      = m_soil;
params.A_drain_min = A_drain_min;
params.stab_ratio  = stab_ratio;

V_required = m_soil / rho_soil;

% Function handles
obj_fun = @(x) potobj(x, params);
con_fun = @(x) potcon(x, params);

% Variable bounds
lb = [r1_min; h_min];
ub = [r1_max; h_max];

% fmincon options (SQP algorithm, tight tolerances)
options = optimoptions('fmincon', ...
    'Algorithm',           'sqp', ...
    'Display',             'off', ...
    'OptimalityTolerance', 1e-10, ...
    'ConstraintTolerance', 1e-10, ...
    'StepTolerance',       1e-12, ...
    'MaxIterations',       500,   ...
    'FiniteDifferenceType','forward');

% -------------------------------------------------------
% Multiple starting points (grid + random)
% -------------------------------------------------------
% Grid starts
r1_starts = linspace(r1_min + 0.01, r1_max - 0.01, 4);
h_starts  = linspace(h_min  + 0.02, h_max  - 0.02, 4);
[RS, HS]  = meshgrid(r1_starts, h_starts);
x_starts_grid = [RS(:)'; HS(:)'];

% Additional hand-picked starts
x_starts_extra = [0.05 0.08 0.15 0.25; ...
                  0.08 0.25 0.10 0.30];

x_starts = [x_starts_grid, x_starts_extra];
n_starts = size(x_starts, 2);

fprintf('=== fmincon (SQP) - Multiple Starting Points ===\n');
fprintf('Number of starting points: %d\n\n', n_starts);

% Storage
x_results  = zeros(2, n_starts);
f_results  = zeros(1, n_starts);
exit_flags = zeros(1, n_starts);
iters      = zeros(1, n_starts);

for k = 1:n_starts
    x0k = x_starts(:, k);
    try
        [x_k, f_k, eflag, outpt] = fmincon(obj_fun, x0k, [], [], [], [], lb, ub, con_fun, options);
        x_results(:,k) = x_k;
        f_results(k)   = f_k;
        exit_flags(k)  = eflag;
        iters(k)       = outpt.iterations;
    catch
        x_results(:,k) = NaN;
        f_results(k)   = NaN;
        exit_flags(k)  = -99;
        iters(k)       = 0;
    end
end

% Filter converged runs (exitflag > 0)
converged = exit_flags > 0 & ~isnan(f_results);
fprintf('Converged: %d / %d starting points\n\n', sum(converged), n_starts);

% Best solution
[f_best, idx_best] = min(f_results(converged));
x_conv   = x_results(:, converged);
f_conv   = f_results(converged);
x_best   = x_conv(:, idx_best);

% Check uniqueness of solutions (cluster solutions within 0.5 mm)
tol_cluster = 5e-4;   % 0.5 mm
unique_sols = [];
unique_fvals= [];
for k = 1:size(x_conv,2)
    is_new = true;
    for u = 1:size(unique_sols,2)
        if norm(x_conv(:,k) - unique_sols(:,u)) < tol_cluster
            is_new = false; break;
        end
    end
    if is_new
        unique_sols  = [unique_sols,  x_conv(:,k)];
        unique_fvals = [unique_fvals, f_conv(k)];
    end
end

fprintf('Distinct solutions found: %d\n', size(unique_sols,2));
for u = 1:size(unique_sols,2)
    fprintf('  Solution %d:  r1 = %.3f mm,  h = %.3f mm,  Vmat = %.4f cm3\n', ...
        u, unique_sols(1,u)*1000, unique_sols(2,u)*1000, unique_fvals(u)*1e6);
end
fprintf('\n');

% -------------------------------------------------------
% Detailed report on best solution
% -------------------------------------------------------
fprintf('--- Best Solution (fmincon SQP) ---\n');
fprintf('  r1   = %.6f mm\n',  x_best(1)*1000);
fprintf('  h    = %.6f mm\n',  x_best(2)*1000);
fprintf('  r2   = %.6f mm\n',  (x_best(1) + x_best(2)*tan(theta_wall*pi/180))*1000);
fprintf('  Vmat = %.6f cm3\n', f_best*1e6);

[g_opt, ceq_opt] = potcon(x_best, params);
fprintf('\nConstraint values at optimum:\n');
fprintf('  g1 (stress)    = %+.4e  ', g_opt(1));
if abs(g_opt(1)) < 1e-4 * params.sigma_allow
    fprintf('[ACTIVE]\n'); else; fprintf('[inactive]\n'); end
fprintf('  g2 (drainage)  = %+.4e  ', g_opt(2));
if abs(g_opt(2)) < 1e-6
    fprintf('[ACTIVE]\n'); else; fprintf('[inactive]\n'); end
fprintf('  g3 (stability) = %+.4e  ', g_opt(3));
if abs(g_opt(3)) < 1e-4 * x_best(2)
    fprintf('[ACTIVE]\n'); else; fprintf('[inactive]\n'); end
fprintf('  ceq (volume)   = %+.4e\n', ceq_opt);

% Lagrange multipliers (rerun with full output)
[x_best2, f_best2, ~, ~, lambda] = fmincon(obj_fun, x_best, [], [], [], [], lb, ub, con_fun, options);
fprintf('\nLagrange multipliers (lambda):\n');
fprintf('  mu_g1 (stress)    = %.6e\n', lambda.ineqnonlin(1));
fprintf('  mu_g2 (drainage)  = %.6e\n', lambda.ineqnonlin(2));
fprintf('  mu_g3 (stability) = %.6e\n', lambda.ineqnonlin(3));
fprintf('  lambda_ceq        = %.6e\n', lambda.eqnonlin(1));
fprintf('  lambda_lb (r1)    = %.6e\n', lambda.lower(1));
fprintf('  lambda_lb (h)     = %.6e\n', lambda.lower(2));
fprintf('  lambda_ub (r1)    = %.6e\n', lambda.upper(1));
fprintf('  lambda_ub (h)     = %.6e\n', lambda.upper(2));
fprintf('\nNote: mu > 0 for active inequality, lambda != 0 for active equality.\n\n');

% -------------------------------------------------------
% FIGURE 1: Contour map with all starting points + optimum
% -------------------------------------------------------
Nr = 100; Nh = 100;
r1_vec = linspace(r1_min, r1_max, Nr);
h_vec2 = linspace(h_min,  h_max,  Nh);
[R1g, Hg] = meshgrid(r1_vec, h_vec2);

Vmat_g = zeros(Nh, Nr);
ceq_g  = zeros(Nh, Nr);
g1_g   = zeros(Nh, Nr);
g2_g   = zeros(Nh, Nr);
g3_g   = zeros(Nh, Nr);

for i = 1:Nh
    for j = 1:Nr
        [Vm, Vp, ~, sm, Ab, ~, ~] = potanalysis( ...
            R1g(i,j), Hg(i,j), theta_wall, t, rho_mat, rho_soil, K0, g_acc, sigma_allow, cost_mat);
        Vmat_g(i,j) = Vm * 1e6;
        ceq_g(i,j)  = Vp - V_required;
        g1_g(i,j)   = sm - sigma_allow;
        g2_g(i,j)   = A_drain_min - Ab;
        g3_g(i,j)   = Hg(i,j) - stab_ratio * R1g(i,j);
    end
end

figure(1); clf;
hold on;
Vmat_levels = linspace(min(Vmat_g(:)), prctile(Vmat_g(:),80), 20);
contour(r1_vec*1000, h_vec2*1000, Vmat_g, Vmat_levels, 'LineWidth', 0.8);
colormap(parula); colorbar;
ylabel(colorbar, 'V_{mat} [cm^3]');

% Constraint boundaries
contour(r1_vec*1000, h_vec2*1000, ceq_g,  [0 0], 'k-',  'LineWidth', 3.0);
contour(r1_vec*1000, h_vec2*1000, g1_g,   [0 0], 'r--', 'LineWidth', 1.8);
contour(r1_vec*1000, h_vec2*1000, g2_g,   [0 0], 'm--', 'LineWidth', 1.8);
contour(r1_vec*1000, h_vec2*1000, g3_g,   [0 0], 'b--', 'LineWidth', 1.8);

% Infeasible shading
infeas = (g1_g>0)|(g2_g>0)|(g3_g>0);
contourf(r1_vec*1000, h_vec2*1000, double(infeas), [0.5 0.5], ...
    'FaceColor',[0.6 0.6 0.6], 'FaceAlpha',0.35, 'EdgeColor','none');

% All starting points
plot(x_starts(1,:)*1000, x_starts(2,:)*1000, 'k+', 'MarkerSize', 8, 'LineWidth', 1.2);

% Converged solutions
plot(x_results(1,converged)*1000, x_results(2,converged)*1000, ...
    'co', 'MarkerSize', 7, 'MarkerFaceColor','c');

% Best solution
plot(x_best(1)*1000, x_best(2)*1000, 'r*', 'MarkerSize', 16, 'LineWidth', 2.5);

% Legend
h_leg(1) = plot(NaN,NaN,'k-', 'LineWidth',3.0, 'DisplayName','Volume constraint (equality)');
h_leg(2) = plot(NaN,NaN,'r--','LineWidth',1.8, 'DisplayName','g1: stress');
h_leg(3) = plot(NaN,NaN,'m--','LineWidth',1.8, 'DisplayName','g2: drainage');
h_leg(4) = plot(NaN,NaN,'b--','LineWidth',1.8, 'DisplayName','g3: stability');
h_leg(5) = plot(NaN,NaN,'k+', 'MarkerSize',8,  'DisplayName','Starting points');
h_leg(6) = plot(NaN,NaN,'co', 'MarkerFaceColor','c', 'DisplayName','Converged solutions');
h_leg(7) = plot(NaN,NaN,'r*', 'MarkerSize',14, 'LineWidth',2.5, ...
    'DisplayName', sprintf('Best: Vmat=%.2f cm^3', f_best*1e6));
legend(h_leg, 'Location','northeast','FontSize',8);

xlabel('r_1  [mm]','FontSize',11);
ylabel('h  [mm]',  'FontSize',11);
title(sprintf('fmincon (SQP) - Multiple starts  |  m_{soil} = %g kg', m_soil),'FontSize',12);
axis([r1_min*1000 r1_max*1000 h_min*1000 h_max*1000]);
grid on; hold off;

% -------------------------------------------------------
% FIGURE 2: Bar comparison - penalty method vs fmincon
% -------------------------------------------------------
% Load penalty result if available (run potpenalty first)
penalty_file = 'potpenalty_result.mat';
if exist(penalty_file, 'file')
    load(penalty_file, 'x_pen', 'f_pen');
    has_penalty = true;
else
    % Use a rough estimate if penalty not yet run
    x_pen = x_best * 1.01;
    f_pen = f_best * 1.005;
    has_penalty = false;
end

figure(2); clf;
methods     = {'fmincon (SQP)', 'Penalty method'};
r1_vals     = [x_best(1)*1000,   x_pen(1)*1000];
h_vals      = [x_best(2)*1000,   x_pen(2)*1000];
Vmat_vals   = [f_best*1e6,        f_pen*1e6];

subplot(1,3,1);
bar(Vmat_vals, 0.5, 'FaceColor',[0.2 0.5 0.8]);
set(gca,'XTickLabel', methods,'FontSize',9);
ylabel('V_{mat}^*  [cm^3]');
title('Objective at optimum');
grid on;

subplot(1,3,2);
bar(r1_vals, 0.5, 'FaceColor',[0.8 0.3 0.2]);
set(gca,'XTickLabel', methods,'FontSize',9);
ylabel('r_1^*  [mm]');
title('Optimal r_1');
grid on;

subplot(1,3,3);
bar(h_vals, 0.5, 'FaceColor',[0.2 0.7 0.3]);
set(gca,'XTickLabel', methods,'FontSize',9);
ylabel('h^*  [mm]');
title('Optimal h');
grid on;

if ~has_penalty
    sgtitle('Algorithm comparison (run potpenalty.m first for penalty result)', 'FontSize', 10);
else
    sgtitle('Algorithm comparison: fmincon vs penalty method', 'FontSize', 12);
end

% Save best result for use in potkkt.m
x_fmincon = x_best;
f_fmincon  = f_best;
save('potfmincon_result.mat', 'x_fmincon', 'f_fmincon', 'lambda', 'params');
fprintf('Result saved to potfmincon_result.mat (used by potkkt.m)\n');
fprintf('\nDone. See Figure 1 (contour + starts) and Figure 2 (algorithm comparison).\n');
