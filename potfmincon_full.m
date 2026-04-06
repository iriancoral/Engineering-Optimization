% potfmincon_full.m
% Extended optimization: 3-variable and 4-variable flower pot models.
%
% 3-var: x = [r1, h, r2],     t fixed at t_min
% 4-var: x = [r1, h, r2, t],  t is free
%
% Runs fmincon (SQP) with multiple starting points for both models.
% Compares results with each other and with the 2-variable baseline.
%
% Produces:
%   Figure 1 - 3-var vs 2-var comparison bar chart
%   Figure 2 - 4-var vs 3-var vs 2-var comparison bar chart
%   Figure 3 - Wall angle and thickness at optimum for both models
%   Console  - Detailed results, constraint activity, multipliers
% -------------------------------------------------------

clear;
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
params.f_drain     = f_drain;
params.stab_ratio  = stab_ratio;
params.t_min       = t_min;

V_required = m_soil / rho_soil;

% fmincon options
options = optimoptions('fmincon', ...
    'Algorithm',           'sqp', ...
    'Display',             'off', ...
    'OptimalityTolerance', 1e-10, ...
    'ConstraintTolerance', 1e-10, ...
    'StepTolerance',       1e-12, ...
    'MaxIterations',       500, ...
    'FiniteDifferenceType','forward');

% -------------------------------------------------------
% Load 2-variable baseline for comparison
% -------------------------------------------------------
if exist('potfmincon_result.mat', 'file')
    base = load('potfmincon_result.mat');
    x_2var = base.x_fmincon;
    f_2var = base.f_fmincon;
    r2_2var = x_2var(1) + x_2var(2)*tan(theta_wall*pi/180);
    fprintf('=== 2-Variable Baseline (from potfmincon_result.mat) ===\n');
    fprintf('  r1 = %.4f mm,  h = %.4f mm,  r2 = %.4f mm,  t = %.1f mm\n', ...
        x_2var(1)*1000, x_2var(2)*1000, r2_2var*1000, t*1000);
    fprintf('  Vmat = %.6f cm3\n\n', f_2var*1e6);
    has_2var = true;
else
    fprintf('Warning: potfmincon_result.mat not found. Run potfmincon.m first.\n');
    fprintf('Proceeding without 2-variable baseline.\n\n');
    has_2var = false;
    f_2var = NaN;
end

% =======================================================
% PART 1: 3-VARIABLE MODEL  x = [r1, h, r2], t fixed
% =======================================================
fprintf('============================================\n');
fprintf('  3-VARIABLE MODEL: x = [r1, h, r2]\n');
fprintf('  Wall thickness fixed: t = %.1f mm\n', t*1000);
fprintf('============================================\n\n');

obj3 = @(x) potobj_full(x, params);
con3 = @(x) potcon_full(x, params);

% Bounds for extended model (not in potparams.m — only needed for 3/4-var)
r2_max = 0.40;       % max top radius [m] (400 mm)
t_max  = 0.010;      % max wall thickness [m] (10 mm)

lb3 = [r1_min; h_min; r1_min];
ub3 = [r1_max; h_max; r2_max];

% Generate starting points scaled to soil mass.
% Estimate characteristic radius from volume: r_est ~ (V_required/pi)^(1/3)
r_est = (V_required / pi)^(1/3);
r_est = max(r1_min + 0.01, min(r1_max - 0.01, r_est));
h_est = r_est;   % roughly cubic proportions as starting guess

% Spread starts around the estimate with varying r1, h, and r2/r1 ratios
r1_frac = [0.5 0.7 0.9 1.1 1.3 0.6 1.5 0.8];
h_frac  = [0.7 1.3 1.0 1.5 0.7 1.8 0.9 1.0];
alpha   = [1.2 1.4 1.3 1.5 1.2 1.6 1.3 1.4];   % r2/r1 ratio

x_starts3 = zeros(3, length(r1_frac));
for si = 1:length(r1_frac)
    r1s = max(lb3(1), min(ub3(1), r_est * r1_frac(si)));
    hs  = max(lb3(2), min(ub3(2), h_est * h_frac(si)));
    r2s = max(lb3(3), min(ub3(3), r1s * alpha(si)));
    x_starts3(:, si) = [r1s; hs; r2s];
end

n_starts3 = size(x_starts3, 2);

x_results3  = zeros(3, n_starts3);
f_results3  = zeros(1, n_starts3);
exit_flags3 = zeros(1, n_starts3);

for k = 1:n_starts3
    try
        [xk, fk, eflag] = fmincon(obj3, x_starts3(:,k), [], [], [], [], lb3, ub3, con3, options);
        x_results3(:,k) = xk;
        f_results3(k)   = fk;
        exit_flags3(k)  = eflag;
    catch
        x_results3(:,k) = NaN;
        f_results3(k)   = NaN;
        exit_flags3(k)  = -99;
    end
end

conv3 = exit_flags3 > 0 & ~isnan(f_results3);
fprintf('Converged: %d / %d starts\n', sum(conv3), n_starts3);

[f_best3, idx3] = min(f_results3(conv3));
x_conv3 = x_results3(:, conv3);
x_best3 = x_conv3(:, idx3);

% Rerun best for multipliers
[x_best3, f_best3, ~, ~, lam3] = fmincon(obj3, x_best3, [], [], [], [], lb3, ub3, con3, options);
[g3, ceq3] = potcon_full(x_best3, params);

theta_opt3 = atan2d(x_best3(3) - x_best3(1), x_best3(2));  % wall angle [deg]

fprintf('\n--- 3-Variable Optimum ---\n');
fprintf('  r1         = %.4f mm\n', x_best3(1)*1000);
fprintf('  h          = %.4f mm\n', x_best3(2)*1000);
fprintf('  r2         = %.4f mm\n', x_best3(3)*1000);
fprintf('  t          = %.1f mm  (fixed)\n', t*1000);
fprintf('  theta_wall = %.2f deg  (computed)\n', theta_opt3);
fprintf('  Vmat       = %.6f cm3\n', f_best3*1e6);
fprintf('  Mass       = %.4f g\n', f_best3*rho_mat*1000);
fprintf('  Cost       = %.4f euro\n', f_best3*rho_mat*cost_mat);

fprintf('\nConstraints at optimum:\n');
cnames = {'g1 (stress)', 'g2 (stability)', 'g3 (r2>=r1)'};
tol_active = 1e-4;
for i = 1:3
    status = 'inactive';
    if abs(g3(i)) < tol_active; status = 'ACTIVE'; end
    fprintf('  %s = %+.4e  [%s]\n', cnames{i}, g3(i), status);
end
fprintf('  ceq (volume) = %+.4e\n', ceq3);

fprintf('\nLagrange multipliers:\n');
for i = 1:3
    fprintf('  mu_%d (%s) = %.6e\n', i, cnames{i}, lam3.ineqnonlin(i));
end
fprintf('  lambda_ceq = %.6e\n', lam3.eqnonlin(1));

if has_2var
    savings3 = (f_2var - f_best3) / f_2var * 100;
    fprintf('\nMaterial savings vs 2-var: %.2f%%\n', savings3);
end
fprintf('\n');

% =======================================================
% PART 2: 4-VARIABLE MODEL  x = [r1, h, r2, t]
% =======================================================
fprintf('============================================\n');
fprintf('  4-VARIABLE MODEL: x = [r1, h, r2, t]\n');
fprintf('============================================\n\n');

obj4 = @(x) potobj_full(x, params);
con4 = @(x) potcon_full(x, params);

% Bounds for 4-var: [r1, h, r2, t]
lb4 = [r1_min; h_min; r1_min; t_min];
ub4 = [r1_max; h_max; r2_max;  t_max];

% Starting points: reuse 3-var starts, add t variations
t_frac = [1.0 1.5 1.2 1.0 2.0 1.3 1.5 1.0];   % multiples of t_min
x_starts4 = zeros(4, size(x_starts3, 2));
for si = 1:size(x_starts3, 2)
    x_starts4(1:3, si) = x_starts3(:, si);
    x_starts4(4, si)   = max(lb4(4), min(ub4(4), t_min * t_frac(si)));
end

n_starts4 = size(x_starts4, 2);

x_results4  = zeros(4, n_starts4);
f_results4  = zeros(1, n_starts4);
exit_flags4 = zeros(1, n_starts4);

for k = 1:n_starts4
    try
        [xk, fk, eflag] = fmincon(obj4, x_starts4(:,k), [], [], [], [], lb4, ub4, con4, options);
        x_results4(:,k) = xk;
        f_results4(k)   = fk;
        exit_flags4(k)  = eflag;
    catch
        x_results4(:,k) = NaN;
        f_results4(k)   = NaN;
        exit_flags4(k)  = -99;
    end
end

conv4 = exit_flags4 > 0 & ~isnan(f_results4);
fprintf('Converged: %d / %d starts\n', sum(conv4), n_starts4);

[f_best4, idx4] = min(f_results4(conv4));
x_conv4 = x_results4(:, conv4);
x_best4 = x_conv4(:, idx4);

% Rerun best for multipliers
[x_best4, f_best4, ~, ~, lam4] = fmincon(obj4, x_best4, [], [], [], [], lb4, ub4, con4, options);
[g4, ceq4] = potcon_full(x_best4, params);

theta_opt4 = atan2d(x_best4(3) - x_best4(1), x_best4(2));

fprintf('\n--- 4-Variable Optimum ---\n');
fprintf('  r1         = %.4f mm\n', x_best4(1)*1000);
fprintf('  h          = %.4f mm\n', x_best4(2)*1000);
fprintf('  r2         = %.4f mm\n', x_best4(3)*1000);
fprintf('  t          = %.4f mm\n', x_best4(4)*1000);
fprintf('  theta_wall = %.2f deg  (computed)\n', theta_opt4);
fprintf('  Vmat       = %.6f cm3\n', f_best4*1e6);
fprintf('  Mass       = %.4f g\n', f_best4*rho_mat*1000);
fprintf('  Cost       = %.4f euro\n', f_best4*rho_mat*cost_mat);

fprintf('\nConstraints at optimum:\n');
for i = 1:3
    status = 'inactive';
    if abs(g4(i)) < tol_active; status = 'ACTIVE'; end
    fprintf('  %s = %+.4e  [%s]\n', cnames{i}, g4(i), status);
end
fprintf('  ceq (volume) = %+.4e\n', ceq4);

fprintf('\nLagrange multipliers:\n');
for i = 1:3
    fprintf('  mu_%d (%s) = %.6e\n', i, cnames{i}, lam4.ineqnonlin(i));
end
fprintf('  lambda_ceq = %.6e\n', lam4.eqnonlin(1));

fprintf('\nBound activity (4-var):\n');
bnd_names = {'r1','h','r2','t'};
for i = 1:4
    lb_act = abs(x_best4(i) - lb4(i)) < 1e-6;
    ub_act = abs(x_best4(i) - ub4(i)) < 1e-6;
    if lb_act
        fprintf('  %s = %.4f mm  [AT LOWER BOUND]\n', bnd_names{i}, x_best4(i)*1000);
    elseif ub_act
        fprintf('  %s = %.4f mm  [AT UPPER BOUND]\n', bnd_names{i}, x_best4(i)*1000);
    else
        fprintf('  %s = %.4f mm  [interior]\n', bnd_names{i}, x_best4(i)*1000);
    end
end

if has_2var
    savings4 = (f_2var - f_best4) / f_2var * 100;
    fprintf('\nMaterial savings vs 2-var: %.2f%%\n', savings4);
end
savings4v3 = (f_best3 - f_best4) / f_best3 * 100;
fprintf('Material savings vs 3-var: %.2f%%\n\n', savings4v3);

% =======================================================
% COMPARISON SUMMARY
% =======================================================
fprintf('============================================\n');
fprintf('  MODEL COMPARISON SUMMARY\n');
fprintf('============================================\n');
fprintf('%-12s  %10s  %8s  %8s  %8s  %8s\n', 'Model', 'Vmat[cm3]', 'r1[mm]', 'h[mm]', 'r2[mm]', 't[mm]');
fprintf('%-12s  %10s  %8s  %8s  %8s  %8s\n', '-----', '---------', '------', '-----', '------', '-----');
if has_2var
    fprintf('%-12s  %10.4f  %8.2f  %8.2f  %8.2f  %8.1f\n', '2-var', f_2var*1e6, ...
        x_2var(1)*1000, x_2var(2)*1000, r2_2var*1000, t*1000);
end
fprintf('%-12s  %10.4f  %8.2f  %8.2f  %8.2f  %8.1f\n', '3-var', f_best3*1e6, ...
    x_best3(1)*1000, x_best3(2)*1000, x_best3(3)*1000, t*1000);
fprintf('%-12s  %10.4f  %8.2f  %8.2f  %8.2f  %8.2f\n', '4-var', f_best4*1e6, ...
    x_best4(1)*1000, x_best4(2)*1000, x_best4(3)*1000, x_best4(4)*1000);
fprintf('\n');

% =======================================================
% FIGURES
% =======================================================

% --- Figure 1: Bar comparison across all models ---
figure(1); clf;

if has_2var
    model_names = {'2-var (r_1, h)', '3-var (r_1, h, r_2)', '4-var (r_1, h, r_2, t)'};
    Vmat_vals = [f_2var*1e6, f_best3*1e6, f_best4*1e6];
    mass_vals = [f_2var*rho_mat*1000, f_best3*rho_mat*1000, f_best4*rho_mat*1000];
else
    model_names = {'3-var (r_1, h, r_2)', '4-var (r_1, h, r_2, t)'};
    Vmat_vals = [f_best3*1e6, f_best4*1e6];
    mass_vals = [f_best3*rho_mat*1000, f_best4*rho_mat*1000];
end

subplot(1,2,1);
bar(Vmat_vals, 0.5, 'FaceColor', [0.2 0.5 0.8]);
set(gca, 'XTickLabel', model_names, 'FontSize', 9);
ylabel('V_{mat}^*  [cm^3]');
title('Material volume at optimum');
grid on;

subplot(1,2,2);
bar(mass_vals, 0.5, 'FaceColor', [0.8 0.3 0.2]);
set(gca, 'XTickLabel', model_names, 'FontSize', 9);
ylabel('Mass  [g]');
title('Pot mass at optimum');
grid on;

sgtitle(sprintf('Model comparison  |  m_{soil} = %g kg', m_soil), 'FontSize', 12);
set(gcf, 'Position', [100 100 1000 400]);
saveas(gcf, 'Pictures/potfmincon_full fig(1).png');

% --- Figure 2: Design variable comparison ---
figure(2); clf;

if has_2var
    r1_vals = [x_2var(1)*1000, x_best3(1)*1000, x_best4(1)*1000];
    h_vals  = [x_2var(2)*1000, x_best3(2)*1000, x_best4(2)*1000];
    r2_vals = [r2_2var*1000,   x_best3(3)*1000, x_best4(3)*1000];
    t_vals  = [t*1000,         t*1000,           x_best4(4)*1000];
    theta_vals = [theta_wall, theta_opt3, theta_opt4];
else
    r1_vals = [x_best3(1)*1000, x_best4(1)*1000];
    h_vals  = [x_best3(2)*1000, x_best4(2)*1000];
    r2_vals = [x_best3(3)*1000, x_best4(3)*1000];
    t_vals  = [t*1000,          x_best4(4)*1000];
    theta_vals = [theta_opt3, theta_opt4];
end

subplot(2,2,1);
bar(r1_vals, 0.5, 'FaceColor', [0.2 0.6 0.4]);
set(gca, 'XTickLabel', model_names, 'FontSize', 8);
ylabel('r_1^* [mm]'); title('Bottom radius'); grid on;

subplot(2,2,2);
bar(h_vals, 0.5, 'FaceColor', [0.6 0.3 0.6]);
set(gca, 'XTickLabel', model_names, 'FontSize', 8);
ylabel('h^* [mm]'); title('Height'); grid on;

subplot(2,2,3);
bar(r2_vals, 0.5, 'FaceColor', [0.8 0.6 0.2]);
set(gca, 'XTickLabel', model_names, 'FontSize', 8);
ylabel('r_2^* [mm]'); title('Top radius'); grid on;

subplot(2,2,4);
bar(t_vals, 0.5, 'FaceColor', [0.3 0.5 0.8]);
set(gca, 'XTickLabel', model_names, 'FontSize', 8);
ylabel('t^* [mm]'); title('Wall thickness'); grid on;

sgtitle(sprintf('Optimal design variables  |  m_{soil} = %g kg', m_soil), 'FontSize', 12);
set(gcf, 'Position', [100 100 900 600]);
saveas(gcf, 'Pictures/potfmincon_full fig(2).png');

% --- Figure 3: Wall angle comparison ---
figure(3); clf;
bar(theta_vals, 0.5, 'FaceColor', [0.5 0.7 0.3]);
set(gca, 'XTickLabel', model_names, 'FontSize', 9);
ylabel('\theta_{wall}^* [deg]');
title(sprintf('Optimal wall angle  |  m_{soil} = %g kg', m_soil), 'FontSize', 12);
grid on;
set(gcf, 'Position', [100 100 600 400]);
saveas(gcf, 'Pictures/potfmincon_full fig(3).png');

% Save results
x_3var = x_best3; f_3var = f_best3;
x_4var = x_best4; f_4var = f_best4;
save('potfmincon_full_result.mat', 'x_3var', 'f_3var', 'lam3', ...
     'x_4var', 'f_4var', 'lam4', 'params');
fprintf('Results saved to potfmincon_full_result.mat\n');
fprintf('Done. See Figures 1-3.\n');
