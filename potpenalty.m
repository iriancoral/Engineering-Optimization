% potpenalty.m
% Step 7: Self-implemented exterior penalty method for 2-variable pot problem.
%
% Converts the constrained problem into a sequence of unconstrained problems:
%   P(x, rp) = f(x) + rp * [ sum(max(0,gi)^2) + ceq^2 ]
%
% Strategy:
%   - Start with low penalty rp, solve unconstrained problem with fminsearch
%   - Increase rp by factor each iteration
%   - Track convergence of x* and f* as rp grows
%
% Produces:
%   Figure 1 - Penalty function landscape at selected rp values
%   Figure 2 - Convergence of x* and f* vs penalty parameter rp
%   Figure 3 - Optimization path overlaid on contour map
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

V_required = m_soil / rho_soil;

% -------------------------------------------------------
% Penalty function definition
% P(x, rp) = f(x) + rp * [ sum_i max(0,gi)^2 + ceq^2 ]
% -------------------------------------------------------
penalty_fun = @(x, rp) penalty_eval(x, rp, params);

% -------------------------------------------------------
% Penalty method settings
% -------------------------------------------------------
rp_start  = 1e2;       % initial penalty parameter
rp_factor = 10;        % multiply rp by this each iteration
n_iter    = 10;         % number of penalty iterations

% Starting point (feasible-ish interior point)
x_start = [0.12; 0.20];

% Bounds for fminsearch region (passed as barrier, not hard bounds)
lb = [r1_min; h_min];
ub = [r1_max; h_max];

fprintf('=== Penalty Method ===\n');
fprintf('Starting point: r1 = %.1f mm,  h = %.1f mm\n', x_start(1)*1000, x_start(2)*1000);
fprintf('rp start = %.0e,  factor = %g,  iterations = %d\n\n', rp_start, rp_factor, n_iter);

% Storage
rp_history   = zeros(1, n_iter);
x_history    = zeros(2, n_iter);
f_history    = zeros(1, n_iter);
fpen_history = zeros(1, n_iter);
viol_history = zeros(1, n_iter);

options = optimset('fminsearch');
options = optimset(options, 'MaxFunEvals', 5000, 'MaxIter', 2000, 'TolFun', 1e-10, 'TolX', 1e-10, 'Display', 'off');

x_cur = x_start;
rp    = rp_start;

for iter = 1:n_iter
    % Solve unconstrained penalty problem
    obj_pen = @(x) penalty_fun(x, rp);
    [x_opt, fpen_opt] = fminsearch(obj_pen, x_cur, options);

    % Clamp to bounds (fminsearch has no bounds, add soft barrier instead)
    x_opt = max(lb, min(ub, x_opt));

    % Evaluate true objective and constraint violation
    f_true       = potobj(x_opt, params);
    [g_val, ceq_val] = potcon(x_opt, params);
    viol = sum(max(0, g_val).^2) + ceq_val^2;

    rp_history(iter)   = rp;
    x_history(:, iter) = x_opt;
    f_history(iter)    = f_true;
    fpen_history(iter) = fpen_opt;
    viol_history(iter) = viol;

    fprintf('iter %d:  rp = %8.1e  |  r1 = %6.2f mm  h = %6.2f mm  |  Vmat = %.4f cm3  |  viol = %.2e\n', ...
        iter, rp, x_opt(1)*1000, x_opt(2)*1000, f_true*1e6, viol);

    % Use current solution as warm start for next iteration
    x_cur = x_opt;
    rp    = rp * rp_factor;
end

x_pen = x_history(:, end);
f_pen = f_history(end);

fprintf('\n--- Penalty Method Result ---\n');
fprintf('  r1   = %.4f mm\n',  x_pen(1)*1000);
fprintf('  h    = %.4f mm\n',  x_pen(2)*1000);
fprintf('  r2   = %.4f mm\n',  (x_pen(1) + x_pen(2)*tan(theta_wall*pi/180))*1000);
fprintf('  Vmat = %.6f cm3\n', f_pen*1e6);
[g_final, ceq_final] = potcon(x_pen, params);
fprintf('  g1 (stress)    = %.4e  (<=0 ok)\n', g_final(1));
fprintf('  g2 (stability) = %.4e  (<=0 ok)\n', g_final(2));
fprintf('  ceq (volume)   = %.4e\n\n',          ceq_final);

% -------------------------------------------------------
% FIGURE 1: Penalty landscape at two rp values
% -------------------------------------------------------
Nr = 80; Nh = 80;
r1_vec = linspace(r1_min, r1_max, Nr);
h_vec2 = linspace(h_min,  h_max,  Nh);
[R1g, Hg] = meshgrid(r1_vec, h_vec2);

rp_plot = [rp_history(2), rp_history(end)];
fig_titles = {sprintf('Penalty landscape  rp = %.0e', rp_plot(1)), ...
              sprintf('Penalty landscape  rp = %.0e', rp_plot(2))};


%Shows how the penalty function landscape evolves as rp increases, with the optimization path marked. 
% The first plot shows a relatively flat landscape where the optimizer can explore more freely, 
% while the second plot shows a much steeper landscape where constraint violations are heavily penalized, 
% guiding the optimizer towards the feasible region and optimal solution.
%The colors represent the value of P(x, rp) — the penalized objective — at each (r1, h) point on the grid.


figure(1); clf;
for p = 1:2
    Pgrid = zeros(Nh, Nr);
    for i = 1:Nh
        for j = 1:Nr
            Pgrid(i,j) = penalty_fun([R1g(i,j); Hg(i,j)], rp_plot(p));
        end
    end
    subplot(1,2,p);
    contourf(r1_vec*1000, h_vec2*1000, Pgrid, 25, 'EdgeColor', 'none');
    colorbar;
    hold on;
    % Mark optimum found at this rp
    plot(x_history(1,p+1)*1000, x_history(2,p+1)*1000, 'r*', 'MarkerSize', 12, 'LineWidth', 2);
    xlabel('r_1 [mm]'); ylabel('h [mm]');
    title(fig_titles{p}, 'FontSize', 10);
    colormap(parula);
    grid on;
end
sgtitle('Penalty function P(x, rp) at two penalty values', 'FontSize', 12);

% -------------------------------------------------------
% FIGURE 2: Convergence of x* and Vmat vs rp
% -------------------------------------------------------
figure(2); clf;
subplot(3,1,1);
semilogx(rp_history, x_history(1,:)*1000, 'b-o', 'LineWidth', 1.5);
ylabel('r_1^*  [mm]', 'FontSize', 10);
title('Penalty method convergence', 'FontSize', 12);
grid on;

subplot(3,1,2);
semilogx(rp_history, x_history(2,:)*1000, 'r-s', 'LineWidth', 1.5);
ylabel('h^*  [mm]', 'FontSize', 10);
grid on;

subplot(3,1,3);
semilogx(rp_history, f_history*1e6, 'k-^', 'LineWidth', 1.5);
ylabel('V_{mat}^*  [cm^3]', 'FontSize', 10);
xlabel('Penalty parameter  r_p', 'FontSize', 10);
grid on;

% -------------------------------------------------------
% FIGURE 3: Optimization path on contour map
% -------------------------------------------------------
% Recompute objective grid
Vmat_grid = zeros(Nh, Nr);
ceq_grid  = zeros(Nh, Nr);
for i = 1:Nh
    for j = 1:Nr
        [Vm, Vp, ~, ~, ~, ~, ~] = potanalysis( ...
            R1g(i,j), Hg(i,j), theta_wall, t, rho_mat, rho_soil, K0, g_acc, sigma_allow, cost_mat, f_drain);
        Vmat_grid(i,j) = Vm * 1e6;
        ceq_grid(i,j)  = Vp - V_required;
    end
end

figure(3); clf;
hold on;
Vmat_levels = linspace(min(Vmat_grid(:)), max(Vmat_grid(:)), 20);
contour(r1_vec*1000, h_vec2*1000, Vmat_grid, Vmat_levels, 'LineWidth', 0.8);
colormap(parula); colorbar;
contour(r1_vec*1000, h_vec2*1000, ceq_grid, [0 0], 'k-', 'LineWidth', 2.5);

% Optimization path with start point
plot(x_start(1)*1000, x_start(2)*1000, 'gs', 'MarkerSize', 12, 'MarkerFaceColor', 'g', ...
    'DisplayName', 'Start');
plot([x_start(1); x_history(1,:)']*1000, [x_start(2); x_history(2,:)']*1000, ...
    'r-', 'LineWidth', 1.5, 'DisplayName', 'Penalty path');

% Number each iteration
for k = 1:n_iter
    plot(x_history(1,k)*1000, x_history(2,k)*1000, 'ro', 'MarkerSize', 8, ...
        'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
    text(x_history(1,k)*1000 + 3, x_history(2,k)*1000 + 3, sprintf('%d', k), ...
        'FontSize', 9, 'FontWeight', 'bold', 'Color', 'r');
end

% Mark final optimum
plot(x_pen(1)*1000, x_pen(2)*1000, 'p', 'MarkerSize', 14, 'Color', [0 0.8 0.8], 'MarkerFaceColor', [0 0.8 0.8], ...
    'DisplayName', sprintf('Penalty optimum (%.2f cm^3)', f_pen*1e6));

legend('Location', 'northeast', 'FontSize', 9);
xlabel('r_1  [mm]', 'FontSize', 11);
ylabel('h  [mm]',   'FontSize', 11);
title('Penalty method path on V_{mat} contours', 'FontSize', 12);
axis([r1_min*1000 r1_max*1000 h_min*1000 h_max*1000]);
grid on;
hold off;

fprintf('Done. See Figures 1-3.\n');

% -------------------------------------------------------
% LOCAL FUNCTION: penalty evaluation
% -------------------------------------------------------
function P = penalty_eval(x, rp, params)
    f           = potobj(x, params);
    [g, ceq]    = potcon(x, params);

    % Exterior penalty: only violated constraints contribute
    pen_ineq = sum(max(0, g).^2);
    pen_eq   = ceq^2;

    P = f + rp * (pen_ineq + pen_eq);
end

% Let's use concrete numbers.

% Say fminsearch is trying two points:

% Point A: r1 = 80 mm, h = 100 mm

% Vmat = 150 cm^3 (low material — looks great!)
% But Vpot is too small for the soil → ceq = -0.002
% Stress is fine → g1 = -5e6 (negative = satisfied)
% Stability fine → g2 = -0.14 (satisfied)
% Point B: r1 = 100 mm, h = 150 mm

% Vmat = 210 cm^3 (more material — looks worse)
% Vpot matches the soil → ceq = 0
% All constraints satisfied
% Without penalties, fminsearch picks A because 150 < 210. But A is useless — the pot 
% can't hold the soil.

% Now add the penalty with rp = 100:
% P(A) = 150e-6 + 100 * [0 + (-0.002)^2] = 150e-6 + 100 * 4e-6 = 150e-6 + 4e-4 = 5.5e-4
% P(B) = 210e-6 + 100 * [0 + 0] = 210e-6
% A is still cheaper! rp = 100 isn't enough to punish the violation.
% With rp = 10^6:
% P(A) = 150e-6 + 10^6 * 4e-6 = 150e-6 + 4.0 = 4.0
% P(B) = 210e-6 + 0 = 210e-6

% Now B wins massively. The penalty makes A's violation so expensive that the optimizer avoids it.
% That's the whole idea: the penalty term rp * ceq^2 makes infeasible points a
% rtificially expensive. As rp grows, any point that violates constraints becomes worse than any 
% feasible point, so the optimizer is forced toward feasible solutions.