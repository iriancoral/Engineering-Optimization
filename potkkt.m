% potkkt.m
% Step 9: KKT condition verification at the optimum found by fmincon.
%
% Checks all four KKT conditions:
%   1. Stationarity:         grad_f + sum(mu_i * grad_gi) + lambda * grad_ceq = 0
%   2. Primal feasibility:   g(x*) <= 0,  ceq(x*) = 0
%   3. Dual feasibility:     mu_i >= 0  for all i
%   4. Complementary slack:  mu_i * g_i(x*) = 0  for all i
%
% Also verifies KKT at a NON-optimal feasible point for contrast.
%
% Requires potfmincon_result.mat (produced by potfmincon.m).
% -------------------------------------------------------

%clc;
potparams;

% Requires potfmincon.m to be run first (variables in workspace)
x_fmincon = x_best;
f_fmincon = f_best;

x_opt = x_fmincon;
h_fd  = 1e-8;    % FD step size (from sensitivity analysis)

fprintf('=== KKT Condition Verification ===\n');
fprintf('Point: r1 = %.4f mm,  h = %.4f mm\n\n', x_opt(1)*1000, x_opt(2)*1000);

% -------------------------------------------------------
% Helper: compute FD gradients of f, g, ceq at a point x
% -------------------------------------------------------
function [gf, gg, gceq] = fd_gradients(x, params, h_step)
    f0         = potobj(x, params);
    [g0, ceq0] = potcon(x, params);
    n          = length(x);
    ng         = length(g0);

    gf   = zeros(n, 1);
    gg   = zeros(n, ng);
    gceq = zeros(n, 1);

    for k = 1:n
        xp      = x;
        xp(k)   = x(k) + h_step;
        fp           = potobj(xp, params);
        [gp, ceqp]   = potcon(xp, params);
        gf(k)        = (fp      - f0)    / h_step;
        gg(k,:)      = (gp      - g0)    / h_step;
        gceq(k)      = (ceqp    - ceq0)  / h_step;
    end
end

% -------------------------------------------------------
% Evaluate at optimum x*
% -------------------------------------------------------
[grad_f, grad_g, grad_ceq] = fd_gradients(x_opt, params, h_fd);
[g_opt, ceq_opt]           = potcon(x_opt, params);

fprintf('Gradients at x* (forward FD, h = %.0e):\n', h_fd);
fprintf('  grad f    = [%+.6e,  %+.6e]\n', grad_f(1),      grad_f(2));
fprintf('  grad g1   = [%+.6e,  %+.6e]\n', grad_g(1,1),    grad_g(2,1));
fprintf('  grad g2   = [%+.6e,  %+.6e]\n', grad_g(1,2),    grad_g(2,2));
fprintf('  grad ceq  = [%+.6e,  %+.6e]\n', grad_ceq(1),    grad_ceq(2));
fprintf('\n');

% -------------------------------------------------------
% Identify active inequality constraints (|g| < tol)
% -------------------------------------------------------
tol_active = 1e-4;
active_ineq = abs(g_opt) < tol_active;

fprintf('Constraint values at x*:\n');
cnames = {'g1 (stress)', 'g2 (stability)'};
for i = 1:2
    status = 'inactive';
    if active_ineq(i); status = 'ACTIVE'; end
    fprintf('  %s = %+.4e  [%s]\n', cnames{i}, g_opt(i), status);
end
fprintf('  ceq (volume) = %+.4e\n\n', ceq_opt);

% -------------------------------------------------------
% KKT CONDITION 1: Stationarity
% Solve for multipliers from active constraints:
%   grad_f + G_active * mu_active + grad_ceq * lambda_ceq = 0
%
% System: [grad_g_active | grad_ceq] * [mu_active; lambda_ceq] = -grad_f
% -------------------------------------------------------
fprintf('--- KKT Check 1: Stationarity ---\n');

% Build system with all active inequality constraints + equality
idx_active = find(active_ineq);
A_sys = [grad_g(:, idx_active), grad_ceq];   % 2 x (n_active+1)
b_sys = -grad_f;                              % 2 x 1

if size(A_sys, 2) == size(A_sys, 1)
    % Square system: direct solve
    mu_lam = A_sys \ b_sys;
    solve_method = 'direct solve';
else
    % Over/underdetermined: use pseudo-inverse
    mu_lam = pinv(A_sys) * b_sys;
    solve_method = 'pseudo-inverse (system not square)';
end

fprintf('Solve method: %s\n', solve_method);

% Separate mu (for active g) and lambda_ceq
mu_active  = mu_lam(1:end-1);
lambda_ceq = mu_lam(end);

% Full mu vector (inactive constraints get mu=0)
mu_full = zeros(2,1);
mu_full(idx_active) = mu_active;

fprintf('Multipliers computed from stationarity:\n');
for i = 1:2
    fprintf('  mu_%d (%s) = %+.6e\n', i, cnames{i}, mu_full(i));
end
fprintf('  lambda_ceq          = %+.6e\n', lambda_ceq);

% Check residual: grad_f + grad_g * mu + grad_ceq * lambda = 0
residual = grad_f + grad_g * mu_full + grad_ceq * lambda_ceq;
fprintf('Stationarity residual = [%+.4e,  %+.4e]\n', residual(1), residual(2));
res_norm = norm(residual);
fprintf('Residual norm = %.4e  ', res_norm);
if res_norm < 1e-4
    fprintf('[SATISFIED]\n\n');
else
    fprintf('[VIOLATED - check active set]\n\n');
end

% Compare with fmincon multipliers
fprintf('fmincon multipliers for comparison:\n');
fprintf('  mu_1 (fmincon) = %+.6e\n', lambda.ineqnonlin(1));
fprintf('  mu_2 (fmincon) = %+.6e\n', lambda.ineqnonlin(2));
fprintf('  lam_ceq        = %+.6e\n\n', lambda.eqnonlin(1));

% -------------------------------------------------------
% KKT CONDITION 2: Primal feasibility
% -------------------------------------------------------
fprintf('--- KKT Check 2: Primal Feasibility ---\n');
pf_ineq = all(g_opt <= tol_active);
pf_eq   = abs(ceq_opt) < 1e-6;
fprintf('  g(x*) <= 0:   %s\n', yn(pf_ineq));
fprintf('  ceq(x*) = 0:  %s  (|ceq| = %.2e)\n\n', yn(pf_eq), abs(ceq_opt));

% -------------------------------------------------------
% KKT CONDITION 3: Dual feasibility
% -------------------------------------------------------
fprintf('--- KKT Check 3: Dual Feasibility (mu >= 0) ---\n');
for i = 1:2
    ok = mu_full(i) >= -1e-8;
    fprintf('  mu_%d = %+.6e  [%s]\n', i, mu_full(i), yn(ok));
end
df_ok = all(mu_full >= -1e-8);
fprintf('  Dual feasibility: %s\n\n', yn(df_ok));

% -------------------------------------------------------
% KKT CONDITION 4: Complementary slackness
% -------------------------------------------------------

fprintf('--- KKT Check 4: Complementary Slackness (mu_i * g_i = 0) ---\n');
cs_vals = mu_full .* g_opt(:);
for i = 1:2
    ok = abs(cs_vals(i)) < 1e-8;
    fprintf('  mu_%d * g_%d = %+.4e  [%s]\n', i, i, cs_vals(i), yn(ok));
end
cs_ok = all(abs(cs_vals) < 1e-6);
fprintf('  Complementary slackness: %s\n\n', yn(cs_ok));

% -------------------------------------------------------
% Summary
% -------------------------------------------------------
all_ok = all([(res_norm < 1e-4), all(pf_ineq), all(pf_eq), all(df_ok), all(cs_ok)]);
fprintf('=== KKT Summary at x* ===\n');
fprintf('  1. Stationarity:          %s\n', yn(res_norm < 1e-4));
fprintf('  2. Primal feasibility:    %s\n', yn(all(pf_ineq) && all(pf_eq)));
fprintf('  3. Dual feasibility:      %s\n', yn(df_ok));
fprintf('  4. Complementary slack:   %s\n', yn(cs_ok));
fprintf('  --> x* is a KKT point:    %s\n\n', yn(all_ok));



% OVERBODIG HIERONDER


% -------------------------------------------------------
% Verify KKT at a NON-optimal feasible point (for contrast)
% -------------------------------------------------------
fprintf('=== KKT Check at a Non-Optimal Feasible Point ===\n');

% Find a point on the volume curve at a different r1 that satisfies g <= 0
r1_test = x_opt(1) * 1.4;    % larger r1 -> different h from volume eq.
th      = params.theta_wall * pi / 180;
V_req   = params.m_soil / params.rho_soil;
% Solve cubic for h_test given r1_test and fixed wall angle
a   = tan(th);
p3  = pi*a^2/3; p2 = pi*a*r1_test; p1 = pi*r1_test^2; p0 = -V_req;
h_roots = roots([p3, p2, p1, p0]);
h_real  = h_roots(imag(h_roots)==0 & real(h_roots)>0);
h_test  = real(h_real(1));
r2_test = r1_test + h_test * tan(th);
A_coef  = (pi/3) * (r1_test^2 + r1_test*r2_test + r2_test^2);
V_req   = params.m_soil / params.rho_soil;
h_test  = V_req / A_coef;

x_test = [r1_test; h_test];
fprintf('Test point: r1 = %.2f mm,  h = %.2f mm\n', x_test(1)*1000, x_test(2)*1000);

[g_test, ceq_test] = potcon(x_test, params);
fprintf('Constraint values:\n');
for i = 1:2
    fprintf('  %s = %+.4e\n', cnames{i}, g_test(i));
end
fprintf('  ceq = %+.4e\n\n', ceq_test);

if any(g_test > 1e-4) || abs(ceq_test) > 1e-4
    fprintf('Point is infeasible - skipping KKT check.\n\n');
else
    [gf_t, gg_t, gc_t] = fd_gradients(x_test, params, h_fd);
    active_t = abs(g_test) < tol_active;
    idx_t    = find(active_t);

    if isempty(idx_t)
        A_t = gc_t;
    else
        A_t = [gg_t(:, idx_t), gc_t];
    end
    mu_lam_t   = pinv(A_t) * (-gf_t);
    lambda_t   = mu_lam_t(end);
    mu_t_full  = zeros(2,1);
    if ~isempty(idx_t)
        mu_t_full(idx_t) = mu_lam_t(1:end-1);
    end

    res_t = gf_t + gg_t * mu_t_full + gc_t * lambda_t;
    fprintf('Stationarity residual = [%+.4e,  %+.4e]\n', res_t(1), res_t(2));
    fprintf('Residual norm = %.4e\n', norm(res_t));
    df_t  = all(mu_t_full >= -1e-8);
    fprintf('Dual feasibility: %s\n', yn(df_t));
    fprintf('--> This is a KKT point: %s  (expected: NO)\n\n', yn(norm(res_t)<1e-4 && df_t));
end





% -------------------------------------------------------
% FIGURE: Gradient vectors at optimum on contour map
% -------------------------------------------------------
Nr = 80; Nh = 80;
r1_vec = linspace(r1_min, r1_max, Nr);
h_vec2 = linspace(h_min,  h_max,  Nh);
[R1g, Hg] = meshgrid(r1_vec, h_vec2);

V_required = params.m_soil / params.rho_soil;
Vmat_g = zeros(Nh, Nr);
ceq_g  = zeros(Nh, Nr);
g2_g   = zeros(Nh, Nr);

for i = 1:Nh
    for j = 1:Nr
        [Vm, Vp, ~, ~, ~, ~, ~] = potanalysis( ...
            R1g(i,j), Hg(i,j), params.theta_wall, params.t, params.rho_mat, ...
            params.rho_soil, params.K0, params.g_acc, params.sigma_allow, params.cost_mat, params.f_drain);
        Vmat_g(i,j) = Vm * 1e6;
        ceq_g(i,j)  = Vp - V_required;
        g2_g(i,j)   = Hg(i,j) - params.stab_ratio * R1g(i,j);
    end
end

figure(1); clf;
hold on;
Vmat_levels = linspace(min(Vmat_g(:)), max(Vmat_g(:)), 20);
contour(r1_vec*1000, h_vec2*1000, Vmat_g, Vmat_levels, 'LineWidth', 0.8);
colormap(parula); colorbar;
contour(r1_vec*1000, h_vec2*1000, ceq_g, [0 0], 'k-', 'LineWidth', 2.5);
contour(r1_vec*1000, h_vec2*1000, g2_g,  [0 0], 'b--', 'LineWidth', 1.8);

plot(x_opt(1)*1000, x_opt(2)*1000, 'r*', 'MarkerSize', 14, 'LineWidth', 2, ...
    'DisplayName', 'x* (optimum)');

h1 = plot(NaN,NaN,'k-', 'LineWidth',2.5);
h2 = plot(NaN,NaN,'b--','LineWidth',1.8);
h3 = plot(NaN,NaN,'r*', 'MarkerSize',14, 'LineWidth',2);
legend([h1 h2 h3], {'Volume constraint (equality)','g2: stability','x* (optimum)'}, ...
    'Location','northeast','FontSize',9);
xlabel('r_1  [mm]', 'FontSize', 11);
ylabel('h  [mm]',   'FontSize', 11);
title('KKT geometry: gradient vectors at x*', 'FontSize', 12);
axis([r1_min*1000 r1_max*1000 h_min*1000 h_max*1000]);
grid on; hold off;

fprintf('Done. See Figure 1 (KKT gradient geometry).\n');

% -------------------------------------------------------
% Local helper function
% -------------------------------------------------------
function s = yn(condition)
    if condition; s = 'YES'; else; s = 'NO'; end
end
