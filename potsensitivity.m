% potsensitivity.m
% Step 6: Sensitivity analysis for the 2-variable flower pot model.
%
% Computes:
%   1. Forward FD gradients of f and all g at a chosen base point x0
%   2. Step size study: gradient error vs h (like Exercise 5.1)
%   3. Logarithmic (relative) sensitivities of f w.r.t. x1=r1, x2=h
%   4. Hessian approximation (central differences) -> convexity check
%
% Run potparams first (called internally here).
% -------------------------------------------------------

clear; clc;
potparams;

% Pack fixed parameters into struct for clean function calls
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

% -------------------------------------------------------
% Base point x0 for sensitivity evaluation
% (interior point, away from constraints)
% -------------------------------------------------------
r1_0 = 0.10;    % 100 mm
h_0  = 0.15;    % 150 mm
x0   = [r1_0; h_0];

fprintf('=== Sensitivity Analysis ===\n');
fprintf('Base point: r1 = %.1f mm,  h = %.1f mm\n\n', x0(1)*1000, x0(2)*1000);

% -------------------------------------------------------
% Helper: evaluate f and g at point x
% -------------------------------------------------------
    function [f, g, ceq] = evaluate(x, p)
        [f, g_vals, ceq_vals] = deal(potobj(x, p), [], []);
        [g, ceq] = potcon(x, p);
    end

% Function handles for cleaner code
f_handle   = @(x) potobj(x, params);
con_handle = @(x) potcon(x, params);

% Base function values
f0         = f_handle(x0);
[g0, ceq0] = con_handle(x0);

fprintf('Function values at x0:\n');
fprintf('  f  (Vmat)     = %.6e m3  (%.4f cm3)\n', f0, f0*1e6);
fprintf('  g1 (stress)   = %.6e  (<= 0 means feasible)\n', g0(1));
fprintf('  g2 (drainage) = %.6e\n', g0(2));
fprintf('  g3 (stability)= %.6e\n', g0(3));
fprintf('  ceq (volume)  = %.6e\n\n', ceq0);

% -------------------------------------------------------
% 1. STEP SIZE STUDY: gradient error vs step size h
% -------------------------------------------------------
fprintf('--- Step size study ---\n');

% Reference gradient using very small complex step (near-exact)
% Complex step: df/dx ~ Im(f(x + i*h)) / h  (no cancellation error)
h_cs  = 1e-20;
grad_ref = zeros(2,1);
for k = 1:2
    x_cs      = x0;
    x_cs(k)   = x0(k) + 1i * h_cs;
    % Evaluate potanalysis with complex input
    [Vmat_cs, ~, ~, ~, ~, ~, ~] = potanalysis( ...
        x_cs(1), x_cs(2), theta_wall, t, rho_mat, rho_soil, K0, g_acc, sigma_allow, cost_mat);
    grad_ref(k) = imag(Vmat_cs) / h_cs;
end
fprintf('Reference gradient (complex step, h=1e-20):\n');
fprintf('  df/dr1 = %.10e\n', grad_ref(1));
fprintf('  df/dh  = %.10e\n\n', grad_ref(2));

% Step sizes to test
h_vec  = 10.^(-1:-1:-14);
n_h    = length(h_vec);

err_f_r1  = zeros(1, n_h);   % error in df/dr1 (forward FD)
err_f_h   = zeros(1, n_h);   % error in df/dh  (forward FD)
err_g1_r1 = zeros(1, n_h);   % error in dg1/dr1

% Compute reference for g1
g1_handle = @(x) potcon(x, params);
grad_g1_ref = zeros(2,1);
for k = 1:2
    x_cs    = x0;
    x_cs(k) = x0(k) + 1i * h_cs;
    [g_cs, ~] = potcon(x_cs, params);
    grad_g1_ref(k) = imag(g_cs(1)) / h_cs;
end

for i = 1:n_h
    hs = h_vec(i);

    % Forward FD for f
    x_p1      = x0; x_p1(1) = x0(1) + hs;
    x_p2      = x0; x_p2(2) = x0(2) + hs;
    grad_f_fd = [(f_handle(x_p1) - f0)/hs;
                 (f_handle(x_p2) - f0)/hs];

    err_f_r1(i) = abs(grad_f_fd(1) - grad_ref(1)) / abs(grad_ref(1));
    err_f_h(i)  = abs(grad_f_fd(2) - grad_ref(2)) / abs(grad_ref(2));

    % Forward FD for g1
    [g_p1, ~] = con_handle(x_p1);
    [g_p2, ~] = con_handle(x_p2);
    grad_g1_fd_r1 = (g_p1(1) - g0(1)) / hs;
    err_g1_r1(i)  = abs(grad_g1_fd_r1 - grad_g1_ref(1)) / abs(grad_g1_ref(1) + eps);
end

figure(1); clf;
loglog(h_vec, err_f_r1,  'b-o', 'LineWidth', 1.5, 'DisplayName', 'df/dr1 (objective)');
hold on;
loglog(h_vec, err_f_h,   'r-s', 'LineWidth', 1.5, 'DisplayName', 'df/dh  (objective)');
loglog(h_vec, err_g1_r1, 'm-^', 'LineWidth', 1.5, 'DisplayName', 'dg1/dr1 (stress constraint)');
% Reference lines for O(h) and O(h^2)
loglog(h_vec, h_vec/h_vec(4),      'k:',  'DisplayName', 'O(h) reference');
loglog(h_vec, h_vec.^2/h_vec(4)^2,'k--', 'DisplayName', 'O(h^2) reference');
xline(1e-8, 'g--', 'LineWidth', 1.0);
text(1.2e-8, 1e-1, 'h = 1e-8', 'FontSize', 8, 'Color', [0 0.6 0]);
xlabel('Step size  h', 'FontSize', 11);
ylabel('Relative gradient error  |FD - ref| / |ref|', 'FontSize', 11);
title('Forward FD step size study', 'FontSize', 12);
legend('Location', 'southwest', 'FontSize', 9);
grid on;
hold off;

% Choose optimal step size
h_opt = 1e-8;
fprintf('Chosen step size for gradient computation: h = %.0e\n\n', h_opt);

% -------------------------------------------------------
% 2. FORWARD FD GRADIENTS at x0 with h_opt
% -------------------------------------------------------
fprintf('--- FD Gradients at x0 (h = %.0e) ---\n', h_opt);

grad_f  = zeros(2,1);
grad_g1 = zeros(2,1);
grad_g2 = zeros(2,1);
grad_g3 = zeros(2,1);
grad_ceq= zeros(2,1);

for k = 1:2
    x_p = x0;  x_p(k) = x0(k) + h_opt;
    f_p          = f_handle(x_p);
    [g_p, ceq_p] = con_handle(x_p);

    grad_f(k)   = (f_p       - f0)    / h_opt;
    grad_g1(k)  = (g_p(1)    - g0(1)) / h_opt;
    grad_g2(k)  = (g_p(2)    - g0(2)) / h_opt;
    grad_g3(k)  = (g_p(3)    - g0(3)) / h_opt;
    grad_ceq(k) = (ceq_p     - ceq0)  / h_opt;
end

fprintf('  grad f   = [%+.6e,  %+.6e]\n', grad_f(1),   grad_f(2));
fprintf('  grad g1  = [%+.6e,  %+.6e]\n', grad_g1(1),  grad_g1(2));
fprintf('  grad g2  = [%+.6e,  %+.6e]\n', grad_g2(1),  grad_g2(2));
fprintf('  grad g3  = [%+.6e,  %+.6e]\n', grad_g3(1),  grad_g3(2));
fprintf('  grad ceq = [%+.6e,  %+.6e]\n', grad_ceq(1), grad_ceq(2));
fprintf('\n');

% -------------------------------------------------------
% 3. LOGARITHMIC (RELATIVE) SENSITIVITIES
% dLf/dLx = (x/f) * (df/dx)
% -------------------------------------------------------
fprintf('--- Logarithmic sensitivities at x0 ---\n');
logsen_r1 = (x0(1) / f0) * grad_f(1);
logsen_h  = (x0(2) / f0) * grad_f(2);

fprintf('  dLf/dL(r1) = %.4f  (|value|>1 means influential)\n', logsen_r1);
fprintf('  dLf/dL(h)  = %.4f\n\n', logsen_h);

if abs(logsen_r1) > abs(logsen_h)
    fprintf('  -> r1 is the more influential design variable on Vmat\n\n');
else
    fprintf('  -> h is the more influential design variable on Vmat\n\n');
end

% Sensitivity of ceq w.r.t. x (shows how volume changes with design vars)
logsen_ceq_r1 = (x0(1) / (abs(ceq0)+eps)) * grad_ceq(1);
logsen_ceq_h  = (x0(2) / (abs(ceq0)+eps)) * grad_ceq(2);
fprintf('  dL(ceq)/dL(r1) = %.4f\n', logsen_ceq_r1);
fprintf('  dL(ceq)/dL(h)  = %.4f\n\n', logsen_ceq_h);

% Bar chart of log sensitivities
figure(2); clf;
bar_vals = [logsen_r1, logsen_h];
bar(bar_vals, 0.4, 'FaceColor', [0.2 0.5 0.8]);
set(gca, 'XTickLabel', {'r_1', 'h'}, 'FontSize', 11);
ylabel('dLf / dLx  (logarithmic sensitivity)', 'FontSize', 11);
title('Logarithmic sensitivities of V_{mat} at base point x_0', 'FontSize', 12);
yline(0, 'k-');
yline( 1, 'r:', 'LineWidth', 1.2);
yline(-1, 'r:', 'LineWidth', 1.2);
text(2.3, 1.05, '|value|=1 threshold', 'FontSize', 8, 'Color', 'r');
grid on;

% -------------------------------------------------------
% 4. HESSIAN (central differences) -> convexity check
% -------------------------------------------------------
fprintf('--- Hessian at x0 (central FD, h = %.0e) ---\n', h_opt);

H_mat = zeros(2,2);
for i = 1:2
    for j = 1:2
        xpp = x0; xpp(i) = xpp(i) + h_opt; xpp(j) = xpp(j) + h_opt;
        xpm = x0; xpm(i) = xpm(i) + h_opt; xpm(j) = xpm(j) - h_opt;
        xmp = x0; xmp(i) = xmp(i) - h_opt; xmp(j) = xmp(j) + h_opt;
        xmm = x0; xmm(i) = xmm(i) - h_opt; xmm(j) = xmm(j) - h_opt;
        H_mat(i,j) = (f_handle(xpp) - f_handle(xpm) - f_handle(xmp) + f_handle(xmm)) ...
                     / (4 * h_opt^2);
    end
end

eig_vals = eig(H_mat);
fprintf('  H = [%.4e  %.4e]\n', H_mat(1,1), H_mat(1,2));
fprintf('      [%.4e  %.4e]\n', H_mat(2,1), H_mat(2,2));
fprintf('  Eigenvalues: lambda1 = %.4e,  lambda2 = %.4e\n', eig_vals(1), eig_vals(2));

if all(eig_vals > 0)
    fprintf('  -> H is positive definite at x0: f is locally CONVEX here.\n');
elseif all(eig_vals >= 0)
    fprintf('  -> H is positive semi-definite at x0: f is locally convex (not strictly).\n');
else
    fprintf('  -> H has a negative eigenvalue at x0: f is NOT convex here.\n');
end

fprintf('\nDone. See Figure 1 (step size study) and Figure 2 (log sensitivities).\n');
