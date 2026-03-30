% potcontour.m
% Step 3: Contour plots for initial problem investigation.
% Simplified 2-variable model: x = [r1, h], r2 = r1 + h*tan(theta_wall), t = t_min fixed.
%
% Produces two figures:
%   Figure 1 - Vmat contours + all constraint boundaries over (r1, h) grid
%   Figure 2 - Vmat along the equality constraint curve (1D view of problem)
%
% Run potparams first, or run this script directly (it calls potparams internally).
% -------------------------------------------------------

clear; clc;
potparams;   % load all fixed parameters into workspace

% Required internal volume
V_required = m_soil / rho_soil;

fprintf('=== Flower Pot Optimization - Contour Analysis ===\n');
fprintf('Soil mass      : %.1f kg\n', m_soil);
fprintf('Required volume: %.4f L (%.6f m3)\n', V_required*1000, V_required);
fprintf('Wall angle     : %.1f deg\n', theta_wall);
fprintf('Wall thickness : %.1f mm\n', t*1000);
fprintf('\n');

% -------------------------------------------------------
% Build evaluation grid
% -------------------------------------------------------
Nr = 120;
Nh = 120;

r1_vec = linspace(r1_min, r1_max, Nr);   % [m]
h_vec  = linspace(h_min,  h_max,  Nh);   % [m]

[R1, H] = meshgrid(r1_vec, h_vec);

Vmat_grid = zeros(Nh, Nr);
ceq_grid  = zeros(Nh, Nr);
g1_grid   = zeros(Nh, Nr);   % wall stress
g2_grid   = zeros(Nh, Nr);   % stability

for i = 1:Nh
    for j = 1:Nr
        r1_ij = R1(i,j);
        h_ij  = H(i,j);

        [Vmat, Vpot, ~, sigma_max, ~, ~, ~] = potanalysis( ...
            r1_ij, h_ij, theta_wall, t, rho_mat, rho_soil, K0, g_acc, sigma_allow, cost_mat, f_drain);

        Vmat_grid(i,j) = Vmat;
        ceq_grid(i,j)  = Vpot - V_required;
        g1_grid(i,j)   = sigma_max - sigma_allow;
        g2_grid(i,j)   = h_ij - stab_ratio * r1_ij;
    end
end

% Unit conversions for plotting
r1_mm    = r1_vec * 1000;           % m  -> mm
h_mm     = h_vec  * 1000;           % m  -> mm
Vmat_cm3 = Vmat_grid * 1e6;         % m3 -> cm3

% -------------------------------------------------------
% Infeasibility mask (inequality constraints only)
% -------------------------------------------------------
infeas = (g1_grid > 0) | (g2_grid > 0);

% -------------------------------------------------------
% Compute Vmat along equality constraint (needed for optimum point)
% For each r1, solve for h from volume equation, then evaluate.
% -------------------------------------------------------

N_line = 400;
r1_line = linspace(r1_min, r1_max, N_line);

Vmat_line  = NaN(1, N_line);
h_line_sol = NaN(1, N_line);
g1_line    = NaN(1, N_line);
g2_line    = NaN(1, N_line);
feas_line  = false(1, N_line);

for j = 1:N_line
    r1j = r1_line(j);
    th  = theta_wall * pi / 180;

    % Vpot is cubic in h: (pi*tan^2/3)*h^3 + (pi*tan*r1)*h^2 + (pi*r1^2)*h - V_req = 0
    a   = tan(th);
    p3  = pi * a^2 / 3;
    p2  = pi * a * r1j;
    p1  = pi * r1j^2;
    p0  = -V_required;
    h_roots = roots([p3, p2, p1, p0]);
    h_real = h_roots(imag(h_roots)==0 & real(h_roots)>0);
    if isempty(h_real); continue; end
    h_sol = real(h_real(1));

    if h_sol < h_min || h_sol > h_max
        continue;
    end

    h_line_sol(j) = h_sol;

    [Vmat, ~, ~, sigma_max, ~, ~, ~] = potanalysis( ...
        r1j, h_sol, theta_wall, t, rho_mat, rho_soil, K0, g_acc, sigma_allow, cost_mat, f_drain);

    Vmat_line(j) = Vmat * 1e6;   % cm3

    g1_line(j) = sigma_max - sigma_allow;
    g2_line(j) = h_sol - stab_ratio * r1j;

    feas_line(j) = (g1_line(j) <= 0) && (g2_line(j) <= 0) && ...
                   (r1j >= r1_min) && (r1j <= r1_max);
end

% Find graphical optimum (min Vmat in feasible segment)
Vmat_feas = Vmat_line;
Vmat_feas(~feas_line) = NaN;
[Vmat_min_est, idx_min] = min(Vmat_feas);

if ~isnan(Vmat_min_est)
    r1_opt_est = r1_line(idx_min) * 1000;   % mm
    h_opt_est  = h_line_sol(idx_min) * 1000; % mm

    fprintf('=== Graphical Optimum Estimate ===\n');
    fprintf('  r1   = %.1f mm\n',   r1_opt_est);
    fprintf('  h    = %.1f mm\n',   h_opt_est);
    fprintf('  r2   = %.1f mm\n',   r1_opt_est + h_opt_est*tan(theta_wall*pi/180));
    fprintf('  Vmat = %.4f cm3\n',  Vmat_min_est);
    fprintf('  g1 (stress)    = %.4e  (<=0 ok)\n', g1_line(idx_min));
    fprintf('  g2 (stability) = %.4e  (<=0 ok)\n', g2_line(idx_min));
    fprintf('\n');
end

% -------------------------------------------------------
% FIGURE 1: Full contour map
% -------------------------------------------------------
figure(1); clf;
hold on;

% --- Objective contours (background) ---
% Choose contour levels to be informative around expected optimum
Vmat_levels = linspace(min(Vmat_cm3(:)), prctile(Vmat_cm3(:), 80), 18);
[C_obj, h_obj] = contour(r1_mm, h_mm, Vmat_cm3, Vmat_levels, ...
    'LineWidth', 1.0, 'EdgeColor', 'flat');
colormap(parula);
cb = colorbar;
ylabel(cb, 'V_{mat} [cm^3]');
clabel(C_obj, h_obj, 'FontSize', 7, 'LabelSpacing', 200);

% --- Shade infeasible region (gray) ---
infeas_double = double(infeas);
contourf(r1_mm, h_mm, infeas_double, [0.5 0.5], ...
    'FaceColor', [0.65 0.65 0.65], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
% Re-draw objective contours on top of shading
contour(r1_mm, h_mm, Vmat_cm3, Vmat_levels, 'LineWidth', 1.0);

% --- Constraint boundaries (g = 0 lines) ---
contour(r1_mm, h_mm, g1_grid, [0 0], 'r--', 'LineWidth', 2.0);   % stress
contour(r1_mm, h_mm, g2_grid, [0 0], 'b--', 'LineWidth', 2.0);   % stability

% --- Equality constraint (volume curve) - thick black ---
contour(r1_mm, h_mm, ceq_grid, [0 0], 'k-', 'LineWidth', 3.0);

% --- Bound lines ---
xline(r1_min*1000, 'g:', 'LineWidth', 1.5);
xline(r1_max*1000, 'g:', 'LineWidth', 1.5);
yline(h_min*1000,  'c:', 'LineWidth', 1.5);
yline(h_max*1000,  'c:', 'LineWidth', 1.5);

% --- Optimum point ---
if ~isnan(Vmat_min_est)
    plot(r1_opt_est, h_opt_est, 'r*', 'MarkerSize', 14, 'LineWidth', 2.0);
end

% --- Manual legend entries ---
h_leg(1) = plot(NaN, NaN, 'k-',  'LineWidth', 3.0, 'DisplayName', 'h_1: volume = V_{req} (equality)');
h_leg(2) = plot(NaN, NaN, 'r--', 'LineWidth', 2.0, 'DisplayName', 'g_1: wall stress limit');
h_leg(3) = plot(NaN, NaN, 'b--', 'LineWidth', 2.0, 'DisplayName', 'g_2: stability limit');
h_leg(4) = plot(NaN, NaN, 'r*',  'MarkerSize', 14, 'LineWidth', 2.0, 'DisplayName', ...
    sprintf('Graphical optimum (%.1f, %.1f) mm', r1_opt_est, h_opt_est));
h_leg(5) = patch(NaN, NaN, [0.65 0.65 0.65], 'FaceAlpha', 0.5, 'DisplayName', 'Infeasible region');
legend(h_leg, 'Location', 'northeast', 'FontSize', 9);

xlabel('r_1  [mm]', 'FontSize', 11);
ylabel('h  [mm]',   'FontSize', 11);
title(sprintf('V_{mat} [cm^3] contours  |  m_{soil} = %g kg, \\theta = %.0f deg, t = %.0f mm', ...
    m_soil, theta_wall, t*1000), 'FontSize', 12);
axis([r1_min*1000 r1_max*1000 h_min*1000 h_max*1000]);
grid on;
hold off;

% -------------------------------------------------------
% FIGURE 2: Vmat along equality constraint (1D view)
% (data already computed above)
% -------------------------------------------------------

figure(2); clf;
hold on;

% Full curve (including infeasible parts)
plot(r1_line*1000, Vmat_line, 'k--', 'LineWidth', 1.5, 'DisplayName', 'V_{mat} (all r_1, h from eq. constraint)');

% Feasible segment (thick blue)
plot(r1_line(feas_line)*1000, Vmat_feas(feas_line), 'b-', 'LineWidth', 3.0, 'DisplayName', 'Feasible segment');

% Estimated optimum
if ~isnan(Vmat_min_est)
    plot(r1_opt_est, Vmat_min_est, 'r*', 'MarkerSize', 14, 'LineWidth', 2.0, ...
        'DisplayName', sprintf('Graphical min: V_{mat} = %.2f cm^3', Vmat_min_est));
end

% Constraint boundary markers on x-axis
% Mark where each g_i = 0 crosses the curve
for j = 2:N_line
    if ~isnan(g2_line(j)) && ~isnan(g2_line(j-1))
        if sign(g2_line(j)) ~= sign(g2_line(j-1))
            xline(r1_line(j)*1000, 'b:', 'LineWidth', 1.5);
        end
    end
end

xlabel('r_1  [mm]', 'FontSize', 11);
ylabel('V_{mat}  [cm^3]', 'FontSize', 11);
title(sprintf('V_{mat} along volume equality constraint  |  m_{soil} = %g kg', m_soil), 'FontSize', 12);
legend('Location', 'best', 'FontSize', 9);
grid on;
hold off;

fprintf('\nDone. See Figure 1 (contour map) and Figure 2 (1D view along volume curve).\n');
