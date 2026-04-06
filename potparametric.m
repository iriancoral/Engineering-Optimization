% potparametric.m
% Parametric study: sweep soil masses for 2-var, 3-var, and 4-var models.
% Generates design relationship data and figures for the report.
%
% Produces:
%   Figure 1 - Optimal dimensions vs soil mass (log-log)
%   Figure 2 - Material volume and mass vs soil mass (log-log)
%   Figure 3 - Wall angle and thickness vs soil mass
%   Figure 4 - Constraint activity vs soil mass
%   Console  - Full results table
%   Saves    - potparametric_result.mat
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
params.f_drain     = f_drain;
params.stab_ratio  = stab_ratio;
params.t_min       = t_min;

% fmincon options
options = optimoptions('fmincon', ...
    'Algorithm',           'sqp', ...
    'Display',             'off', ...
    'OptimalityTolerance', 1e-10, ...
    'ConstraintTolerance', 1e-10, ...
    'StepTolerance',       1e-12, ...
    'MaxIterations',       500, ...
    'FiniteDifferenceType','forward');

% Extended model bounds (same as potfmincon_full.m)
r2_max = 0.40;
t_max  = 0.010;

% Soil masses to sweep
masses = [1, 5, 10, 25];
n_mass = length(masses);

% Storage
res2 = struct('r1', zeros(1,n_mass), 'h', zeros(1,n_mass), 'r2', zeros(1,n_mass), ...
              'Vmat', zeros(1,n_mass), 'g1', zeros(1,n_mass), 'g2', zeros(1,n_mass), ...
              'conv', false(1,n_mass));
res3 = struct('r1', zeros(1,n_mass), 'h', zeros(1,n_mass), 'r2', zeros(1,n_mass), ...
              'theta', zeros(1,n_mass), 'Vmat', zeros(1,n_mass), ...
              'g1', zeros(1,n_mass), 'g2', zeros(1,n_mass), 'g3', zeros(1,n_mass), ...
              'conv', false(1,n_mass));
res4 = struct('r1', zeros(1,n_mass), 'h', zeros(1,n_mass), 'r2', zeros(1,n_mass), ...
              't', zeros(1,n_mass), 'theta', zeros(1,n_mass), 'Vmat', zeros(1,n_mass), ...
              'g1', zeros(1,n_mass), 'g2', zeros(1,n_mass), 'g3', zeros(1,n_mass), ...
              'conv', false(1,n_mass), 't_on_lb', false(1,n_mass));

fprintf('=== Parametric Study: Soil Mass Sweep ===\n\n');

for mi = 1:n_mass
    ms = masses(mi);
    params.m_soil = ms;
    V_req = ms / rho_soil;

    % Characteristic scale for starting points
    r_est = max(r1_min + 0.01, min(r1_max - 0.01, (V_req/pi)^(1/3)));
    h_est = r_est;

    fprintf('--- m_soil = %g kg  (V_req = %.4f L) ---\n', ms, V_req*1000);

    % ===================================================
    % 2-VARIABLE MODEL
    % ===================================================
    obj2 = @(x) potobj(x, params);
    con2 = @(x) potcon(x, params);
    lb2 = [r1_min; h_min];
    ub2 = [r1_max; h_max];

    % Multiple starts scaled to soil mass
    r1f = [0.6 0.8 1.0 1.2 1.5 0.7 1.3 0.9];
    hf  = [0.8 1.2 1.0 1.5 0.7 1.4 0.9 1.1];
    best_f2 = Inf; best_x2 = [r_est; h_est];
    for si = 1:length(r1f)
        x0 = [max(lb2(1), min(ub2(1), r_est*r1f(si))); ...
              max(lb2(2), min(ub2(2), h_est*hf(si)))];
        try
            [xk, fk, ef] = fmincon(obj2, x0, [],[],[],[], lb2, ub2, con2, options);
            if ef > 0 && fk < best_f2
                best_f2 = fk; best_x2 = xk;
            end
        catch; end
    end
    r2_2 = best_x2(1) + best_x2(2)*tan(theta_wall*pi/180);
    [g2v, ~] = potcon(best_x2, params);
    res2.r1(mi) = best_x2(1); res2.h(mi) = best_x2(2); res2.r2(mi) = r2_2;
    res2.Vmat(mi) = best_f2; res2.g1(mi) = g2v(1); res2.g2(mi) = g2v(2);
    res2.conv(mi) = best_f2 < Inf;

    fprintf('  2-var: r1=%6.1fmm h=%6.1fmm r2=%6.1fmm Vmat=%8.2fcm3\n', ...
        best_x2(1)*1000, best_x2(2)*1000, r2_2*1000, best_f2*1e6);

    % ===================================================
    % 3-VARIABLE MODEL
    % ===================================================
    obj3 = @(x) potobj_full(x, params);
    con3 = @(x) potcon_full(x, params);
    lb3 = [r1_min; h_min; r1_min];
    ub3 = [r1_max; h_max; r2_max];

    alpha_vals = [1.2 1.4 1.3 1.5 1.2 1.6 1.3 1.4];
    best_f3 = Inf; best_x3 = [r_est; h_est; r_est*1.3];
    for si = 1:length(r1f)
        r1s = max(lb3(1), min(ub3(1), r_est*r1f(si)));
        hs  = max(lb3(2), min(ub3(2), h_est*hf(si)));
        r2s = max(lb3(3), min(ub3(3), r1s*alpha_vals(si)));
        try
            [xk, fk, ef] = fmincon(obj3, [r1s;hs;r2s], [],[],[],[], lb3, ub3, con3, options);
            if ef > 0 && fk < best_f3
                best_f3 = fk; best_x3 = xk;
            end
        catch; end
    end
    [g3v, ~] = potcon_full(best_x3, params);
    th3 = atan2d(best_x3(3)-best_x3(1), best_x3(2));
    res3.r1(mi) = best_x3(1); res3.h(mi) = best_x3(2); res3.r2(mi) = best_x3(3);
    res3.theta(mi) = th3; res3.Vmat(mi) = best_f3;
    res3.g1(mi) = g3v(1); res3.g2(mi) = g3v(2); res3.g3(mi) = g3v(3);
    res3.conv(mi) = best_f3 < Inf;

    fprintf('  3-var: r1=%6.1fmm h=%6.1fmm r2=%6.1fmm theta=%5.1fdeg Vmat=%8.2fcm3\n', ...
        best_x3(1)*1000, best_x3(2)*1000, best_x3(3)*1000, th3, best_f3*1e6);

    % ===================================================
    % 4-VARIABLE MODEL
    % ===================================================
    obj4 = @(x) potobj_full(x, params);
    con4 = @(x) potcon_full(x, params);
    lb4 = [r1_min; h_min; r1_min; t_min];
    ub4 = [r1_max; h_max; r2_max; t_max];

    t_frac = [1.0 1.5 1.2 1.0 2.0 1.3 1.5 1.0];
    best_f4 = Inf; best_x4 = [r_est; h_est; r_est*1.3; t_min];
    for si = 1:length(r1f)
        r1s = max(lb4(1), min(ub4(1), r_est*r1f(si)));
        hs  = max(lb4(2), min(ub4(2), h_est*hf(si)));
        r2s = max(lb4(3), min(ub4(3), r1s*alpha_vals(si)));
        ts  = max(lb4(4), min(ub4(4), t_min*t_frac(si)));
        try
            [xk, fk, ef] = fmincon(obj4, [r1s;hs;r2s;ts], [],[],[],[], lb4, ub4, con4, options);
            if ef > 0 && fk < best_f4
                best_f4 = fk; best_x4 = xk;
            end
        catch; end
    end
    [g4v, ~] = potcon_full(best_x4, params);
    th4 = atan2d(best_x4(3)-best_x4(1), best_x4(2));
    res4.r1(mi) = best_x4(1); res4.h(mi) = best_x4(2); res4.r2(mi) = best_x4(3);
    res4.t(mi) = best_x4(4); res4.theta(mi) = th4; res4.Vmat(mi) = best_f4;
    res4.g1(mi) = g4v(1); res4.g2(mi) = g4v(2); res4.g3(mi) = g4v(3);
    res4.conv(mi) = best_f4 < Inf;
    res4.t_on_lb(mi) = abs(best_x4(4) - t_min) < 1e-6;

    fprintf('  4-var: r1=%6.1fmm h=%6.1fmm r2=%6.1fmm t=%5.2fmm theta=%5.1fdeg Vmat=%8.2fcm3\n', ...
        best_x4(1)*1000, best_x4(2)*1000, best_x4(3)*1000, best_x4(4)*1000, th4, best_f4*1e6);

    % Material savings
    sav3 = (res2.Vmat(mi) - res3.Vmat(mi)) / res2.Vmat(mi) * 100;
    sav4 = (res2.Vmat(mi) - res4.Vmat(mi)) / res2.Vmat(mi) * 100;
    fprintf('  Savings: 3-var vs 2-var = %.2f%%,  4-var vs 2-var = %.2f%%\n\n', sav3, sav4);
end

% =======================================================
% SUMMARY TABLE
% =======================================================
fprintf('============================================\n');
fprintf('  PARAMETRIC STUDY SUMMARY\n');
fprintf('============================================\n\n');

fprintf('%-6s | %-8s %-8s %-8s %-8s %-10s | %-8s %-8s %-8s %-8s %-10s | %-8s %-8s %-8s %-8s %-8s %-10s\n', ...
    'm[kg]', 'r1_2', 'h_2', 'r2_2', 't_2', 'Vmat_2', ...
    'r1_3', 'h_3', 'r2_3', 'th_3', 'Vmat_3', ...
    'r1_4', 'h_4', 'r2_4', 't_4', 'th_4', 'Vmat_4');
for mi = 1:n_mass
    fprintf('%-6.0f | %7.1f  %7.1f  %7.1f  %7.1f  %9.2f | %7.1f  %7.1f  %7.1f  %7.1f  %9.2f | %7.1f  %7.1f  %7.1f  %7.2f  %7.1f  %9.2f\n', ...
        masses(mi), ...
        res2.r1(mi)*1000, res2.h(mi)*1000, res2.r2(mi)*1000, t*1000, res2.Vmat(mi)*1e6, ...
        res3.r1(mi)*1000, res3.h(mi)*1000, res3.r2(mi)*1000, res3.theta(mi), res3.Vmat(mi)*1e6, ...
        res4.r1(mi)*1000, res4.h(mi)*1000, res4.r2(mi)*1000, res4.t(mi)*1000, res4.theta(mi), res4.Vmat(mi)*1e6);
end

% =======================================================
% SCALING LAW FIT (log-log)
% =======================================================
fprintf('\n--- Scaling Law Fit (3-var model) ---\n');
log_m = log(masses);
log_r1 = log(res3.r1);
log_h  = log(res3.h);
log_Vm = log(res3.Vmat);

% Linear regression: log(y) = a*log(m) + b  =>  y = exp(b) * m^a
p_r1 = polyfit(log_m, log_r1, 1);
p_h  = polyfit(log_m, log_h, 1);
p_Vm = polyfit(log_m, log_Vm, 1);

fprintf('  r1  ~ m^%.3f  (theory: 1/3 = 0.333)\n', p_r1(1));
fprintf('  h   ~ m^%.3f  (theory: 1/3 = 0.333)\n', p_h(1));
fprintf('  Vmat ~ m^%.3f  (theory: 2/3 = 0.667)\n', p_Vm(1));

% Check wall angle constancy
fprintf('\n--- Wall Angle Across Masses (3-var) ---\n');
for mi = 1:n_mass
    fprintf('  m = %2d kg:  theta = %.2f deg\n', masses(mi), res3.theta(mi));
end
fprintf('  -> Self-similar scaling: angle stays constant if exponents match 1/3.\n');

% =======================================================
% FIGURE 1: Optimal dimensions vs soil mass (log-log)
% =======================================================
figure(1); clf;
ms_plot = masses;

subplot(1,3,1);
loglog(ms_plot, res2.r1*1000, 'bo-', 'LineWidth', 1.5, 'DisplayName', '2-var'); hold on;
loglog(ms_plot, res3.r1*1000, 'rs-', 'LineWidth', 1.5, 'DisplayName', '3-var');
loglog(ms_plot, res4.r1*1000, 'g^-', 'LineWidth', 1.5, 'DisplayName', '4-var');
% Reference line m^(1/3)
m_ref = linspace(masses(1), masses(end), 50);
loglog(m_ref, res3.r1(2)*1000 * (m_ref/masses(2)).^(1/3), 'k:', 'LineWidth', 1.0, 'DisplayName', 'm^{1/3} ref');
ylabel('r_1^* [mm]'); xlabel('m_{soil} [kg]');
title('Bottom radius'); legend('Location','northwest','FontSize',8); grid on;

subplot(1,3,2);
loglog(ms_plot, res2.h*1000, 'bo-', 'LineWidth', 1.5); hold on;
loglog(ms_plot, res3.h*1000, 'rs-', 'LineWidth', 1.5);
loglog(ms_plot, res4.h*1000, 'g^-', 'LineWidth', 1.5);
loglog(m_ref, res3.h(2)*1000 * (m_ref/masses(2)).^(1/3), 'k:', 'LineWidth', 1.0);
ylabel('h^* [mm]'); xlabel('m_{soil} [kg]');
title('Height'); grid on;

subplot(1,3,3);
loglog(ms_plot, res2.r2*1000, 'bo-', 'LineWidth', 1.5); hold on;
loglog(ms_plot, res3.r2*1000, 'rs-', 'LineWidth', 1.5);
loglog(ms_plot, res4.r2*1000, 'g^-', 'LineWidth', 1.5);
loglog(m_ref, res3.r2(2)*1000 * (m_ref/masses(2)).^(1/3), 'k:', 'LineWidth', 1.0);
ylabel('r_2^* [mm]'); xlabel('m_{soil} [kg]');
title('Top radius'); grid on;

sgtitle('Optimal dimensions vs soil mass (log-log)', 'FontSize', 12);
set(gcf, 'Position', [100 100 1200 400]);
saveas(gcf, 'Pictures/potparametric fig(1).png');

% =======================================================
% FIGURE 2: Material volume and mass vs soil mass
% =======================================================
figure(2); clf;

subplot(1,2,1);
loglog(ms_plot, res2.Vmat*1e6, 'bo-', 'LineWidth', 1.5, 'DisplayName', '2-var'); hold on;
loglog(ms_plot, res3.Vmat*1e6, 'rs-', 'LineWidth', 1.5, 'DisplayName', '3-var');
loglog(ms_plot, res4.Vmat*1e6, 'g^-', 'LineWidth', 1.5, 'DisplayName', '4-var');
loglog(m_ref, res3.Vmat(2)*1e6 * (m_ref/masses(2)).^(2/3), 'k:', 'LineWidth', 1.0, 'DisplayName', 'm^{2/3} ref');
ylabel('V_{mat}^* [cm^3]'); xlabel('m_{soil} [kg]');
title('Material volume'); legend('Location','northwest','FontSize',8); grid on;

subplot(1,2,2);
loglog(ms_plot, res2.Vmat*rho_mat*1000, 'bo-', 'LineWidth', 1.5, 'DisplayName', '2-var'); hold on;
loglog(ms_plot, res3.Vmat*rho_mat*1000, 'rs-', 'LineWidth', 1.5, 'DisplayName', '3-var');
loglog(ms_plot, res4.Vmat*rho_mat*1000, 'g^-', 'LineWidth', 1.5, 'DisplayName', '4-var');
ylabel('Pot mass [g]'); xlabel('m_{soil} [kg]');
title('Pot mass'); legend('Location','northwest','FontSize',8); grid on;

sgtitle('Material volume and pot mass vs soil mass', 'FontSize', 12);
set(gcf, 'Position', [100 100 1000 400]);
saveas(gcf, 'Pictures/potparametric fig(2).png');

% =======================================================
% FIGURE 3: Wall angle and thickness vs soil mass
% =======================================================
figure(3); clf;

subplot(1,2,1);
plot(masses, res3.theta, 'rs-', 'LineWidth', 1.5, 'DisplayName', '3-var'); hold on;
plot(masses, res4.theta, 'g^-', 'LineWidth', 1.5, 'DisplayName', '4-var');
yline(theta_wall, 'b--', 'LineWidth', 1.2, 'DisplayName', '2-var (fixed 10 deg)');
ylabel('\theta_{wall}^* [deg]'); xlabel('m_{soil} [kg]');
title('Optimal wall angle'); legend('Location','best','FontSize',8); grid on;

subplot(1,2,2);
plot(masses, res4.t*1000, 'g^-', 'LineWidth', 1.5, 'DisplayName', '4-var'); hold on;
yline(t*1000, 'b--', 'LineWidth', 1.2, 'DisplayName', '2/3-var (fixed 3 mm)');
ylabel('t^* [mm]'); xlabel('m_{soil} [kg]');
title('Optimal wall thickness'); legend('Location','best','FontSize',8); grid on;

sgtitle('Wall angle and thickness vs soil mass', 'FontSize', 12);
set(gcf, 'Position', [100 100 1000 400]);
saveas(gcf, 'Pictures/potparametric fig(3).png');

% =======================================================
% FIGURE 4: Material savings (3-var and 4-var vs 2-var)
% =======================================================
figure(4); clf;

sav3_pct = (res2.Vmat - res3.Vmat) ./ res2.Vmat * 100;
sav4_pct = (res2.Vmat - res4.Vmat) ./ res2.Vmat * 100;

bar_data = [sav3_pct; sav4_pct]';
b = bar(bar_data, 'grouped');
b(1).FaceColor = [0.8 0.3 0.2];
b(2).FaceColor = [0.2 0.6 0.4];
set(gca, 'XTickLabel', arrayfun(@(m) sprintf('%d kg', m), masses, 'UniformOutput', false));
ylabel('Material savings vs 2-var [%]');
legend('3-var (free r_2)', '4-var (free r_2 + t)', 'Location', 'best', 'FontSize', 9);
title('Material savings from extended models', 'FontSize', 12);
grid on;
set(gcf, 'Position', [100 100 600 400]);
saveas(gcf, 'Pictures/potparametric fig(4).png');

% =======================================================
% SAVE
% =======================================================
save('potparametric_result.mat', 'masses', 'res2', 'res3', 'res4', ...
     'p_r1', 'p_h', 'p_Vm', 'params');
fprintf('\nResults saved to potparametric_result.mat\n');
fprintf('Done. See Figures 1-4.\n');
