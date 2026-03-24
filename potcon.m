function [g, ceq] = potcon(x, params)
% Constraint function wrapper for flower pot optimization (fmincon format).
%
% INPUT:
%   x      : design variable vector [r1; h]   [m]
%   params : struct with all fixed parameters (from potparams.m)
%
% OUTPUT:
%   g      : inequality constraints g <= 0
%   ceq    : equality constraints   ceq = 0
%
% CONSTRAINTS:
%   g(1): sigma_max <= sigma_allow   (wall stress)
%   g(2): h - stab_ratio*r1 <= 0    (stability: pot not too tall/narrow)
%   ceq : Vpot - V_required = 0     (soil volume must fit)
% Note: drainage handled directly in potanalysis (15% of base area removed)

r1 = x(1);
h  = x(2);

[~, Vpot, ~, sigma_max, ~, ~, ~] = potanalysis( ...
    r1, h, ...
    params.theta_wall, params.t, params.rho_mat, params.rho_soil, ...
    params.K0, params.g_acc, params.sigma_allow, params.cost_mat, params.f_drain);

% Required internal volume
V_required = params.m_soil / params.rho_soil;

% --- Inequality constraints: g <= 0 -------------------

% g1: wall stress must not exceed allowable stress
g(1) = sigma_max - params.sigma_allow;

% g2: stability - pot must not be too tall relative to base
g(2) = h - params.stab_ratio * r1;

% --- Equality constraint: ceq = 0 ---------------------

% Pot internal volume must match required soil volume
ceq = Vpot - V_required;

% end
