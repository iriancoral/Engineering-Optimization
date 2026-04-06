function [g, ceq] = potcon_full(x, params)
% Constraint function for extended flower pot optimization (fmincon format).
% Supports 3-variable (r1, h, r2) and 4-variable (r1, h, r2, t) models.
%
% INPUT:
%   x      : [r1; h; r2] or [r1; h; r2; t]   [m]
%   params : struct with fixed parameters
%
% OUTPUT:
%   g      : inequality constraints g <= 0
%   ceq    : equality constraint    ceq = 0
%
% CONSTRAINTS:
%   g(1): sigma_max <= sigma_allow        (wall stress)
%   g(2): h - stab_ratio*r1 <= 0          (stability)
%   g(3): r1 - r2 <= 0                    (pot must taper outward: r2 >= r1)
%   ceq : Vpot - V_required = 0           (soil volume)

r1 = x(1);
h  = x(2);
r2 = x(3);

if length(x) >= 4
    t = x(4);
else
    t = params.t;
end

[~, Vpot, ~, sigma_max, ~, ~, ~] = potanalysis_full( ...
    r1, h, r2, t, params.rho_mat, params.rho_soil, ...
    params.K0, params.g_acc, params.sigma_allow, params.cost_mat, params.f_drain);

V_required = params.m_soil / params.rho_soil;

% Inequality constraints: g <= 0
g(1) = sigma_max - params.sigma_allow;          % wall stress
g(2) = h - params.stab_ratio * r1;              % stability
g(3) = r1 - r2;                                 % r2 >= r1 (outward taper)

% Equality constraint: ceq = 0
ceq = Vpot - V_required;
end
