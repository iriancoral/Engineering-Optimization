function f = potobj_full(x, params)
% Objective function for extended flower pot optimization.
% Supports 3-variable (r1, h, r2) and 4-variable (r1, h, r2, t) models.
%
% INPUT:
%   x      : [r1; h; r2] or [r1; h; r2; t]   [m]
%   params : struct with fixed parameters
%
% OUTPUT:
%   f      : material volume Vmat              [m³]

r1 = x(1);
h  = x(2);
r2 = x(3);

if length(x) >= 4
    t = x(4);
else
    t = params.t;   % fixed wall thickness (3-variable model)
end

[Vmat, ~, ~, ~, ~, ~, ~] = potanalysis_full( ...
    r1, h, r2, t, params.rho_mat, params.rho_soil, ...
    params.K0, params.g_acc, params.sigma_allow, params.cost_mat, params.f_drain);

f = Vmat;
end
