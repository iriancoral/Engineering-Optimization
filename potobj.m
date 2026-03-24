function f = potobj(x, params)
% Objective function wrapper for flower pot optimization.
% Returns material volume Vmat to be minimized.
%
% INPUT:
%   x      : design variable vector [r1; h]   [m]
%   params : struct with all fixed parameters (from potparams.m)
%
% OUTPUT:
%   f      : material volume Vmat             [m³]

r1 = x(1);
h  = x(2);

[Vmat, ~, ~, ~, ~, ~, ~] = potanalysis( ...
    r1, h, ...
    params.theta_wall, params.t, params.rho_mat, params.rho_soil, ...
    params.K0, params.g_acc, params.sigma_allow, params.cost_mat);

f = Vmat;

% end
