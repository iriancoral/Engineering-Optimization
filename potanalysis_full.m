function [Vmat, Vpot, s, sigma_max, A_base, Vmass, Vcost] = ...
    potanalysis_full(r1, h, r2, t, rho_mat, rho_soil, K0, g_acc, sigma_allow, cost_mat, f_drain)
% Core analysis function for the EXTENDED flower pot optimization.
% 3-4 variable model: r2 and t are direct inputs (not derived from theta_wall).
%
% INPUT:
%   r1, h, r2, t  : design variables [m]
%   (remaining)    : fixed parameters from potparams.m
%
% OUTPUT:
%   Vmat      : material volume (wall + base)   [m³]  <- OBJECTIVE
%   Vpot      : internal volume (truncated cone) [m³]
%   s         : slant height                     [m]
%   sigma_max : max hoop stress at bottom        [Pa]
%   A_base    : internal base area               [m²]
%   Vmass     : pot mass                         [kg]
%   Vcost     : material cost                    [euro]

% Slant height of lateral surface
s = sqrt(h^2 + (r2 - r1)^2);

% Material volume: lateral shell + base disc (minus drainage holes)
Vmat = t * (pi*(r1 + r2)*s + pi*r1^2*(1 - f_drain));

% Internal volume of truncated cone
Vpot = (pi*h/3) * (r1^2 + r1*r2 + r2^2);

% Max hoop stress along wall height.
% sigma(z) = K0*rho_soil*g*z * r(z) / t, where r(z) = r2 - z*(r2-r1)/h
% This is quadratic in z with maximum at z* = r2*h / (2*(r2-r1)).
% If z* < h (i.e. r2 > 2*r1), the max is at z*; otherwise at z = h (bottom).
if r2 > 2*r1
    z_crit    = r2 * h / (2*(r2 - r1));
    r_crit    = r2 - z_crit * (r2 - r1) / h;   % = r2/2
    sigma_max = K0 * rho_soil * g_acc * z_crit * r_crit / t;
else
    % Maximum at bottom (z = h, r = r1)
    sigma_max = K0 * rho_soil * g_acc * h * r1 / t;
end

% Base area
A_base = pi * r1^2;

% Mass and cost
Vmass = Vmat * rho_mat;
Vcost = Vmass * cost_mat;
end
