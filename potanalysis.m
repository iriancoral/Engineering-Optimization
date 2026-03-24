function [Vmat, Vpot, s, sigma_max, A_base, Vmass, Vcost] = ...
    potanalysis(r1, h, theta_wall, t, rho_mat, rho_soil, K0, g_acc, sigma_allow, cost_mat, f_drain)
% Core analysis function for flower pot optimization.
% Simplified model: r2 = r1 + h*tan(theta_wall) (fixed wall angle), t = fixed wall thickness.
%
% -------------------------------------------------------
% SI units: m, N, kg, Pa, euro
% -------------------------------------------------------
%
% INPUT PARAMETERS:
%   Design variables:
%     r1          : Bottom (inner) radius of pot             [m]
%     h           : Height of pot                            [m]
%
%   Fixed parameters (from potparams.m):
%     theta_wall  : Wall taper angle from vertical           [deg]
%     t           : Wall thickness                           [m]
%     rho_mat     : Pot wall material density                [kg/m³]
%     rho_soil    : Moist potting soil density               [kg/m³]
%     K0          : Lateral earth pressure coefficient       [-]
%     g_acc       : Gravitational acceleration               [m/s²]
%     sigma_allow : Allowable wall stress                    [Pa]
%     cost_mat    : Material cost per kg                     [euro/kg]
%     f_drain     : Fraction of base area occupied by holes  [-]  (5 holes, 15%)
%
% OUTPUT PARAMETERS:
%   Vmat      : Material volume of pot wall + base           [m³]  <- OBJECTIVE
%   Vpot      : Internal volume of pot (truncated cone)      [m³]
%   s         : Slant height of lateral wall                 [m]
%   sigma_max : Maximum hoop stress in wall (at bottom)      [Pa]
%   A_base    : Internal base area                           [m²]
%   Vmass     : Mass of pot wall material                    [kg]
%   Vcost     : Material cost of pot                         [euro]

% --- Derived geometry ----------------------------------

% Top (inner) radius from fixed wall angle
% r2 = r1 + h * tan(theta_wall): pot widens by h*tan(theta) from bottom to top
r2 = r1 + h * tan(theta_wall * pi / 180);

% Slant height of lateral surface
s = sqrt(h^2 + (r2 - r1)^2);

% --- Objective: material volume ------------------------
% Lateral surface area of truncated cone + base disc minus drainage holes
% Holes occupy f_drain fraction of base area (5 holes, scaling with r1)
% (outer surface approximated as inner, since t << r1)
Vmat = t * (pi*(r1 + r2)*s + pi*r1^2*(1 - f_drain));

% --- Internal volume of truncated cone -----------------
Vpot = (pi*h/3) * (r1^2 + r1*r2 + r2^2);

% --- Structural: max hoop stress at bottom of wall -----
% Lateral soil pressure at depth h (max at base): p = K0 * rho_soil * g * h
% Hoop stress from thin-shell membrane theory: sigma = p * r / t
p_max     = K0 * rho_soil * g_acc * h;
sigma_max = p_max * r1 / t;

% --- Base area -----------------------------------------
A_base = pi * r1^2;

% --- Mass and cost -------------------------------------
Vmass = Vmat * rho_mat;
Vcost = Vmass * cost_mat;

% end
