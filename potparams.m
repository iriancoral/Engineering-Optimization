% Flower pot optimization - fixed parameters
% Engineering Optimization ME46060
%
% SI units: m, N, kg, Pa, euro
% -------------------------------------------------------

% --- Material properties (polypropylene plastic) CES EDUPACK PP (homopolymer, high flow)-------
rho_mat     = 910;          % Pot wall material density       [kg/m³]
sigma_allow = 20e6;         % Allowable wall stress           [Pa]
t_min       = 0.003;        % Minimum wall thickness          [m]  (3 mm)
cost_mat    = 2.0;          % Material cost                   [euro/kg]

% --- Soil properties -----------------------------------
rho_soil    = 1200;         % Moist potting soil density      [kg/m³]
K0          = 0.5;          % Lateral earth pressure coeff.   [-]  (Jaky: K0 = 1-sin(30°))
phi_soil    = 30;           % Internal friction angle         [deg]
g_acc       = 9.81;         % Gravitational acceleration      [m/s²]

% --- Target soil mass (varies in parametric study) -----
m_soil      = 5;            % Target soil mass                [kg]

% --- Simplified model: fixed wall angle ----------------
theta_wall  = 10;           % Wall taper angle from vertical   [deg]
% r2 = r1 + h*tan(theta_wall): angle is fixed, alpha = r2/r1 varies with h and r1

% --- Wall thickness (simplified model: t fixed at t_min)
t           = t_min;        % Wall thickness used in analysis [m]

% --- Geometric bounds (design variable bounds) ---------
r1_min      = 0.03;         % Minimum bottom radius           [m]  (30 mm)
r1_max      = 0.30;         % Maximum bottom radius           [m]  (300 mm)
h_min       = 0.05;         % Minimum height                  [m]  (50 mm, min root depth)
h_max       = 0.40;         % Maximum height                  [m]  (400 mm)

% --- Stability constraint ------------------------------
stab_ratio  = 3.0;          % h <= stab_ratio * r1 (prevents tipping)

% --- Drainage holes ------------------------------------
n_holes     = 5;            % Number of drainage holes         [-]
f_drain     = 0.15;         % Holes occupy 15% of base area   [-]
% Hole diameter scales with pot: d_hole = 2*r1*sqrt(f_drain/n_holes)

% end
