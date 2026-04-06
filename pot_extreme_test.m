% pot_extreme_test.m
% Find where the 4-var model diverges from 3-var (t > t_min).
% Tests: (a) very large soil masses with PP, (b) weak material at normal masses.

clear; potparams;
params.theta_wall=theta_wall; params.t=t; params.rho_mat=rho_mat; params.rho_soil=rho_soil;
params.K0=K0; params.g_acc=g_acc; params.sigma_allow=sigma_allow; params.cost_mat=cost_mat;
params.f_drain=f_drain; params.stab_ratio=stab_ratio; params.t_min=t_min;

options = optimoptions('fmincon','Algorithm','sqp','Display','off', ...
    'OptimalityTolerance',1e-10,'ConstraintTolerance',1e-10, ...
    'StepTolerance',1e-12,'MaxIterations',500,'FiniteDifferenceType','forward');

r2_max=0.50; t_max_val=0.050;
lb4=[r1_min;h_min;r1_min;t_min];
ub4=[r1_max;h_max;r2_max;t_max_val];

fprintf('=== Extreme cases: where does t deviate from t_min? ===\n\n');

% Helper: run 4-var optimization for given params
r1f=[0.6 0.8 1.0 1.2 1.5 0.7 1.3 0.9];
hf =[0.8 1.2 1.0 1.5 0.7 1.4 0.9 1.1];
af =[1.2 1.4 1.3 1.5 1.2 1.6 1.3 1.4];
tf =[1.0 2.0 3.0 5.0 8.0 1.5 4.0 1.0];

% --- (a) Polypropylene at increasing masses ---
fprintf('--- PP (sigma_allow = 20 MPa), increasing mass ---\n');
masses_ext = [25 50 100 200 500 1000];
for mi = 1:length(masses_ext)
    ms = masses_ext(mi);
    params.m_soil = ms;
    params.sigma_allow = 20e6;
    V_req = ms/rho_soil;
    r_est = max(r1_min+0.01, min(ub4(1)-0.01, (V_req/pi)^(1/3)));

    best_f=Inf; best_x=[r_est;r_est;r_est*1.3;t_min];
    for si=1:length(r1f)
        r1s=max(lb4(1),min(ub4(1),r_est*r1f(si)));
        hs =max(lb4(2),min(ub4(2),r_est*hf(si)));
        r2s=max(lb4(3),min(ub4(3),r1s*af(si)));
        ts =max(lb4(4),min(ub4(4),t_min*tf(si)));
        try
            [xk,fk,ef]=fmincon(@(x)potobj_full(x,params),[r1s;hs;r2s;ts], ...
                [],[],[],[],lb4,ub4,@(x)potcon_full(x,params),options);
            if ef>0 && fk<best_f; best_f=fk; best_x=xk; end
        catch; end
    end
    [gv,~]=potcon_full(best_x,params);
    t_lb = abs(best_x(4)-t_min)<1e-6;
    fprintf('  m=%5dkg: r1=%6.1f h=%6.1f r2=%6.1f t=%6.3fmm g1=%+.2e t_on_lb=%d\n', ...
        ms, best_x(1)*1000, best_x(2)*1000, best_x(3)*1000, best_x(4)*1000, gv(1), t_lb);
end

% --- (b) Weak material (sigma_allow = 0.5 MPa) at normal masses ---
fprintf('\n--- Weak material (sigma_allow = 0.5 MPa) ---\n');
params.sigma_allow = 0.5e6;
masses_weak = [1 5 10 25];
for mi = 1:length(masses_weak)
    ms = masses_weak(mi);
    params.m_soil = ms;
    V_req = ms/rho_soil;
    r_est = max(r1_min+0.01, min(ub4(1)-0.01, (V_req/pi)^(1/3)));

    best_f=Inf; best_x=[r_est;r_est;r_est*1.3;t_min*3];
    for si=1:length(r1f)
        r1s=max(lb4(1),min(ub4(1),r_est*r1f(si)));
        hs =max(lb4(2),min(ub4(2),r_est*hf(si)));
        r2s=max(lb4(3),min(ub4(3),r1s*af(si)));
        ts =max(lb4(4),min(ub4(4),t_min*tf(si)));
        try
            [xk,fk,ef]=fmincon(@(x)potobj_full(x,params),[r1s;hs;r2s;ts], ...
                [],[],[],[],lb4,ub4,@(x)potcon_full(x,params),options);
            if ef>0 && fk<best_f; best_f=fk; best_x=xk; end
        catch; end
    end
    [gv,~]=potcon_full(best_x,params);
    t_lb = abs(best_x(4)-t_min)<1e-6;
    g1_act = abs(gv(1)) < 0.01*params.sigma_allow;
    th = atan2d(best_x(3)-best_x(1), best_x(2));
    fprintf('  m=%2dkg: r1=%6.1f h=%6.1f r2=%6.1f t=%6.3fmm theta=%5.1fdeg g1=%+.2e g1_active=%d Vmat=%.2fcm3\n', ...
        ms, best_x(1)*1000, best_x(2)*1000, best_x(3)*1000, best_x(4)*1000, th, gv(1), g1_act, best_f*1e6);
end

% --- (c) Reduced sigma_allow with PP to find transition point ---
fprintf('\n--- PP with reduced sigma_allow, m_soil=5kg, find transition ---\n');
params.m_soil = 5;
sig_vals = [20e6, 10e6, 5e6, 2e6, 1e6, 0.5e6, 0.2e6, 0.1e6];
for si_idx = 1:length(sig_vals)
    params.sigma_allow = sig_vals(si_idx);
    V_req = 5/rho_soil;
    r_est = (V_req/pi)^(1/3);

    best_f=Inf; best_x=[r_est;r_est;r_est*1.3;t_min*3];
    for si=1:length(r1f)
        r1s=max(lb4(1),min(ub4(1),r_est*r1f(si)));
        hs =max(lb4(2),min(ub4(2),r_est*hf(si)));
        r2s=max(lb4(3),min(ub4(3),r1s*af(si)));
        ts =max(lb4(4),min(ub4(4),t_min*tf(si)));
        try
            [xk,fk,ef]=fmincon(@(x)potobj_full(x,params),[r1s;hs;r2s;ts], ...
                [],[],[],[],lb4,ub4,@(x)potcon_full(x,params),options);
            if ef>0 && fk<best_f; best_f=fk; best_x=xk; end
        catch; end
    end
    [gv,~]=potcon_full(best_x,params);
    t_lb = abs(best_x(4)-t_min)<1e-6;
    g1_act = abs(gv(1)) < 0.01*params.sigma_allow;
    th = atan2d(best_x(3)-best_x(1), best_x(2));
    fprintf('  sig_allow=%8.1e: t=%6.3fmm theta=%5.1fdeg g1=%+.2e g1_active=%d t_on_lb=%d Vmat=%.2fcm3\n', ...
        sig_vals(si_idx), best_x(4)*1000, th, gv(1), g1_act, t_lb, best_f*1e6);
end

fprintf('\nDone.\n');
