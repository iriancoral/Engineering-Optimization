# Final Report Progress Tracker

## Report Phrases to Use

**Stability constraint justification:**
> "The stability constraint (h <= 3*r1) is not a physical tipping analysis but a geometric bound
> that prevents the optimizer from finding unrealistic solutions such as r1 = 10 mm with h = 400 mm,
> which would satisfy all other constraints but clearly does not represent a practical flower pot."

**Drainage constraint justification:**
> "The drainage constraint (pi*r1^2 >= 3 cm^2) is never active at the optimum — with r1_min = 30 mm
> the minimum base area is already 28.3 cm^2, far exceeding the 3 cm^2 requirement."

**Upper bounds justification:**
> "Upper bounds r1_max = 300 mm and h_max = 400 mm represent practical limits for a hand-held pot.
> These bounds were found to be non-active at the optimum for all tested soil masses."

**Material justification:**
> "PP homopolymer high flow (CES EduPack) selected for its widespread use in commercial pot
> manufacturing, low density, low cost, and sufficient strength for applied soil pressure loads."

**Thin-shell approximation:**
> "A thin-shell approximation (t << r1, h) is used throughout the model. Small geometric corrections
> of order t are neglected (base area, soil height, inner vs outer radius, cylinder vs cone hoop stress).
> This is justified since the focus of this project is on the optimization methodology, not on detailed
> structural analysis."

**Hoop stress formula — cylinder vs cone:**
> "The hoop stress formula sigma = p*r/t (cylinder) is also correct for horizontal soil pressure on a
> conical wall. The conical shell formula sigma = p_normal * r / (t * cos(theta)) requires the pressure
> perpendicular to the wall. The horizontal soil pressure decomposes as p_normal = p_horizontal * cos(theta).
> Substituting: sigma = p_horizontal * cos(theta) * r / (t * cos(theta)) = p_horizontal * r / t.
> The cos(theta) cancels — no correction needed."

**Step size study — why error curves don't follow O(h) exactly:**
> "The O(h) reference line assumes the function has significant higher-order terms (h^2, h^3, ...) in
> the Taylor expansion. If the function is nearly linear — like ours — those higher-order terms are tiny,
> so the truncation error is already very small even at large h. That is why the error curve sits below
> the O(h) line and does not follow it closely. O(h) represents the worst case for forward FD; our
> functions are well-behaved, so we get better-than-worst-case accuracy. The reference line is still
> useful — it shows the maximum expected slope, and the error curves never go steeper than O(h)."

**Wall angle model (alpha -> theta_wall):**
> "Fixed wall angle theta_wall = 10 deg gives r2 = r1 + h*tan(theta_wall), physically more meaningful
> than a fixed ratio r2/r1. Consistent with injection moulding practice. theta_wall < phi_soil = 30 deg
> so soil does not slide along the wall and at-rest K0 pressure model remains valid."

---
**Course:** ME46060 Engineering Optimization
**Project:** Optimal Flower Pot Design
**Student IDs:** 4894413 & 5170540
**Last updated:** 2026-04-04
**GitHub:** https://github.com/iriancoral/Engineering-Optimization (branch: master)
- Push: git add . -> git commit -m "msg" -> git push
- Pull: git pull
- Partner clones with: git clone https://github.com/iriancoral/Engineering-Optimization.git

---

## Work Split
- **Irian:** Full 2-variable model + all analysis scripts (COMPLETE)
- **Partner:** 3-4 variable extension, parametric study, commercial comparison, report writing

---

## 2-Variable MATLAB Code — COMPLETE

All files in: `Final report/`

| File | Purpose | Status |
|------|---------|--------|
| potparams.m | Fixed parameters (material, soil, geometry, bounds) | Done |
| potanalysis.m | Core analysis: (r1,h) -> (Vmat, Vpot, s, sigma_max, A_base, Vmass, Vcost) | Done |
| potobj.m | Objective wrapper for optimizers | Done |
| potcon.m | Constraint wrapper: g(1..3) inequalities + ceq equality | Done |
| potcontour.m | Step 3: contour plots over (r1,h) grid + 1D view along volume curve | Done |
| potsensitivity.m | Step 6: FD step size study, log sensitivities, Hessian/convexity check | Done |
| potpenalty.m | Step 7: self-implemented exterior penalty method (course requirement) | Done |
| potfmincon.m | Step 8: fmincon SQP, multiple starts, saves potfmincon_result.mat | Done |
| potkkt.m | Step 9: KKT verification + contour plot (no gradient arrows — too small to see) | Done |

### Run Order
```
potcontour  ->  potsensitivity  ->  potpenalty  ->  potfmincon  ->  potkkt
```
**Note:** potpenalty saves potpenalty_result.mat (used by potfmincon for comparison)
**Note:** potfmincon saves potfmincon_result.mat (required by potkkt)

---

## Partner Still Needs to Write

| File | Purpose |
|------|---------|
| potanalysis_full.m | Extend potanalysis with r2 and/or t as free variables |
| potcon_full.m | Updated constraints for 3-4 variable model |
| potfmincon_full.m | fmincon on 3-4 variable problem |
| potparametric.m | Soil mass sweep: 1, 5, 10, 25 kg |
| potcommercial.m | Compare optimum with real commercial pots |

---

## Recent Changes (2026-04-04)

- potkkt.m: Removed `clc` (was clearing terminal output when figures loaded in -nodesktop mode)
- potkkt.m: Removed gradient arrow code (arrows invisible due to tiny gradient values vs axis scale)
- potkkt.m: Added g2 stability constraint line + legend to contour plot
- potkkt.m: Changed contour levels from prctile(80) to full max range
- potfmincon.m: Removed g1 stress contour + legend entry (not visible in plot range)

## KKT Results Summary

- No active inequality constraints (both g1, g2 inactive) → optimum determined by volume equality only
- KKT reduces to Lagrange conditions (mu1=mu2=0, lambda=-0.0497)
- Stationarity residual ≈ 1.5e-11 (satisfied)
- All 4 KKT conditions satisfied
- Insight: stress constraint has large margin → wall thickness `t` could be a design variable in the extended model

---

## Steps Still To Do (Both)

- Steps 4 & 5: Written analysis of boundedness, monotonicity, convexity — no new code needed, based on contour + sensitivity results
- Step 14: Commercial pot comparison data (partner)
- Step 15: Write the actual report document

---

## Key Design Decisions

| Parameter | Value | Reason |
|-----------|-------|--------|
| Material | Polypropylene plastic | Lightweight, cheap, realistic |
| rho_mat | 910 kg/m3 | PP density |
| sigma_allow | 20 MPa | PP yield stress |
| t (wall thickness) | 3 mm fixed in 2-var model | Minimum practical thickness |
| cost_mat | 2 euro/kg | PP material cost |
| rho_soil | 1200 kg/m3 | Typical potting soil |
| K0 | 0.5 | Jaky formula, phi = 30 deg |
| theta_wall | 10 deg | Fixed wall angle, r2 = r1 + h*tan(10deg), replaces fixed alpha |
| Stability constraint | h <= 3 * r1 | Geometric bound, prevents unrealistic tall/narrow pots |
| Drainage holes | 5 holes, f_drain = 0.15 | 15% of base area removed from Vmat, scales with r1 |
| n_holes | 5 | Fixed number of drainage holes |
| f_drain | 0.15 | Holes occupy 15% of base area, d_hole = 2*r1*sqrt(0.15/5) |

---

## Report Structure (target 16-20 pages)

| Section | Owner | Status |
|---------|-------|--------|
| Statement of Contribution | Both | Not started |
| 1. Introduction | Partner | Not started |
| 2. Problem Formulation | Irian | Not started |
| 3. Initial Investigation (contours, boundedness, convexity, sensitivity) | Irian | Not started |
| 4. Simplified Optimization (2-var, penalty + fmincon, KKT) | Irian | Not started |
| 5. Full Optimization (3-4 var, parametric study, algorithm comparison) | Partner | Not started |
| 6. Commercial comparison | Partner | Not started |
| 7. Conclusions & Recommendations | Both | Not started |
| References | Both | Not started |
| Appendix (MATLAB code) | Both | Not started |
