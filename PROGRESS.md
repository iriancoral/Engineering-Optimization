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

**Wall angle model (alpha -> theta_wall):**
> "Fixed wall angle theta_wall = 10 deg gives r2 = r1 + h*tan(theta_wall), physically more meaningful
> than a fixed ratio r2/r1. Consistent with injection moulding practice. theta_wall < phi_soil = 30 deg
> so soil does not slide along the wall and at-rest K0 pressure model remains valid."

---
**Course:** ME46060 Engineering Optimization
**Project:** Optimal Flower Pot Design
**Student IDs:** 4894413 & 5170540
**Last updated:** 2026-03-24

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
| potkkt.m | Step 9: KKT verification + gradient geometry plot | Done |

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
| alpha (taper ratio) | 1.3 (r2 = 1.3 * r1) | Fixed in 2-var, free in 3-4 var |
| Stability constraint | h <= 3 * r1 | Pot must not tip over |
| Drainage constraint | pi * r1^2 >= 3e-4 m2 | Minimum drainage hole area |

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
