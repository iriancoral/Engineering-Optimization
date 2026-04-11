# Final Report Progress Tracker

**Course:** ME46060 Engineering Optimization
**Project:** Optimal Flower Pot Design
**Student IDs:** 4894413 (partner) & 5170540 (Mehrag)
**Last updated:** 2026-04-11
**Final report deadline:** Fri Apr 17, 23:00 (~6 days remaining)

**Repos:**
- Code: `Engineering-Optimization/` (GitHub: https://github.com/iriancoral/Engineering-Optimization, branch master)
- Report: `FlowerPotOptimization/` (separate repo, LaTeX)

---

## Executive summary

The 2-variable AND extended (3/4-variable + parametric) MATLAB pipelines are
**all implemented and the report is ~80% drafted**. What's left is mostly
finishing/editing — not new code or new analysis.

| Area | Status |
|------|--------|
| 2-var MATLAB code | DONE |
| 3-/4-var MATLAB code | DONE |
| Parametric soil-mass sweep | DONE |
| Figures (16 PNGs) | DONE — all in `FlowerPotOptimization/figures/` |
| Section 1 Introduction | Drafted (1 leftover Dutch TODO) |
| Section 2 Problem Formulation | Done, polished |
| Section 3 Initial Investigation | Done, 1 Dutch TODO in monotonicity |
| Section 4 Simplified Optimization | Done, 2 Dutch TODOs to resolve |
| Section 5 Full Optimization | Done |
| Section 6 Conclusions | **STUB — to write** |
| Statement of Contribution | **STUB — table empty** |
| Appendix (code listings) | **STUB — code/ folder empty, listings commented out** |
| Commercial pot comparison | **DROPPED from report scope** (was in old plan, no section/code exists) |

---

## MATLAB code (Engineering-Optimization/) — STATE ON DISK

### 2-variable pipeline (DONE)
| File | Purpose |
|------|---------|
| `potparams.m` | Fixed parameters (PP material, soil, geometry, bounds) |
| `potanalysis.m` | Core analysis: (r1,h) → (Vmat, Vpot, sigma, area, mass, cost) |
| `potobj.m` | Objective wrapper |
| `potcon.m` | Constraint wrapper (g1 stress, g2 stability, ceq volume) |
| `potcontour.m` | Step 3 — contour map + 1D Vmat-vs-r1 along volume curve |
| `potsensitivity.m` | Step 6 — FD step study, log-sensitivities, Hessian |
| `potpenalty.m` | Step 7 — self-implemented exterior penalty (course requirement) |
| `potfmincon.m` | Step 8 — fmincon SQP, 20 multi-starts → `potfmincon_result.mat` |
| `potkkt.m` | Step 9 — KKT verification |

**Run order:** `potcontour → potsensitivity → potpenalty → potfmincon → potkkt`

### 3-/4-variable extended pipeline (DONE)
| File | Purpose |
|------|---------|
| `potanalysis_full.m` | Extended analysis with r2 (and t) free; quadratic stress in z |
| `potobj_full.m` | Extended objective wrapper |
| `potcon_full.m` | Extended constraints (adds g3: r1 ≤ r2) |
| `potfmincon_full.m` | fmincon SQP for 3-var and 4-var models → `potfmincon_full_result.mat` |
| `potparametric.m` | Soil-mass sweep ms ∈ {1, 5, 10, 25} kg → `potparametric_result.mat` |
| `pot_extreme_test.m` | Sanity check at design-domain extremes |

**Result files saved:**
`potfmincon_result.mat`, `potfmincon_full_result.mat`, `potparametric_result.mat`

### Key numerical results
- 2-var optimum (5 kg): r1* = 104.51 mm, h* = 102.63 mm, **Vmat* = 310.5747 cm³**
- 3-var optimum (5 kg): r1* = 89.1 mm, h* = 104.8 mm, r2* = 134.3 mm, theta* = 23.4°, **Vmat* = 304.04 cm³** (2.1% saving)
- 4-var: identical to 3-var (t lands on lower bound t_min = 3 mm at every soil mass)
- Scaling law: r1, h, r2 ∝ ms^(1/3); Vmat ∝ ms^(2/3) — exact across 1–25 kg
- KKT: both inequalities inactive at all optima → optimum determined by volume equality only; lambda_ceq ≈ −4.97e-2

---

## Report (FlowerPotOptimization/) — STATE ON DISK

### Build setup
- `main.tex`: 11pt article, biblatex (sorting=none), TikZ, listings, all 7 inputs wired
- `references.bib`: 8 entries (Das, Timoshenko, Mellouli, Papalambros, Rao, Langelaar lec5, Langelaar lec10, Wikipedia num diff)
- `figures/`: 16 PNGs all present and used in section files (potcontour ×2, potsensitivity ×3, potpenalty ×3, potfmincon ×1, potfmincon_full ×3, potparametric ×4)
- `code/`: empty except `.gitkeep`

### Section status

**`00_contribution.tex` — STUB (~20 lines)**
- Table skeleton present, all rows blank except "Problem definition: (together)"
- TODO: fill out with actual split before submission

**`01_introduction.tex` — DRAFT (1 paragraph)**
- Frames the problem, names methods, references Papalambros/Das/Timoshenko
- TODO: contains literal Dutch comment `!!!!MOETEN WE DAN WEL NOG DOEN!!!!` next to the parametric-study sentence — the parametric study IS done, so just rephrase to remove the TODO.

**`02_problem_formulation.tex` — COMPLETE**
- Geometry + design variables, full TikZ frustum figure with r1/r2/h/t/theta_w
- Objective eq., 3 constraints (ceq, g1 stress, g2 stability), boxed math problem
- Drainage-area constraint discussed and explicitly omitted (always inactive)

**`03_initial_investigation.tex` — COMPLETE w/ 1 TODO**
- Boundedness via cubic-in-h derivation, contour figure + 1D Vmat curve
- Monotonicity subsection title still says `VOOR MIER NOG CHECKEN` (Dutch: "still to check by Mier") → strip TODO marker after partner review
- Convexity: numerical Hessian at x0=(100,150)mm, indefinite eigenvalues, justified via 1D feasibility curve + multi-start
- Numerical noise / FD step-size study: V-shape, h*=1e-8 chosen
- Sensitivity: gradients, log-sensitivities (r1 dominates h), discussion of fixed-parameter sensitivities

**`04_simplified_optimization.tex` — COMPLETE w/ 2 TODOs**
- Reduced formulation, graphical solution, exterior penalty (own implementation), fmincon SQP, comparison table, full 4-step KKT verification
- TODO 1 (line ~126): `MISSCHIEN OVERBODIG HIER:` — decide whether the trailing paragraph about t being over-designed should stay or be cut (it's a useful bridge to Section 5, recommend KEEP and remove the marker)
- TODO 2 (line ~135): `\subsection{Observations WSS ONNODIG KAN MOGELIJK GESKIPT WORDEN}` (Dutch: "probably unnecessary, can possibly be skipped") — decide whether to keep/cut. Content is short and adds aspect-ratio and material-share commentary; consider folding into Section 5 observations or dropping.

**`05_full_optimization.tex` — COMPLETE**
- Extended formulation (3- and 4-var), generalized stress eq with critical-depth analysis, results table for 5 kg, parametric study table for 1/5/10/25 kg, log-log scaling figures, scaling-law derivation, role of t analysis, observations
- No TODOs

**`06_conclusions.tex` — STUB (commented-outline only, ~16 lines)**
- TODO: write the section. Outline already in file:
  - summary of optimal pot geometry
  - material savings vs. commercial pots (will need to be reframed since commercial comparison was dropped — frame as savings vs. 2-var baseline / vs. fixed-angle assumption)
  - which algorithms/formulations worked best
  - recommendations for designers
  - recommendations for further work
  - reflection on optimization techniques + limitations

**`appendix_code.tex` — STUB**
- Header `\section{MATLAB Code}` only; all `\lstinputlisting{}` calls commented out
- Two preconditions before this can be uncommented:
  1. Copy the .m files into `FlowerPotOptimization/code/` (currently only has .gitkeep)
  2. Decide which files to include — likely all 12 (potparams, potanalysis(_full), potobj(_full), potcon(_full), potcontour, potsensitivity, potpenalty, potfmincon(_full), potkkt, potparametric)

---

## What needs to be done before Apr 17

### Required (blocks submission)
1. **Write Section 6 Conclusions** — fill out the outline already in `06_conclusions.tex`. ~1 page.
2. **Fill out Statement of Contribution table** in `00_contribution.tex`. Coordinate with partner.
3. **Resolve the 4 Dutch TODOs** still in section files:
   - `01_introduction.tex` line 6: `!!!!MOETEN WE DAN WEL NOG DOEN!!!!`
   - `03_initial_investigation.tex` §3.2 title: `VOOR MIER NOG CHECKEN`
   - `04_simplified_optimization.tex` line ~126: `MISSCHIEN OVERBODIG HIER:`
   - `04_simplified_optimization.tex` §4.7 title: `WSS ONNODIG KAN MOGELIJK GESKIPT WORDEN`
4. **Appendix code listings** — copy .m files into `FlowerPotOptimization/code/` and uncomment the `\lstinputlisting` calls in `appendix_code.tex`.

### Recommended polish
5. **Compile and check page count** — target 16–20 pages excluding appendices. If over, prime candidates to trim are the redundant FD step-size paragraph in §3.4 (currently has a long commented-out version AND a written version saying the same thing).
6. **Cite AI tools** if used (course requirement — plagiarism rules).
7. **Final proofread** for English consistency (some sections were drafted by partner in mixed Dutch/English).

### Already decided NOT to do (so PROGRESS reflects reality)
- Commercial pot comparison (was in original plan and old PROGRESS, but no section exists in `main.tex` and no `potcommercial.m` was written). Section 5 instead frames savings vs. the 2-var baseline. If desired, a paragraph comparing the optimal 2.1% improvement to typical commercial wall angles (5–10° vs. our 23.4°) is already in §5.8 — that essentially fulfils the spirit of the comparison without needing a separate section.

---

## Reusable report phrases (kept from previous PROGRESS for consistency)

**Stability constraint justification:**
> "The stability constraint (h ≤ 3·r1) is not a physical tipping analysis but a geometric bound
> that prevents the optimizer from finding unrealistic solutions such as r1 = 10 mm with h = 400 mm,
> which would satisfy all other constraints but clearly does not represent a practical flower pot."

**Drainage constraint justification:**
> "The drainage constraint (π·r1² ≥ 3 cm²) is never active at the optimum — with r1_min = 30 mm
> the minimum base area is already 28.3 cm², far exceeding the 3 cm² requirement."

**Material justification:**
> "PP homopolymer high flow (CES EduPack) selected for its widespread use in commercial pot
> manufacturing, low density, low cost, and sufficient strength for applied soil pressure loads."

**Hoop stress formula — cylinder vs cone:**
> "The hoop stress formula sigma = p·r/t (cylinder) is also correct for horizontal soil pressure on a
> conical wall. The conical shell formula requires the pressure perpendicular to the wall;
> p_normal = p_horizontal·cos(theta), and substituting cancels the cos(theta). No correction needed."

**Wall angle model (alpha → theta_wall):**
> "Fixed wall angle theta_wall = 10° gives r2 = r1 + h·tan(theta_wall), physically more meaningful
> than a fixed ratio r2/r1. Consistent with injection moulding practice. theta_wall < phi_soil = 30°
> so soil does not slide along the wall and the at-rest K0 pressure model remains valid."

---

## Key design parameters (reference table)

| Parameter | Value | Reason |
|-----------|-------|--------|
| Material | Polypropylene (PP homopolymer, high flow) | Lightweight, cheap, realistic |
| rho_mat | 910 kg/m³ | PP density |
| sigma_allow | 20 MPa | PP yield stress |
| t (2-var) | 3 mm fixed | Minimum practical injection-mould thickness |
| t (4-var bounds) | [3, 10] mm | Lower bound = manufacturing minimum |
| cost_mat | 2 €/kg | PP material cost |
| rho_soil | 1200 kg/m³ | Typical potting soil |
| K0 | 0.5 | Jaky's formula, phi = 30° |
| theta_wall (2-var) | 10° fixed | r2 = r1 + h·tan(10°), injection-moulding-realistic |
| g2 (stability) | h ≤ 3·r1 | Geometric bound vs unrealistic tall/narrow |
| f_drain | 0.15 | 15% of base area removed from Vmat |
| n_holes | 5 | Fixed number of drainage holes |
| Bounds | r1 ∈ [30,300], h ∈ [50,400], r2 ∈ [30,400] mm | Hand-held pot range |
