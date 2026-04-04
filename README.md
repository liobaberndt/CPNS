# 22q11.2DS sleep thalamocortical modelling (2026)

## Abstract

22q11.2 deletion syndrome (22q11.2DS) is a strong genetic risk factor for neuropsychiatric conditions, including schizophrenia, yet the underlying synaptic mechanisms remain unclear. Sleep EEG suggests thalamocortical dysfunction, but scalp data alone lack mechanistic resolution. Computational modelling can bridge this gap by inferring receptor-level dynamics from EEG.

We applied a conductance-based thalamocortical Dynamic Causal Model (DCM) to sleep-wake EEG from children with 22q11.2DS (*n* = 28) and their neurotypical siblings (*n* = 17), estimating contributions of AMPA, NMDA, GABA<sub>A</sub>, and GABA<sub>B</sub> conductances. Building on these estimates, we investigated which receptor systems, if perturbed, could shift circuit dynamics toward sibling patterns. To address this, we implemented *in silico* pharmacology by systematically scaling receptor-mediated conductances.

Increasing NMDA receptor (NMDA-R) efficacy consistently produced the strongest improvements in alignment with sibling spectra (effect size = 0.32 in light NREM, 0.45 in deep NREM), whereas AMPA- or GABA-based manipulations were weaker. The most influential pathways were recurrent NMDA-R excitation among superficial pyramidal cells and NMDA-R excitation from spiny stellates to superficial pyramidal populations. Exploratory regressions linked greater thalamocortical delay during deep sleep to greater sleep problems (*β* = 0.34, *p*<sub>FDR</sub> = 0.006), and AMPA-mediated excitation of interneurons during wakefulness to anxiety symptoms (*β* = −0.32, *p*<sub>FDR</sub> = 0.032).

These findings implicate NMDA receptor hypofunction as a key mechanism in 22q11.2DS and suggest it may serve as a treatment target. We further show that DCM-based virtual pharmacology can simulate drug-level interventions, and we are now testing whether NMDA receptor modulation restores network activity in preclinical models (e.g., mouse models of 22q11.2DS).

---

*Funding (manuscript version): This research was funded in whole, or in part, by the Wellcome Trust (226709/Z/22/Z). For the purpose of Open Access, the author has applied a Creative Commons Attribution (CC BY) public copyright licence to any Author Accepted Manuscript version arising from this submission.*

---

## Code in this repository

Paths are under `code/`.

| Script | Role |
|--------|------|
| `RunTCM_Script_transfun_22qSleep_April2024_fixed.m` | **Primary DCM pipeline:** specifies the thalamocortical model (forward model `@atcm.tc_hilge2`, spectral observer `@atcm.fun.alex_tf`, priors/preparation via `atcm.*`), loads EEG, and runs inversion. |
| `PEB_22q.m` | Second-level **PEB** (22q vs sibling) per sleep stage and receptor family; saves `GCM` / `PEB` / `BMA` outputs. |
| `lasso_regr.m` | Nested CV LASSO linking DCM parameters to behavioural outcomes (with FDR). |
| `Simulation_22q_combined.m` | Combined **in silico** receptor sweep / simulation across stages. |
| `sim_analysis.m` | Summarises and compares simulation outputs. |

### `code/helperfunctions/`

Shared MATLAB helpers used by the pipelines above. Add `code` and `code/helperfunctions` to the MATLAB path as in each script’s header block.

### External modelling libraries: **atcm** and **aoptim**

- **atcm** - thalamocortical modelling for DCM (integrators, spectral functions, parameter handling). Upstream: [github.com/alexandershaw4/atcm](https://github.com/alexandershaw4/atcm).  
- **aoptim** - optimisation routines used in the inversion workflow. Upstream: [github.com/alexandershaw4/aoptim](https://github.com/alexandershaw4/aoptim).
