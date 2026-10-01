# Research directions

[中文](research-directions.zh-CN.md) · [README](../README.md)

Literature checked on 2026-10-01.

The central opportunity is to connect prescribed QBO targets, realized model winds and East Asian responses under two background climates. The existing factorial design and boundary operators provide a concrete starting point.

## Background dependence of the regional response

Keep the two background periods, W/E targets and primary endpoints fixed. Estimate the interaction

\[
I=(Y_{W,late}-Y_{E,late})-(Y_{W,early}-Y_{E,early}).
\]

The primary outcomes remain July–August ascent and wind-driven moisture-flux convergence evaluated with a common reference humidity. Rainfall dipole variability is secondary. Observational work has already established relevant questions about East Asian drought/flood prediction and QBO vertical structure; the controlled background interaction is the next question for this project. [Zhang et al. (2024)](https://www.nature.com/articles/s41467-023-44445-y), [Luo et al. (2023)](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2023GL105863).

Archive the actual boundary settings in all four cells. This repository currently builds SST and sea-ice backgrounds; greenhouse gases, ozone and aerosols must be documented from the model configurations. Attribution follows the boundary conditions that actually change. The historical SST/ice contrast combines forced change and internal variability.

Pair W/E initial conditions within each background. Resample independent member pairs, using cross-background blocks only when the design supplies a justified pairing. Grid cells and daily samples are not independent ensemble members. Estimate the variance and runtime in a bounded pilot, then freeze the ensemble size before examining production outcomes. Report within-background effects, their interaction and confidence intervals under the existing endpoint rules.

## Realized QBO control

Identical target files can produce different realized winds. Here the boundary weights depend on tropopause geometry, and the zonal minimum can change with the background. This creates a practical identification question: how much of a regional interaction accompanies a difference in the realized QBO contrast?

An offline evaluation of six archived initialization/short-run states gives a concrete starting point. Within 5° of the equator, v5 retains about 49.33% of the original relaxation rate at the layer centered near 61.52 hPa and zero at 73.75 hPa; the hard uniform-support bound is also zero at the latter level. Levels are labeled by midpoint pressure, while the operator uses layer-bottom pressure. The next physical test is whether control aloft can realize the intended lower-level QBO structure through model dynamics. [Archived values and source hash](archived-support.json).

Recent studies show that QBO nudging improves some teleconnections while other biases persist, and that lower-stratospheric amplitude and analysis choices affect teleconnection estimates. [Andrews et al. (2026)](https://wcd.copernicus.org/articles/7/1797/2026/), [Garfinkel et al. (2026)](https://wcd.copernicus.org/articles/7/1133/2026/).

For fixed latitudes, levels and time windows, save target/realized winds, bias and RMSE; the realized W−E profile in each background; 30–50–70 hPa structure and descent timing; effective weights and relaxation rates; accumulated increments; and zonal/eddy components of the applied tendency. Compare these diagnostics with prespecified tolerances and retain member distributions.

Report the QBO contrast alongside the regional response. Dividing by a near-zero contrast is unstable. Adjusting a regression for realized winds, a post-treatment variable, changes the estimand and does not automatically recover a controlled background effect.

A methods extension can quantify feasible control under the frozen boundary constraint. If each column permits at most \(w_i\), a longitude-uniform weight must satisfy \(a\leq w_i\) for every column, so its largest feasible value is \(\min_i w_i\). The scientific work is to measure how that bound varies across seasons and backgrounds on the fixed grid, and how it relates to realized wind error. Count missing-diagnostic restrictions separately from physical geometry. Real model states and wind responses are needed alongside numerical kernel checks. The QBOi phase-2 nudging protocol offers a useful diagnostic reference and remains a preprint. [Anstey et al. (2026)](https://egusphere.copernicus.org/preprints/2026/egusphere-2026-1165/).

## Predictive value across decades

A separate forecasting study could compare a fixed baseline, baseline plus March QBO70, and baseline plus a prespecified low-dimensional QBO profile. A physically motivated background interaction can then be tested as one additional model. Select any regularization or dimension reduction inside the training data.

Use rolling-origin outer validation with contiguous held-out years. Fit climatologies, trends, scaling and EOFs within each training window, respecting data-release delays. The full-period target composites used for the prescribed-forcing experiment have a different role from training-window transformations in a forecast evaluation.

Publish year-level predictions, held-out correlation and RMSE, improvement over the baseline, and time-block uncertainty. Probability forecasts additionally need Brier scores and calibration. Historical data already used during exploration should be distinguished from a genuinely untouched validation period.

A recent study finds a changing QBO–typhoon-track relationship in the western North Pacific. It motivates a transportability check, while providing no direct result about the stability of this project's rainfall relationship. [Kim et al. (2026)](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2026GL125034/).

## External comparisons and order of work

QBOi perpetual-ENSO experiments and LESFMIP large ensembles can help assess background dependence, model differences and internal variability. ENSO-related changes in QBO period are more consistent across the QBOi models than amplitude changes. [Kawatani et al. (2025)](https://wcd.copernicus.org/articles/6/1045/2025/). Inventory available winds and regional endpoint variables before requesting regional/level subsets of larger archives. The LESFMIP data used in the cited study are available through ESGF; QBOi archive access requires contacting the coordinators and applying through JASMIN. [Archive access](https://wcd.copernicus.org/articles/7/1797/2026/#section8).

The immediate sequence is full-driver equivalence, realized-QBO fidelity, then the factorial ensemble. A small annual forecasting dataset and a fixed validation protocol can be prepared independently. New endpoints or mechanisms require a separate scientific design.

The focused search covered East Asian precipitation, QBO vertical structure, nudging, background dependence and nonstationarity. Records were deduplicated by DOI, with abstracts and publication status checked on publisher pages. Six linked works are journal articles; Anstey et al. is a preprint.
