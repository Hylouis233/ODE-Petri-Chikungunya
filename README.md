# Dual-Model Framework for CHIKV Transmission Modeling: ODE and Petri Net Analysis of the 2025 Foshan Outbreak

<div align="center">

[English](README.md) | [中文](README_CN.md)

</div>

<div align="center">

![arXiv](https://img.shields.io/badge/arXiv-XXXX.XXXXX-red.svg)
![Zenodo](https://img.shields.io/badge/Zenodo-xx.xxxx%2Fzenodo.XXXXXXX-blue.svg)
![ORCID](https://img.shields.io/badge/ORCID-0000--0002--7423--0090-green.svg)
![GitHub](https://img.shields.io/github/stars/Hylouis233/ODE-Petri-Chikungunya?style=social)

</div>

This repository contains the complete codebase and data for the research paper "Dual-Model Framework for CHIKV Transmission Modeling: ODE and Petri Net Analysis of the 2025 Foshan Outbreak". This study provides a comprehensive comparative analysis of Ordinary Differential Equation (ODE) and Petri Net models for Chikungunya virus transmission dynamics.

## Epidemic Background

### 2025 Foshan Chikungunya Outbreak Overview

In 2025, a Chikungunya outbreak occurred in Foshan City, Guangdong Province, China. This study constructs a dual-model framework based on ODE and Petri Nets to conduct a comparative analysis of Chikungunya virus transmission dynamics and basic reproduction number estimation methods. The research employs the SEICR (Susceptible-Exposed-Infectious-Chronic-Recovered) compartmental model, specifically distinguishing between symptomatic and asymptomatic infected individuals with different transmission characteristics, and evaluates the effectiveness of intervention measures through three-phase fitting.

![Foshan Cases Barplot](figure/foshan_cases.png)

The above figure shows the epidemiological data of the 2025 Chikungunya outbreak in Foshan City, including time series of daily incidence and cumulative cases, which were used for model fitting and validation.

### Epidemic Development Overview

In July 2025, a Chikungunya outbreak occurred in Foshan City, Guangdong Province, originating from imported cases discovered through active surveillance in Shunde District on July 8. The epidemic spread rapidly and was effectively controlled through prevention and control measures:

<div align="center">

| Time Point | Confirmed Cases | Main Distribution Areas | Epidemic Characteristics | Prevention and Control Measures |
|------------|-----------------|------------------------|-------------------------|--------------------------------|
| July 19 | 1,873 cases | Shunde District: 1,790 cases, Chancheng District: 49 cases, Nanhai District: 34 cases | Mainly concentrated in Lecong, Beijiao, Chencun towns in Shunde District | Early surveillance and case detection |
| July 22 | 2,934 cases | Shunde District only | Continued growth, spreading to Macau Special Administrative Region | Enhanced surveillance |
| July 30 | 5,158 cases | Guangdong Province cumulative, Foshan City accounting for 93% | Level III response activated, covering 53 towns and streets | Level III response activated, designated 53 designated hospitals, providing 3,696 mosquito-proof isolation beds |
| August 25 | <50 cases/day for 9 consecutive days | Multi-regional distribution | Epidemic entered low-level sporadic state | Continuous surveillance, maintaining prevention and control measures |
| August 26 | Level III response terminated | City-wide | Prevention and control normalized, achieving community clearance | Level III response terminated, prevention strategy normalized |
| September 14-20 | 0 new cases | City-wide | No new cases for consecutive days, epidemic effectively controlled | Mosquito vector density reduced to safe levels, social activities resumed normal |

</div>

### Three-Phase Intervention Measures

Based on the official timeline and implemented control measures, this study divides the epidemic into four phases for modeling:

<div align="center">

| Phase | Time Range | Transmission Rate | Main Intervention Measures |
|-------|------------|-------------------|---------------------------|
| Initial Phase | 2025-06-16 to 2025-07-09 | β₀ | No control measures, natural transmission |
| Phase 1 | 2025-07-09 to 2025-07-14 | β₁ | Early surveillance and case detection |
| Phase 2 | 2025-07-14 to 2025-07-18 | β₂ | Medical institutions initiated nucleic acid testing |
| Phase 3 | After 2025-07-18 | β₃ | Centralized mosquito breeding site cleanup |

</div>

## Research Overview

### Research Abstract

This study addresses the 2025 Chikungunya outbreak in Foshan City, Guangdong Province, China, by constructing a dual-model framework based on Ordinary Differential Equations (ODE) and Petri Nets (PN) for comparative analysis of Chikungunya transmission dynamics and reproduction number estimation methods. The research employs SEICR (Susceptible-Exposed-Infectious-Chronic-Recovered) compartmental modeling to establish comparative analysis of two distinct modeling approaches under identical epidemiological scenarios, and evaluates the timeliness of epidemic control measures through a three-phase intervention fitting protocol.

Model validation results demonstrate excellent accuracy and consistency for both ODE and Petri Net models, with Mean Absolute Error (MAE) of 18.77 and 18.91 cases, respectively, and Root Mean Square Error (RMSE) of 36.52 and 36.54 cases. Both models predicted epidemic peaks at day 32 (406 cases), occurring 3 days earlier than the observed peak at day 35 (432 cases), with a peak value error of 6.0%.

Reproduction number analysis revealed that the initial basic reproduction number (R₀) was approximately 14.67 (ODE)/13.90 (PN), with effective reproduction numbers (R_eff) progressively decreasing through three intervention phases: 7.85/7.86 after Phase 1, 7.59/7.56 after Phase 2, and finally reaching 0.059 in Phase 3, successfully achieving epidemic transmission blockade.

Sensitivity analysis demonstrated that recovery rate (γ) represents the most sensitive parameter affecting R₀, with a Sobol index of 0.9672, explaining 96.72% of total R₀ variation, indicating that optimizing clinical treatment and accelerating patient recovery constitutes the most effective epidemic control strategy. This study presents the first systematic numerical comparison between ODE and Petri Net models, providing a novel dual-model analytical framework for vector-borne disease modeling such as Chikungunya, with significant theoretical guidance and practical value for epidemic prevention and control strategy formulation.

### Three-Phase Intervention Fitting

![Three-Phase Intervention Fitting](figure/combined_three_phase_fitting.png)

### Key Highlights
- **Dual-Model Comparison**: First systematic numerical comparison between ODE and Petri Net models.
- **Fitting Performance**: Both models achieve excellent fitting performance (MAE: 18.77/18.91 cases, RMSE: 36.52/36.54 cases).
- **Effective Intervention**: Three-phase control measures successfully reduced effective reproduction number from 14.7 to 0.059.
- **Comprehensive Analysis**: Including sensitivity analysis, confidence intervals, and uncertainty quantification.
- **Data Extraction**: Automatic case data extraction from official charts using OpenCV image processing technology, achieving open-source data sharing.

### Epidemic Control Effectiveness

Through the implementation of three-phase intervention measures, the epidemic was effectively controlled:

- **Phase 1 Effect**: Transmission coefficient reduced from β₀ to β₁, approximately 47% reduction
- **Phase 2 Effect**: Transmission coefficient further reduced to β₂, cumulative 48% reduction  
- **Phase 3 Effect**: Transmission coefficient dramatically reduced to β₃, cumulative 99% reduction, effective reproduction number dropped to 0.059

### Model Performance

- **Peak Prediction**: Both models predicted peak at day 32 (406 cases), while observed peak was at day 35 (432 cases), with peak prediction error of 6.0%.
- **Control Timeline**: Models predicted epidemic control achieved at day 65 (daily cases <0.1).

### Reproduction Numbers

- **Initial R₀**: ODE 14.67 vs Petri Net 13.90
- **Phase 1 R_eff**: 7.85 vs 7.86
- **Phase 2 R_eff**: 7.59 vs 7.56
- **Phase 3 R_eff**: 0.059 vs 0.059

### Sensitivity Analysis

- **Most Sensitive Parameter**: Recovery rate (γ) is the most influential parameter with Sobol index of 0.9672.
- **Variance Explanation**: γ explains 96.72% of R₀ variance.
- **Key Insight**: Accelerating patient recovery is the most effective control strategy.

## Model Validation and Sensitivity Analysis

![Model Validation Diagnostics](figure/model_validation_diagnostics.png)

The above figure shows the comprehensive validation and diagnostic analysis results of the Chikungunya epidemiological model, comparing the performance of Ordinary Differential Equation (ODE) and Petri Net models in residual analysis, prediction intervals, and phased reproduction number (R₀) comparison. Subplots include: (a) ODE model residual analysis, (b) Petri Net model residual analysis, (c) prediction interval analysis, (d) phased R₀ value comparison.

Subplots (a) and (b) show the residual distributions (observed - predicted) for ODE and Petri Net models, respectively. The residual statistics for both models are as follows: ODE model residual mean of -5.841, standard deviation of 36.049; Petri Net model residual mean of -6.005, standard deviation of 36.046. Both models exhibit similar systematic patterns in residuals rather than ideal random white noise characteristics. In the early epidemic phase (approximately first 20 days), residuals are close to zero, indicating good model fitting. However, during the epidemic peak period (approximately days 20-35), residuals systematically become negative, followed by large positive fluctuations. This non-random pattern suggests that the models may not fully capture dynamic changes caused by intervention measures or other external factors, particularly in the epidemic decline phase, where model predictions show systematic bias, underestimating the excess effectiveness achieved by government epidemic control measures.

Subplot (c) compares the prediction trajectories of both models and their 95% prediction intervals. The vast majority of observed values fall within the similar 95% prediction intervals of both models, validating the models' predictive capability. Subplot (d) compares reproduction number R₀ estimates across four different phases (initial phase and three intervention phases). Both models consistently reveal the phased decline trend of R₀ values with the implementation of intervention measures. The initial phase R₀ value is significantly higher than 1.0, indicating rapid epidemic transmission. Subsequent intervention measures significantly reduced R_eff values, finally achieving R_eff < 1.0 in the third intervention phase, marking effective epidemic transmission control.

![Sensitivity Analysis](figure/sensitivity_analysis.png)

The above figure shows the visualization of parameter sensitivity analysis results. (a) Single-parameter sensitivity analysis shows the effects of recovery rate (γ), symptomatic probability (p), and relative asymptomatic infectivity (ε) on the basic reproduction number (R₀). (b) Distribution and confidence intervals of R₀ values under different calculation methods. (c) Three-dimensional sensitivity surface of recovery rate (γ) and symptomatic probability (p) relative to R₀. (d-f) Scatter plots and trend lines for γ, ε, and p, respectively.

The above figure intuitively validates quantitative results through six subplots. Subplot (a) shows that the γ curve is steepest and exhibits non-linear decline, indicating high sensitivity of R₀ to γ. In contrast, p and ε show linear moderate growth with lower sensitivity, while the p curve exhibits steeper linear gradient than the ε curve. Subplot (b) presents violin plots showing that three methods produce similar R₀ means but reveal different distribution patterns and confidence intervals. Subplot (c) shows a three-dimensional surface where R₀ values are maximum when γ is low and p is high, rapidly declining as γ increases or p decreases, providing a multi-dimensional perspective on parameter interactions. Subplots (d-f) combine single-parameter sensitivity analysis with multi-parameter Monte Carlo sampling results: curves represent systematic R₀ calculations using Petri Net next-generation matrix method at 50 parameter value points, while scatter points represent experimental data from 1000 multi-parameter Monte Carlo sampling runs with randomly assigned parameter values. These scatter plots confirm the quantitative analysis conclusions that γ has strong negative correlation with R₀, while ε and p have moderate positive correlation with R₀.

In conclusion, recovery rate γ represents the core parameter for controlling Chikungunya transmission. Therefore, intervention strategies aimed at accelerating patient recovery will produce the most significant epidemic control effects.

## Data Preview

### Data Extraction Methods

This study uses the `data/data_extract.py` script to automatically extract daily incidence data from official epidemiological charts published by Foshan health authorities through OpenCV image processing technology. This method employs HSV color space detection to accurately identify blue and red bar charts in the images and automatically calculate their geometric properties.

### Extracted Data Overview

The epidemic data extracted from the `foshan_chikungunya_outbreak_epidemiological_data_2025.csv` file contains the following key information:

<div align="center">

| Date | Onset Cases | Reported Cases | Data Type | Notes |
|------|-------------|----------------|-----------|-------|
| 2025-06-16 | 1 | - | Onset Date | First case |
| 2025-06-28 | 1 | - | Onset Date | - |
| 2025-06-30 | 3 | - | Onset Date | - |
| 2025-07-01 | 1 | - | Onset Date | - |
| 2025-07-02 | 8 | - | Onset Date | - |
| 2025-07-08 | 54 | 5 | Onset/Report | Epidemic notification |
| 2025-07-09 | 100 | 0 | Onset/Report | Enhanced case management and vector control measures |
| 2025-07-14 | 169 | 190 | Onset/Report | Medical institutions initiated nucleic acid testing |
| 2025-07-18 | 383 | 360 | Onset/Report | Three-day centralized mosquito breeding site cleanup |
| 2025-07-19 | 388 | 674 | Onset/Report | Reported cases reached peak |
| 2025-07-21 | 432 | 371 | Onset/Report | Onset cases reached peak |
| 2025-07-26 | 95 | 387 | Onset/Report | Onset numbers began to decline |

</div>

### Open-Source Data Files

- **Main Data File**: `foshan_chikungunya_outbreak_epidemiological_data_2025.csv`
- **Extraction Script**: `data/data_extract.py`
- **Data Format**: CSV format, containing fields such as date, case numbers, data type

## Results Reproduction

- Python 3.8+
- Required packages listed in `requirements.txt`:

```bash
pip install -r requirements.txt
```

## SEICR Model

### ODE Model

We adopt the SEICR model with states $(S,E,I_s,I_a,C,R)$. The force of infection is $\lambda = \beta S (I_s + \epsilon I_a)/N$, where $\epsilon$ is the relative infectiousness of asymptomatic infections. The dynamics follow the balance equations described in our documentation, with demographic inflow/outflow and transition rates $\sigma,\gamma,\omega,\mu, p, q$.

The system of ordinary differential equations governing the SEICR model dynamics is given by:

$$
\begin{align}
\frac{dS}{dt} &= \Lambda - \lambda S - \mu S \\
\frac{dE}{dt} &= \lambda S - \sigma E - \mu E \\
\frac{dI_s}{dt} &= p\sigma E - \gamma I_s - \mu I_s \\
\frac{dI_a}{dt} &= (1-p)\sigma E - \gamma I_a - \mu I_a \\
\frac{dC}{dt} &= q\gamma I_s - \omega C - \mu C \\
\frac{dR}{dt} &= (1-q)\gamma I_s + \gamma I_a + \omega C - \mu R
\end{align}
$$

where the force of infection is defined as:

$$
\lambda = \frac{\beta S (I_s + \epsilon I_a)}{N}
$$

The total population satisfies:

$$
N = S + E + I_s + I_a + C + R
$$

Model parameters are defined as follows: $\Lambda$ is the population birth/immigration rate, $\mu$ is the natural mortality rate, $\beta$ is the transmission rate, $\epsilon$ is the relative infectiousness of asymptomatic infections, $\sigma$ is the exposure transition rate, $p$ is the proportion of exposed individuals who become symptomatic, $\gamma$ is the acute recovery rate, $q$ is the proportion of symptomatic individuals who progress to chronic stage, $\omega$ is the chronic recovery rate.

### Petri Net Model Construction

Based on the aforementioned advantages of Petri Nets in epidemiological modeling, we construct a Petri Net implementation comparable to the SEICR ODE model. This PN model contains places corresponding to SEICR compartments and transitions representing infection, progression, recovery, mortality, and birth processes. Variable-rate transitions implement the same rate functions as the ODE terms, ensuring structural similarity.

The Petri Net model can be formally defined as a six-tuple $PN = (P, T, F, W, M_0, \Lambda)$, with symbol explanations as follows:

<div align="center">

| Symbol | Definition | Description |
|--------|------------|-------------|
| $P$ | $\{S, E, I_s, I_a, C, R\}$ | Set of places, corresponding to the six compartments of the SEICR model |
| $T$ | $\{t_{inf}, t_{exp\_s}, t_{exp\_a}, t_{rec\_s}, t_{rec\_a}, t_{chr}, t_{rec\_c}, t_{birth}, t_{death}\}$ | Set of transitions, representing all state transition processes |
| $F$ | $\subseteq (P \times T) \cup (T \times P)$ | Flow relation, defining connections between places and transitions |
| $W$ | $F \rightarrow \mathbb{R}^+$ | Weight function, defining arc weights |
| $M_0$ | $P \rightarrow \mathbb{N}$ | Initial marking, representing initial population in each compartment |
| $\Lambda$ | $T \rightarrow \mathbb{R}^+$ | Transition rate function, defining firing rates for each transition |

</div>

The rate functions for each transition remain consistent with the corresponding terms in the ODE model:

<div align="center">

| Transition | Rate Function | Description |
|------------|---------------|-------------|
| $t_{inf}$ | $\lambda_{inf}(t) = \frac{\beta(t) S(t) (I_s(t) + \epsilon I_a(t))}{N(t)}$ | Infection transition |
| $t_{exp\_s}$ | $\lambda_{exp\_s}(t) = p\sigma E(t)$ | Exposed to symptomatic infectious transition |
| $t_{exp\_a}$ | $\lambda_{exp\_a}(t) = (1-p)\sigma E(t)$ | Exposed to asymptomatic infectious transition |
| $t_{rec\_s}$ | $\lambda_{rec\_s}(t) = (1-q)\gamma I_s(t)$ | Symptomatic infectious to chronic transition |
| $t_{chr}$ | $\lambda_{chr}(t) = q\gamma I_s(t)$ | Symptomatic infectious direct recovery transition |
| $t_{rec\_a}$ | $\lambda_{rec\_a}(t) = \gamma I_a(t)$ | Asymptomatic infectious recovery transition |
| $t_{rec\_c}$ | $\lambda_{rec\_c}(t) = \omega C(t)$ | Chronic patient recovery transition |
| $t_{birth}$ | $\lambda_{birth}(t) = \Lambda$ | Population birth transition |
| $t_{death}$ | $\lambda_{death}(t) = \mu N(t)$ | Natural death transition |

</div>

$$
\begin{align}
\lambda_{inf}(t) &= \frac{\beta(t) S(t) (I_s(t) + \epsilon I_a(t))}{N(t)} \\
\lambda_{exp\_s}(t) &= p\sigma E(t) \\
\lambda_{exp\_a}(t) &= (1-p)\sigma E(t) \\
\lambda_{rec\_s}(t) &= (1-q)\gamma I_s(t) \\
\lambda_{chr}(t) &= q\gamma I_s(t) \\
\lambda_{rec\_a}(t) &= \gamma I_a(t) \\
\lambda_{rec\_c}(t) &= \omega C(t) \\
\lambda_{birth}(t) &= \Lambda \\
\lambda_{death}(t) &= \mu N(t)
\end{align}
$$

### Basic Reproduction Number (R₀) Calculation

The basic reproduction number (R₀) is a key indicator for understanding transmission potential. We use two methods for cross-validation calculation.

#### Definition Method

This method is based on the definition R₀ = β × (effective infectious period).

$$
R_0 = \beta \times T_{eff}
$$

where the effective infectious period $T_{eff}$ is:

$$
T_{eff} = p \times \frac{1}{\gamma_s} + (1-p) \times \epsilon \times \frac{1}{\gamma_a}
$$

According to the definition method, only $I_s$ (symptomatic infectious individuals) and $I_a$ (asymptomatic infectious individuals) contribute to subsequent transmission in the model. The expected number of secondary cases generated by a typical initial individual can be weighted by pathway proportions and relative infectiousness. When natural mortality $\mu$ is present, the definition method yields:

$$
R_{\mathrm{def}}(\beta) = \beta\,\frac{\sigma}{\sigma+\mu}\left( \frac{p}{\gamma+\mu} + \frac{(1-p)\,\epsilon}{\gamma+\mu} \right)
$$

where the factor $\sigma/(\sigma+\mu)$ reflects the survival proportion from the latent period to the infectious period, and $1/(\gamma+\mu)$ represents the average infectious period for $I_s$ and $I_a$. For short outbreak durations or when $\mu\ll \sigma,\gamma$, equation above approximates to the commonly used expression:

$$
R_{\mathrm{def}}(\beta) \approx \beta\left( \frac{p}{\gamma} + \frac{(1-p)\,\epsilon}{\gamma} \right)
$$

#### Next-Generation Matrix Method

Based on the spectral radius of the next-generation matrix K = FV⁻¹, where F is the new infection matrix and V is the transition matrix.

##### ODE Next-Generation Matrix Method

For the next-generation matrix (NGM) method in ODEs, we first linearize the system at the disease-free equilibrium (DFE) where $(S^*,E^*,I_s^*,I_a^*,C^*,R^*) = (N,0,0,0,0,0)$. The disease-free equilibrium represents the state where no infection is present in the population, where the susceptible population equals the total population and all infectious compartments are empty. This equilibrium serves as the reference point for linearization and reproduction number calculation.

Let the infectious subsystem be $(E,I_s,I_a)$ and define the new infection matrix and transition matrix as:

$$
F = \begin{pmatrix}
0 & \beta & \beta\,\epsilon\\
0 & 0 & 0\\
0 & 0 & 0
\end{pmatrix},\qquad
V = \begin{pmatrix}
\sigma+\mu & 0 & 0\\
-p\,\sigma & \gamma+\mu & 0\\
-(1-p)\,\sigma & 0 & \gamma+\mu
\end{pmatrix}
$$

The next-generation matrix is then $K=FV^{-1}$, and the reproduction number is its spectral radius $R_{\mathrm{NGM}}=\rho(K)$.

##### Petri Net Next-Generation Matrix Method

In the Petri Net framework, the construction of the next-generation matrix method is based on the same infectious subsystem $(E,I_s,I_a)$. The disease-free equilibrium (DFE) in the Petri Net framework corresponds to the marking vector $m^*=(N,0,0,0,0,0)$, representing that the susceptible place contains the entire population while all infection-related places are empty.

For the infectious subsystem $(E,I_s,I_a)$, the Petri Net's Jacobian matrix $J$ is decomposed at the DFE as $J=F-V$, where:

$$
F^{\mathrm{PN}} = \begin{pmatrix}
0 & \beta & \beta\,\epsilon\\
0 & 0 & 0\\
0 & 0 & 0
\end{pmatrix},\qquad
V^{\mathrm{PN}} = \begin{pmatrix}
\sigma+\mu & 0 & 0\\
-p\,\sigma & \gamma+\mu & 0\\
-(1-p)\,\sigma & 0 & \gamma+\mu
\end{pmatrix}
$$

The construction of matrix $F^{\mathrm{PN}}$ is based on the rate function of the new infection transition in the Petri Net: the infection transition $t_{\mathrm{inf}}$ from susceptible $S$ to exposed $E$, with rate function $\lambda_{\mathrm{inf}}(t) = \beta(t) S(t) (I_s(t) + \epsilon I_a(t))/N(t)$. At the DFE where $S \approx N$, the non-zero elements of $F^{\mathrm{PN}}$ correspond to the partial derivatives of the infection transition with respect to $I_s$ and $I_a$.

Therefore, the next-generation matrix in the Petri Net framework is $K^{\mathrm{PN}} = (F^{\mathrm{PN}})(V^{\mathrm{PN}})^{-1}$, and the basic reproduction number is $R_{\mathrm{NGM}}^{\mathrm{PN}} = \rho(K^{\mathrm{PN}})$.

### 95% Confidence Interval Results

<div align="center">

| Phase | ODE (Definition) | ODE (NGM) | PN (NGM) | Transmission Reduction |
|-------|------------------|-----------|----------|----------------------|
| Initial | 14.67 [11.30, 20.82] | 14.67 [11.29, 20.81] | 13.90 [12.59, 15.30] | - |
| Phase 1 | 7.85 [6.07, 10.93] | 7.84 [6.07, 10.92] | 7.86 [7.07, 8.65] | 46.5% |
| Phase 2 | 7.59 [5.86, 10.45] | 7.59 [5.85, 10.45] | 7.56 [6.85, 8.35] | 48.3% |
| Phase 3 | 0.059 [0.046, 0.082] | 0.059 [0.046, 0.082] | 0.059 [0.053, 0.065] | 99.6% |

</div>

### Model Parameters

#### Demographic Parameters

<div align="center">

| Parameter | Symbol | Value | Unit | Meaning |
|-----------|--------|-------|------|---------|
| Total Population | $N$ | 9,615,400 | people | Foshan resident population (2023) |
| Birth Rate | $\Lambda$ | 317 | people/day | Human birth or net immigration rate |
| Death Rate | $\mu$ | 1.504×10⁻⁵ | day⁻¹ | Natural human mortality rate |

</div>

Lambda (daily birth/net immigration) and mu (daily mortality rate) calculations are based on Foshan City 2024 population statistics.

#### Disease Progression Parameters

<div align="center">

| Parameter | Symbol | Value | Unit | Meaning |
|-----------|--------|-------|------|---------|
| Incubation Rate | $\sigma$ | 0.2 | day⁻¹ | Average incubation period 5 days |
| Symptomatic Proportion | $p$ | 0.75 | dimensionless | 75% of infected individuals are symptomatic |
| Symptomatic Recovery Rate | $\gamma_s$ | 0.143 | day⁻¹ | Average symptomatic infectious period 7 days |
| Asymptomatic Recovery Rate | $\gamma_a$ | 0.143 | day⁻¹ | Average asymptomatic infectious period 7 days |
| Chronic Proportion | $q$ | 0.43 | dimensionless | 43% of symptomatic patients develop chronic stage |
| Chronic Recovery Rate | $\omega$ | 0.0026 | day⁻¹ | Average chronic period 385 days |

</div>

#### Transmission Parameters

<div align="center">

| Parameter | Symbol | Value | Unit | Meaning |
|-----------|--------|-------|------|---------|
| Effective Transmission Rate | $\beta$ | 0.9470 | day⁻¹ | Optimal value based on three-phase intervention fitting |
| Asymptomatic Relative Transmissibility | $\epsilon$ | 0.375 | dimensionless | Asymptomatic transmissibility is 37.5% of symptomatic |

</div>

The `epsilon` parameter (relative transmissibility of asymptomatic infections) was updated from 0.1 to 0.375, supported by scientific evidence and computational processes. This change affected the model's `beta` value and three-phase intervention fitting results, achieving R²=0.941. The calculation is based on biological data such as viral load differences and mosquito infection probabilities, and this value has been validated through other studies and three-phase fitting models.

## Contributing

We welcome contributions! Whether for model extensions and improvements, additional validation studies, or applications to other epidemics.

If you have questions about the article or markdown documentation or find errors, please raise them through [GitHub Issues](https://github.com/Hylouis233/ODE-Petri-Chikungunya/issues).

## Citation

If you use this work, please cite:

```bibtex
@article{liu2025dual,
  title={Dual-Model Framework for CHIKV Transmission Modeling: ODE and Petri Net Analysis of the 2025 Foshan Outbreak},
  author={Liu, Hong and Tian, Jingjing and Tang, Sheng and Yang, Zifeng},
  journal={arXiv preprint arXiv:XXXX.XXXXX},
  year={2025},
  note={Under review}
}
```

### Data Citation

If you only use the open-source data from this project, please cite:

```bibtex
@dataset{liu2025foshan,
  title={Foshan Chikungunya Outbreak Epidemiological Data 2025},
  author={Liu, Hong and Tian, Jingjing and Tang, Sheng and Yang, Zifeng},
  year={2025},
  publisher={Zenodo},
  doi={xx.xxxx/zenodo.XXXXXXX},
  url={https://doi.org/xx.xxxx/zenodo.XXXXXXX}
}
```
