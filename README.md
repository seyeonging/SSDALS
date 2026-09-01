# Semi-supervised Discriminant Analysis under Label Shift

This repository contains research code developed for an ongoing project on **semi-supervised linear discriminant analysis (LDA) under label shift**.

The project focuses on classification problems where labeled observations are available from a source population, while the target population consists mainly of unlabeled observations and may have a substantially different class distribution.

The goal of this work is to investigate how information from both populations can be incorporated into discriminant analysis so that the resulting classifier is better adapted to the target population.

> **Note:** This repository corresponds to ongoing and unpublished research.
> Detailed methodological derivations, theoretical results, and complete experimental results are intentionally omitted.

---

## Background

In many classification problems, a model is trained on data collected from one population and subsequently applied to another population.

A standard supervised learning approach implicitly assumes that the training and target distributions are sufficiently similar. In practice, however, the distribution observed at deployment may differ from the distribution used for model training.

This project considers **label shift**, where the class proportions change between the source and target populations while the class-conditional feature distributions are assumed to remain stable.

In notation,

* the **source population** provides labeled observations,
* the **target population** provides unlabeled observations,
* the marginal distribution of the class label may differ between the two populations.

Under such a setting, directly applying a classifier trained only on the source population may lead to suboptimal prediction for the target population.

At the same time, simply ignoring the unlabeled target observations discards potentially useful information about the population in which prediction will actually be performed.

This motivates a semi-supervised approach that combines labeled source observations with unlabeled target observations.

---

## Research Objective

The main objective of this project is to study **target-adapted Linear Discriminant Analysis under label shift**.

The current work investigates an estimation framework that uses an LDA model together with information from the unlabeled target population to estimate parameters relevant to the target distribution.

The proposed procedure is developed using an **Efficient Influence Function (EIF)-based framework**.

At a high level, the research pipeline consists of:

1. fitting an initial LDA model using the labeled source data,
2. extracting probabilistic predictions for source and target observations,
3. estimating quantities associated with the distribution shift,
4. incorporating information from both labeled and unlabeled observations,
5. estimating target-population model parameters,
6. evaluating both parameter recovery and downstream classification performance.

The detailed estimating equations and theoretical derivations are not included in this repository description because the project is currently under preparation.

---

## Baseline Methods

To evaluate the proposed approach, the current implementation includes two established approaches for adaptation under label shift:

* **BBSE**
* **RLLS**

These methods are used as primary baselines because they estimate changes in class proportions under label shift and provide a natural comparison with the proposed target-adaptation procedure.

### BBSE

**Black Box Shift Estimation (BBSE)** estimates the change in label distribution using the predictions of a classifier trained on the source population.

The basic idea is to treat the source classifier as a black-box predictor and use its prediction behavior to infer how the class proportions have changed in the target population.

In the current implementation, BBSE uses:

* predictions from the source-trained LDA classifier,
* a confusion-matrix-based estimate obtained from the labeled source observations,
* the predicted class distribution of the unlabeled target observations.

These quantities are combined to estimate the target class proportions and corresponding label-shift weights.

The estimated weights are then incorporated into the estimating procedure so that source observations are reweighted according to the estimated difference between source and target label distributions.

BBSE provides a particularly useful baseline because it directly addresses label shift without requiring target labels.

---

### RLLS

**Regularized Learning under Label Shift (RLLS)** is another label-shift adaptation method used as a baseline in this project.

Like BBSE, RLLS aims to estimate how the label distribution changes between the source and target populations. However, it introduces a **regularized estimation procedure** to improve stability when estimating the shift.

This is particularly relevant when the amount of labeled source data is limited or when the confusion-matrix-based estimation problem becomes numerically unstable.

In the current implementation, RLLS uses:

* the prediction distribution from the source population,
* the prediction distribution from the target population,
* a joint confusion-matrix estimate based on source observations,
* a regularized optimization problem for estimating the label-shift adjustment.

The resulting estimate is transformed into class-specific importance weights and incorporated into the corresponding estimating equation.

Including both BBSE and RLLS allows the proposed method to be compared against both a standard confusion-matrix-based label-shift estimator and a regularized alternative designed to improve estimation stability.

---

## Proposed Approach

The main method investigated in this repository is an **EIF-based semi-supervised estimation procedure for LDA under label shift**.

Rather than only estimating class-proportion weights and applying them to the labeled source observations, the proposed framework is designed to incorporate information from both:

* labeled source data, and
* unlabeled target data.

The implementation includes components for:

* fitting the LDA working classifier,
* calculating posterior class probabilities,
* estimating source-to-target label-shift quantities,
* constructing the necessary weighting terms,
* calculating augmentation-related quantities,
* evaluating the resulting estimating function,
* numerically solving for target-population parameters.

Only a high-level description is provided here because the exact construction constitutes part of the ongoing methodological work.

---

## Simulation Study

Simulation experiments are used to evaluate the statistical behavior of the proposed estimator and compare it with the baseline approaches.

The simulations generate separate source and target populations under controlled label-shift settings.

The source population contains labeled observations, whereas target labels are treated as unavailable during model estimation.

This setup allows the experiments to reproduce the main setting considered in the research:

```text
Labeled Source Data
        |
        v
   Initial LDA
        |
        +------------------+
        |                  |
        v                  v
Source Information   Unlabeled Target Data
        |                  |
        +--------+---------+
                 |
                 v
       Label-Shift Adaptation
                 |
                 v
     Target Parameter Estimation
                 |
                 v
     Classification Evaluation
```

Multiple simulation replicates are conducted to assess the stability of each method.

Parallel computation is used to make repeated simulation experiments computationally feasible.

---

## Methods Compared

The current simulation framework compares the following estimators:

| Method                   | Description                                                                         |
| ------------------------ | ----------------------------------------------------------------------------------- |
| **Proposed / Efficient** | EIF-based semi-supervised target adaptation procedure                               |
| **BBSE**                 | Confusion-matrix-based label-shift estimation                                       |
| **RLLS**                 | Regularized estimation of label-shift weights                                       |
| **Naive weighting**      | Simplified importance-weighted estimating procedure used as an additional reference |

For downstream classification experiments, source-based and oracle classifiers may also be included as reference points.

These reference classifiers are not treated as directly comparable adaptation methods. Instead, they provide useful lower/upper reference points for understanding the effect of adaptation.

---

## Evaluation

The experiments evaluate the methods from two perspectives.

### 1. Parameter Estimation

The first objective is to determine how accurately each approach estimates the parameters associated with the target population.

The current simulation pipeline summarizes repeated estimates using quantities including:

* empirical mean of the estimates,
* bias,
* standard error,
* mean squared error (MSE),
* root mean squared error (RMSE).

This allows the comparison to focus not only on final prediction accuracy but also on the statistical quality and stability of the parameter estimates.

---

### 2. Classification Performance

The second objective is to evaluate whether improvements in parameter estimation translate into better target-population classification.

An independent target test sample is generated for this purpose.

Classification performance is summarized using:

* **Accuracy**
* **Matthews Correlation Coefficient (MCC)**
* True Positives
* True Negatives
* False Positives
* False Negatives

MCC is considered together with Accuracy because the target class distribution can be highly unbalanced under strong label shift.

The evaluation pipeline therefore examines both overall prediction accuracy and the balance of classification performance across the two classes.

---

## Implementation

The project is implemented in **R**.

The main packages currently used include:

```r
MASS
mvtnorm
snowfall
dplyr
tidyr
```

`MASS` is used primarily for LDA and multivariate data generation, while `mvtnorm` is used for multivariate normal density calculations.

The simulation code also uses parallel processing to support repeated experiments.

---

## Repository Structure

```text
.
├── function_all_02.13.r
├── simulation_all_outlier_02.25.r
└── simulation_all_outlier_eval_ver.r
```

### `function_all_02.13.r`

Contains the main functions used throughout the experiments.

The current implementation includes functions related to:

* LDA model fitting,
* posterior probability calculation,
* score-function evaluation,
* label-shift weight estimation,
* EIF-based estimation components,
* BBSE,
* RLLS,
* additional comparison estimators.

---

### `simulation_all_outlier_02.25.r`

Contains the main simulation pipeline for parameter estimation.

The script includes:

* source and target data generation,
* repeated simulation runs,
* parameter optimization,
* comparison across methods,
* numerical error handling,
* screening of unstable estimates,
* parallel computation,
* summary of estimation performance.

---

### `simulation_all_outlier_eval_ver.r`

Extends the simulation framework to evaluate downstream classification performance.

In addition to parameter-estimation summaries, this script includes:

* generation of an independent target test set,
* classification using estimated parameters,
* Accuracy and MCC calculation,
* confusion-matrix-based performance summaries,
* comparison with reference source and oracle classifiers.

---

## Current Research Workflow

The overall research workflow can be summarized as:

```text
1. Define the label-shift setting
            |
            v
2. Fit a source LDA working model
            |
            v
3. Estimate target-distribution information
            |
            v
4. Apply each adaptation method
   ├── Proposed EIF-based approach
   ├── BBSE
   └── RLLS
            |
            v
5. Estimate target LDA parameters
            |
            v
6. Repeat simulation experiments
            |
            v
7. Evaluate estimation error
            |
            v
8. Evaluate target classification
```

The code in this repository represents the implementation and empirical evaluation stages of this workflow.

---

## Research Status

This repository is associated with **ongoing research**.

The current codebase is primarily intended to document the implementation and experimental development of the project.

Some aspects of the code, experiment design, and repository structure may therefore change as the study progresses.

In particular, the following are intentionally not provided in detail at this stage:

* full theoretical derivations,
* complete estimating equations,
* proofs,
* detailed algorithmic derivations,
* full simulation settings,
* complete numerical results,
* manuscript-specific implementation details.

These components will be documented separately after the corresponding research work is finalized.

---

## Disclaimer

This repository should be considered **research-stage code rather than a production software package**.

The implementation is being actively developed for methodological validation and simulation experiments. Interfaces, function names, and experiment settings may change in future versions.
