# SDALS : Semi-supervised Discriminant Analysis under Label Shift

This repository contains research code for an ongoing project on **semi-supervised linear discriminant analysis (LDA) under label shift**.

## Overview

The project considers a classification setting where labeled observations are available from a source population, while only unlabeled observations are available from a target population with a different class distribution.

Under label shift, the class proportions may change between the source and target populations while the class-conditional feature distributions remain stable. The main goal of this project is to investigate how both labeled source data and unlabeled target data can be used to construct an LDA model better adapted to the target population.

The current implementation studies an estimation procedure based on the **Efficient Influence Function (EIF)** framework. Detailed methodological derivations are omitted because the work is still ongoing and unpublished.

## Baseline Methods

The proposed approach is mainly compared with two established label-shift adaptation methods.

### BBSE

**Black Box Shift Estimation (BBSE)** estimates changes in class proportions using the prediction behavior of a classifier trained on the source population.

In this implementation, BBSE uses the source confusion matrix together with the predicted class distribution on the unlabeled target data to estimate label-shift weights. These weights are then used to adjust the estimation procedure for the target population.

### RLLS

**Regularized Learning under Label Shift (RLLS)** also estimates the change in class proportions, but introduces regularization to improve numerical stability.

The implementation uses source and target prediction distributions together with a joint confusion-matrix estimate, and solves a regularized optimization problem to obtain class-specific adjustment weights.

Together, BBSE and RLLS provide natural baselines for evaluating the proposed method under label shift.

## Simulation Study

Simulation experiments are conducted to compare the proposed estimator with BBSE, RLLS, and a simple weighted reference method.

The experiments evaluate both:

* target-parameter estimation
* downstream target classification performance

Parameter estimation is summarized using metrics such as **Bias, SE, MSE, and RMSE**.

Classification performance is evaluated using:

* Accuracy
* Matthews Correlation Coefficient (MCC)
* TP / TN / FP / FN

Repeated simulation runs are parallelized to improve computational efficiency.

## Repository Structure

```text
function_all_02.13.r
    Core functions for LDA fitting, proposed estimation,
    BBSE, RLLS, and related components.

simulation_all_outlier_02.25.r
    Simulation experiments for parameter estimation
    and method comparison.

simulation_all_outlier_eval_ver.r
    Extended simulation code for downstream
    classification performance evaluation.
```

## Environment

The project is implemented in **R** and mainly uses:

```r
MASS
mvtnorm
snowfall
dplyr
tidyr
```

## Research Status

This repository contains **ongoing and unpublished research**.

The code is provided to document the implementation and experimental development of the project. Full theoretical derivations, detailed algorithms, simulation settings, and complete numerical results are intentionally omitted at this stage.

## Acknowledgement

This repository was developed with support from the 서울시립대학교 데이터 사이언스 플러스 차세대 융합인재 양성사업단 - http://dsplus.uos.ac.kr/
