# DSMGA2
Individual Study by Franky and Sammy

## Branch

### master
- mutual information

### simple
Replace mutual information with following:
```
p00 * p11 + p01 * p10
```

### kai
Replace mutual information with chi-square test.

### optBM
Prevent rare populations from disappearing during backmixing.

## Experiment

- kai_100_3: kai, ell-100, step-3
- master_100_3: mi, ell-100, step-3
- optBM_100_3: optBM, ell-100, step-3
- optBM_100_3_fixed_1: optBM, ell-100, step-3, threshold set to 1
- optBM_1_100_3: optBM, threshold set to one standard deviation, ell-100, step-3
- simple_100_3: simple, ell-100, step-3
