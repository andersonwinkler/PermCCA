# PermCCA

This repository contains a reference implementation for statistically valid permutation tests with canonical correlation analysis.

The paper reference is:
* Winkler AM, Renaud O, Smith SM, Nichols TE. [**Permutation Inference for Canonical Correlation Analysis**](https://doi.org/10.1016/j.neuroimage.2020.117065). *NeuroImage*. 2020; 117065 (Open Access, in press).

The implementation of both Huh-Jhun and Theil are slightly more involved in the code than in the paper. Both implementations are expected to lead to identical results. Here in the repository, the implementation of Huh-Jhun uses null(Z'), which appears to be numerically more stable in Matlab/Octave. For Theil, there are more options for selecting observations to be dropped (via selection matrix, numerical indices, logical indices, or random selection). These can be helpful in practice. In the paper only the generalised case of using a selection matrix is shown.
