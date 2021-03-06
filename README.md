
Inference for time-varying signals using locally stationary processes
=====================================================================

Code to reproduce the simulations presented in the research paper

> Rachele Anderson, Maria Sandsten, "Inference for time-varying signals using Locally Stationary Processes", Journal of Computational and Applied Mathematics, Volume 347, Pages 24-35, 2019.

available online at: https://doi.org/10.1016/j.cam.2018.07.046

Keywords
--------
Locally Stationary Process, Time-varying signals, Time-series modelling, Statistical inference, Covariance estimation, EEG signals

Highlights
----------
- A novel inference method for Locally Stationary Processes is proposed.
- The algorithm divides the inference problem into two simpler sub-problems. 
- Convergence, accuracy, robustness improve respect to traditional approaches.

Abstract of the paper
----------------------
Locally Stationary Processes (LSPs) in Silverman’s sense, deﬁned by the modulation in time of a stationary covariance function, are valuable in stochastic modelling of time-varying signals. However, for practical applications, methods to conduct reliable parameter inference from measured data are required. In this paper, we address the lack of suitable methods for estimating the parameters of the LSP model, by proposing a novel inference method. The proposed method is based on the separation of the two factors defining the LSP covariance function, in order to take advantage of their individual structure and divide the inference problem into two simpler sub-problems. The method’s performance is tested in a simulation study and compared with traditional sample covariance based estimation. An illustrative example of parameter estimation from EEG data, measured during a memory encoding task, is provided.

How to
------
The code is written in Matlab (R2017b). Optimization Toolbox and Curve Fitting toolbox are required for running the code.

Run main.m to reproduce the simulation study presented in the paper.

The code is divided in sections:

1. Set parameters to simulate the data, simulation and optimization settings
2. Plot true covariance, model function and example of realizations from this parametric setting
3. Inference: For different number of realizations (num_real_vec), repeat num_sim simulations where the parameters and the covariance matrix are estimated with HATS and SCM
4-5-6. Compute distance from true covariance, confidence intervals for estimates, and plot figures 



