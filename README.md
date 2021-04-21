# Granger Causality (Wiener-Granger Causality)

## Definition:

If the autoregressive prediction of the first time series at present time could be improved by including the past information of the second time series, we say that the second time series has a causal influence on the first.

In the frequency domain, the total spectral power of the effect variable is decomposed into its intrinsic power and the causal contribution from Y and the ratio of the total power to the intrinsic power indicates the presence of causal influence.

Total power = intrinsic power + causal power

## History:

**1956**, Norbert Wiener proposed the notion that one variable or time series could be called causal to another if the ability to predict the second variable is improved by incorporating information about the first.

_Wiener, N., 1956. The theory of prediction. In: Becknbach. E. (Ed.), Modern Mathematics for Engineers. McGraw-Hill, New York._

**1969**, the econometrician Clive Granger implement this idea in the context of linear autoregressive models of stochastic processes.

_Granger, C. W. (1969). Investigating causal relations by econometric models and cross-spectral methods.** Econometrica: Journal of the Econometric Society**, 424-438._

**1982/4**, John Geweke developed a spectral decomposition of WGC, where the total time-domain WGC is equal to the sum of spectral WGC components over all frequencies from zero to the Nyquist frequency.

_Geweke, J. (1982). Measurement of linear dependence and feedback between multiple time series. **Journal of the American statistical association**, 77(378), 304-313._

_Geweke, J. (1984). Measures of conditional linear dependence and feedback between time series. **Journal of the American Statistical Association**, 79(388), 907-915._

Autoregressive modeling, the basis of the current parametric Granger causality techniques, has proven effective for data modeled by low-order AR processes. However, AR methods sometimes fail to capture complex spectral features in data that require higher order AR models. Additionally, the proper determination of model order remains a concern, although this concern may be mitigated by the recently proposed Bayesian framework.

**2008**, the author proposed **a nonparametric approach** based on widely used Fourier and wavelet transforms to estimate both pairwise and conditional measures of Granger causality, eliminating the need of explicit autoregressive data modeling.

_Dhamala, M., Rangarajan, G., & Ding, M. (2008). Analyzing information flow in brain networks with nonparametric Granger causality. **Neuroimage**, 41(2), 354-362._

_Dhamala, M., Liang, H., Bressler, S. L., & Ding, M. (2018). Granger-Geweke causality: Estimation and interpretation. **NeuroImage**, 175, 460-463. **[MATLAB CODE]**_

**2008**, a Matlab toolbox was developed for granger causality at time domain.

_Cui, J., Xu, L., Bressler, S. L., Ding, M., & Liang, H. (2008). BSMART: a Matlab/C toolbox for analysis of multichannel neural time series. **Neural Networks**, 21(8), 1094-1104. [function: armorf]_

**2010**, another MATLAB toolbox for GC analysis was released. [GCCA]

_Seth, A. K. (2010). A MATLAB toolbox for Granger causal connectivity analysis. **Journal of neuroscience methods**, 186(2), 262-273._

**2014**, the 2010 version toolbox was updated to MVGC.

_Barnett, L., & Seth, A. K. (2014). The MVGC multivariate Granger causality toolbox: a new approach to Granger-causal inference. **Journal of neuroscience methods**, 223, 50-68._

**Several reviews about GC:**

_Bressler, S. L., & Seth, A. K. (2011). Wiener–Granger causality: a well established methodology. **Neuroimage**, 58(2), 323-329._

_Friston, K., Moran, R., & Seth, A. K. (2013). Analysing connectivity with Granger causality and dynamic causal modelling. **Current opinion in neurobiology**, 23(2), 172-178._

_Seth, A. K., Barrett, A. B., & Barnett, L. (2015). Granger causality analysis in neuroscience and neuroimaging.** Journal of Neuroscience**, 35(8), 3293-3297._

## Analytic methods:

In general, sampling rates between 250 Hz and 1000 Hz will be appropriate.

### (1) Preprocessing:

Subtract the ERP and ensemble mean; detrend and z-score the data

Check the data, make sure it is stationary. KPSS test or ADF test

### Parametric analysis:

### (2) Model order

Use the ‘_armorf _’ function in BSMART toolbox to estimate the best order.

Use BIC (recommended, ‘LWR’ or ‘OLS’ mode) and AIC to determine the suitable order.

_Models with a small order have fewer parameters to estimate, which means that the estimation of those parameters (e.g., via the Matlab function armorf.m) will be more robust to noise. On the other hand, models with a small order are insensitive to longer time lags and thus may fail to detect true interactions that have a long temporal lag._

_Models with a larger order are sensitive to longer time lags and will allow you to extract low-frequency interactions but have more parameters to estimate and thus require more data. Thus, with larger model orders, more trials and longer time segments are helpful._

**Attentions: how to choose the model order?**

**_You should treat the “ optimal ” order as a statistically guided recommendation rather than an absolute rule. When you select a model order, you should apply that order to all conditions, time segments, electrode pairs, and subjects._**

### (3) Main procedure

Confirm the order is the best by executing these three tests: whiteness, consistency and stability tests.

**Time domain:**

Use the ‘armorf’ function to estimate the errors and compute the GC using ‘ln’ function.

Time segments ranging from as little as 50 ms to 100ms to 300-400ms up to 1s or 2 s (P. 377, Mike 2014).

**Frequency domain:**

Use the ‘pwcausal’ function to estimate the GC values.


### Nonparametric analysis: 
use ´compute_allnpCGCvar3´ function (it has alreday including all the following procedures)

### (2) Compute Spectral matrix

Use mutitaper analysis (‘sig2mTspect_nv’ function), frequency domain

Wavelet transform, time-frequency domain

### (3) Spectral matrix factorization------KEY

Use Wilson’s algorithm (‘wilson_sf’ function)

### (4) Compute the GC values

Use the ‘hz2cgcAll’ function to estimate the GC values.

