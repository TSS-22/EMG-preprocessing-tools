# EMG-denoising
Functions to denoise High Density EMG signals.

## Table of Contents
* [EMG Contamination](#EMG-Contamination)
* [Movement artifact, PLI, WGI contamination](#Movement-artifact,-PLI,-WGI-contamination)
* [Baseline wander filter](#Baseline-wander-filter)
* [Technologies](#Technologies)
* [Setup](#Setup)
* [How to cite](#How-to-cite)
* [Sources](#Sources)
* [Disclaimer](#Disclaimer)
* [Contact](#Contact)

## EMG Contamination
This ECG-artifact filter is based on the paper Mak, J. N., Hu, Y., & Luk, K. D. (2010). An automated ECG-artifact removal method for trunk muscle surface EMG recordings. Medical engineering & physics, 32(8), 840-848. DOI: 10.1016/j.medengphy.2010.05.007 and on "Rectification and non-linear pre-processing of EMG signals for cortico-muscular analysis." from Myers et al. DOI: 10.1016/s0165-0270(03)00004-9 for the pre-processing of the HD EMG.

The following function remove ECG contamination from High Density EMG signals. This is the matlab implementation of the work described in "Mak, J. N., Hu, Y., & Luk, K. D. (2010). An automated ECG-artifact removal method for trunk muscle surface EMG recordings. Medical engineering & physics, 32(8), 840-848." With the following changes:
- A moving average filter has been added to the pre processing step in order to smooth out potential artifacts on the ECG ICA component that could appear in field conditions. Those artifact would be detected as potential ECG spike and would teherefor make the component fail the RR interval check.
- To improve reliability of the RR interval check, it had been added a 10% margin of error, meaning at least 90% of the peak of the tested component should meet the RR interval check requirement.

### Example
Before:
![Before signals](./Miscellaneous/img/rawSig.png)

After:
![After signals](./Miscellaneous/img/ecgFilt.png)

### Principle
- Hilbert transform
- Median filter (order 50)
- Moving average filter (5% of the acquisiton frequency) /!\ this is an empirical addon in order to make the filter work in field condition /!\
- FastICA to take out the different signal sources
- Detection of the ECG sources
- Separation of the EMG and ECG sources and reconstruction of the signals

## Movement artifact, PLI, WGI contamination
This filter remove movement artifact, power line interference and white gaussian noise. It is based on the paper "Al Harrach, M., Boudaoud, S., Hassan, M., Ayachi, F. S., Gamet, D., Grosset, J. F., & Marin, F. (2017). Denoising of HD-sEMG signals using canonical correlation analysis. Medical & biological engineering & computing, 55(3), 375-388."

IMPORTANT: The algorithm needs a 0.5 seconds of signal without any muscle activity. So either add it a posteriori or be carefull to have enough room in your signal. The selection and the PNR algorithms depends on it.

/!\ WARNING /!\\

If your HD EMG signals are contaminated by ECG artifacts, it is vividly recommended to filter them out first before using this filter. Indeed the chance of successful use dramatically lower if there is ECG contamination of the signals.

### Example
Before:
![Before signals](./Miscellaneous/img/ecgFilt.png)

After:
![After signals](./Miscellaneous/img/ccaFilt.png)

### Principle
- Canonical Correlation Analysis to extract linear combination of the HD EMG signals (CCA components)
- Selection of the noisy components to filter via the use of an intensity ratio
- Correlation analysis between filtered and unfiltered signals to keep the maximum EMG information
- Selective CCA step, where filtered and unfiltered channels PNR is compared. The one with highest PNR is choosen to be kept for the final filtered EMG matrix. This is done in order to prevent noise contamination of the high PNR channel due to signals reconstruction

## Baseline wander fluctuation
This filter remove baseline wander also called baseline fluctuation. It is based on the paper "Fasano, A., & Villani, V. (2014). Baseline wander removal for bioelectrical signals by quadratic variation reduction. Signal Processing, 99, 48-57." for the filter and is using sparsity measure from "Hurley, N., & Rickard, S. (2009). Comparing measures of sparsity. IEEE Transactions on Information Theory, 55(10), 4723-4741."
Two part of the algorithm are personal implementation and have not been peer reviewed:
- This algorithm is using an inversion of a very large matrix, therefore mathematical optimization. I will remind that the optimization done here is of my own doing as the optimization process is not detailed by the authors.
- A gradient descent assessing the signal's sparsity is used to determine to optimal lambda from eq (16) of Fasano et al., 2014, as the process of finding the optimal lambda is not detailed or explicitly stated by the author. This way is therefore of my own doing

### Example
Before:
![Before signals](./Miscellaneous/img/bwFilt1.png)

After:
![After signals](./Miscellaneous/img/bwFiltafter2.png)

### Principle
The filter works by minimizing the variability of the signal by minimizing it quadratic variation, a measure of the signal variability. Baseline wander is estimated solving a constrained convex optimization problem where quadratic variation enters as a constraint.

Optimization:
This filter uses a matrix inversion of a matrix of the size N\*N, thus being extremely costly memory and computationally wise. Because the matrix is a tridiagonal, symmetric, positive-definite system it is possible to “linearize” the calculus. To do so the equation (16) have been developed, and the linear system solver “\\” has been used to “get rid of” the standard inverse calculation. This is possible due to the aforementioned specific properties of the system. More details and explanation can be found in the books from Golub matrix computation. Nonetheless, clear optimization process not being detailed in Fasano et al, 2014, the optimization process implemented here might not be the one used by the author of the article.

Gradient descent:
The process to find the optimal lambda for the filter is not detailed in Fasano et al., 2014. Therefore, the process implemented here is of my own design. The logic behind it is, we are facing a signal away from the baseline it should be on, which is 0. The filter is trying to get the signal back on this baseline. As the signal get filtered closer and closer to the baseline, the number of 0 value increase.
We therefore can formulate the problem as a sparsity problem. The gradient descent work on the following logic: 
if sparsity of sFilt(n-1) - sparsity of sFilt(n) < 0 then stop --> we found the optimal lambda.
The formula used have been chosen from Hurley, N., & Rickard, S. (2009). Comparing measures of sparsity. IEEE Transactions on Information Theory, 55(10), 4723-4741, and selected for their behavior that seemed to get close to the optimal lambda found experimentally.

## Technologies
The project has been developped on MATLAB 2018a 9.4.0.813654, on windows 10 64 bits.

## Setup
The filters are matlab function, to run them, just download them then add them to your matlab path.

## How to cite
Robinault Lucien (2021). EMG preprocessing tools V1.0.0 (https://github.com/TSS-22/EMG-preprocessing-tools), GitHub. DOI: 10.5281/zenodo.4441389

## Sources
- Al Harrach, M., Boudaoud, S., Hassan, M., Ayachi, F. S., Gamet, D., Grosset, J. F., & Marin, F. (2017). Denoising of HD-sEMG signals using canonical correlation analysis. Medical & biological engineering & computing, 55(3), 375-388. DOI: 10.1007/s11517-016-1521-x
- Mak, J. N., Hu, Y., & Luk, K. D. (2010). An automated ECG-artifact removal method for trunk muscle surface EMG recordings. Medical engineering & physics, 32(8), 840-848. DOI: 10.1016/j.medengphy.2010.05.007
- Myers, L. J., Lowery, M., O'malley, M., Vaughan, C. L., Heneghan, C., Gibson, A. S. C., ... & Sreenivasan, R. (2003). Rectification and non-linear pre-processing of EMG signals for cortico-muscular analysis. Journal of neuroscience methods, 124(2), 157-165. DOI: 10.1016/s0165-0270(03)00004-9
- Fasano, A., & Villani, V. (2014). Baseline wander removal for bioelectrical signals by quadratic variation reduction. Signal Processing, 99, 48-57.
- Hurley, N., & Rickard, S. (2009). Comparing measures of sparsity. IEEE Transactions on Information Theory, 55(10), 4723-4741.

## Disclaimer
The information and tools contained in this repository are provided in good faith and no warranty, representation, statement or undertaking is given regarding any information or tool connected with this repository and any warranty, representation, statement or undertaking whatsoever that may be expressed or implied by statute, custom or otherwise is hereby expressly excluded.
The use of the tools in this repository and any information in this repository is entirely at the risk of the user.
Under no other circumstances the author should be liable for any costs, losses, expenses or damages (whether direct or indirect, consequential, special, economic or financial including any loss of profits) whatsoever that may be incurred through the use of any information or tools contained in this repository. This repository may contain inaccurate information. Which the author is under no responsibility to update or correct any such information or to even maintain this repository. Which the author reserves its right to change any information or any part of this repository without notice.

## Contact
lucien.robinault@protonmail.com
