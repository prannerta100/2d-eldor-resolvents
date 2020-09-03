# 2d-eldor-resolvents
### Paper 1 :2D ELDOR improved simulation package, very slow motional, works up to the rigid limit of protein motions

*P. Gupta et al, Microsecond dynamics in proteins by two-dimensional ESR: Predictions, J. Chem. Phys. 152, 214112 (2020)*

### Abstract:

Two-dimensional electron–electron double resonance (2D-ELDOR) provides extensive insight into molecular motions. Recent developments permitting experiments at higher frequencies (95 GHz) provide molecular orientational resolution, enabling a clearer description of the nature of the motions. In this work, simulations are provided for the example of domain motions within proteins that are themselves slowly tumbling in solution. These show the nature of the exchange cross-peaks that are predicted to develop in real time from such domain motions. However, we find that the existing theoretical methods for computing 2D-ELDOR experiments over a wide motional range begin to fail seriously when applied to very slow motions characteristic of proteins in solution. One reason is the failure to obtain accurate eigenvectors and eigenvalues of the complex symmetric stochastic Liouville matrices describing the experiment when computed by the efficient Lanczos algorithm in the range of very slow motion. Another, perhaps more serious, issue is that these matrices are “non-normal,” such that for the very slow motional range even rigorous diagonalization algorithms do not yield the correct eigenvalues and eigenvectors. We have employed algorithms that overcome both these issues and lead to valid 2D-ELDOR predictions even for motions approaching the rigid limit. They are utilized to describe the development of cross-peaks in 2D-ELDOR at 95 GHz for a particular case of domain motion.

<p align="center">
  <img src='https://aip.scitation.org/na101/home/literatum/publisher/aip/journals/content/jcp/2020/jcp.2020.152.issue-21/5.0008094/20200604/images/medium/5.0008094.figures.online.f1.jpg'/>
</p>

### Dependencies:
1. MATLAB
2. `expokit` package, from https://www.maths.uq.edu.au/expokit/matlab/
3. `nlspmc` package for generating Stochastic Liouville Matrices, from https://www.acert.cornell.edu/index_files/acert_ftp_links.php

### Running the simulations:
run the Paper1/Exch/plot_all_exchange_spectra.m MATLAB script to generate exchange 2D-ELDOR spectra


Send feedback to prannerta100@gmail.com (Pranav Gupta).
