# Iemi2017_JNeurosci
Matlab script for simulating the temporal contamination of a wavelet-transformed 10~Hz signal into pre-signal period

Time-frequency analysis was carried out using a wavelet transform (Morlet wavelets, frequency range: 1-30 Hz, number of cycles increasing linearly from 1 to 12). Thus, a wavelet at 10~Hz was 4.4 cycles long and had a temporal resolution σ_t of 0.14 s and a spectral resolution σ_f of 4.53~Hz. Frequencies and time points were sampled every 2 Hz and 40 ms, respectively. Since we were primarily interested in the prestimulus time range, no baseline correction was applied.

Since wavelet analysis is computed by convolution of the data with a function that is extended in time, it is conceivable that prestimulus effects close to stimulus onset are actually affected by post-stimulus data. To determine the extent of this contamination, we applied a time-frequency transform using the same settings as in the main analysis to synthetic data, i.e. a sinusoidal oscillation of exactly 10 Hz (Figure 3 in Iemi et al, 2016). The extent of temporal contamination caused by a wavelet is determined by the wavelet's temporal resolution σ_t, which is defined as twice the standard deviation of the Gaussian envelope (Tallon-Baudry et al, 1996). This simulation illustrates that prestimulus time points are indeed contaminated by post-stimulus data points. However, the magnitude of this contamination is virtually null at time points earlier than onset - σ_t. Thus, we consider effects as truly "prestimulus" only if they occur prior to this limit, which is indicated by a red line in Figures 4, 5, 6 and 7 in Iemi et al (2017).

To run this script, eeglab is required.

Tallon-Baudry, C., Bertrand, O., Delpuech, C., & Pernier, J. (1996). Stimulus specificity of phase-locked and non-phase-locked 40 Hz visual responses in human. The Journal of Neuroscience, 16(13), 4240–9. doi:10.1016/j.neuropsychologia.2011.02.038
