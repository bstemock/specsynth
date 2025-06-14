import numpy as np
import code.specsynth as ss
import matplotlib.pyplot as plt

# Basic setup
ATOM = ss.get_atomic("data/atoms.dat")      # atomic data used by the Absorber class
CON = ss.get_constants("data/const.dek")    # constants used by the Absorber class
HIRES = ss.Instrument(CON, R=45000, presel=3.0, rdnoise=3.0, slit=1.0, n=3.0, resfac=3.0)   # Instrument
MgII = ss.Absorber(CON, ATOM, trans="MgII2796", vels=[-100.0, 100.0], Inst=HIRES, seed=42,
                   snr=70, zabs=1.000000, v=[-22.2, 25.2, 52.0], logN=[12.1, 13.7, 11.7], b=[5.0, 7.2, 2.9])
print(type(MgII.atom), MgII.atom)

# Visualize absorber at different steps in the process
fig, axs = plt.subplots(1, 3, sharey="all")
fig.suptitle("Basic Spectrum-Building Functionality")

for i in range(3):
    axs[i].set_xlabel(r"Wavelength ($\mathrm{\AA}$)")
    axs[i].plot(MgII.waves, np.ones(MgII.waves.shape[0]), c="b", lw=0.5, ls="--")
    axs[i].plot(MgII.waves, np.zeros(MgII.waves.shape[0]), c="b", lw=0.5, ls="--")

axs[0].set_title("Pre-Convolution (Oversampled)")
axs[0].set_ylabel("Normalized Flux")
axs[0].step(MgII.waves_os, MgII.f_norm_preconv, c="k", lw=0.5, where="mid")

axs[1].set_title("Post-Convolution, Noiseless")
axs[1].step(MgII.waves, MgII.f_norm_noiseless, c="k", lw=0.5, where="mid")

axs[2].set_title("Noisy Spectrum")
axs[2].step(MgII.waves, MgII.f_norm, c="k", lw=0.5, where="mid", label="Flux")
axs[2].step(MgII.waves, MgII.I_sig, c="g", lw=0.5, where="mid", label="Uncertainty")
plt.show()
plt.close(fig)

# Line detection functionality
MgII_ = ss.Absorber(CON, ATOM, trans="MgII2803", vels=[-100.0, 100.0], Inst=HIRES, seed=42,
                    snr=70, zabs=1.000000, v=[25.2, 52.0], logN=[13.7, 11.7], b=[7.2, 2.9])

# Check for 5 sigma detections of dominant transition
ew_spec2796, ew_sig2796 = ss.get_ew_spec(CON, ATOM, HIRES, "MgII2796", 1.000000,
                                         MgII.vels, MgII.waves, MgII.f_norm, MgII.I_sig)
abs_vels2796 = ss.get_abs_regions(MgII.vels, ew_spec2796, ew_sig2796, sigma_threshold=5.0)
print("Absorbing regions before doublet check:", abs_vels2796)

# Check for aligned 3 sigma detections in non-dominant transitions
ew_spec2803, ew_sig2803 = ss.get_ew_spec(CON, ATOM, HIRES, "MgII2803", 1.000000,
                                         MgII_.vels, MgII_.waves, MgII_.f_norm, MgII_.I_sig)
abs_flags2803 = ss.get_abs_regions(MgII_.vels, ew_spec2803, ew_sig2803, sigma_threshold=3.0,
                                   dominant=False, region_vels=abs_vels2796)
abs_vels2796 = ss.dblt_checker(abs_vels2796, abs_flags2803)
print("Absorbing regions after doublet check:", abs_vels2796)
