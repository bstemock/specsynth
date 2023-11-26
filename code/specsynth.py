# functions used by specgen.py to generate spectra, locate absorption regions, and clean spectra
# DO NOT place absorption lines within half of the ISF kernel number of pixels of the edge
# Bryson Stemock, bstemock@nmsu.edu

import numpy as np
import pandas as pd

# constants for the function "voigt"
b1 = 0.850
t = [0.314240376, 0.947788391, 1.59768264, 2.27950708, 3.02063703, 3.8897249]
c = [1.01172805, -0.75197147, 1.2557727e-2, 1.00220082e-2, -2.42068135e-4, 5.00848061e-7]
s = [1.393237, 0.231152406, -0.155351466, 6.21836624e-3, 9.19082986e-5, -6.27525958e-7]


# read constants from <path to data/const.dek>
def get_constants(path):
    df = pd.read_table(path, index_col="constants", comment="#", delim_whitespace=True)
    hbar = 0.5 * df.loc["h"] / df.loc["pi"]
    sig = df.loc["pi"] * df.loc["pi"] * df.loc["kerg"] * df.loc["kerg"] * df.loc["kerg"] * df.loc["kerg"] /\
        (60.0 * hbar * hbar * hbar * df.loc["c"] * df.loc["c"])
    constants = ["deg2rad", "rad2deg", "hbar", "rbohr", "fine", "sig", "asol", "weinlam", "weinfre", "pc"]
    values = [
        df.loc["pi"] / 180.0,
        180.0 / df.loc["pi"],
        hbar,
        hbar * hbar / (df.loc["me"] * df.loc["e"] * df.loc["e"]),
        df.loc["e"] * df.loc["e"] / (hbar * df.loc["c"]),
        sig,
        4.0 * sig / df.loc["c"],
        df.loc["h"] * df.loc["c"] / (df.loc["kerg"] * 4.965114232),
        2.821439372 * df.loc["kerg"] / df.loc["h"],
        3.261633 * df.loc["lyr"]
    ]
    df_ = pd.DataFrame(data=values, index=constants, columns=df.columns)
    df = pd.concat((df, df_))
    return pd.Series(df["values"])


# read atomic data from <path to data/atoms.dat>
def get_atomic(path):
    df = pd.read_table(path, index_col="trans", delim_whitespace=True)
    return df


# convert wavelengths to rest-frame velocities
def wave_to_vel(CON, waves, wave_cen, zabs):
    return CON["ckms"] * (waves / (1.0 + zabs) - wave_cen) / wave_cen


# convert rest-frame velocities to wavelengths
def vel_to_wave(CON, vels, wave_cen, zabs):
    return (1.0 + zabs) * ((vels * wave_cen) / CON["ckms"] + wave_cen)


# compute pixel velocities oversampled by a factor of Inst.resfac to be used for convolution
def get_vels_os(vel_range, Inst):
    p_, p = np.floor(vel_range[0] / Inst.dv), np.ceil(vel_range[1] / Inst.dv)
    return np.arange(Inst.resfac * p_, Inst.resfac * p, 1.0) * Inst.dv / Inst.resfac


# voigt profile generation functions
def voigt_reg1(CON, x, dlamD, y1, y2):
    u = np.zeros(len(x))
    for i in range(6):
        r, r_ = x - t[i], x + t[i]
        d, d_ = 1.0 / (r * r + y2), 1.0 / (r_ * r_ + y2)
        d1, d2, d3, d4 = d * y1, d * r, d_ * y1, d_ * r_
        u += c[i] * (d1 + d3) - s[i] * (d2 - d4)
    return u / (np.sqrt(CON["pi"]) * dlamD)


def voigt_reg2(CON, x, dlamD, y, y1, y2):
    y3 = y + 3.0
    u = np.where(np.abs(x) < 12.0, np.exp(-x * x), 0.0)
    for i in range(6):
        r, r_ = x - t[i], x + t[i]
        r2, r2_ = r * r, r_ * r_
        d, d_ = 1.0 / (r2 + y2), 1.0 / (r2_ + y2)
        d1, d2, d3, d4 = d * y1, d * r, d_ * y1, d_ * r_
        u += y * (c[i] * (r * d2 - 1.50 * d1) + s[i] * y3 * d2) / (r2 + 2.250)
        u += y * (c[i] * (r_ * d4 - 1.50 * d3) - s[i] * y3 * d4) / (r2_ + 2.250)
    return u / (np.sqrt(CON["pi"]) * dlamD)


# detect absorption regions
def get_ew_spec(CON, Abs, Inst):
    J = int(2.0 * Inst.presel)
    N = int(2 * J + 1)

    # pad J pixels onto either end of necessary arrays
    down_pad = np.arange(Abs.vels[0] - (J + 1) * Inst.dv, Abs.vels[0] - Inst.dv, Inst.dv)       # velocity down_pad
    up_pad = np.arange(Abs.vels[-1] + Inst.dv, Abs.vels[-1] + (J + 1) * Inst.dv, Inst.dv)       # velocity up_pad
    vels = np.concatenate((down_pad, Abs.vels, up_pad))
    pad_num = int((vels.shape[0] - Abs.vels.shape[0]) / 2)                      # number of pixels padded onto one end
    down_pad_ = np.mean(Abs.I_sig[:J]) * np.ones(pad_num)                       # sigma spec down_pad
    up_pad_ = np.mean(Abs.I_sig[-J:]) * np.ones(pad_num)                        # sigma spec up_pad
    Dsig = np.concatenate((down_pad_, Abs.I_sig, up_pad_))                      # sigma_D = sigma_{f_norm} = I_sig
    f_norm = np.concatenate((np.ones(pad_num), Abs.f_norm, np.ones(pad_num)))

    # compute additional required arrays and construct empty arrays
    waves = (1.0 + Abs.zabs) * (vels * Abs.atom["wave_cen"] / CON["ckms"] + Abs.atom["wave_cen"])
    D = 1.0 - f_norm                                                                            # flux decrement
    sig_ISF = Abs.waves / (Inst.R * 2.35482)                                                    # sigma_{ISF}
    x_ki = np.empty((Abs.waves.shape[0], N))                                                    # x_ki
    P_ni = np.empty((Abs.waves.shape[0], N))                                                    # P_n
    dlambda = np.empty((Abs.waves.shape[0]))                                    # lambda_{i + 1} - lambda_i
    prod1 = np.empty((Abs.waves.shape[0], N))                                   # P_n * D_k, used in sum to compute w_i
    prod2 = np.empty((Abs.waves.shape[0], N))           # P_n**2 * sigma_{D_k}**2, used in sum to compute sigma_{w_i}

    # populate empty arrays and compute P**2
    for i in range(Abs.waves.shape[0]):
        for n in range(1, N + 1):
            k = i + n - 1 - J
            x_ki[i, n - 1] = (waves[i + pad_num] - waves[k + pad_num]) / sig_ISF[i]
        for n in range(1, N + 1):
            k = i + n - 1 - J
            denom = np.exp(-x_ki[i] * x_ki[i])
            P_ni[i, n - 1] = np.exp(-x_ki[i, n - 1] * x_ki[i, n - 1]) / np.sum(denom)
            prod1[i, n - 1] = P_ni[i, n - 1] * D[k + pad_num]
            prod2[i, n - 1] = np.sqrt(P_ni[i, n - 1] * P_ni[i, n - 1] * Dsig[k + pad_num] * Dsig[k + pad_num])
        dlambda[i] = waves[i + pad_num + 1] - waves[i + pad_num]
    P2 = np.sum(P_ni * P_ni, axis=1)

    # compute ew spectrum and ew uncertainty spectrum
    ew_spec = dlambda * np.sum(prod1, axis=1) / P2
    ew_sig = dlambda * np.sum(prod2, axis=1) / P2
    return ew_spec, ew_sig


# detect absorption lines above the given sigma threshold and return absorption region limits
def get_abs_regions(Abs, ew_spec, ew_sig, sigma_threshold=3.0, dominant=True, region_vels=None):
    detect1 = np.where(ew_spec > ew_sig, 1, 0)                              # True where ew_spec above 1 sigma threshold
    detectN = np.where(ew_spec > sigma_threshold * ew_sig, True, False)     # True where ew_spec above sigma_threshold

    # if dominant transition, scan entire spectrum for abs regions
    if dominant:
        detection_vels = []                                                 # list to contain absorption region limits

        # scan through entire spectrum to find abs regions
        i = 0
        while i < ew_spec.shape[0]:
            if detectN[i]:                                                  # if detection at 5 sigma
                try:
                    down_pix = i - np.where(detect1[i::-1] == 0)[0][0] + 1
                except IndexError:                                          # if abs region doesn't recover to 1 sigma
                    down_pix = 0
                try:
                    up_pix = i + np.where(detect1[i:] == 0)[0][0] - 1
                except IndexError:                                          # if abs region doesn't recover to 1 sigma
                    up_pix = ew_spec.shape[0] - 1
                detection_vels.append([Abs.vels[down_pix], Abs.vels[up_pix]])
                i = up_pix + 1
            else:
                i += 1
        return detection_vels

    # if non-dominant transition, check for sigma_threshold detection within dominant transition abs regions
    elif not dominant:
        detection_flags = []                                                # list to contain True/False for each region

        # scan through each dominant transition abs regions to find abs regions
        for region in region_vels:
            ind = np.where((Abs.vels > region[0]) & (Abs.vels < region[1]))[0]
            if True in detectN[ind[0] - 1:ind[-1] + 2]:
                detection_flags.append(True)
            else:
                detection_flags.append(False)
        return detection_flags


def dblt_checker(region_vels, region_flags):
    return [region for i, region in enumerate(region_vels) if region_flags[i]]


def cleanspec(Abs, region_vels):
    abs_pix = []
    for region in region_vels:
        ind = np.where((Abs.vels > region[0]) & (Abs.vels < region[1]))[0]
        if ind[0] != 0 and Abs.vels[ind[0]] != region[0]:
            ind = np.concatenate(([ind[0] - 1], ind))
        if ind[-1] != len(Abs.vels) - 1 and Abs.vels[ind[-1]] != region[1]:
            ind = np.concatenate((ind, [ind[-1] + 1]))
        abs_pix.append(ind)
    clean = np.ones(len(Abs.f_norm))
    for i in abs_pix:
        clean[i] = Abs.f_norm[i]
    Abs.f_norm = clean
    return Abs


# absorber class
class Absorber:
    def __init__(self, CON, ATOM, trans, vel_range, Inst, seed, snr, zabs, v, logN, b):
        self.trans = trans                                                                  # transition name (string)
        self.atom = ATOM.loc[trans]                                                         # atomic data
        self.zabs = zabs                                                                    # absorber redshift
        self.vels_os = get_vels_os(vel_range, Inst)                                         # oversampled velocity grid
        self.waves_os = vel_to_wave(CON, self.vels_os, self.atom["wave_cen"], self.zabs)    # oversampled wavelengths
        self.tau = self.get_tau(CON, self.zabs, v, logN, b)                                 # optical depth
        self.f_norm_preconv = np.exp(-self.tau)                                     # normalized flux before convolution
        self.vels, self.f_norm_noiseless = self.convolve_with_isf(Inst)             # properly-pixelized velocity, flux
        self.waves = vel_to_wave(CON, self.vels, self.atom["wave_cen"], self.zabs)  # redshifted wavelengths
        self.f_norm, self.I_sig = self.add_noise(Inst, snr, seed)                   # finalized flux and sigma spectrum

    # compute y from atomic parameters and dlamD
    def get_y(self, CON, dlamD):
        y = self.atom["damp"] * self.atom["wave_cen"] * self.atom["wave_cen"] * 1e-8 / (4.0 * CON["pi"] *
                                                                                        CON["c"] * dlamD)
        return y

    # compute voigt function for a given x, y where x is an array
    def voigt(self, CON, x, dlamD):
        y = self.get_y(CON, dlamD)                                          # ratio of Lorentzian HWHM to Doppler width
        y1 = y + 1.50                                                       # constants for function "voigt"
        y2 = y1 * y1
        b2 = 18.10 * y + 1.650
        if y > b1:
            u = voigt_reg1(CON, x, dlamD, y1, y2)
            return u
        elif y <= b1:
            u = np.where(np.abs(x) < b2, voigt_reg1(CON, x, dlamD, y1, y2), voigt_reg2(CON, x, dlamD, y, y1, y2))
            return u

    # compute optical depth and loop over each component
    def get_tau(self, CON, zabs, v, logN, b):
        tau = np.zeros(self.waves_os.shape)
        for i in range(len(v)):
            dlamD = b[i] * self.atom["wave_cen"] / CON["ckms"]                                          # Doppler width
            shifted_wave_cen = self.atom["wave_cen"] * (1.0 + v[i] / CON["ckms"]) * (1.0 + zabs)
            x = (self.waves_os - shifted_wave_cen) / (dlamD * (1.0 + zabs))  # wave number scale in Doppler width units
            u = self.voigt(CON, x, dlamD)                                                               # voigt function
            tau += (10 ** logN[i] * CON["pi"] * CON["e"] * CON["e"] * self.atom["wave_cen"] * self.atom["wave_cen"]
                    * 1e-8 * self.atom["osc_str"] * u) / (CON["me"] * CON["c"] * CON["c"])
        return tau

    # convolve f_norm with instrument kernel and then average over every "resfac" number of pixels
    def convolve_with_isf(self, Inst):
        pad = np.ones(int(Inst.num_pix))
        conv = np.concatenate((pad, np.convolve(self.f_norm_preconv, Inst.gauss_kernel, mode="valid"), pad))
        vels = np.mean(self.vels_os.reshape(-1, int(Inst.resfac)), axis=1)
        f_norm = np.mean(conv.reshape(-1, int(Inst.resfac)), axis=1)
        return vels, f_norm

    # add noise to spectrum
    def add_noise(self, Inst, snr, seed):
        np.random.seed(seed)
        I_c = 0.5 * snr * snr * (1.0 + np.sqrt(1.0 + 4.0 * Inst.rdnoise * Inst.rdnoise / (snr * snr)))  # continuum flux
        I_sig = np.sqrt(I_c * self.f_norm_noiseless + Inst.rdnoise * Inst.rdnoise) / I_c                # sigma spectrum
        noise = I_sig * np.random.normal(size=I_sig.shape)                                              # actual noise
        f_norm = self.f_norm_noiseless + noise                                              # spectrum with noise added
        return f_norm, I_sig


# instrument class
# HIRES, UVES - R=45000, p = 3
# COS - ISF isn't a simple gaussian, has weird wings, will have to download the kernel
class Instrument:
    def __init__(self, CON, R, presel, rdnoise, slit, n, resfac):
        self.slit = slit
        self.R = R
        self.presel = presel
        self.rdnoise = rdnoise
        self.resfac = resfac
        self.vresel = CON["ckms"] / R                                           # resolution element (FWHM) in velocity
        self.sigma_pix = presel * resfac / 2.35482                              # std dev of gaussian in oversampled pix
        self.dv = self.vresel / presel                                          # dv per pixel
        self.num_pix = np.ceil(n * self.sigma_pix)                              # pixels in half of the gaussian ISF
        self.pix_isf = np.arange(-self.num_pix, self.num_pix + 1.0, 1.0)
        self.gauss_kernel = self.get_gauss_kernel(CON)

    # compute normalized gaussian ISF kernel
    def get_gauss_kernel(self, CON):
        return np.exp(-self.pix_isf * self.pix_isf / (2.0 * self.sigma_pix * self.sigma_pix)
                      ) / np.sqrt(2.0 * CON["pi"] * self.sigma_pix * self.sigma_pix)
