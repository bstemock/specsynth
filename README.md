# Specsynth

Specsynth produces continuum-normalized quasar absorption line spectra with random noise. The advantage of this package is the user's ability to visualize the spectrum throughout its generation as demonstrated in example.py. Added functionality includes line detection and a check to ensure components are aligned for doublets. This code is modeled after fortran code written by Churchill (1997).

## Getting Started

Class and function descriptions are provided further on in this readme and the script example.py provides a helpful tutorial with the proper usage of relevent classes and functions.

### Prerequisites

Specsynth is written in Python 3 and uses the common libraries numpy and pandas and example.py uses matplotlib.pyplot to provide the user with visualizations of an example spectrum.

## Classes

### Absorber(CON, ATOM, trans, vels, Inst, seed, snr, zabs, v, logN, b)

Contains atomic data for the absorber as well as rest-frame velocity and wavelength grids and the normalized flux of the spectrum before convolution with the instrument ISF, after convolution with the instrument ISF, and after the addition of random noise.

**Parameters:**

* **CON : *dict***  
	The dictionary produced using data/const.dek.
* **ATOM *dict***  
        The dictionary produced using data/atoms.dat.
* **trans : *str***  
        Transition name. Must match a transition listed in the first column of atoms.dat.
* **vels : *array\_like***  
        Either the beginning and end values for the desired rest-frame velocity window or the desired velocity grid.
* **Inst : *instrument class***  
	The instrument class, initialized separately.
* **seed : *int***  
	Random seed.
* **snr : *float***  
	Signal-to-noise ratio.
* **zabs : *int***  
	Absorber redshift.
* **v : *list***  
	List of velocities in km s<sup>-1</sup> for each absortion line in the spectrum. Must have the same length as **logN** and **b**.
* **logN : *list***  
        List of log (column densities / cm<sup>-2</sup>) for each absortion line in the spectrum. Must have the same length as **v** and **b**.
* **b : *list***  
        List of Doppler b parameters in km s<sup>-1</sup> for each absortion line in the spectrum. Must have the same length as **v** and **logN**.

**Attributes:**

* **trans : *str***  
	Atomic transition of the absorber.
* **atom : *pandas.core.series.Series***  
	Atomic data for the absorber transition. Includes rest-frame wavelength in angstroms, oscilator strength, damping constant, atomic mass in amu, and ionization potential in eV.
* **zabs : *float***  
	Absorber redshift.
* **vels\_os : *1darray***  
	Oversampled rest-frame velocity grid to ensure smooth convolution with the intrument's ISF. Oversampling factor is **Inst.resfac**.
* **waves\_os : *1darray***  
	Oversampled wavelength grid to ensure smooth convolution with the instrument's ISF. Oversampling factor is **Inst.resfac**.
* **tau : *1darray***  
	Optical depth at each pixel.
* **f\_norm\_preconv : *1darray***  
	Normalized flux before convolution with the instrument ISF.
* **vels : *1darray***  
	Rest-frame velocity grid.
* **waves : *1darray***  
	Wavelength grid.
* **f\_norm\_noiseless : *1darray***  
	Normalized flux after convolution, but before the addition of random noise.
* **f\_norm : *1darray***  
	Normalized flux.
* **I\_sig : *1darray***  
	Uncertainty spectrum.


### Instrument(CON, R, presel, rdnoise, slit, n, resfac)

Contains specified instrument parameters and the ISF used for convolution by the **Absorber** class. Currently generates only a gaussian ISF.

**Parameters:**

* **CON : *dict***  
        The dictionary produced using data/const.dek.
* **R : *float***  
	Resolving power.
* **presel : *float***  
	Pixels per resolution element.
* **rdnoise : *float***  
	Read noise.
* **slit : *float***  
	Slit width in cm. Currently not operational.
* **n : *float***  
	Number of standard deviations to include in the gaussian ISF.
* **resfac : *int***  
	Oversampling factor. The ISF and Absorber pixel grids will be oversampled by this value to provide a smoother convolution.

**Attributes:**

* **slit : *float***  
	Slit width.
* **R : *float***  
	Resolving power.
* **presel : *float***  
	Pixels per resolution element.
* **vresel : *float***  
	Resolution element in km s<sup>-1</sup>.
* **rdnoise : *float***  
        Read noise.
* **resfac : *int***  
        Oversampling factor. The ISF and Absorber pixel grids will be oversampled by this value to provide a smoother convolution.
* **sigma_pix : *float***  
	Standard deviation of the gaussian kernel on the oversampled grid.
* **dv : *float***  
	Velocity width of each pixel.
* **num\_pix : *float***  
	Number of pixels in half of the gaussian ISF.
* **pix\_isf : *1darray***  
	Oversampled pixel grid for the gaussian ISF.
* **gauss\_kernel : *1darray***  
	Oversampled gaussian ISF. Used for convolution by the **Absorber** class.

## Functions

Many functions in these files are used internally by the Absorber and Instrument classes and are not intended for use by users. These functions are absent from the following list, but more information is available upon request.

### get\_constants(path)

Retrives constants from const.dek.

**Parameters:**

* **path : *str***  
	Path to const.dek.

**Returns:**

* **CON : *dict***  
	Dictionary containing all constants in const.dek.

### get\_atomic(path)

Retrieves atomic data from atoms.dat.

**Parameters:**

* **path : *str***  
	Path to atoms.dat.

**Returns:**

* **ATOM : *dict***  
	Dictionary containing atomic data for all transitions in atoms.dat.

### wave\_to\_vel(CON, waves, wave\_cen, zabs)

Converts wavelengths in angstroms to rest-frame velocities in km s<sup>-1</sup>.

**Parameters:**

* **CON : *dict***  
        The dictionary produced using data/const.dek.
* **waves : *float* or *array-like***  
	Wavelength(s) in angstroms to be converted.
* **wave\_cen : *float***  
	Rest-frame wavelength of the absorber transition.
* **zabs : *float***  
	Absorber redshift.

**Returns:**

* **vels : *float* or *array-like***  
	Rest-frame velocity/velocities, in km s<sup>-1</sup>.

### vel\_to\_wave(CON, vels, wave\_cen, zabs)

Converts rest-frame velocities in km s<sup>-1</sup> to wavelengths in angstroms.

**Parameters:**

* **CON : *dict***  
        The dictionary produced using data/const.dek.
* **vels : *float* or *array-like***  
        Rest-frame velocity/velocities in km s<sup>-1</sup> to be converted.
* **wave\_cen : *float***  
        Rest-frame wavelength of the absorber transition.
* **zabs : *float***  
        Absorber redshift.

**Returns:**

* **waves : *float* or *array-like***  
        Wavelength(s) in angstroms.

### get\_ew\_spec(CON, Abs, Inst)

Calculates the equivalent width spectrum of an absorber.

**Parameters:**

* **CON : *dict***  
        The dictionary produced using data/const.dek.
* **Abs : *Absorber***  
	Absorber class object.
* **Inst : *Instrument***  
	Instrument class object.

**Returns**

* **ew\_spec : *1darray***  
	Equivalent width spectrum.
* **ew\_sig : *1darray***  
	Equivalent width uncertainty spectrum.

### get\_abs\_regions(Abs, ew\_spec, ew\_sig, sigma\_threshold=3.0, dominant=True, region\_vels=None)

Uses the equivalent width and equivalent width uncertainty spectra of an abosrber to detect absorption lines in the spectrum produced by an Absorber class object.

**Parameters:**

* **Abs : *Absorber***  
        Absorber class object.
* **ew\_spec : *1darray***  
        Equivalent width spectrum.
* **ew\_sig : *1darray***  
        Equivalent width uncertainty spectrum.
* **sigma\_threshold : *1darray***  
	Detection threshold. Pixels in the equivalent width spectrum that exceed this number times the corresponding equivalent width uncertainty spectrum value will be flagged as detections.
* **dominant : *bool***  
	If *True*, the provided absorber is treated as the dominant member of the doublet and the entire spectrum is checked for detections. If *False*, only the detection regions provided to **region\_vels** will be checked for detections.
* **region\_vels : *list***  
	List of detection regions using the format of this function's output.

**Returns:**

* **detection\_vels : *list*** or **detection\_flags : *list***  
	If **dominant=True**, returns a list of detection regions. Each detection region is a list containing rest-frame velocities in km s<sup>-1</sup> for the beginning and end pixels for the given detection region. If **dominant=False**, returns a list of booleans indicating whether or not a detection occurred in each of the provided detection regions.

### dblt\_checker(region\_vels, region\_flags)

Returns only detection regions in which both doublet members yielded a detection.

**Parameters:**

* **region\_vels : *list***  
	A list of detection region beginning and end velocities for the dominant member of the doublet.
* **region\_flags : *list***  
	A list of booleans indicating whether or not the non-dominant member of the doublet yielded a detection in each detection region.

**Returns:**

* **detection\_regions : *list***  
	The provided **region\_vels** list, but without all elements that lacked a detection in that region as indicated by **region\_flags**.

### cleanspec(Abs, region\_vels)

Sets all pixel flux values to unity outside of absorbing regions.

**Parameters:**

* **Abs : *Absorber***  
        Absorber class object.
* **region\_vels : *list***  
	List of detection regions. Each detection region is a list containing rest-frame velocities in km s<sup>-1</sup> for the beginning and end pixels for the given detection region.

**Returns:**

* **Abs : *Absorber***  
	The modified absorber class object. In **Abs.f\_norm**, all flux values outside of the provided detection regions have been set to unity.

### get_Naod(ATOM, trans, flux, sigma)

Returns the apparent optical depth (AOD) column density spectrum (cm<sup>-2</sup> pix<sup>-1</sup>) of a given flux spectrum.

**Parameters:**

* **ATOM : *dict***  
        Dictionary containing atomic data for all transitions in atoms.dat.
* **trans : *str***
        Transition name. Must match a transition listed in the first column of atoms.dat.
* **flux : *1darray***
	Normalized flux.
* **sigma : *1darray***
	Uncertainty spectrum.

**Returns:**

* **Naod : *1darray***  
        The AOD column density spectrum in units of cm<sup>-2</sup> pix<sup>-1</sup>.

## Authors

Bryson Stemock

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE.txt](LICENSE.txt) file for details.

## Acknowledgements

This project is adapted from code written by Christopher W. Churchill (1997).
