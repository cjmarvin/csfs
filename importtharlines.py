import numpy as np

def import_thar_lines():
    print "Importing ThAr spectra..."

    # ThAr SPECTRA FILES
    kerber = "MODELS/THAR/ThAr_Kerber_2008.dat"
    lovis = "MODELS/THAR/table1.dat"

    # import lovis (2007) [Ang] -> [mm]
    lovis_data = np.loadtxt(lovis, usecols=(0,3))
    lovis_wavelengths = lovis_data[:, 0] * 1.0e-7
    lovis_intensities = 0.5 * lovis_data[:, 1] / lovis_data[:, 1].max()

    # import kerber (2008) [nm] -> mm
    kerber_data = np.loadtxt(kerber, usecols=(0,5))
    kerber_wavelengths = kerber_data[:, 0] * 1.0e-6
    kerber_intensities = kerber_data[:, 1] / kerber_data[:, 1].max()

    # CONCATENATE LOVIS AND KERBER WAVELENGTHS
    wavelengths = np.concatenate([lovis_wavelengths, kerber_wavelengths])
    intensities = np.concatenate([lovis_intensities, kerber_intensities])
    indices = np.argsort(wavelengths)
    intensities = intensities[indices]
    wavelengths = np.sort(wavelengths)

    return wavelengths, intensities

if __name__ == "__main__":
    wavelengths, intensities = import_thar_lines()
    np.savetxt("thar_lines.txt", zip(wavelengths, intensities))