from hcifits import HCIFits
from astropy.table import Table
from astropy.io import fits
import numpy as np
import pandas as pd

def explode_hcifits(path, name):

    # read and save hcifits file
    h = HCIFits.read(path + f"/{name}")

    # read SNR MAP
    snr_map = h.snr_map
    # Check if h.snr_map is an ImageHDU object and save it to a FITS file
    if isinstance(snr_map, fits.ImageHDU):
        snr_map.writeto(path +'/SNR_MAP.fits', overwrite=True)
    else:
        print("h.snr_map is not an ImageHDU object")

    # Read and save SENSITIVITY MAP
    sensitivity_map = h.sensitivity_map
    # Check if h.sensitivity_map is an ImageHDU object and save it to a FITS file
    if isinstance(sensitivity_map, fits.ImageHDU):
        sensitivity_map.writeto(path +'/SENSITIVITY_MAP.fits', overwrite=True)
    else:
        print("h.sensitivity_map is not an ImageHDU object")

    # Read and save REDUCED image
    reduced_data = h.reduced_data
    # Check if h.reduced_data is an ImageHDU object and save it to a FITS file
    if isinstance(reduced_data, fits.ImageHDU):
        reduced_data.writeto(path + '/REDUCED_DATA.fits', overwrite=True)
    else:
        print("h.reduced_data is not an ImageHDU object")

    # Read and Save data_information table
    data_info_table = Table.read(h.data_information)
    # Extract Wavelength column from data_info_table
    wavelengths = data_info_table['Wavelength']
    nbr_wavelength = np.size(data_info_table['Wavelength'])
    data_info_table.write(path + '/DATA_INFORMATION.csv', format='csv', overwrite=True)

    # Read and transform source_detection table
    source_detection_table = Table.read(h.source_detection)

    # lists for the columns of the DataFrame
    candidates = []
    wavelength_list = []
    snr_list = []
    dra_list = []
    err_dra_list = []
    ddec_list = []
    err_ddec_list = []
    sep_list = []
    err_sep_list = []
    pa_list = []
    err_pa_list = []
    flux_cs_list = []
    err_flux_cs_list = []
    flux_mag_list = []
    err_flux_mag_list = []
    flux_jy_list = []
    err_flux_jy_list = []
    flux_erg_list = []
    err_flux_erg_list = []
    contrast_list = []
    err_contrast_list = []
    wp3_id_list = []
    wp3_status_list = []

    # Transform multidimensional to unidimensional columns
    for row in source_detection_table:
        for i in range(nbr_wavelength):
            candidates.append(row['Candidate'])
            wavelength_list.append(wavelengths[i])
            snr_list.append(row['SNR'][i])
            dra_list.append(row['dRA'][i])
            err_dra_list.append(row['err_dRA'][i])
            ddec_list.append(row['dDEC'][i])
            err_ddec_list.append(row['err_dDEC'][i])
            sep_list.append(row['Sep'][i])
            err_sep_list.append(row['err_Sep'][i])
            pa_list.append(row['PA'][i])
            err_pa_list.append(row['err_PA'][i])
            flux_cs_list.append(row['Flux_cs'][i])
            err_flux_cs_list.append(row['err_Flux_cs'][i])
            flux_mag_list.append(row['Flux_mag'][i])
            err_flux_mag_list.append(row['err_Flux_mag'][i])
            flux_jy_list.append(row['Flux_Jy'][i])
            err_flux_jy_list.append(row['err_Flux_Jy'][i])
            flux_erg_list.append(row['Flux_erg'][i])
            err_flux_erg_list.append(row['err_Flux_erg'][i])
            contrast_list.append(row['Contrast'][i])
            err_contrast_list.append(row['err_Contrast'][i])
            wp3_id_list.append(row['WP3_ID'])
            wp3_status_list.append(row['WP3_Status'])

    # Creat a pandas DataFrame
    df = pd.DataFrame({
        'Candidate': candidates,
        'Wavelength': wavelength_list,
        'SNR': snr_list,
        'dRA': dra_list,
        'err_dRA': err_dra_list,
        'dDEC': ddec_list,
        'err_dDEC': err_ddec_list,
        'Sep': sep_list,
        'err_Sep': err_sep_list,
        'PA': pa_list,
        'err_PA': err_pa_list,
        'Flux_cs': flux_cs_list,
        'err_Flux_cs': err_flux_cs_list,
        'Flux_mag': flux_mag_list,
        'err_Flux_mag': err_flux_mag_list,
        'Flux_Jy': flux_jy_list,
        'err_Flux_Jy': err_flux_jy_list,
        'Flux_erg': flux_erg_list,
        'err_Flux_erg': err_flux_erg_list,
        'Contrast': contrast_list,
        'err_Contrast': err_contrast_list,
        'WP3_ID': wp3_id_list,
        'WP3_Status': wp3_status_list
    })
    # Save the source_detection table
    df.to_csv(path + '/SOURCE_DETECTION.csv', index=False)

    # Read DETECTION_LIMIT table
    detection_limit_table = Table.read(h.detection_limit)
    nbr_radius = len(detection_limit_table['Radius'][0])

    # lists for the columns of the DataFrame
    radius = []
    wavelength_list = []
    detection_Limit = []

    # Transform multidimensional to unidimensional columns
    for row in detection_limit_table:
        for i in range(nbr_wavelength):
            for j in range(nbr_radius):
                radius.append(row['Radius'][j])
                wavelength_list.append(wavelengths[i])
                detection_Limit.append(row['Detection_Limit'][j])

    # Creat a pandas DataFrame
    df = pd.DataFrame({
        'Radius': radius,
        'Wavelength': wavelength_list,
        'Detection_Limit': detection_Limit,
    })
    # Save the detection limit table
    df.to_csv(path + '/DETECTION_LIMIT.csv', index=False)