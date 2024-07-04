import os
from astropy.io import fits
import numpy as np
import pandas as pd

def save_fits_image(hdu, path, name):
    image_path = os.path.join(path, f"{name}.fits")
    hdu.writeto(image_path, overwrite=True)
    print(f"Saved image: {image_path}")

def save_fits_table(hdu, path, name):
    table_path = os.path.join(path, f"{name}.csv")
    data = hdu.data
    columns = hdu.columns.names

    # Dictionary to store data columns
    df_dict = {}

    for col in columns:
        col_data = data[col]
        if isinstance(col_data[0], np.ndarray):
            # Merge multidimensional data into unique columns
            for idx in range(col_data.shape[1]):
                col_name = f"{col}_{idx}"
                if col in df_dict:
                    df_dict[col].extend([row[idx] for row in col_data])
                else:
                    df_dict[col_name] = [row[idx] for row in col_data]
        else:
            df_dict[col] = list(col_data)

    # Find the maximum length of all columns
    max_length = max(len(v) for v in df_dict.values())

    # Ensure all columns have the same length
    for col, values in df_dict.items():
        if len(values) < max_length:
            df_dict[col].extend([None] * (max_length - len(values)))

    # Create the DataFrame from the dictionary
    df_final = pd.DataFrame(df_dict)
    df_final.to_csv(table_path, index=False)
    print(f"Saved table: {table_path}")

def inspect_and_extract_fits_file(file_path, save_path):
    # Open the FITS file
    with fits.open(file_path) as hdul:
        # Iterate through each HDU (Header/Data Unit) to process each type
        for hdu in hdul:
            hdu_type = type(hdu).__name__
            hdu_name = hdu.name

            try:
                if hdu_type == 'ImageHDU':
                    save_fits_image(hdu, save_path, hdu_name)
                elif hdu_type == 'BinTableHDU':
                    save_fits_table(hdu, save_path, hdu_name)
                elif hdu_type == 'PrimaryHDU':
                    continue  # Skip PRIMARY HDU
                else:
                    print(f"Other HDU found: Name = {hdu_name}, Type = {hdu_type}")
            except Exception as e:
                print(f"WARNING: {hdu_name} extraction failed with error: {e}")

        # Check for expected elements and print a warning if any are missing
        expected_elements = ['DATA_INFORMATION', 'REDUCED_DATA', 'SNR_MAP', 'SENSITIVITY_MAP', 'DETECTION_LIMIT', 'SOURCE_DETECTION']
        for element in expected_elements:
            if element not in [hdu.name for hdu in hdul]:
                print(f"WARNING: {element} is missing")


# path to your FITS file and the save directory
PATH = '/Users/herve/Documents/recherche/tmp/'
file_path = PATH + 'HIP_114189_DC_id_151321.fits'
save_path = PATH
inspect_and_extract_fits_file(file_path, save_path)
