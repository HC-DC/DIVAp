# Should we keep the headers of the cube ?
# IMPORTANT: Saving a read FITS file will not work because of headers.
import logging
from collections import defaultdict
from datetime import datetime

from astropy import units as u
from astropy.io import fits
from astropy.io.fits import BinTableHDU, HDUList, ImageHDU, PrimaryHDU
from astropy.table import Table

INFO_TABLE_COLNAMES = [
    "Image_Number", "Orientation", "Combined_rotation_angle",
    "Number_of_Exposures", "Exposure_Time", "Observation_Start",
    "Observation_End", "UT_Midpoint_Date_of_Observation",
    "UT_Midpoint_Time_of_Observation", "Wavelength", "Bandwidth",
    "Polarization"]
INFO_TABLE_DTYPES = [
    "i2", "f8", "f8", "i2", "f8", "f8", "f8", "U160", "U160", "f8", "f8",
    "U160"]
DET_TABLE_COLNAMES = [
    "Candidate", "SNR", "dRA", "err_dRA", "dDEC", "err_dDEC", "Sep", "err_Sep",
    "PA", "err_PA", "Flux_cs", "err_Flux_cs", "Flux_mag", "err_Flux_mag",
    "Flux_Jy", "err_Flux_Jy", "Flux_erg", "err_Flux_erg", "Contrast",
    "err_Contrast"]
DET_TABLE_DTYPES = [
    "i2", "f8", "f8", "f8", "f8", "f8", "f8", "f8", "f8", "f8", "f8", "f8",
    "f8", "f8", "f8", "f8", "f8", "f8", "f8", "f8"]

log = logging.getLogger(__name__)


def create_empty_info_table():
    """Create an empty table to store image information."""
    table = Table(
        names=INFO_TABLE_COLNAMES,
        dtype=INFO_TABLE_DTYPES
    )

    table['Orientation'].unit = u.deg
    table['Combined_rotation_angle'].unit = u.deg
    table['Exposure_Time'].unit = u.s
    table['Wavelength'].unit = u.Angstrom
    table['Bandwidth'].unit = u.Angstrom

    return table


def create_empty_det_table(nb_images):
    """Create and empty source detections table.

    Parameters
    ----------
    nb_images: int
        Number of images.

    """
    table = Table(
        names=DET_TABLE_COLNAMES,
        dtype=DET_TABLE_DTYPES
    )

    # The information on each candidate is a vector which length is the number
    # of images.
    for col in DET_TABLE_COLNAMES[1:]:  # No for Candidate column
        table[col].shape = (0, nb_images)

    # Set the units
    for col in DET_TABLE_COLNAMES:
        if col in ["dRA", "err_dRA", "dDEC", "err_dDEC", "Sep", "err_Sep"]:
            table[col].unit = u.arcsec
        elif col in ["PA", "err_PA"]:
            table[col].unit = u.deg
        elif col.endswith("_cs"):
            table[col].unit = u.count / u.s
        elif col.endswith("mag"):
            table[col].unit = u.mag
        elif col.endswith("_Jy"):
            table[col].unit = u.Jy
        elif col.endswith("_erg"):
            table[col].unit = u.erg / u.cm**2 / u.s / u.Angstrom
        elif col.endswith("Contrast"):
            table[col].unit = u.dmag

    return table


def create_detection_limit_table(radiuses, limits):
    """Creates a dection limit table."""


def nim_nx_ny(image_hdu):
    """Return the number of maps and their size."""
    if image_hdu.header['NAXIS'] == 2:
        return 1, image_hdu.header['NAXIS2'], image_hdu.header['NAXIS1']

    if image_hdu.header['NAXIS'] == 3:
        return (image_hdu.header['NAXIS3'], image_hdu.header['NAXIS2'],
                image_hdu.header['NAXIS1'])

    raise ValueError("Wrong number of axis.")


def hdu_type(hdu):
    """String indicating the type of HDU for the headers."""
    if isinstance(hdu, ImageHDU):
        return "IMAGE"

    if isinstance(hdu, BinTableHDU):
        return "BINTABLE"

    raise ValueError("Unknown HDU type.")

class HCIFits:
    """Class representing an HCI-FITS.

    For a description of the parameters, please read the HCI-FITS
    documentation.

    To produce a well organised main header as in the DIVA file, the `headers`
    attribute contains a dictionary of dictionaries. The first key is the
    (pseudo) section in the header, the second key is the keyword and the
    second level value is a tuple (keyword value, commentary).

    The `extra_hdus` attribute contains a list of HDUs that will be added to
    the FITS file.

    Parameters
    ----------
    data_information: `astropy.table.Table` or `astropy.io.fits.BinTableHDU`
        Table containing the information about the reduced data maps.
    reduced_data: `astropy.io.fits.ImageHDU`
        ImageHDU containing the reduced data: a Nx×Ny image or a Nim×Nx×Ny
        cube.
    snr_map: `astropy.io.fits.ImageHDU`
        ImageHDU containing the SNR maps corresponding to the reduced images.
    sensitivity_map: `astropy.io.fits.ImageHDU`
        ImageHDU containing the sensitivity maps of the reduced images. It must
        contain the NSIGMA keyword to indicate the confidence level of the
        sensitivity map.
    detection_limit: `astropy.io.fits.BinTableHDU`
        Detection limit table.
    source_detection: `astropy.table.Table` or `astropy.io.fits.BinTableHDU`
        Source detection table.

    """

    def __init__(self, data_information, reduced_data, snr_map,
                 sensitivity_map, detection_limit, source_detection=None,
                 extra_hdus=None):
        if isinstance(data_information, BinTableHDU):
            self.data_information = data_information
        else:
            self.data_information = BinTableHDU(data_information)
        self.data_information.name = "DATA_INFORMATION"
        self.reduced_data = reduced_data
        self.reduced_data.name = "REDUCED_DATA"
        self.snr_map = snr_map
        self.snr_map.name = "SNR_MAP"
        self.sensitivity_map = sensitivity_map
        self.sensitivity_map.name = "SENSITIVITY_MAP"
        self.detection_limit = detection_limit
        self.detection_limit.name = "DETECTION_LIMIT"
        if source_detection is not None:
            if isinstance(source_detection, BinTableHDU):
                self.source_detection = source_detection
            else:
                self.source_detection = BinTableHDU(source_detection)
            self.source_detection.name = "SOURCE_DETECTION"
        else:
            self.source_detection = None
        self.extra_hdus = []
        if extra_hdus is not None:
            self.extra_hdus = extra_hdus
        self.headers = defaultdict(dict)

    def add_extra_hdu(self, hdu):
        """Add an extra HDU."""
        if hdu.name is None:
            raise ValueError("The HDU must have an EXTNAME.")
        self.extra_hdus.append(hdu)

    def is_valid(self, strict=False):
        """Check if the content conforms to the standard.

        This method check the content of the HCIFits object against the
        standard.  It logs major problems as errors and minor problems as
        warnings.

        By default, if there is no major problem the object is valid; if strict
        is True, the object is valid if there is no major nor minor problem.

        Parameters
        ----------
        strict: bool
            If False (default) the object is valid if there is no major
            problem; if true, the object is valid if there is no major nor
            minor problem.

        """
        minor_pb, major_pb = False, False

        # Sizes

        nim, nx, ny = nim_nx_ny(self.reduced_data)
        try:
            if nim_nx_ny(self.snr_map) != (nim, nx, ny):
                major_pb = True
                log.error("SNR map has wrong size.")
        except ValueError:
            major_pb = True
            log.error("SNR map has wrong axis number.")
        try:
            if nim_nx_ny(self.sensitivity_map) != (nim, nx, ny):
                major_pb = True
                log.error("Sensitivity map has wrong size.")
        except ValueError:
            major_pb = True
            log.error("Sensitivity map has wrong axis number.")

        if len(self.data_information.data) != nim:
            major_pb = True
            log.error("Data information table has wrong length.")

        # Table columns

        for tpl_col in BinTableHDU(create_empty_info_table()).columns:
            try:
                col = self.data_information.columns[tpl_col.name]
                if tpl_col.format != col.format:
                    major_pb = True
                    log.error(
                        "Data information column %s has wrong format: %s , "
                        "expected %s.", col.name, col.format, tpl_col.format)
                if tpl_col.unit is not None and col.unit is None:
                    major_pb = True
                    log.error(
                        "Data information column %s has not unit.", col.name)
                if tpl_col.unit != col.unit:
                    minor_pb = True
                    log.warning(
                        "Data information column %s has an unexpected unit: "
                        "%s, expected %s", col.name, col.unit, tpl_col.unit)
            except KeyError:
                major_pb = True
                log.error(
                    "Data information column %s is missing", col.name)

        if self.source_detection is not None:
            for tpl_col in BinTableHDU(create_empty_det_table(nim)).columns:
                try:
                    col = self.source_detection.columns[tpl_col.name]
                    if tpl_col.format != col.format:
                        major_pb = True
                        log.error(
                            "Source detection column %s has wrong format: "
                            "%s, expected %s.", col.name, col.format,
                            tpl_col.format)
                    if tpl_col.unit is not None and col.unit is None:
                        major_pb = True
                        log.error(
                            "Source detection column %s has not unit.",
                            col.name)
                    if tpl_col.unit != col.unit:
                        minor_pb = True
                        log.warning(
                            "Source detection column %s has a unexpected "
                            "unit: %s, expected %s.", col.name, col.unit,
                            tpl_col.unit)
                except KeyError:
                    major_pb = True
                    log.error(
                        "Source detection column %s is missing", col.name)

        if strict:
            return not (major_pb and minor_pb)
        return not major_pb

    def save(self, filename, force=False):
        """Save to a HCI-FITS file.

        Parameters
        ----------
        filename: string
            Name of the FITS file.
        force: bool
            By default a non-valid object won't be saved. Set force to True to
            save to a non-valid HCI-FITS.

        """
        if not self.is_valid() and not force:
            raise ValueError("The object is not valid.")

        hdu_list = HDUList()

        primary_hdu = PrimaryHDU()
        header = primary_hdu.header
        hdu_list.append(primary_hdu)

        # Set generic headers
        header.append(('', '', ''), end=True)
        header.append(('', 'FILE INFORMATION', ''), end=True)
        header.append(('', '', ''), end=True)
        header.append(
            ('FILETYPE',
             'HLSP',
             'type of data found in data file'),
            end=True
        )
        header.append(
            ('ORIGIN',
             'DIVA team',
             'FITS file originator'),
            end=True
        )
        header.append(
            ('DATE',
             datetime.utcnow().strftime('%Y-%m-%d'),
             'date this file was written (yyyy-mm-dd)'),
            end=True)
        header.append(('', '', ''), end=True)
        header.append(('', '', ''), end=True)

        # Add HDUs
        header.append(('', '', ''), end=True)
        header.append(('', 'FITS FILE DATA STRUCTURE', ''), end=True)
        header.append(('', '', ''), end=True)

        hdus_to_add = [self.data_information, self.reduced_data,
                       self.snr_map, self.sensitivity_map,
                       self.detection_limit]
        if self.source_detection is not None:
            hdus_to_add.append(self.source_detection)
        hdus_to_add += self.extra_hdus

        header.append(('NEXTEND', len(hdus_to_add),
                       "number of standard extensions"),
                      end=True)
        for idx, hdu in enumerate(hdus_to_add):
            header.append(
                (f"EXT{idx+1}NAME",
                 hdu.name,
                 f"extension {idx+1} name"),
                end=True)
            header.append(
                (f"EXT{idx+1}TYPE",
                 hdu_type(hdu),
                 f"extension {idx+1} type"),
                end=True)
            hdu_list.append(hdu)

        # Supplementary headers
        for section in self.headers:
            header.append(("", "", ""), end=True)
            header.append(("", section, ""), end=True)
            header.append(("", "", ""), end=True)
            for keyword, (value, comment) in self.headers[section].items():
                header.append((keyword, value, comment), end=True)

        hdu_list.writeto(filename)

    @classmethod
    def read(cls, filename):
        """Read a HCI-FITS file."""
        main_extensions = ['DATA_INFORMATION', 'REDUCED_DATA', 'SNR_MAP',
                           'SENSITIVITY_MAP', 'DETECTION_LIMIT',
                           'SOURCE_DETECTION', 'PRIMARY']
        with fits.open(filename, memmap=False) as hdu_list:
            data_information = hdu_list['DATA_INFORMATION'].copy()
            reduced_data = hdu_list['REDUCED_DATA'].copy()
            snr_map = hdu_list['SNR_MAP'].copy()
            sensitivity_map = hdu_list['SENSITIVITY_MAP'].copy()
            detection_limit = hdu_list['DETECTION_LIMIT'].copy()
            try:
                source_detection = hdu_list['SOURCE_DETECTION'].copy()
            except KeyError:
                source_detection = None
            extra_hdus = [hdu.copy() for hdu in hdu_list if hdu.name not in
                          main_extensions]
            primary_header = hdu_list['PRIMARY'].header.copy()

        new_hcifits = cls(
            data_information=data_information,
            reduced_data=reduced_data,
            snr_map=snr_map,
            sensitivity_map=sensitivity_map,
            detection_limit=detection_limit,
            source_detection=source_detection,
            extra_hdus=extra_hdus
        )
        new_hcifits.primary_header = primary_header

        return new_hcifits
