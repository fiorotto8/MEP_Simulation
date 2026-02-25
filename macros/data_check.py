from astropy.io import fits

# Replace 'your_file.fits' with the actual file path
file_path = 'output.fits'

# Open the FITS file
hdulist = fits.open(file_path)

# Loop through HDUs and print header and data names
for i, hdu in enumerate(hdulist):
    print(f"HDU {i + 1}:")
    print("Header Names:")
    print(hdu.header.keys())
    print("Data Names:")
    print(hdu.data.names if hdu.data is not None else "No data in this HDU")
    print("-" * 40)

# Close the FITS file
hdulist.close()
