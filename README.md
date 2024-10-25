This code is intended for the stitching of spectra collected from different STIS gratings, although it should function with spectra from other instruments as well.
All data files are assumed to be .txt files, with columns labeled "wave", "flux", and "error", separated by tabs.
All pixels in your data are assumed to have nonzero values for their reported uncertainties.

Before running the code, please create a folder in which the only contained files are the data files to be stitched.

Then run
```
python stischer.py
```
to start.

Please note that overlapping segments of spectra are merged through linear interpolation and weighted means determined via reported uncertainty.

To Do:
- Expand capabilities to enable different file formats, such as differently labeled columns, or even the actual .fits files.
- Allow for zero values of uncertainties.

Inspired by CoAdd.pro by Xinfeng Xu.
