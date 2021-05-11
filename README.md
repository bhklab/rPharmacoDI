# rPharmacoDI

This package is the R compliment to the PharmacoDI Python package. It provides the functions needed
for interfacing with a PSet object to generate a .csv representation of each PSet.

Additional functions for downloading precomputed data which is needed for creating database tables in PharmacoDI
are also included and processed to .csv.

The goal of this package is to eliminate or minimize the need for dependencies external to a given package/module,
as would be the case using `rp2` from Python as well as for using `reticulate` from R.