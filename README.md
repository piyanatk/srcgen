srcgen generates HEALPix map of extragalactic radio sources by randomly fill
the sky with sources based on Poisson distribution, drawing their fluxes and
spectral indexes from CDFs derived from Subrahmanyan et al. 2002.

This program is a striped down version of simsky (original author: Judd Bowman),
which also includes routines to generate diffuse emission (that seems to be
unfinished).

Installation
------------
srcgen requires HEALPix and CFITSIO packages, included in the tar directory.
Note that the more recent version of HEALPic and CFITSIO will not work due to
change in libraries. In general you can:
1. `tar xvf tar/Healpix_2.15a_2010Jun18.tar.gz`
2. `cd Healpix_2.15a`
3. `./configure`
4. Build the HEALPix C++ library with the cfitsio tarbal in tar directory.
5. Edit the Makefile in the base directory to point to HEALPix libraies that
you build
6. run `make`
