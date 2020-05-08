![Logo](/doc/images/logo/lasif_logo.png)
---
[![Build Status](https://travis-ci.org/krischer/LASIF.png?branch=master)](https://travis-ci.org/krischer/LASIF)
[![GPLv3](http://www.gnu.org/graphics/gplv3-88x31.png)](https://github.com/krischer/LASIF/blob/master/LICENSE)


## Specfem3d_fwi flavored LASIF

This fork of LASIF is being developed for the use with specfem3d_fwi  
Also, original Lasif codes have been modified for supporting python3 and Qt5.  

Prerequesites:
You must have installed Python3, and libraries `geos`, `openmpi`, `mpi4py`

I recommend to create a virtual environment (first install virtualenv with `pip install virtualenv`
`virtualenv Lasif-Py3`
Activate the Lasif environment
`source Lasif-Py3/bin/activate`

to install,  
`pip3 install -r requirements.txt`  
then  
`pip3 install .`  


Documentation about can be found here: [LASIF](http://krischer.github.io/LASIF)


### Paper

For more details and the reasoning behind LASIF, please also see our paper:

*Lion Krischer, Andreas Fichtner, Saule Zukauskaite, and Heiner Igel* (2015),
**Large‐Scale Seismic Inversion Framework**, Seismological Research Letters, doi:10.1785/0220140248.


* [DOI Link](http://dx.doi.org/10.1785/0220140248)
