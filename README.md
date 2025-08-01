# "Magglow" - Magnetic Bullet Afterglow 

<img align="left" width="422" height="321" alt="Image" src="https://github.com/user-attachments/assets/d63eb294-e4f5-47cd-b972-8359382bf4f8" />

**Magglow** is a Julia module to calculate Gamma-ray Burst (GRB) afterglow light curves and spectra based on Magnetic Bullet model [Y. Kusafuka & K. Asano (2024)](https://ui.adsabs.harvard.edu/abs/2025MNRAS.536.1822K/abstract). 
The framework implements semi-analytic models for forward and reverse shock dynamics propagating in a stratified CSM, leptonic (synchrotron with self-absorption and self-Compton scattering with Klein-Nishina corrections) and hadronic (pp colisions and photomeson interactions) multimessenger emission mechnisms, and arbitrary viewing angle.

Details of the methods can be found in [Y. Kusafuka, K. Obayashi, K. Asano, & R. Yamazaki (in prep)](). 
<!-- This code is under active development.  -->

<!-- Documentation is available at <https://afterglowpy.readthedocs.io/> -->

<br clear="left"/>

## Attribution

If you use this code in a publication, please refer to the package by name and cite [Y. Kusafuka, K. Obayashi, K. Asano, & R. Yamazaki (in prep)]().

The Magnetic Bullet model used in this code is successfully applied several afterglows, some of which can be found in the following works. 
 - GRB 221009A: [Y. Kusafuka & K. Asano (2025)](https://ui.adsabs.harvard.edu/abs/2025arXiv250201437K/abstract)
 - GRB 110213A: [Y. Kusafuka, K. Obayashi, K. Asano, & R. Yamazaki (in prep)]() ...
 <!-- - GRB 080710A: [K. Obayashi, Y. Kusafuka, Y. Sudo, K. Asano, & R. Yamazaki (2025)]() ... -->

## Features

"Magglow" computes leptonic & hadronic emission from both forward & reverse shocks of a relativistic magnetized jet based on Magnetic Bullet model.  

It includes:
- Relativistic forward and reverse shock evolution 
- Any stratified density medium
- Equal arrival time surface
- Arbitrary viewing angles
- leptonic photon emission (synchrotron/SSC) 
- hadronic neutrino production (pp collision/photo-meson interaction)
- Internal absorption processes (SSA/gg absorption)
- Klien-Nishina effects

It does not include (yet):
- Lateral spreading
- Structured jet
- Sedov-Taylor solution
- Polarization
- EBL absorption
<!-- - Gravitational lensing -->

## Installation


If you are working on a local copy of this repo and would like to install from source, you can the run the following from the top level directory of the project.
```bash
$ git clone https://github.com/yo3-sun/Magglow.git 
```

## Using

In your Julia code, import the library with 
```
include("PATH to Magglow.jl")
using .Magglow  
```

The main function of interest is `MagneticBulletAfterglow!(z,DL,t_in,nu_in,Input,Output)`. 

See `example/LC_sample.jl` for a simple example.

`z` is a source redshift.

`DL` is a source luminosity distance. 

`t_in` is an array of observed time.  

`nu_in` is an array of observed frequency.

`Input` has 20 arguments:
- `1 E0` isotropic equivalent energy in erg
- `2 thetaC` half-width of the jet core in radians (jetType specific)
- `3 thetaW` "wing" truncation angle of the jet, in radians
- `4 b` power for power-law structure, &theta;<sup>-b</sup>
- `5 L0` Fiducial luminosity for energy injection, in erg/s, typically 0.
- `6 q` Temporal power-law index for energy injection, typically 0.
- `7 ts` Fiducial time-scale for energy injection, in seconds, typically 0.
- `8 n0` Number density of ISM, in cm<sup>-3</sup>
- `9 p` Electron distribution power-law index (p>2)
- `10 epsilon_e` Thermal energy fraction in electrons
- `11 epsilon_B` Thermal energy fraction in magnetic field
- `12 xi_N` Fraction of electrons that get accelerated
- `13 d_L` Luminosity distance in cm
- `14 d_L` Luminosity distance in cm
- `15 d_L` Luminosity distance in cm
- `16 d_L` Luminosity distance in cm
- `17 d_L` Luminosity distance in cm
- `18 d_L` Luminosity distance in cm
- `19 d_L` Luminosity distance in cm
- `20 d_L` Luminosity distance in cm

`Output` is an array of output observed flux:
- `1  FS e-synchrotron`    Forward shock
- `2  FS p-synchrotron`    Forward shock
- `3  FS e-SSC`            Forward shock
- `4  RS e-synchrotron`    Reverse shock
- `5  RS p-synchrotron`    Reverse shock
- `6  RS e-SSC`            Reverse shock
- `7  FS pp e-neutrino`    Forward shock
- `8  FS pp mu-neutrino`   Forward shock
- `9  FS pp pi0 gamma`     Forward shock
- `10 FS pg e-neutrino`   Forward shock
- `11 FS pg mu-neutrino`  Forward shock
- `12 FS pg pi0 gamma`    Forward shock
- `13 RS pp e-neutrino`   Reverse shock
- `14 RS pp mu-neutrino`  Reverse shock
- `15 RS pp pi0 gamma`    Reverse shock
- `16 RS pg e-neutrino`   Reverse shock
- `17 RS pg mu-neutrino`  Reverse shock
- `18 RS pg pi0 gamma`    Reverse shock

5 arguments of array `is_calc` expresses the radiation processes: 
- `1 electron synchrotron`
- `2 proton synchrotron`
- `3 electron SSC`
- `4 pp collision`
- `5 photo-meson interaction`
