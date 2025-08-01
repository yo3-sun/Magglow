# "Magglow" - Magnetic Bullet Afterglow 

<img align="left" width="422" height="321" alt="Image" src="https://github.com/user-attachments/assets/d63eb294-e4f5-47cd-b972-8359382bf4f8" />

**Magglow** is a Julia module to calculate Gamma-ray Burst (GRB) afterglow light curves and spectra based on Magnetic Bullet model [Y. Kusafuka & K. Asano (2024)](https://ui.adsabs.harvard.edu/abs/2025MNRAS.536.1822K/abstract). 
The framework implements semi-analytic models for forward and reverse shock dynamics propagating in a stratified CSM and leptonic (synchrotron with self-absorption and self-Compton scattering with Klein-Nishina corrections) and hadronic (pp colisions and photomeson interactions) multimessenger emission mechnisms observed from arbitrary viewing angle.

Details of the methods can be found in [Y. Kusafuka, K. Obayashi, K. Asano, & R. Yamazaki (in prep)](). 
<!-- This code is under active development.  -->

<!-- Documentation is available at <https://afterglowpy.readthedocs.io/> -->

<br clear="left"/>

The Magnetic Bullet model used in this code is successfully applied several afterglows, some of which can be found in the following works. 
 - GRB 221009A: [Y. Kusafuka & K. Asano (2025)](https://ui.adsabs.harvard.edu/abs/2025arXiv250201437K/abstract)
 - GRB 110213A: [Y. Kusafuka, K. Obayashi, K. Asano, & R. Yamazaki (in prep)]() ...
 <!-- - GRB 080710A: [K. Obayashi, Y. Kusafuka, Y. Sudo, K. Asano, & R. Yamazaki (2025)]() ... -->

## Attribution

If you use this code in your research, please cite the relevant papers:
- [Y. Kusafuka & K. Asano (2024)](https://ui.adsabs.harvard.edu/abs/2025MNRAS.536.1822K/abstract)
- [Y. Kusafuka, K. Obayashi, K. Asano, & R. Yamazaki (in prep)]().

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

## Usage

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

`Input` has 20 arguments: <details>
<summary>20 arguments <i>(click to expand/collapse)</i></summary>
<br>

- `1 E0` isotropic equivalent energy in erg
- `2 G0` initial Loretnz factor 
- `3 S0` initial magnetization 
- `4 D0` initial thickness, in cm
- `5 n0` CSM number density, in cm<sup>-3</sup>
- `6 k`  density slope (ISM:0 --- 2:wind)
- `7 ts` Fiducial time-scale for energy injection, in seconds, typically 0.
- `8 n0` Number density of ISM, in cm<sup>-3</sup>
- `9 epsiron_e` thermal energy fraction in electrons for FS
- `10 epsilon_B` thermal energy fraction in magnetic field for FS
- `11 p` particle spectral index for FS (p>2)
- `12 fe` fraction of electrons that get accelerated for FS
- `13 theta_j` opening angle, in rad
- `14 theta_o` viewing angle, in rad
- `15 epsiron_e,RS` thermal energy fraction in electrons for RS
- `16 epsilon_B,RS` thermal energy fraction in magnetic field for RS
- `17 p_RS` particle spectral index for RS (p>2)
- `18 fe_RS` fraction of electrons that get accelerated for RS
- `19 epsiron_p,RS` thermal energy fraction in protons for RS
- `20 fp_RS` fraction of protons that get accelerated for RS

</details>

`Output` is an array of output observed flux:
- `1  FS e-synchrotron`    
- `2  FS p-synchrotron`   
- `3  FS e-SSC`           
- `4  RS e-synchrotron`   
- `5  RS p-synchrotron`  
- `6  RS e-SSC`            
- `7  FS pp e-neutrino`    
- `8  FS pp mu-neutrino` 
- `9  FS pp pi0 gamma`   
- `10 FS pg e-neutrino`  
- `11 FS pg mu-neutrino` 
- `12 FS pg pi0 gamma`    
- `13 RS pp e-neutrino`  
- `14 RS pp mu-neutrino`  
- `15 RS pp pi0 gamma`   
- `16 RS pg e-neutrino`   
- `17 RS pg mu-neutrino`  
- `18 RS pg pi0 gamma`    

5 arguments of array `is_calc` expresses the radiation processes: 
- `1 electron synchrotron`
- `2 proton synchrotron`
- `3 electron SSC`
- `4 pp collision`
- `5 photo-meson interaction`
