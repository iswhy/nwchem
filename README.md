<!---[![Build Status](https://img.shields.io/endpoint.svg?url=https%3A%2F%2Factions-badge.atrox.dev%2Fnwchemgit%2Fnwchem%2Fbadge%3Fref%3Dmaster&style=flat)](https://actions-badge.atrox.dev/nwchemgit/nwchem/goto?ref=master)--->
[![NWChem_CI](https://img.shields.io/github/actions/workflow/status/nwchemgit/nwchem/github_actions.yml?style=plastic)](https://github.com/nwchemgit/nwchem/actions)   
<!---[![Build Travis-CI Status](https://img.shields.io/travis/nwchemgit/nwchem.svg)](https://travis-ci.org/nwchemgit/nwchem/builds)--->
[![License](https://img.shields.io/badge/license-ECL2-blue.svg)](https://raw.githubusercontent.com/nwchemgit/nwchem/master/LICENSE.md)  
![GitHub top language](https://img.shields.io/github/languages/top/nwchemgit/nwchem)
![Github Downloads All Releases](https://img.shields.io/github/downloads/nwchemgit/nwchem/total)   
![Homebrew Downloads](https://img.shields.io/homebrew/installs/dy/nwchem?label=Homebrew%20downloads)![Conda-Forge Downloads](https://img.shields.io/conda/dn/conda-forge/nwchem?label=Conda-forge%20downloads)![Docker Pulls](https://img.shields.io/docker/pulls/nwchemorg/nwchem-700.fedora.sockets?label=Docker%20pulls)    
![Github Downloads Release 6.8](https://img.shields.io/github/downloads/nwchemgit/nwchem/v6.8-release/total)
![Github Downloads Release 6.8.1](https://img.shields.io/github/downloads/nwchemgit/nwchem/6.8.1-release/total)  
![Github Downloads Release 7.0.0](https://img.shields.io/github/downloads/nwchemgit/nwchem/v7.0.0-release/total)
![Github Downloads Release 7.0.2](https://img.shields.io/github/downloads/nwchemgit/nwchem/v7.0.2-release/total)  
![Github Downloads Release 7.2.0](https://img.shields.io/github/downloads/nwchemgit/nwchem/v7.2.0-release/total)
![Github Downloads Release 7.2.1](https://img.shields.io/github/downloads/nwchemgit/nwchem/v7.2.1-release/total)  
![Github Downloads Release 7.2.2](https://img.shields.io/github/downloads/nwchemgit/nwchem/v7.2.2-release/total)
![Github Downloads Release 7.2.3](https://img.shields.io/github/downloads/nwchemgit/nwchem/v7.2.3-release/total)


[![All Releases tags](https://img.shields.io/github/release/nwchemgit/nwchem/all.svg)](https://github.com/nwchemgit/nwchem/releases)
![GitHub Release Date](https://img.shields.io/github/release-date/nwchemgit/nwchem)
[![EMSL](https://rawgit.com/nwchemgit/nwchem/master/contrib/git.nwchem/emsl_logo2.svg)](https://www.emsl.pnl.gov)
 
## NWChem: Open Source High-Performance Computational Chemistry
<img alt="NWChem logo" src="https://raw.githubusercontent.com/nwchemgit/nwchem/master/contrib/git.nwchem/MS3_logo_cropped.png" align=middle width="102pt" height="117pt"/>  
NWChem is actively developed by a consortium of developers and maintained by the Environmental Molecular Sciences Laboratory (EMSL), a US DOE Office of Science User Facility
located at the Pacific Northwest National Laboratory (PNNL) in Washington State

### Scalable computational chemistry tools

Source download instructions for developers

* Master
```
git clone https://github.com/nwchemgit/nwchem.git
```
* 7.2.3 Release
```
git clone  -b hotfix/release-7-2-0 https://github.com/nwchemgit/nwchem.git nwchem-7.2.3
```
### Release 7.2.3

Release files available at the page  
https://github.com/nwchemgit/nwchem/releases/tag/v7.2.3-release    
Tar files for download:  
https://github.com/nwchemgit/nwchem/releases/download/v7.2.3-release/nwchem-7.2.3-release.revision-d690e065-srconly.2024-08-27.tar.bz2  
https://github.com/nwchemgit/nwchem/releases/download/v7.2.3-release/nwchem-7.2.3-release.revision-d690e065-nonsrconly.2024-08-27.tar.bz2  
Release notes available at  
https://github.com/nwchemgit/nwchem/blob/master/release.notes.7.2.3.md  

### Documentation
https://nwchemgit.github.io



## WFAT Module
The current branch of the NWChem quantum chemistry program contains the implementation of the *weak-field asymptotic theory (WFAT) in the integral representation (IR)* as a module. WFAT is a method for the simulation of molecular tunneling ionization due to static fields. A quick tutorial on how to use this module as well as the definition of the input directives are given in the following.

***Note:***
At the moment, the WFAT module can only be used for states or orbitals obtained from a DFT calculation using the `nwxc` DFT library. Extending the implementation to allow use of the `libxc` library and to be independent of DFT, hence allowing more general wave function types, is still in our plan.



### Ionization calculation using one-electron WFAT (OE-WFAT)
OE-WFAT is used to calculate tunneling ionization rate from a one-electron function, such as molecular orbitals (MOs). Since, tunneling ionization from a weakly correlated state can be well approximated with an ionization from a particular molecular orbital, one can use OE-WFAT to simulate tunneling ionization in such a situation. The following input deck is an example of OE-WFAT calculation for the HOMO of methyl bromide ionized by a static electric field of magnitude 0.004 a.u. (atomic unit).
```
start
title "CH3Br OE-WFAT(0) calculation"
echo

######################################################
geometry noautosym
   C   0.000   0.000  -1.538
  Br   0.000   0.000   0.424
   H   1.044   0.000  -1.872
   H  -0.522   0.904  -1.872
   H  -0.522  -0.904  -1.872
end

basis spherical
  * library cc-pvdz
end

set dft:gonwxc .true.
dft
  maxiter 500
  mult 1
  XC new b3lyp
end

task dft
######################################################

######################################################
wfat
  beta   0.0 121 180.0
  gamma  0.0 149 360.0
  paraqn "n0" 0 -1 +1
  print   srfc rate
  movec_id  a 0 -1
  lmax 10
  field
    type dc
    max 0.0040 au
  end
end
task wfat oe
######################################################
```
The line `set dft:gonwxc .true.` instructs NWChem to activate the `nwxc` library. The spefication `new` for the `XC` input further tells NWChem to choose the B3LYP XC functional from the `nwxc` library. The option `noautosym` is invoked in the directive `geometry`. It is desirable to prevent NWChem to automatically reorient the molecule during the simulation, which can happen while it is trying to determine the molecule's point group. If the the geometry has been reoriented, the zero-orientation reference for the orientation-dependent ionization rates will be associated with the reoriented geometry, not the one the user has defined in the input deck, which may have been arranged according to the user's preference.

The `wfat` block is where WFAT-related parameters should be defined. `beta` and `gamma` specify the $(\beta,\gamma)$ pairs of Euler orientation angles. The number of parabolic quantum channels used for calculating the total rate is specified through `paraqn`. The options specified for the `print` directive instruct the module to print the structure factor (`srfc`) and ionization rate (`rate`). Since the HOMO of methyl bromide is doubly denegerate, one specifies `a 0 -1` for `movec_id`. `a`, `0`, and `-1` stand for spin alpha channel, the first HOMO, and the second HOMO, respectively. Note that the MOs to be ionized are identified so that their ID starts from zero for the HOMO, -1 for the next MO with a lower (or the same) orbital energy, and so on. Our implementation of WFAT in the IR provides both the partial-wave and non-partial-wave formulations. The input above specifies a partial-wave OE-WFAT calculation where the maximum angular momentum is specified as `lmax 10`. The external field can specified through the `field` input block. In the example above, the field is static (DC field) and its magnitude is 0.004 a.u. To specify that the simulation should be run with the OE-WFAT algorithm, one uses `task wfat oe`. The more detailed description of each of these input directives can be found [below](wfat-input-directives).

The first example above will run a OE-WFAT calculation in the so-called leading-order approximation, which essentially means that the field-dependence in the method is truncated up to the zeroth order. If more accuracy in terms of field-dependence of the ionization rate is desired, one should specify the approximation level to be the first order. A snippet of the `wfat` block below exemplifies how to use the first-order algorithm.
```
wfat
  wfat1
  wfat1_exclude \
    "a  0 0 -1" \
    "a  0 0 +1" \
    "a -1 0 -1" \
    "a -1 0 +1"
  degenthr 1.0e-3
  beta   0.0 121 180.0
  gamma  0.0 149 360.0
  paraqn "n0" 0 -1 +1
  print   srfc rate
  movec_id  a 0 -1
  lmax 14
  field
    type dc
    max 0.004 au
  end
end
task wfat oe
```


### Ionization calculation using many-electron WFAT (ME-WFAT)

