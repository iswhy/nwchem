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
The current branch of the NWChem quantum chemistry program contains the implementation of the *weak-field asymptotic theory (WFAT) in the integral representation (IR)* as a module. WFAT is a rigorous method for the simulation of molecular tunneling ionization due to static fields. Our current implementation of WFAT as a NWChem module supports one-electron WFAT in the leading-order approximation (OE-WFAT(0)), one-electron WFAT including the first-order correction (OE-WFAT(1)), and many-electron WFAT in the leading-order approximation (ME-WFAT(0)). Typical applications of WFAT include the calculation of orientation-dependent strong-field ionization rates of a molecule and the associated transverse electron momentum distributions in the adiabatic regime, that is, when the laser's oscillation period is considerably longer than the characteristic time scale of the dynamics. A tutorial on how to use this module as well as the definition of the input directives are given in the following.

***Note***:
At the moment, the WFAT module can only be used for states or orbitals obtained from a DFT calculation using the `nwxc` DFT library. Extending the implementation to allow use of the `libxc` library and to be independent of DFT, hence allowing more general wave function types, is in our plan.



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

The `wfat` block is where WFAT-related parameters should be defined. `beta` and `gamma` specify the $(\beta,\gamma)$ pairs of Euler orientation angles. The number of parabolic quantum channels used for calculating the total rate is specified through `paraqn`. The options specified for the `print` directive instruct the module to print the structure factor (`srfc`) and ionization rate (`rate`). Since the HOMO of methyl bromide is doubly denegerate, one specifies `a 0 -1` for `movec_id`. `a`, `0`, and `-1` stand for spin alpha channel, the first HOMO, and the second HOMO, respectively. Note that the MOs to be ionized are identified so that their ID starts from zero for the HOMO, -1 for the next MO with a lower (or the same) orbital energy, and so on. Our implementation of WFAT in the IR provides both the partial-wave and non-partial-wave formulations. The input above specifies a partial-wave OE-WFAT calculation where the maximum angular momentum is specified as `lmax 10`. The external field can specified through the `field` input block. In the example above, the field is static (DC field) and its magnitude is 0.004 a.u. To specify that the simulation should be run with the OE-WFAT algorithm, use the `oe` option for the `task wfat` directive. The more detailed description of each of these input directives can be found [below](wfat-input-directives).

The first input example above will run a OE-WFAT calculation in the so-called leading-order approximation, which essentially means that the field-dependence in the method is truncated up to the zeroth order. If more accuracy in terms of field-dependence of the ionization rate is desired, one should specify the approximation level to be the first order. A snippet of the `wfat` block below exemplifies how to use the first-order algorithm.
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
In addition to the directives present in the first example, three additional directive have been added: `wfat1`, `wfat1_exclude`, and `degenthr`. The directive `wfat1` must be present in order to run the simulation using OE-WFAT with the first-order correction (OE-WFAT(1)). The remaining two additional directives are optional for running OE-WFAT(1) simulations, nevertheless, they almost always need to be set. `wfat1_exclude` controls which channels (defined as the combination of an ionized orbital and a parabolic quantum number) are to be treated using OE-WFAT(0) (the leading-order approximation). In order to obtain calculation results that follow a systematic progression in the order of the field strength, one needs to choose the approperiate level of approximation for a particular parabolic channel given the order of the field in the field dependence of that parabolic channel. In particular, if the dominant channel $(0,0)$ is treated using OE-WFAT(1) and channels $(0,\pm1)$ are to be included, then the latter channels should be treated using OE-WFAT(0), see [trinh2015](trinh2015) and [wahyutama2025](wahyutama2025). The channels to be treated using OE-WFAT(0) (excluded from OE-WFAT(1)) are specified as a string of four characters, such as `a -1 0 1`. The first two of the four-character sequence specify the spin channel and the ID of the ionized orbital, while the last two specify the $n_\xi$ and the $m$ components of the parabolic quantum number. `degenthr` sets a threshold to detect degenerate orbitals. It is sometimes necessary to set `degenthr` to a value different from its default because the energy difference between degenerate states in a molecule can vary by more than one orders of magnitude from the energy difference in another molecule. See the [WFAT Input Directives Section below](wfat-input-directives) for the more detailed definitions of these input directives.



### Ionization calculation using many-electron WFAT (ME-WFAT)
In the previous section, two examples on how to run OE-WFAT(0) and OE-WFAT(1) simulations are presented. Our implementation of WFAT also includes a capability to treat molecular tunneling ionization directly from the exact multi-electron Hamiltonian without approximating the molecule as a single-active-electron system as in OE-WFAT. While ME-WFAT results in more accurate ionization rates from weakly-correlated states than OE-WFAT, it is most suitably used to calculate tunneling ionization rates from highly-correlated states.

***Note***:
Although theoretically, ME-WFAT can treat a vast range of wave function structures, our current implementation is still limited to single-determinantal wave functions, such as Hartree-Fock or DFT wave functions.

The following snippet is an example of input deck to run ME-WFAT(0) calculation on methyl bromide in which the ionization channel starts from the DFT ground state of the neutral to the DFT ground state of the cation.
```
start
title "CH3Br ME-WFAT(0) - the first degenerate rate"
echo

######################################################
geometry noautosym
 Br   0.000   0.000   1.934
  C   0.000   0.000   0.000
  H   1.032   0.000  -0.333
  H  -0.516   0.894  -0.333
  H  -0.516  -0.894  -0.333
end

basis spherical
  * library aug-cc-pvtz
end
######################################################

set dft:gonwxc .true.

######################################################
charge 0
dft
  maxiter 500
  mult 1
  xc new xcampbe96 1.0 cpbe96 1.0         #HFexch 1.0
  cam 0.30 cam_alpha 0.70 cam_beta 0.30
  vectors output "CH3Br.ion0.movecs"
end

task dft
######################################################

#####################################################
charge 1
dft
  maxiter 500
  mult 2
  xc new xcampbe96 1.0 cpbe96 1.0         #HFexch 1.0
  cam 0.30 cam_alpha 0.70 cam_beta 0.30
  vectors input "CH3Br.ion0.movecs" output "CH3Br.ion1.movecs"
end

task dft
######################################################

######################################################
basis "wfat_rsbas" spherical
  x  s
     1.00  1.00
end

wfat
  print srfc yield
  beta   0.0 121 180.0
  gamma  0.0 151 360.0
  paraqn  "n0" 0 -1 +1
  field
    type dc
    max 0.012
  end
  si_type dyson
end
task wfat me
######################################################
```
ME-WFAT calculations require two wave functions, that of the neutral and that of the cation. That is why there are two DFT calculations in the input deck above indicated by the charge settings `charge 0` and `charge 1`. The ME-WFAT code reads the single-determinantal wave functions of the neutral and cation through their respective occupied molecular orbitals. It searches for files named `<prefix>.ion0.movecs` and `<prefix>.ion1.movecs` in the current directory for the occupied MOs of the neutral and cation, respectively. The line `vectors output "CH3Br.ion0.movecs"` ensures that the file containing neutral's MOs can be recognized by the WFAT module. The MO saving specification for the cation, namely, `vectors input "CH3Br.ion0.movecs" output "CH3Br.ion1.movecs"` is slightly different due to the presence of the input (initial guess) MOs specification. While this is optional, it helps in cases where the cation's ground state is degenerate, such as for methyl bromide cation by fixing a particular cation ground state resulting from this simulation. Without the input option, different but still degenerate cation ground state may result when this simulation is rerun. In the example above, we see a new directive, `si_type dyson`. This directive is used to control the type of self-interaction correction in the XC potential, see [wahyutama2022](wahyutama2022).

Note that the neutral and cation must have the same geometry, hence the vertical/Franck-Condon ionization is assumed. This example also illustrates one quirk of the `nwxc` XC library, namely, when a range-separated XC functional is used (indicated by the line `cam ... cam_alpha ... cam_beta ...`), the Hartree-Fock exact exchange component (`HFexch`) must not be specified. Furthermore, when the range-separated functional is used, the WFAT module requires that an additional basis be defined (see the `basis` input block in the input deck above). To instruct the module to run the ME-WFAT(0) simulation, use the `me` option for the `task wfat` directive. 

As mentioned above, the ground state of methyl bromide cation is degenerate, doubly degenerate to be more precise. To calculate tunneling ionization rates for the ionization channel that ends up in the second degenerate state, use this line
```
vectors input "CH3Br.ion0.movecs" swap beta 21 22 output "CH3Br.ion1.movecs"
```
instead of the corresponding line in the previous input deck example. In the previous example, the cation ground state is based on the removal of the 22-nd MO of the neutral. By using the option `swap beta 21 22`, NWChem will remove the 21-st MO to calculate the second degenerate ground state of the cation. The 21-st and 22-nd MOs of the neutral are degenerate.



### Input directives for the WFAT module



## References
