
dnelson-idl
===========

* **NOTE: This codebase is no longer maintained, since 2015, and is only historical.**

Main analysis, plotting, and visualization codebase for working with GADGET/AREPO cosmological and idealized simulations of galaxy formation and related topics. Allows the reproduction of the scientific results of:

* [Nelson et al. 2013](http://ui.adsabs.harvard.edu/abs/2013MNRAS.429.3353N) - Moving mesh cosmology: Tracing cosmological gas accretion
* [Nelson et al. 2015a](http://ui.adsabs.harvard.edu/abs/arXiv:1410.5425) - The impact of feedback on cosmological gas accretion
* [Nelson et al. 2015b](http://ui.adsabs.harvard.edu/abs/arXiv:1503.02665) - Zooming in on accretion - I. The structure of halo gas

More recent papers have transitioned to the `python` codebase.


## Installation

On a typical cluster you can run (and/or change your `.bashrc` to include):

```bash
  module load idl
  export IDL_PATH=$IDL_PATH:+~/idl/
  idl
```

Note that this codebase was last used with IDL 8.1.


## Getting Started

### Snapshot Data

To load a particular field/partType from a snapshot:

```idl
IDL> sP = simParams(run='illustris',res=1820,redshift=2.0)

IDL> gas_masses = loadSnapshotSubset(sP=sP,partType='gas',field='mass')
```

If you want the field only for a particular group, then given its offset index and length:

```idl
IDL> sP = simParams(run='illustris',res=1820,redshift=2.0)

IDL> indRange = [offset_index, offset_index + length]
IDL> gas_masses = loadSnapshotSubset(sP=sP,partType='gas',field='mass',indRange=indRange)
```

### Group Catalogs

To load the complete group+subgroup catalog:

```idl
IDL> sP = simParams(run='illustris',res=910,redshift=2.0)

IDL> gc = loadGroupCat(sP=sP,/skipIDs,/skipOffsets)
IDL> print, 'Total number of groups and subgroups: ', gc.nGroupsTot, gc.nSubgroupsTot
```

### Merger Trees

Currently no code to deal with LHaloTree, SubLink, or Rockstar. 


## Other Functionality

* CalcBoxRemap - "A volume and local structure preserving mapping of periodic boxes" (C extension)
* CalcCoolTime - implementation of the Katz, Weinberg, Hernquist (1996) primordial chemical network (C extension)
* CalcHSML - tree-based, periodic smoothing length estimation in 1D, 2D, and 3D (C extension)
* CalcHSMLds - as above, but estimated at the locations of a different pointset (C extension)
* CalcMatch - multi-threaded implementations of `sort()` and `match()` for integer types/IDs (C extension)
* CalcNN - tree-based, periodic nearest neighbor search (C extension)
* CalcSphMap - SPH kernel-based density projection (C extension)
* CalcTHVal - adaptive-size tophat kernel interpolation of scattered point sets in 1D, 2D, and 3D (C extension)
* CalcTHValWt - as above but weighted on an additional property of each point (C extension)
* ICs - generation of idealized test initial conditions including gas spheres, shocktubes
* spirals - analysis of idealized N-body galactic disk simulations of spiral arm formation
* TI - thermal instability tests and analysis
* tracersMC - utilities and analysis for Monte Carlo tracer particle outputs in AREPO
* tracersVel - utilities and analysis for Velocity Field tracer particle outputs in AREPO

