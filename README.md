**GLoBESfit**
=============
***A `GLoBES`-based fit to the global reactor antineutrino dataset***

Overview
========
`GLoBESfit` is a software suite used to analyze the global reactor antineutrino dataset. As the name suggests, `GLoBESfit` is based on the publicly available `GLoBES` framework – this is all the user needs to run `GLoBESfit`!

By default, `GLoBESfit` is designed to examine the evidence for the existence of one additional, sterile neutrino species. However, the user to free to modify the expressions for the oscillation probabilities to probe any new-physics scenario that they wish to investigate

`GLoBESfit` consists of two separate executables:
  1. `GLoBESfit_rate` executes an analysis of measurements of the inverse beta decay (IBD) *rate measurements* performed at the following experiments: Bugey-3, Bugey-4, Chooz, Daya Bay, Double Chooz, Gösgen, ILL, Krasnoyarsk (1987, 1994, 1999), Nucifer, Palo Verde, RENO, Rovno (1988, 1991) and Savannah River.
  2. `GLoBESfit_spectra` executes an analysis of measurements of IBD *spectral measurements* performed at the following experiments: Bugey-3, DANSS, Daya Bay, Double Chooz, NEOS and RENO.

See the documentation for more details on the experiments we've included and how we perform these analyses!

Documentation
=============

`GLoBESfit` will be extensively documented in an upcoming publication - check back soon for an arXiv number!

Installation
============

`GLoBESfit` can be downloaded from source:
```bash
git clone https://github.com/JMBerryman/GLoBESfit.git
cd GLoBESfit
```

A `Makefile` has been included as a part of `GLoBESfit`. The only user-required input here is the location of their `GLoBES` directory. By default, the relevant line of the `Makefile` reads:
```bash
prefix = /home/jeffb17/Research/globes-3.2.17.0
```
Once the user has changed this line as appropriate, they need simply run
```bash
make
```
At this point, the user can execute the rate analysis with
```bash
./GLoBESfit_rate
```
and the spectral-ratio analysis can be executed with
```bash
./GLoBESfit_spectra
```
It's as easy as that! (The user should review the documentation for explanations on how to use or ignore different experiments within a given analysis, how to use different reactor antineutrino fluxes, etc. This isn't available quite yet, but we're getting there!)

Referencing `GLoBESfit`
=======================

We hope you find `GLoBESfit` useful! If you do use any of these tools as a part of your own research, then we ask that you reference the following papers:

```bash
@article{Berryman:2019hme,
      author         = "Berryman, Jeffrey M. and Huber, Patrick",
      title          = "{Reevaluating Reactor Antineutrino Anomalies with Updated Flux Predictions}",
      journal        = "Phys. Rev.",
      volume         = "D101",
      year           = "2020",
      number         = "1",
      pages          = "015008",
      doi            = "10.1103/PhysRevD.101.015008",
      eprint         = "1909.09267",
      archivePrefix  = "arXiv",
      primaryClass   = "hep-ph",
      SLACcitation   = "%%CITATION = ARXIV:1909.09267;%%"
}
```

(More references will be added soon!)
