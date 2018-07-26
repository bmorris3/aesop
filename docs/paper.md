---
title: 'aesop: ARC Echelle Spectroscopic Observation Pipeline'
tags:
  - Python
  - astronomy
  - spectroscopy
authors:
  - name: Brett M. Morris
    orcid: 0000-0003-2528-3409
    affiliation: 1
  - name: Trevor Dorn-Wallenstein
    orcid: 0000-0003-3601-3180
    affiliation: 1
affiliations:
 - name: Astronomy Department, University of Washington, Seattle, WA, USA
   index: 1
date: 15 July 2018
bibliography: paper.bib
---

# Summary

The ARC Echelle Spectroscopic Observation Pipeline, or ``aesop``, is a high resolution 
spectroscopy software toolkit tailored for observations from the Astrophysics Research 
Consortium (ARC) Echelle Spectrograph mounted on the ARC 3.5 m Telescope at Apache 
Point Observatory. ``aesop`` picks up where the traditional IRAF reduction scripts leave 
off, offering an open development, object-oriented Pythonic analysis framework for echelle
spectra. 

Basic functionality of ``aesop`` includes: (1) blaze function normalization by polynomial 
fits to observations of early-type stars, (2) an additional/alternative robust least-squares 
normalization method, (3) radial velocity measurements (or offset removals) via 
cross-correlation with model spectra, including barycentric radial velocity calculations, 
(4) concatenation of multiple echelle orders into a simple 1D spectrum, and (5) approximate
flux calibration. 

Some handy additional utilities include methods for download PHOENIX model spectra
[@Husser:2013], and methods for measuring the Mount Wilson Observatory-calibrated 
CaII H & K "S" index [@Morris:2017]. aesop was built from the Astropy 
package-template, and thus includes self-building documentation and continuous integration
[@astropy:2018].

# Acknowledgements

We acknowledge guidance from Suzanne L. Hawley and Emily Levesque, and the invaluable 
framework and dev team behind the astropy package-template.

# References
