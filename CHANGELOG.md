# Changelog
All notable changes to this project will be documented in this file.


## [v0.1.0]

First release of the nanobody pre-processing pipeline!! :grin:

## [v0.2.0]

NanoLogix has been updated to include a (beta) reporting module! Currently only works with singularity (docker support coming soon)

## [v0.3.0]

Name changed to alpseq! Full docker support added also. Reporting module has been updated to include cladograms (trees) of the clones

## [v0.4.0]

Updates to clustering scheme (threshold is now 80%), added support for [matchbox](https://github.com/jakob-schuster/matchbox) based annotation of nanobodies (faster)

## [v0.5.0]

Fixed some bugs (broken tag for process, missing docker container for matchbox, error in clustering when the number of rounds does not equal 2)
Added some new features:
- parameter validation
- FWR4 region (for matchbox annotation) can now be set as a parameter
- the error rate used for matchbox annotation can now be set as a parameter
- added in-frame check for matchbox annotated sequences

