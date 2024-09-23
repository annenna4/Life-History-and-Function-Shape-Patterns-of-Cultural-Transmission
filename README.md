# Life History and Functionality in Cultural Transmission: Insight from Mismatches between Neutral Theory and Data

## Paper
This repository contains code used for both versions of the age-structured neutral model described in the manuscript:

'Life History and Functionality in Cultural Transmission: Insight from Mismatches between Neutral Theory and Data' by Anne Kandler,  Rafael D'Andrea, James O'Dwyer

> *Abstract*: While models of neutral evolution have been broadly applied in understanding cultural dynamics, they provide an inadequate description of one of the archetypal examples of cultural transmission: patterns of human first name choice. Guided by the multiple mismatches between neutral theory and first names, we propose several large scale mechanisms that could plausibly explain the deviations, and thus shed light on the mechanisms driving three key patterns. First, we show that introducing a non-trivial life history of cultural variants, induced by an age-constrained transmission process, produces instantaneous variant abundances that deviate from standard neutral predictions but coincide with the observed power law behaviour of the variant abundance distribution in US census data.  Second, we show that thresholded first name records, where all names not exceeding a certain count threshold have been removed, are well-described by an effective neutral model characterising the behaviour of all but the rarest cultural variants. As a result, we can identify an effective innovation rate for a given dataset.  Third, applying this effective neutral theory to multiple first name data sets reveals a novel scaling relationship: the larger the total number of births per unit time (which serves as a proxy for population size) the lower is the fitted per-capita effective innovation rate. Introducing an ‘anti-dominance’ bias through a local-global tradeoff allowed us to replicate this relationship, suggesting the importance of the functionality of a first name choice implying mechanisms avoiding dominance by a single variant may be at play in cultural transmission.

## Code

All code is implemented in Matlab. The script `main_ageSim.m` can be used to
generate populations using the simulation models described in the paper. 


## License
[![CC BY-SA 4.0][cc-by-sa-shield]][cc-by-sa]

This work is licensed under a
[Creative Commons Attribution-ShareAlike 4.0 International License][cc-by-sa].

[![CC BY-SA 4.0][cc-by-sa-image]][cc-by-sa]

[cc-by-sa]: http://creativecommons.org/licenses/by-sa/4.0/
[cc-by-sa-image]: https://licensebuttons.net/l/by-sa/4.0/88x31.png
[cc-by-sa-shield]: https://img.shields.io/badge/License-CC%20BY--SA%204.0-lightgrey.svg

