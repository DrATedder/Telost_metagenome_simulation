# Telost_metagenome_simulation
A series of scripts developed for an MSc project to (semi-)realistically simulate Telost metagenomes under *A. salmonicida* infection.

## Overview
***
Conceptually, this project is realtively straightforward - pull down some likely microbiome OTUs, and then simulate some metagenomes under an infection scenario using `InSilicoSeq` [see here](https://github.com/HadrienG/InSilicoSeq). The reality is a little different (althought I do achieve this simplistic aim), with more hardcoding of OTU names (for downloading genomes from `NCBI`) than I would have liked (it turns out, that there isn't a way to group targets based on likelihood of being part of the telost species microbiome...).

There are two principle scripts:

1) `pull_OTU_genomes.py` - Uses `ncbi.datasets` to pull down OTUs meeting a pre-determined 'likely component' status (i.e. what limited evidence we have suggests these may actually be microbiome components of one of the three telost species).

2) `simulate_metagenomes.py` - This script generates `abundance` files, and then shell scripts (one per telost species) to run `InSilicoSeq` and generate simulated reads.

### Rules & logic
***

This sounds pretty simple, but it requires quite a few biological and logic decisions that probably warrant some consideration:

* Generates 20 control & 20 infected samples per telost species (can be edited; balanced sample design not essential)
* Each sample has endogenous DNA contamination between 0.5-2%
* Pathogen abundance in 'infected' samples
  ```
  13% to 37% Aeromonas salmonicida
  ```
  **Note**. A [triangular distribution](https://en.wikipedia.org/wiki/Triangular_distribution) used for realism.
* Each telost species a defined 'core' microbiome:

```
control samples: ~80% of remaining reads
infected samples: ~55% of remaining reads (reduced due to dysbiosis)
  ```
  **Note**. Core abundances are generated using a [Dirichlet distribution](https://en.wikipedia.org/wiki/Dirichlet_distribution) to introduce biological variability.

* Remaining abundance considered 'background':
  - randomly sampled from non-ore genomes
  - higher diversity in 'control' samples than 'infected' samples.

 ### Prerequisites (likely not exhaustive)

 ```python
random
zipfile
pathlib
ncbi.datasets
subprocess
Bio

```

