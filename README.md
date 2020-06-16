MethICA
========

Methylation signature analysis with independent component analysis.

Description: apply independent component analysis (ICA) to methylation data to extract epigenetic signatures and perform downstream analyses and representations to interpret the biological meaning of methylation components.

## Installation
========
Install from the GitHub repository using devtools:

install.packages("devtools")
library(devtools)
devtools::install_github("FunGeST/MethICA")

## Input files
========
For the extraction of methylation signature _MethICA_ requires one mandatory input file -- a **methylation level** file (bval) describing methylation levels of individual CpG sites or genomic regions (row) in each sample (column). Two additional files are needed for the interpretation of methylation components: a **CpG annotation** file describing the epigenetic context of CpG sites, and a **sample annotation** providing any relevant information about the samples.

**The input files should have the following format. Example input files are provided with the package.**

`1]. bval `: __methylation level data__

* `Rownames` : CpG IDs
* `Colnames` : Sample IDs

**Optional: (the header is required, but the order of the columns can change)**

`2]. CpG annotation`:


* `TargetID` : CpG ID. Any alphanumeric string.
* `state` : Chromatin state (active/inactive)
* `domain` : Chromatin domain (ROADMAP nomenclature)
* `CpG_context` : Methylation domain (HMD/PMD/LMR/UMR)
* `cgi_feature` : CpG-island based feature (Island, Shore, Shelf, out of cgi)
* `gene_feature` : Gene-based feature (Body, TSS500, out)
* `gene_name` : Gene name
* `decile` : Replication timing decile (1: early --> 10: late)
* `CHROM` : Chromosome. Between chr1 and chr22 or the chrX or chrY ('chr' prefix required).
* `MAPINFO` : CpG position (CHROM and MAPINFO are only required to annotate the CpG table with epigenomic features using chromatin.feature function)


`3]. annot_data`: __sample annotation data__

* `Sample`: Sample identifier. Any alphanumeric string.
* add in column the annotations you are interested in testing

## Running MethICA
================

* The [RUNNING\_METHICA\_EXAMPLE](https://github.com/FunGeST/MethICA) folder contains an example dataset and an R script of a typical MethICA analysis. Please try!</br>
* [*Introduction to MethICA*] provides a comprehensive example of the MethICA workflow with detailed  explanations of each function.</br> 


## License
========

Copyright (C) 2020 Léa Meunier

MethICA is a free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
