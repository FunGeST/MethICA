Abstract

Installation
========
Install from the GitHub repository using devtools:

    install.packages("devtools")
    library(devtools)
    devtools::install_github("FunGeST/MethICA")

Input files
========
For the extraction of methylation signature _MethICA_ requires one mandatory input file -- a **methylation level** file (bval) describing methylation level of CpG (row) in the samples series (column). For the representation to help the interpretation of methylation signature, a further two files are required -- a **CpG annotation** file with interest CpG information, and a minimal **sample annotation** file.

**The input files should have the following format. Example input files are provided with the package.**

`1]. bval `: __methylation level data__

* rownames : CpG id
* colnames : Samples id

**Optional: (the header is required, but the order of the columns can change)**

`2]. CpG annotation`:

* `CpG`: Sample CpG. Any alphanumeric string.
* `CHROM`: Chromosome. Between chr1 and chr22 or the chrX or chrY ('chr' prefix required).

`3]. annot_data`: __sample annotation data__

* `Sample`: Sample identifier. Any alphanumeric string.

Running MethICA
================

* The [RUNNING_METHICA_EXAMPLE](https://github.com/FunGeST/) folder contains an example dataset and an R script of a typical MethICA analysis using this data. Please try!</br>
* [*Introduction to MethICA*] provides a comprehensive example of the MethICA workflow with detailed  explanations of each function.</br> 


License
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
