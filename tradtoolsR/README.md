## Installation
![tradtools pipeline](https://img.shields.io/badge/TRAD--Seq-tradtools-brightgreen)
![tradtools version](https://img.shields.io/github/v/release/devin-ak/tradtools?label=version)

Step 1. [Install R](https://cran.r-project.org/)

Step 2. From within R, [install devtools](https://www.r-project.org/nosvn/pandoc/devtools.html):
  
  `install.packages('devtools')`

Step 3. Install the [tradtoolsR package](https://github.com/devin-AK/tradtools) from github. Note that the R package is in a subdirectory of the main tradtools repo.

`devtools::install_github('devin-AK/tradtools', subdir='tradtoolsR')`

Step 4: Load the package!
  
  `library(tradtoolsR)`