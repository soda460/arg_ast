

# arg_ast : R code to visualize phenotype-genotype associations between resistance genes identified by whole genome sequencing and antibiotic susceptibility testing datasets.


<p align="left"><img src="misc/AUG.png" alt="Phenotype Genotype associations with antibiotic resistance data" width="1000"></p>



## Installation

```shell
git clone https://github.com/soda460/arg_ast
```

If need install these R libraries

```R
install.packages('tidyverse')
install.packages('ggplot2')
install.packages('ggrepel')

# optional
install.packages('svglite')
```

## Usage

Try to fit your phenotype and genotype data to the format of the csv files loaded by the main program : arg_ast.R

The script data_wrangling.R is highly specific to our datasets.


## Contributors

  * Jean-Simon Brouard, Ph.D.  
Biologist in bioinformatics, Science and Technology Branch  
Agriculture and Agri-Food Canada / Government of Canada  
jean-simon.brouard@canada.ca


  * Dominic Poulin-Laprade, Ph.D.  
Research Scientist, Science and Technology Branch  
Agriculture and Agri-Food Canada / Government of Canada  
Dominic.Poulin-Laprade@Canada.ca


## Notes

Technically, this script is also useful for me because it make use of quosure objects in R functions.



