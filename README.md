# BayesianPickWinner
The tool calculates power for the Bayesian pick-the-winner design in a randomized phase II clinical trial. It also generates report and statistical plan so clinicians could easily incorporate it into the clinical trial protocol.
(reference: A Bayesian pick-the-winner design in a randomized phase II clinical trial; Chen et al; submitted).

## Features

* The shiny applictaion provides power calculation for the Bayesian pick-the-winner design in a randomized phase II clinical trial and generates a statistical plan to be easily incorporated  into the clinical trial protocol. 

## Reference

* Chen DT, Huang PY, Lin HY, Chiappori AA, Gabrilovich DI, Haura EB, Antonia SJ, Gray JE. A Bayesian pick-the-winner design in a randomized phase II clinical trial. Oncotarget. 2017 Jul 7;8(51):88376-88385. doi: 10.18632/oncotarget.19088. PMID: 29179442; PMCID: PMC5687612.
  
## Installation

Simply run the following from an R console:

```r
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("dungtsa/BayesianPickWinner",force = TRUE)
```

## Getting Started

```r
require("BayesianPickWinner")
bayes.pick.winner()
```
-------------------------------
Snapshot of shiny app: initial 
![snapshot of shiny app: initial](inst/img/shiny1.png)

-------------------------------
Snapshot of shiny app: output
![snapshot of shiny app: output](inst/img/shiny2.png)

-------------------------------
Snapshot of shiny app: output in Word format (through "download" button)
![snapshot of shiny app: output](inst/img/shiny3.png)
