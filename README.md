# Post-holiday outpatient visits, 2024

R code for analysing daily outpatient visits during the first week after public holidays. The script includes descriptive statistics, quasi-Poisson regression, model diagnostics, sensitivity analyses and a two-panel figure.

## Usage

1. Download the data from https://doi.org/10.17632/nd2dtwh42y.2.
2. Place `outpatient_visits_2024.csv` and `analysis.R` in the same folder.
3. Set that folder as the R working directory and run:

```r
source("analysis.R")
plot_main()
```

Results appear in the Console and the figure in the graphics window. Tested with R 4.4.2. No additional packages are required.
