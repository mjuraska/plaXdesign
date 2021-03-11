## plaXdesign: Assessing vaccine durability and harm following placebo crossover in two-arm vaccine efficacy trials
`plaXdesign` is an R package implementing sample size and power calculations presented in Follmann D. et al. (2021), Assessing durability of vaccine effect following blinded crossover in COVID-19 vaccine efficacy trials.

To install the package, first install the package `devtools` by running `install.packages("devtools")` and then run 
```
devtools::install_github("mjuraska/plaXdesign").
```

To reproduce Table 2 in Follmann et al. (2021) using the analytical approach, run
```
library(plaXdesign)
tab <- plaXpower(theta1=c(rep(200, 4), rep(25, 4)),
                 theta2=c(rep(200, 2), rep(100, 2), rep(25, 2), rep(12, 2)),
                 ve1=c(rep(0.9, 4), rep(0.5, 4)), 
                 ve2=c(0.75, 0.9, 0.75, 0.9, -1, -3, -1, -3))
tab

# for 3 decimal places
print(tab, digits=3)
```
For the simulation approach, run
```
plaXpower(theta1=c(rep(200, 4), rep(25, 4)),
          theta2=c(rep(200, 2), rep(100, 2), rep(25, 2), rep(12, 2)),
          ve1=c(rep(0.9, 4), rep(0.5, 4)), 
          ve2=c(0.75, 0.9, 0.75, 0.9, -1, -3, -1, -3),
          method="simulation",
          seed=42)
```
