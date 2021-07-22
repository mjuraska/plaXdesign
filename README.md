## plaXdesign: Assessing vaccine durability and harm following placebo crossover in two-arm vaccine efficacy trials
`plaXdesign` is an R package implementing sample size and power calculations presented in [Follmann, Fintzi, Fay et al. (2021), A deferred-vaccination design to assess durability of COVID-19 vaccine effect after the placebo group is vaccinated](https://www.acpjournals.org/doi/10.7326/M20-8149).

To install the package, first install the package `devtools` by running `install.packages("devtools")` and then run 
```
devtools::install_github("mjuraska/plaXdesign").
```

To reproduce Table 2 in Follmann, Fintzi, Fay et al. (2021), run
```
library(plaXdesign)
tab <- plaXpower(theta1=rep(c(200, 25), each=4),
                 theta2=rep(c(200, 100, 25, 12), each=2),
                 ve1=rep(c(0.9, 0.5), each=4), 
                 ve2=c(0.75, 0.9, 0.75, 0.9, -1, -3, -1, -3),
                 method="simulation",
                 iter=100000,
                 seed=42)
tab

# for 3 decimal places
print(tab, digits=3)
```
To replace the simulation with the analytical approach, run
```
plaXpower(theta1=rep(c(200, 25), each=4),
          theta2=rep(c(200, 100, 25, 12), each=2),
          ve1=rep(c(0.9, 0.5), each=4), 
          ve2=c(0.75, 0.9, 0.75, 0.9, -1, -3, -1, -3))
```
