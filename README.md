# spaceNet
**Latent space models for multivariate networks**

Latent space models for binary multivariate networks (multiplex). The general model assumes that the nodes in the multiplex lie in a low-dimensional latent space. The probability of two nodes being connected is inversely related to their distance in this latent space: nodes close in the space are more likely to be linked, while nodes that are far apart are less likely to be connected. The model is defined in a hierarchical Bayesian framework and estimation is carried out via MCMC algorithm.

<br>
To install the development version from GitHub:

```
# install.packages("devtools")
devtools::install_github("michaelfop/spaceNet")
```

<br>

### References

D'Angelo, S., Murphy, T. B., Alfò, M. (2018).<br>
**Latent space modeling of multidimensional networks with application to the exchange of votes in Eurovision Song Contest**.<br>
*arXiv:1803.07166*.
