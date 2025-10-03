
<table style="width:100%;">

<tr>

<td style="text-align:left; vertical-align:middle;">

<h1 style="margin:0;">

MUSC
</h1>

<!-- <p style="margin:0.25rem 0 0;"> -->

<!--   <strong>Authors:</strong> Chloe Mattila, Elizabeth Hill, Brian Neelon, and Souvik Seal -->

<!-- </p> -->

</td>

<!-- <td style="text-align:right; width:1%; white-space:nowrap;"> -->

<!--   <!-- invisible spacer: keeps layout, not visible -->

–\> <!--   <td style="text-align:right; padding-left:400px;"> -->
<!--   <img src="MUSC_logo_simple.png" width="140" height="150"/> -->
<!-- </td> -->
</tr>

</table>

<!-- README.md is generated from README.Rmd. Please edit that file -->

The *R* package implements the models proposed in the manuscript “MUSC:
MUlti-level variable selection for Spatially indexed Count data.” It
enables Bayesian variable selection for regression of overdispersed
negative binomial spatially indexed count data, with or without feature
grouping. MUSC supports both standard priors, horseshoe and
spike-and-slab, and their multi-level group extensions when feature
grouping is applied. The package is broadly applicable to spatial
datasets where the outcome is spatially indexed count data (e.g.,
state-level cancer incidence counts) and the predictors may be grouped
into meaningful categories (e.g., population/SES, behavior, smoking).

## Install and load MUSC

We install and load the developmental version of MUSC from GitHub.

``` r

# MUSC installation
devtools::install_github('Cmattila/MUSC', quiet = TRUE, force = TRUE)
# package 'magrittr' successfully unpacked and MD5 sums checked
# package 'Rcpp' successfully unpacked and MD5 sums checked
# package 'RcppArmadillo' successfully unpacked and MD5 sums checked
# package 'spdep' successfully unpacked and MD5 sums checked
library(MUSC)


# # List of packages for Import
# library(svMisc); library(MCMCpack);library(Matrix);
#  library(spdep); library(BayesLogit);library(spam); 
# library(coda); library(bayestestR); library(ggplot2);
# library(dplyr); library(grid); library(stats);
# library(graphics); library(utils); library(rlang); library(mvnfast)
```

## Simulate mildly correlated feature data, $p=50$

Then we simulate $(p-1)=199$ features with a within group correlation of
$\rho=0.25$. For each non-null coefficient $\beta_j$, the signal
magnitude was independently drawn from Unif(0.5, 2); all remaining
coefficients, intercept included, were set to zero.We set the group size
to $M_g = 5$, yielding $G = ⌈(p − 1)/M_g$⌉ groups, with all groups of
size five except one, which had four, with the thress selected groups
under Case 2 where within a group some features are active and some are
inactive, have three, three and four active features.

``` r
set.seed(2025*2025)
N <- 100 # number of observations/ spatial units
p_dim <- 50 # number of features (including features)
signal_vec <- rep(2, 10)

Mg = 5 # group size
G = p_dim/Mg # number of groups
group_ind <- generate_group_ind(p_dim, G) # group indicator
# Warning in generate_group_ind(p_dim, G): p - 1 is not divisible by G. Groups
# will be uneven.

#generating the features where nonnulls are under Case 1 where 2 fulls groups are active
nonNull_indices_p50_Case1 <- generate_nonnull_indices(p_dim, group_ind, num_nonnull = 10,
                                                        nonnull_group_structure = "Case1")

Null_indices_p50_Case1 <- setdiff(1:50, nonNull_indices_p50_Case1)  # Complement of nonNull indices

Features <- feature_sim_grouped(N = N, p = p_dim, group_ind = group_ind, 
                            nonNull_indices = nonNull_indices_p50_Case1, signal_vec = signal_vec)
```

## Construct the random effect $\phi$ and neighborhood list

We construct a $10 \times 10$ grid under a queen adjacency structure to
simulate the spatial random effect, $\mathbf{\phi}$ with $\nu = 0.1$
with which $\mathbf{y}$ will be generated and the associated
neighborhood list.

``` r
set.seed(2025*2025)

# Create a 10x10 grid of coordinates since we have an N=100 and assuming no repeated measures
coords100 <- expand.grid(x = 1:10, y = 1:10)

# Create a neighborhood structure (using "queen" to classify neighbors as sharing edges and corners) 
nb100 <- spdep::cell2nb(10, 10, type = "queen") 

plot(nb100,coords100)   # Visualization of a 10x10 lattice graph with queen-style adjacency (neighbors share edges or corners)
```

<img src="README_files/figure-gfm/constructing-phi-1.png" width="50%" />

``` r

n_mat100 <- spdep::nb2mat(nb100, style = "B") #b indicates binary 0/1
Q100 <- diag(rowSums(n_mat100)) - n_mat100    #iCAR precision 

## Simulating phi, iCAR spatial random effects
phi <- MUSC::phi_True_func(nu_vals=0.1, Q=Q100, n_space=100)
# nu = 0.1 phi range = -3.727779 4.306155 actual SD = 1.579542

# plotting phi
df <- data.frame(x = coords100$x,
                       y = coords100$y,
                       phi = phi$phi_nu_0.1)

phi_plot <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_tile(ggplot2::aes(fill = phi)) +
        ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
        ggplot2::labs(title = paste("Spatial Effect: ", expression(phi), " with ", expression(nu), "=0.1"),
             fill = expression(phi)) +
        ggplot2::theme_minimal() +
        ggplot2::coord_fixed()
print(phi_plot)
```

<img src="README_files/figure-gfm/constructing-phi-2.png" width="50%" />

## Simulate the spatially indexed count data

Then we simulate $N=100$ outcomes with spatial index phi, simulated
populations for each spatial unit, and a scale of per 100,000.

``` r
set.seed(2025*2025)

scale = 100000
pop <- sample(10000:50000, N, TRUE)
 
yy <- MUSC::y_gen_fun(N = N, r = 1, K = Features$K, beta_true = Features$beta_true, phi_true = phi$phi_nu_0.1, Scale = scale, pop_col = pop, offset = TRUE)
```

## Fit the SS and HS grouped varible selection models

The main function of “MUSC” fits negative binomial models using either
the standard horseshoe prior (which.prior = “HS”) or spike and slab
prior (which.prior = “SS”), or their grouped extensions (specifying the
vector of group indicators; *group_ind*). Additionally, these models can
be specified to include a spatial random effect (specifying the
neighborhood list, an object of class *nb*; *NeighborhoodList*). $X$
denotes the matrix of spatial-unit-level features (dimension
$N\times p-1$). *niter* denotes the number of MCMC iterations.*pop_col*
is the vector of the population for each spatial unit and *scale* is the
desired scale for the rates (i.e. per 100,000).

``` r
set.seed(2025*2025)

SS_group_offset_space <- MUSC::MUSC(X = scale(Features$X),                 # feature matrix
                              y = yy$y,                        # outcome vector
                              group_ind =  group_ind,         # group indicator vector
                              NeighborhoodList = nb100,       # Neighborhood list
                              pop_col = pop,                  # vector of spatial unit populations
                              scale = scale,                  # rate scale
                              which.prior = "SS",             # method
                              niter = 10000,                  # MCMC iterations; use more in practice, brief here for speed
                              verbose = FALSE)                # not printing MCMC iteration progress
# [1] "500 / 10000"
# [1] "1000 / 10000"
# [1] "1500 / 10000"
# [1] "2000 / 10000"
# [1] "2500 / 10000"
# [1] "3000 / 10000"
# [1] "3500 / 10000"
# [1] "4000 / 10000"
# [1] "4500 / 10000"
# [1] "5000 / 10000"
# [1] "5500 / 10000"
# [1] "6000 / 10000"
# [1] "6500 / 10000"
# [1] "7000 / 10000"
# [1] "7500 / 10000"
# [1] "8000 / 10000"
# [1] "8500 / 10000"
# [1] "9000 / 10000"
# [1] "9500 / 10000"
# [1] "10000 / 10000"

print(str(SS_group_offset_space)) # check the returned list object
# List of 9
#  $ Beta  : num [1:5000, 1:50] 1.63 1.74 1.69 1.55 1.64 ...
#  $ r     : num [1:5000, 1] 0.301 0.259 0.286 0.36 0.284 ...
#  $ Delta : num [1:5000, 1:50] 1 1 1 1 1 1 1 1 1 1 ...
#  $ pDelta: num [1:5000, 1:50] 1 1 1 1 1 1 1 1 1 1 ...
#  $ phi   : num [1:5000, 1:100] 0.01844 0.19825 1.04582 0.07089 -0.00758 ...
#  $ sphiSq: num [1:5000, 1] 0.3 0.322 0.309 0.269 0.278 ...
#  $ Zeta2 : num [1:5000, 1] 0.349 0.403 0.421 0.298 0.252 ...
#  $ Tau2  : num [1:5000, 1:10] 0.1008 0.4203 0.5589 0.6684 0.0776 ...
#  $ Lambda: num [1:5000, 1:49] 0 0 0 0 0 0 0 0 0 0 ...
# NULL


HS_group_offset_space <- MUSC::MUSC(X = scale(Features$X),                 # feature matrix
                              y = yy$y,                       # outcome vector
                              group_ind =  group_ind,         # group indicator vector
                              NeighborhoodList = nb100,       # Neighborhood list
                              pop_col = pop,                  # vector of spatial unit populations
                              scale = scale,                  # rate scale
                              which.prior = "HS",             # method
                              niter = 10000,                  # MCMC iterations; use more in practice, brief here for speed
                              verbose = FALSE)                # not printing MCMC iteration progress
# [1] "500 / 10000"
# [1] "1000 / 10000"
# [1] "1500 / 10000"
# [1] "2000 / 10000"
# [1] "2500 / 10000"
# [1] "3000 / 10000"
# [1] "3500 / 10000"
# [1] "4000 / 10000"
# [1] "4500 / 10000"
# [1] "5000 / 10000"
# [1] "5500 / 10000"
# [1] "6000 / 10000"
# [1] "6500 / 10000"
# [1] "7000 / 10000"
# [1] "7500 / 10000"
# [1] "8000 / 10000"
# [1] "8500 / 10000"
# [1] "9000 / 10000"
# [1] "9500 / 10000"
# [1] "10000 / 10000"

print(str(HS_group_offset_space)) # check the returned list object
# List of 9
#  $ Beta  : num [1:5000, 1:50] 2.39 2.24 1.89 2.21 1.86 ...
#  $ r     : num [1:5000, 1] 0.166 0.251 0.237 0.209 0.194 ...
#  $ Delta : num [1:5000, 1:50] 1 1 1 1 1 1 1 1 1 1 ...
#  $ pDelta: num [1:5000, 1:50] 1 1 1 1 1 1 1 1 1 1 ...
#  $ phi   : num [1:5000, 1:100] -0.1 -0.274 0.24 -0.177 -0.421 ...
#  $ sphiSq: num [1:5000, 1] 0.106 0.118 0.133 0.198 0.19 ...
#  $ Zeta2 : num [1:5000, 1] 0.00966 0.00855 0.00928 0.00954 0.00915 ...
#  $ Tau2  : num [1:5000, 1:10] 0.828 0.3 0.287 0.22 0.211 ...
#  $ Lambda: num [1:5000, 1:49] 8.37 26.77 4.1 3.85 1.31 ...
# NULL
```

## MCMC convergence plots

The above function returns a list of several objects, including a matrix
of the post-burn-in posterior beta estimates, with dimension *niter*
$\times p$. We randomly select four betas and plot their convergence.

``` r
set.seed(2025)

ran_sam <- sample(2:p_dim, 4)

titles <- c(
  paste0("HS - Beta ", ran_sam)
)

# Plot
par(mfrow = c(2, 2), mar = c(3, 4, 2, 1), oma = c(0, 0, 2, 0))
for (i in seq_along(ran_sam)) {
  MUSC::MCMC_plot(
    HS_group_offset_space$Beta[, ran_sam[i]],
    main = titles[i]
  )
}
```

<img src="README_files/figure-gfm/MCMC-diagnostic-plots-1.png" width="100%" />

## Visualize the posterior inclusion probablities (PIP) for SS

Display the posterior inclusion probability (PIP) of each beta with a
selection threshold of $0.75$ for the SS based model. Within the plot,
along the x-axis are the names of the features and above are their
respective groups. Truly non-null features are shown in a deeper color.

``` r
post_means <- colMeans(SS_group_offset_space$pDelta)
selected <- as.numeric(post_means > 0.75)
tp <- sum(selected[nonNull_indices_p50_Case1]) # truly nonnulls detected 
fn <- length(nonNull_indices_p50_Case1) - tp   # truly nonnulls not detected
fp <- sum(selected[Null_indices_p50_Case1])    # truly null detected; nulls include intercept
tn <- length(Null_indices_p50_Case1) - fp      # truly null not detected; nulls include intercept



# predictors correspond to columns 2:p (since col 1 is the intercept)
predictors <- paste0( 1:(p_dim - 1))              # keys for plotting
predictor_labels <- predictors                     # or swap in pretty labels if you have them

# Group labels: map each predictor to its group id (1..G)
# group_ind is length p-1, aligned to predictors
group_map <- setNames(paste0(group_ind), predictors)
group_levels <- paste0( sort(unique(group_ind)))  # order under x-axis

# plotting in uniform color
# res <- plot_pdelta_bars_mcmc( 
#   predictors        = predictors,
#   pdelta_mat        = SS_group_offset_space$pDelta[,-1],      # iterations × (p-1) matrix, removing the intercept column 
#   outcome_name      = "SS (group, spatial, offset)",    # whatever label you want in the legend
#   predictor_labels  = predictor_labels,          # pretty x-axis labels (optional)
#   group_map         = group_map,                 # bands under x-axis (optional)
#   group_levels      = group_levels,              # order of bands (required if using groups)
#   use_groups        = TRUE,
#   x_order           = "group",                   # "group" to sort by bands, else "given"
#   thresholds        = c(0.75),               # reference lines
#   legend_title      = "Model/Outcome",
#   colors            = "steelblue",               # single outcome → one color
#   base_size         = 12,
#   bracket_text_size = 4
# )
# print(res$plot)


#plotting with truly nonnulls in a darker shade
res_highlight<- MUSC::plot_pdelta_bars_mcmc_highlights(
  predictors        = predictors,
  pdelta_mat        = SS_group_offset_space$pDelta[,-1],      # iterations × (p-1) matrix, removing the intercept column 
  outcome_name      = "SS (group, spatial, offset)",    # whatever label you want in the legend
  predictor_labels  = predictor_labels,          # pretty x-axis labels (optional)
  group_map         = group_map,                 # bands under x-axis (optional)
  group_levels      = group_levels,              # order of bands (required if using groups)
  use_groups        = TRUE,
  x_order           = "group",                   # "group" to sort by bands, else "given"
  thresholds        = c(0.75),               # reference lines
  legend_title      = "Model/Outcome",
  base_size         = 12,
  bracket_text_size = 4,
  highlight_idx  = nonNull_indices_p50_Case1-1, #shifting indexes by 1 to accomodate predictors 1:p-1
  highlight_fill = "#1F77B4"  # highlight color for truly nonnulls
)
# Coordinate system already present.
# ℹ Adding new coordinate system, which will replace the existing one.
# Coordinate system already present.
# ℹ Adding new coordinate system, which will replace the existing one.

print(res_highlight$plot)
```

<img src="README_files/figure-gfm/PIP-visualization-1.png" width="100%" />

## Visualize the 95% equal-tail interval based credible intervals for HS

Display the 95% ETI CrI of each beta for the HS based model. Within the
plot, along the x-axis are the names of the features and above are their
respective groups. Truly non-null features are shown in a deeper color.

``` r
inclusion <- MUSC::local_credible(HS_group_offset_space$Beta, local.p.ths = 0.95)  # 0/1 vector
tp <- sum(inclusion[nonNull_indices_p50_Case1]) # truly nonnulls detected 
fn <- length(nonNull_indices_p50_Case1) - tp   # truly nonnulls not detected
fp <- sum(inclusion[Null_indices_p50_Case1])    # truly null detected; nulls include intercept
tn <- length(Null_indices_p50_Case1) - fp      # truly null not detected; nulls include intercept


# predictors correspond to columns 2:p (since col 1 is the intercept)
predictors <- paste0( 1:(p_dim - 1))              # keys for plotting
predictor_labels <- predictors                     # or swap in pretty labels if you have them

# Group labels: map each predictor to its group id (1..G)
# group_ind is length p-1, aligned to predictors
group_map <- setNames(paste0(group_ind), predictors)
group_levels <- paste0( sort(unique(group_ind)))  # order under x-axis

#plotting in uniform color
# res1 <- plot_cri(
#   outcomes        = "beta",
#   betas           = list(beta=HS_group_offset_space$Beta), 
#   method          = "HS",
#   predictors      = predictors,
#   predictor_labels= predictor_labels,
#   group_map       = group_map,
#   group_levels    = group_levels,
#   use_groups      = TRUE,
#   cred_level      = 0.95,
#   star            = "median",          # posterior median (or mean) shown as a star
#   x_order         = "group",
#   legend_title    = "Outcome",
#   show_legend     = FALSE,             # auto-hides with one outcome if desired
#   verbose         = TRUE,
#   base_size       = 8, 
#   colors = "darkgreen",
#   bracket_text_size = 2
# )
# print(res1$plot)


#plotting with truly nonnulls in a darker shade
res_HS_highlight <- MUSC::plot_cri_highlights(
  outcomes = c("HS"),
  predictors = predictors,
  betas = list(
    HS = HS_group_offset_space$Beta[,-1]   # iterations x p-1 (excluding intercept)
  ),
  star = "median",
  group_map       = group_map,
  group_levels    = group_levels,
  x_order = "group",
  highlight_idx = (nonNull_indices_p50_Case1-1),
  verbose         = TRUE,
  base_size       = 12, 
  colors = "darkgreen",
  bracket_text_size = 4
)
#   outcome_key pretty n_iter n_pred
# 1          HS     HS   5000     49
print(res_HS_highlight$plot)
```

<img src="README_files/figure-gfm/CrI-visualization-1.png" width="100%" />
