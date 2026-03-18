# Residualize trial-wise responses against nuisance terms

Fits linear models within optional groups and returns response residuals
in a new column. This provides a general trial-level nuisance-regression
step for behavioral analyses without hard-coding any specific
covariates.

## Usage

``` r
residualize_trials(
  data,
  response,
  terms = NULL,
  by = NULL,
  family = "gaussian",
  type = c("response", "working", "deviance", "pearson"),
  out = NULL
)
```

## Arguments

- data:

  data.frame containing response and nuisance columns.

- response:

  character scalar naming the response column.

- terms:

  optional character vector of nuisance column names. If empty or
  `NULL`, the response is copied unchanged.

- by:

  optional character vector of grouping columns. Separate models are fit
  within each group.

- family:

  model family. May be either a character scalar such as `"gaussian"`,
  `"binomial"`, or `"inverse.gaussian"`, or a family object returned by
  [`family`](https://rdrr.io/r/stats/family.html).

- type:

  residual type passed to
  [`residuals.glm`](https://rdrr.io/r/stats/glm.summaries.html).

- out:

  optional output column name. Defaults to `paste0(response, "_resid")`.

## Value

data.frame with an added residual column.
