SIR with disease-induced mortality and density-dependent transmission
================
Bill Behrman
2020-04-11

  - [Model](#model)
  - [Plot](#plot)

``` r
# Libraries
library(tidyverse)
library(diffeqr)

#===============================================================================

# Setup ODE solver
diffeq_setup() %>% 
  invisible()
```

## Model

``` r
# Derivative function
f <- JuliaCall::julia_eval("
function f(du, u, p, t)
  du[1] = p[4] - p[1] * u[1] * u[2] - p[3] * u[1]
  du[2] = p[1] * u[1] * u[2] - (p[2] + p[3]) / (1 - p[5]) * u[2]
  du[3] = p[2] * u[2] - p[3] * u[3]
  return nothing
end
")

# Model parameters
p <- 
  c(
    beta = 520 / 365,
    gamma = 1 / 7,
    mu = 1 / (70 * 365),
    nu = 1 / (70 * 365),
    rho = 0.5
  )

# Initial conditions
u0 <- 
  c(
    x0 = 0.2,
    y0 = 1e-6,
    z0 = 1 - 0.2 - 1e-6
  )

# Time range for solution and time increment
t0 <- 0
t_end <- 300 * 365
t_inc <- t_end / 1000

tspan <- c(t0, t_end)
saveat <- seq(t0, t_end, t_inc)

# Tolerances for solution
reltol <- 1e-8
abstol <- 1e-8

# Solve model ODE
sol <- 
  ode.solve(
    f = "f", 
    u0 = u0,
    tspan = tspan,
    p = p,
    reltol = reltol,
    abstol = abstol,
    saveat = saveat
  ) %>%
  map_dfc(as_tibble) %>% 
  rename(x = V1, y = V2, z = V3, t = value)
```

## Plot

``` r
sol %>% 
  pivot_longer(cols = -t, names_to = "class", values_to = "population") %>% 
  ggplot(aes(t, population, color = class)) +
  geom_line() +
  scale_x_continuous(
    breaks = scales::breaks_width(50 * 365),
    labels = scales::label_number(scale = 1 / 365)
  ) +
  scale_color_discrete(
    name = NULL,
    breaks = c("x", "y", "z"),
    labels = c("Susceptible", "Infectious", "Recovered"),
  ) +
  theme(legend.position = "bottom") +
  labs(
    title = 
      "SIR with disease-induced mortality and density-dependent transmission",
    x = "Time (years)",
    y = "Population in each class"
  )
```

![](model_2.3_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->
