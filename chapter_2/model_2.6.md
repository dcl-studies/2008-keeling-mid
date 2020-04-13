SEIR with births and deaths
================
Bill Behrman
2020-04-12

  - [Model](#model)
  - [Plots](#plots)

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
  du[1] = p[3] - p[1] * u[1] * u[3] - p[3] * u[1]
  du[2] = p[1] * u[1] * u[3] - p[4] * u[2] - p[3] * u[2]
  du[3] = p[4] * u[2] - p[2] * u[3] - p[3] * u[3]
  du[4] = p[2] * u[3] - p[3] * u[4]
  return nothing
end
")

# Model parameters
p <- 
  c(
    beta = 520 / 365,
    gamma = 1 / 7,
    mu = 1 / (70 * 365),
    sigma = 1 / 14
  )

# Initial conditions
u0 <- 
  c(
    s0 = 0.1,
    e0 = 1e-4,
    i0 = 1e-4,
    r0 = 1 - 0.1 - 1e-4 - 1e-4
  )

# Time range for solution and time increment
t0 <- 0
t_end <- 60 * 365
t_inc <- t_end / 200

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
  rename(s = V1, e = V2, i = V3, r = V4, t = value)
```

## Plots

``` r
plot <- 
  sol %>% 
  pivot_longer(cols = -t, names_to = "class", values_to = "prop") %>% 
  ggplot(aes(t, prop, color = class)) +
  geom_line() +
  scale_x_continuous(
    breaks = scales::breaks_width(10 * 365),
    labels = scales::label_number(accuracy = 1, scale = 1 / 365)
  ) +
  scale_color_discrete(
    name = NULL,
    breaks = c("s", "e", "i", "r"),
    labels = c("Susceptible", "Exposed", "Infectious", "Recovered"),
  ) +
  theme(legend.position = "bottom") +
  labs(
    title = "SEIR with births and deaths",
    x = "Time (years)",
    y = "Percentage of the population in each class"
  )

plot +
  scale_y_continuous(
    breaks = scales::breaks_width(0.1),
    minor_breaks = NULL,
    labels = scales::label_percent(accuracy = 1)
  )
```

![](model_2.6_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
plot +
  scale_y_continuous(
    breaks = scales::breaks_width(2e-4),
    labels = scales::label_percent(accuracy = 0.01)
  ) +
  coord_cartesian(ylim = c(0, 12e-4))
```

![](model_2.6_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->
