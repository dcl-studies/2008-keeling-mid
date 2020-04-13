SIR with carrier state
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
  du[1] = p[5] - p[1] * u[1] * u[2] - p[4] * p[1] * u[1] * u[3] - p[5] * u[1]
  du[2] = p[1] * u[1] * u[2] + p[4] * p[1] * u[1] * u[3] - p[2] * u[2] - 
    p[5] * u[2]
  du[3] = p[2] * p[6] * u[2] - p[3] * u[3] - p[5] * u[3]
  du[4] = p[2] * (1 - p[6]) * u[2] + p[3] * u[3] - p[5] * u[4]
  return nothing
end
")

# Model parameters
p <- 
  c(
    beta = 0.2,
    gamma_i = 1 / 100,
    gamma_c = 1 / 1000,
    epsilon = 0.1,
    mu = 1 / (50 * 365),
    q = 0.4
  )

# Initial conditions
u0 <- 
  c(
    s0 = 0.1,
    i0 = 1e-4,
    c0 = 1e-3,
    r0 = 1 - 0.1 - 1e-4 - 1e-3
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
  rename(s = V1, i = V2, c = V3, r = V4, t = value)
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
    breaks = c("s", "i", "c", "r"),
    labels = c("Susceptible", "Infectious", "Carrier", "Recovered"),
  ) +
  theme(legend.position = "bottom") +
  labs(
    title = "SIR with carrier state",
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

![](model_2.7_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
plot +
  scale_y_continuous(
    breaks = scales::breaks_width(0.02),
    labels = scales::label_percent(accuracy = 1)
  ) +
  coord_cartesian(ylim = c(0, 0.11))
```

![](model_2.7_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->
