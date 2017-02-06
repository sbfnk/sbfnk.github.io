# Collection of LibBi models

## Deterministic SIR model, observations of prevalence

```c
model SIR_deterministic {
  const N = 1000; // population size
  const d_infection = 14; // duration of infection: 2 weeks
  
  state S, I, R;  // susceptible, infectious, recovered

  obs Prevalence;  // observations

  param R0; // basic reproduction number

  sub parameter {
    R0 ~ uniform(1, 3)
  }

  sub initial {
    S <- N - 1
    I <- 1
    R <- 0
  }

  sub transition { // daily time step

    inline beta = R0 / d_infection
    inline gamma = 1 / d_infection

    ode {
      dS/dt = - beta * S * I / N
      dI/dt = beta * S * I / N - gamma * I
      dR/dt = gamma * I
    }
  }

  sub observation {
    Prevalence ~ gaussian(mean = I, std = sqrt(I))
  }
}
```


## Deterministic SIR model, observations of incidence

```c
model SIR_deterministic {
  const h = 7; // incidence time step: 1 week
  const N = 1000; // population size
  const d_infection = 14; // duration of infection: 2 weeks
  
  state S, I, R, Z;  // susceptible, infectious, recovered, incident

  obs Incidence;  // observations

  param rep; //reporting rate
  param R0; // basic reproduction number

  sub parameter {
    rep ~ uniform(0,1)
    R0 ~ uniform(1,3)
  }

  sub initial {
    S <- N - 1
    I <- 1
    R <- 0
    Z <- 1
  }

  sub transition {

    n_transmission ~ wiener() // noise terms
    n_recovery ~ wiener() // noise terms

    Z <- (t_now % 7 == 0 ? 0 : Z) // reset incidence

    inline beta = R0 / d_infection * exp(n_transmission)
    inline gamma = 1 / d_infection * exp(n_recovery)

    ode {
      dS/dt = - beta * S * I / N
      dI/dt = beta * S * I / N - gamma * I
      dR/dt = gamma * I
      dZ/dt = beta * S * I / N
    }
  }

  sub observation {
    Incidence ~ truncated_gaussian(mean = rep * Z, std = sqrt(rep * Z + 1), lower = 0)
  }

}
```

## Stochastic SIR model, observations of incidence

```c
model SIR {
  const h = 7; // incidence time step: 1 week
  const N = 1000; // population size
  const d_infection = 14; // duration of infection: 2 weeks
  
  noise n_transmission;  // noise term
  noise n_recovery;  // noise term

  state S, I, R, Z;  // susceptible, infectious, recovered

  obs Incidence;  // observations

  param rep; //reporting rate
  param R0; // basic reproduction number

  sub parameter {
    rep ~ uniform(0,1)
    R0 ~ uniform(1,3)
  }

  sub initial {
    S <- N - 1
    I <- 1
    R <- 0
    Z <- 1
  }

  sub transition {

    n_transmission ~ wiener() // noise terms
    n_recovery ~ wiener() // noise terms

    Z <- (t_now % 7 == 0 ? 0 : Z) // reset incidence

    inline i_beta = R0 / d_infection * exp(n_transmission)
    inline i_gamma = 1 / d_infection * exp(n_recovery)

    ode (alg='RK4(3)', h=1e-1, atoler=1e-2, rtoler=1e-5) {
      dS/dt = - i_beta * S * I / N
      dI/dt = i_beta * S * I / N - i_gamma * I
      dR/dt = i_gamma * I
      dZ/dt = i_beta * S * I / N
    }
  }

  sub observation {
    Incidence ~ truncated_gaussian(mean = rep * Z, std = sqrt(rep * (1 - rep) * Z + 1), lower = 0)
  }

```

# Example observation data frame

```r
obs <- data.frame(value = c(1,6,2,26,99,57,78,57,15,9,4,1,1,1,0,2,0)
                  time = c(0,7,14,21,28,35,42,49,56,63,70,77,84,91,98,105,112))
```
