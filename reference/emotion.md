# EEG and EMG recordings in a computerised gambling study

To analyse the functional relationship between electroencephalography
(EEG) and facial electromyography (EMG), Gentsch et al. (2014)
simultaneously recorded EEG and EMG signals from 24 participants while
they were playing a computerised gambling task. The given subset
contains aggregated observations of 23 participants. Curves were
averaged over each subject and each of the 8 study settings, resulting
in 23 times 8 curves.

## Usage

``` r
data("emotion")
```

## Format

A list with the following 10 variables.

- `power`:

  factor variable with levels *high* and *low*

- `game_outcome`:

  factor variable with levels *gain* and *loss*

- `control`:

  factor variable with levels *high* and *low*

- `subject`:

  factor variable with 23 levels

- `EEG`:

  matrix; EEG signal in wide format

- `EMG`:

  matrix; EMG signal in wide format

- `s`:

  time points for the functional covariate

- `t`:

  time points for the functional response

## Details

The aim is to explain potentials in the EMG signal by study settings as
well as the EEG signal (see Ruegamer et al., 2018).

## Source

Gentsch, K., Grandjean, D. and Scherer, K. R. (2014) Coherence explored
between emotion components: Evidence from event-related potentials and
facial electromyography. Biological Psychology, 98, 70-81.

Ruegamer D., Brockhaus, S., Gentsch K., Scherer, K., Greven, S. (2018).
Boosting factor-specific functional historical models for the detection
of synchronization in bioelectrical signals. Journal of the Royal
Statistical Society: Series C (Applied Statistics), 67, 621-642.

## Examples

``` r
data("emotion", package = "FDboost")

# fit function-on-scalar model with random effect and power effect
fos_random_power <- FDboost(EMG ~ 1 + brandomc(subject, df = 2)
                            + bolsc(power, df = 2),
                            timeformula = ~ bbs(t, df = 3),
                            data = emotion)
#> Use a smooth offset.
if (FALSE) { # \dontrun{                            
# fit function-on-function model with intercept and historical EEG effect
# where limits specifies the used lag between EMG and EEG signal
fof_historical <- FDboost(EMG ~ 1 + bhist(EEG, s = s, time = t,
                          limits = function(s,t) s < t - 3),
                          timeformula = ~ bbs(t, df = 3), data = emotion,
                          control = boost_control(mstop = 200))                            
} # }
```
