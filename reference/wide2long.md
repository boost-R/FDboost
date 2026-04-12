# Transform id and time of wide format into long format

Transform id and time from wide format into long format, i.e., time and
id are repeated accordingly so that two vectors of the same length are
returned.

## Usage

``` r
wide2long(time, id)
```

## Arguments

- time:

  the observation points

- id:

  the id for the curve

## Value

a list with `time` and `id`
