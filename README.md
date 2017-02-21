# Statistical catch-at-age model

This is a stock assessment model that includes two methods to model the catches:
  - uncensored (catch is modelled with an observation error)
  - censored (catch is estimated to be between an upper and lower limit)
 
Other characteristics are:
- random effects: logFy and logN (annual fishing mortality and abundance)
- separable fishing mortality (and logFy as random walk)
- several recruitment options (Ricker, Beverton-Holt, random walk, etc.)
- process error with AR1 structure over years and ages
- catch-at-age modelled using the continuation-ratio logit transform
- annual SSB index
- estimation of several common reference points (Fmsy, SSBmsy, F20-30-40%, F01, Fmed, Fmax, Fcol)
 
Details can be found in Van Beveren et al. (submitted to ICES Journal of Marine Science)
 
# Usage

Requires:
- R
- the TMB package (Template Model Builder): https://github.com/kaskr/adcomp

# References

This code is part of a submitted article.

The censored catch, structured process error and crl transformed catch-at-age chunks are from:






