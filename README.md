# Censored Catch Assessment Model

This is a Statistical catch-at-age stock assessment model that includes two methods to model catch:
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
 
Details can be found in:

Van Beveren, E., Castonguay, M., Doniol-Valcroze, T., Cadigan, N., Plourde, S., Duplisea, D. (2017). How catch underreporting can bias stock assessment and advice in northwest Atlantic mackerel and a possible resolution using censored catch. Fish. Res. http://dx.doi.org/10.1016/j.fishres.2017.05.015 
 
# Usage

Requires:
- R
- the TMB package (Template Model Builder): https://github.com/kaskr/adcomp

# References

This code is part of:
- Van Beveren, E., Duplisea,D., Castonguay,C., Doniol-Valcroze,T., Plourde, S., Cadigan, N., (2017), How catch underreporting can bias stock assessment of and advice for northwest Atlantic mackerel and a possible resolution using censored catch.  Fish. Res.

The censored catch, structured process error and crl transformed catch-at-age cpp code chunks are from:

- Cadigan, N. 2016a. A state-space stock assessment model for northern cod, including under-reported catches and variable natural mortality rates. Can. J. Fish. Aquat. Sci., 73: 296–308. http://www.nrcresearchpress.com/doi/pdfplus/10.1139/cjfas-2015-0047

- Cadigan, N. 2016b. Updates to a Northern Cod (Gadus morhua) State-Space Integrated Assessment Model. DFO Can. Sci. Advis. Sec. Res. Doc., 2016/022. Centre for Fisheries Ecosystem Research, St. John’s, NL. http://publications.gc.ca/collections/collection_2016/mpo-dfo/Fs70-5-2016-022-eng.pdf

For more information on censored catch:

- Hammond, T. R., and Trenkel, V. M. 2005. Censored catch data in fisheries stock assessment. ICES Journal of Marine Science, 62: 1118–1130. https://academic.oup.com/icesjms/article/62/6/1118/617407/Censored-catch-data-in-fisheries-stock-assessment

- Bousquet, N., Cadigan, N., Duchesne, T., and Rivest, L.-P. 2010. Detecting and correcting underreported catches in fish stock assessment: trial of a new method. Can. J. Fish. Aquat. Sci., 67: 1247–1261.http://www.nrcresearchpress.com/doi/abs/10.1139/F10-051#.WKyhadcrKM8




