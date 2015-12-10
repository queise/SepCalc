# SepCalc
Calculation of asteroseismic parameters (separations, ratios and slopes), and chi-squared comparision with observed data.

Needs some inputs, all present as example files in the [example folder](example):
- the observed acoustic frequencies of a star (star_observed.freqs)
- other observed stellar properties (input.dat and star_observed.*)
- the oscillation frequencies of a stellar model, calulated with [ADIPLS](http://astro.phys.au.dk/~jcd/adipack.n/) (star_model.freqs)
- other properties of the stellar model (star_model.fin)

Creates several output files (chose which ones in input.dat):
- star_model.Ls012
- star_model.r010
- star_model.ech
- star_model.Chi2
- StelChars.dat
- Seismparams.dat
- Chi.dat

Please read the info in [SepCalc.f](src/SepCalc.f) for more information on how to use the code.
