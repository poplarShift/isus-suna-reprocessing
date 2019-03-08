# Spectrophotometric Nitrate Detection using the ISUS/SUNA detector family

[For](https://poplarshift.github.io/papers/randelhoff2017vertical.pdf) [my](https://poplarshift.github.io/papers/randelhoff2015seasonal.pdf) [own](https://poplarshift.github.io/papers/randelhoff2016regional.pdf) [research](https://poplarshift.github.io/papers/randelhoff2016vertical.pdf) I wrote a software package that handles reading, integration and reprocessing data from the [ISUS/SUNA family of nitrate sensors](https://www.seabird.com/nutrient-sensors/suna-v2-nitrate-sensor/family?productCategoryId=54627869922) with CTD temperature and salinity data. Much of it is really just data handling (such as e.g. aligning temperature/salinity profiles with nitrate profiles). This isn’t strictly relevant to the reprocessing itself - especially if you have exactly synchronous timestamps thanks due to e.g. telemetry, in which case you only need to consider accounting for whether sensors are pumped or not. In fact, synchronicity within a second or two is usually enough for vertical profiles; for horizontal transects (such as on ROVs), or on moored timeseries, one could probably allow even more. It all depends on how steep the gradients are.

I entirely switched to Python (and never looked back) soon after completing my PhD, and hence was never motivated to turn this piece of Matlab code into robust software beyond development status that could be used by others without knowing about the software details. But due to repeated requests I'm making the essential bits of the processing code available here. If you have a lot of vertical profiles of ISUS/SUNA data and don't want to start from scratch (and for some reason feel like using Matlab...), send me a mail and I may be able to help you.

The centerpiece of the reprocessing, i.e. Sakamoto et al.'s algorithm (2009), is implemented in the file [ISUS_REPROCESSOR_v2_2.m](ISUS_REPROCESSOR_v2_2.m). You can also find an example of a ISUS/SUNA calibration file, and how to translate it into a .mat file that the reprocessing algorithm understands. The input to the `ISUS_REPROCESSOR_v2_2.m` function are:
- `CHANNEL`: the absorption data, one value per channel per sample,
- `CHANNEL_DF`: the dark frames as measured in the same ISUS data file,
- in-situ temperature `T`,
- salinity `S`, and
- a string `calfile` that supplies the name of the calibration `.mat` file.
The last argument `progressbarflag` is not important when using this function as a standalone.

You can find more detail on the accuracy of the processing and algorithm in our 2016 paper:

Randelhoff, A., Fer, I., Sundfjord, A., Tremblay, J.-É., & Reigstad, M. (2016). Vertical Fluxes of Nitrate in the Seasonal Nitracline of the Atlantic Sector of the Arctic Ocean. Journal of Geophysical Research: Oceans, 121(7), 5282–5295. https://doi.org/10.1002/2016JC011779,

(also available [under this link](https://poplarshift.github.io/papers/randelhoff2016vertical.pdf)) and on [my webpage](https://poplarshift.github.io).
