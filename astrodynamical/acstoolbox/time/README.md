# Time

This README requires a [MathJax](https://www.mathjax.org/) browser plugin for the mathematicsto be rendered in between then `$...$`. Many are widely available for all the major browsers.

## Introduction

Time is an important fundamental dimension, especially in astrodynamics where high speed objects travel great distances in short amounts of time. ACS design and GNC algorithms rely on 3 main time standards: UT1, Terrestrial Time and Barycentric Dynamical Time. The most reliable and precise time interval accepted by the astronomical community is the SI second which is based on the quantum transition of elecrons in a Cs-133 atom.

### Universal Time (UT0, UT1, UTC)

Universal time (UT) is evaluated using the Earth's rotation, which varies over time. UT approximates mean solar time over Greenwich, so that the mean sun is directly on the meridian at 12h00. To enable this, each second in this standard shrinks and grows so that each mean solar day is equal to 86400 UT seconds. Thus each UT second is **not** necessarily equal to a SI second. 

There are a number of realizations to UT. UT1 accounts for variations in the Earth's rotation due to polar motion, while UT0 does not. UTC which is the basis of civil time and found on everyday clocks. It is a convenient approximation to UT1 and their offset $\Delta t = t_{ut1} - t_{utc}$ is at most offset by $\pm0.9$ s.

### Terrestrial Time (TT)

The Terrestrial Time (TT) standard is a continuosly running time scale, unaffected by the irregularities of the Earth's rotation and orbit around the Sun. 1 day in TT has precisely 86400 SI seconds.

###  Barycentric Dynamical Time (TDB)

Barycentric Dynamical Time (TDB) is a time standard whose reference system is at rest with respect to solar system's barycenter. It differs from TT by periodic terms related to the Earth's orbit.

For geocentric missions orbiting around Earth, Julian centuries evaluated from J2000 in TT are a valid approximation to the same quantity evaluated in TDB (i.e., $c_{tdb} \approx c_{tt}$).

<!-- TODO Additional Information-->
<!-- TODO API -->