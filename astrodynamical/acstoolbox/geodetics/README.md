# Earth

This README requires a [MathJax](https://www.mathjax.org/) browser plugin for the mathematicsto be rendered in between then `$...$`. Many are widely available for all the major browsers.

## Site Components

Site locations, such as those found on standard maps, are usually given in geodetic coordinates. But we require geocentric values for calculations that involve gravitational potential operations. Therefore, we must be able to represent position vectors in geocentric values.

First, we need to determine two auxiliary quantities obtained from the geometrical properties of an ellipse.

$$\mathrm{C}\_{\bigoplus}= \frac{{R\_{\bigoplus}}}{\sqrt{1 - e\_{\bigoplus}^2sin^2(\phi\_{gd})}}$$

$$\mathrm{S}\_{\bigoplus}= \frac{{R\_{\bigoplus}(1 - e\_{\bigoplus})}}{\sqrt{1 - e\_{\bigoplus}^2sin^2(\phi\_{gd})}}$$

Where  the mean equatorial radius of the Earth $R_{\bigoplus}$ = 6378.1363 km and the  Earth's eccentricity $e_{\bigoplus}$ = 0.081 819 221 456. $C_{\bigoplus}$ is commonly known as the radius of curvature of the meridian.

Next, we can determine the site vector at a location on Earth. We can assume that the height above sea level ($h_{ellp}$) is equal to the height above the ellipsoid.

$$\vec{r_{IJK}} = \begin{bmatrix} (C_{\bigoplus} + h_{ellp})\cdot cos(\phi_{gd}) \cdot cos(\lambda) \\\\(C_{\bigoplus} + h_{ellp})\cdot cos(\phi_{gd}) \cdot sin(\lambda)\\\\(S_{\bigoplus}+h_{ellp})\cdot sin(\phi_{gd})\end{bmatrix} $$