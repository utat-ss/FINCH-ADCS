# Ephemerides

This README requires a [MathJax](https://www.mathjax.org/) browser plugin for the mathematicsto be rendered in between then `$...$`. Many are widely available for all the major browsers.



## Sun

Determining the position of the Sun is of particular use to satellites for power generation analysis (using solar cells) and for Attitude Determination and Control Subsystems which utilize sun sensors. The Astronomical Almanac lists an algorithm which evaluates the vector from the Earth to the Sun in a Mean-Equator of Date (MOD) frame with the following properties:

1. Accuracy: 0.01º
1. Validity: 1950 to 2050

First, the mean longitude of the sun $\lambda_{\odot}$ is calculated in the MOD frame using the Julian Centuries in UT1 ($T_{ut1}$):

$$\lambda_\odot = 280.460 + 36000.771 \cdot T_\textrm{ut1}$$

Next the mean anomaly for the Sun $\mathrm{M}\_\odot$ is found using $T_{ut1}$:

$$M_\odot =  357.5277233 + 35999.05034 \cdot T_ut1$$

The ecliptic latitude ($\phi\_\xi$) and longitude ($\lambda\_\xi$) of the Sun are determined with the aforementioned quantities. The ecliptic latitude is often approximated as $0º$, while the ecliptic longitude can be expressed as:

$$\lambda_\xi = \lambda_\odot + 1.914666471 \cdot \sin M_\odot + 0.019994643 \cdot \sin\left(2M_\odot \right)$$

Before expressing the Sun's coordinates with respect to the Earth, we also determine the obliquity of the ecliptic:

$$\epsilon = 23.439291 - 0.0130042\cdot T_\textrm{ut1}$$

By considering the Earth's semi-major axis $a\_\oplus$ and Earth's eccentricity $e\_\oplus$ we can evaluate the magnitude of the Sun's position:

$$\begin{aligned} r_\odot = a_\oplus &\left[ 1 + \frac{e_\oplus^2}{2}  + \left( -e_\oplus + \frac{3e_\oplus^2}{8} -\frac{5e_\oplus^5}{192} + \cdots \right)\right] \cos M_\odot   \\
&+ \left(-\frac{e_\odot^2}{2} + \frac{e_\odot^4}{3} = \frac{e_\odot^6}{16} + \cdots \right) \cos \left( 2M_\odot \right) + \mathcal{O}(e_\oplus^3)
\end{aligned}$$

Finally, we can combine all of the above to evaluate the Sun's position in the MOD frame:

$$\hat{\mathbf{r}}\_\odot = \begin{bmatrix} \cos\lambda_\xi & \cos\epsilon \cdot \sin\lambda_\xi & \sin\epsilon \cdot \sin \lambda_\xi \end{bmatrix}^T$$

<!-- TODO: Validated using Astronomical Almanac-->
<!-- TODO: Mathbf underscore must be \_-->
<!-- TODO: matrix newline must be \\\\, since \\ is rendered as \ in markdown-->


## Moon

<!--TODO-->

Knowledge of the Moon's location is of particular use for optical tracking satellites with software and hardware that presents direct observation of the Moon.

### Moon Position Vector

We often need to know the position vector from the Earth to the Moon. The ephemerides of the Jet Propulsion Laboratory (JPL) are the most accurate, but the Astronomical Almanac provided us with a less precise formula.

First, we need to determine the number of centuries elapsed from epoch J2000 using Julian Date in Barycentric dynamic time ($\mathcal{J}_{tdb}$)

$$\mathrm{T_{tdb}} = \frac{\mathcal{J}_{tdb} - 2\,451\,545.0} {36\,525}$$

Next, we need to find an expression for the mean anomalies ($M_{☾}$)

$$\mathrm{M_{☾}} = 134.9 + 477\,198.85\cdot T_{tdb}$$

The mean argument of latitude ($u_{M_{☾}}$) can be found using ($T_{tdb}$)

$$u_{M_{☾}} = 93.3 + 483\,202.03  \cdot T_{tdb}$$

The elongation ($D_\odot$) can also be found using $T_{tdb}$:

$$D_\odot = 117.85 + 445\,267.115\cdot T_{tdb}$$
  
As with the sun, we can find the ecliptic longitude ($\lambda_\xi$), ecliptic latitude ($\phi_\xi$) and the parallax ($\vartheta$) using the following equations.

$$\lambda_\xi = \lambda_{☾} + 6.29 \cdot sin(M_☾)- 1.27\cdot \sin\left(M_☾-2D_\odot\right) + 0.66\cdot sin(2D_{\odot})+ 0.21\cdot sin(2M_☾) - 0.19\cdot sin(M_\odot) - 0.11\cdot sin(u_{M_{☾}})$$

The ecliptic latitude ($\phi_\xi$) isn't zero, so we can find it as

$$\phi_\xi = 5.13\cdot sin(u_{M_{☾}}) + 0.28\cdot sin(M_☾ + u_{M_{☾}}) - 0.28\cdot sin(u_{M_{☾}} - M_☾) - 0.17\cdot sin(u_{M_{☾}} - 2{D}_{\odot})$$

The horizontal parallex ($\mathrm{\vartheta}$) is

$$\mathrm{\vartheta} = 0.9508 + 0.0518\cdot cos(M_{☾}) + 0.0095\cdot cos(M_{☾} - 2D_\odot) + 0.0078\cdot cos(2D_{\odot}) + 0.0028\cdot cos(2M_{☾})$$

The obliquity of the ecliptic can be calculated using $T_{tdb}$

$$\epsilon = 23.439291 - 0.0130042 \cdot T_{tdb}$$

The next step is to find the magnitude of the moon position vector.

$$\mathrm{r}_{☾} = \frac{1}{sin(\vartheta)}$$

Finally, we can combine all of the above to evaluate the Moon Position Vector expressed in the mean frame:

$$\mathbf{r}\_\phi = r_{☾} \cdot \begin{bmatrix} cos(\phi_\xi)cos(\lambda_\xi) \\\\ cos(\epsilon)cos(\phi_\xi)cos(\lambda_\xi) - sin(\epsilon)sin(\phi_\xi) \\\\ sin(\epsilon)cos(\phi_\xi)sin(\lambda_\xi) + cos(\epsilon)sin(\phi_\xi) \end{bmatrix}$$

## Phase of the Moon

The percentage of the Moon's surface that reflects light back to Earth is constantly change. There are four **phases of the Moon** (new, full, last, and first quarter phases). The Moon's phase can be calculated using the Sun's and Moon's ecliptic longitudes.

$$\mathrm{phase}\_{☾}=\lambda\_{ecliptic\_{\odot}} - \lambda\_{ecliptic\_{☾}}$$

The values for each phase are:

 - $0^{\circ}$: new moon
 -  $90^{\circ}$: first quarter
 -  $180^{\circ}$: full
 -  $270^{\circ}$: last quarter

The phase is close to the angle between the Sun and the Moon. The percentage of the Moon's surface that is illuminated can be approximated as

$$ \\%disk= \frac{100\\%}{2}\cdot(1 - cos(phase_☾))$$

The Earth blocks solar light from reaching the Moon (lunar eclipse) does occur on occasion. The lunar phases are due to the relative position of the Moon with respect to an Earth observer.

 The phases of the Moon also have other names. **Waxing** and **waning** Moons describe if the apparent Moon illumination is increasing or decreasing, respectively. **Gibbous** is when more than 50% of the Moon appears illuminated.
 
  Each season (winter, etc.) has 3 full Moons. The blue moon occurs whenever there are 4 full Moons in a season, the third full moon is the blue Moon.
