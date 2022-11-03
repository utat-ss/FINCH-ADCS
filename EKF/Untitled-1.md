# Extended Kalman Filter

We will use the EKF it to compute the attitude as a quaternion with the observations of tri-axial sensors.

The **state** is the physical state, which can be desribed by dynamic variables. The noise in the measurements means that there is a certain degree of uncertainty in them. 

A dynamical system is the physical statw whose state evolve over time, so diffrential equations are normally used to model them. There is also noise in the dynamics of the system, **process noise**, which means we cannot be entirely deterministic, but we can get indirect noisy measurements. 

**Filtering** refers to the process of filtering out the noise in the measurements to provide an optimal estimate for a state, given the observed measurement. 

The instantaneous state of the system is represented with a vector update through discrete time increments to generate the next state. The simplest of the state space models are linear models which can be expressed with equations of the following form

$$X_{t} = F_{x_{t-1}}+W_{x}$$
$$Z_{t} = H_{x_{t}}+ W_{z}$$

$X_{t}$: is the state of the system describing the condition of the n elements of time t.

$Z_{t}$: are the measurements at time t. 

$W_{x}$:is the process noise at time t.

$F$ is called either the State Transition Matrix or the Fundamental Matrix ot the Fundamental Matrix, and sometimes is represented with $\phi$

$H$: is the measurement model matrix.

Many linear models are also described with continuous time state equations of the form:

$$\dot{X_{t}} = \frac{dx_t}{dt} = Ax_{t} + Lw_{t}$$

Where A and L are constant matrices characterizing the bahaviour of the model, and $w_{t}$ is a white noise with a power spectal density $\sigma_{w}^2$

THe main diffrence is that A models as a set of linear differential equation, and is continuous. F is discrete, and represents a set of linear equations (not differential equations) which transitions $X_{t-1}$ to $X_{t}$ over a discrete time step $\Delta t$

A commmon way to optain F uses the matrix exponential, which can be expanded with a Taylor series. 

The main goal is to find equation that recursively finds the value of $X_t$ in term of $X_{t-1}$

## Extended Kalman Filter

The fuctions have been Gaussian and linear so far and, thus, the output was always another Gaussian, but Gaussians are not closed under nonlinear fuction.

The EKF handles nonlinearity by forming a Gaussian approximation to the joint distribution of state x and measurements z using Taylor Series based transformation. 