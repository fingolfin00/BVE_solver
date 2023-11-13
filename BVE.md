# Numerical Integration of the Barotropic Vorticity Equation

We integrated the system which describes the Barotropic Vorticity Z over the rectangular domain of coordinates, corresponding roughly to North America, every half an hour, for 24 hours.
The provided initial conditions consisted of the field of Z at 500 hPa derived from observations for the 5th of January 1949, representing t0.
The Leapfrog method and a linear extrapolation of the boundaries for the Laplacian were used to achieve the numerical integration, resulting in a time-dependent forecast.

The numerical integration proceeds as follow. 
The initial-conditions file, consisting of data for geopotential eight in a 19x16 grid, is read as the starting point of the forecast. From this intial values of Z we compute the Laplacin with the "make_Laplacian" function that reads the Z field and return a 17x15, thus exclunding boundaries. Boundaries are then computed with the "extrapolate" function that extrapolate linearly the boundaries according to:
\[X_{i,0}=2*X_{i,1}-X_{i,2} , X_{i,M}=2*X_{i,M-1}-X_{i,M-2} ;\]
and analogous for intial and finale row. The four corners are then updated.


The following plot shows contour lines for the forecast at t= t0 + 24h, in red, compared to t0, in black. 

![](forecast.png)

![](analysis.png)

![](tendency.png)
# Prova
