# Numerical Integration of the Barotropic Vorticity Equation

We integrated the system which describes the Barotropic Vorticity Z over the rectangular domain of coordinates, corresponding roughly to North America, every half an hour for 24 hours.
The provided initial conditions consisted of the  field of Z at 500 hPa derived from observations for the 5th of January 1949, representing t0.
The Leapfrog method and a linear extrapolation of the boundaries for the Laplacian were used to achieve the numerical integration, resulting in a time-dependent forecast.
The following plot shows contour lines for the forecast at t= t0 + 24h, in red, compared to t0, in black. 
![](forecast.png)

![](analysis.png)

![](tendency.png)
#Prova
