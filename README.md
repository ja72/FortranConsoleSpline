# FortranConsoleSpline

Utility module for cubic spline interpolation, with example program.

## Usage

Enable the module in your code with

```fortran
use mod_splines
```

and then define a `type(spline)` object in your variable 
declaration section. For example:

```fortran
type(spline) :: sp
```

Initialize the spline with two arrays, one for `x` values and one
for `y` values. Both arrays must have the same number of elements

```fortran
! Create a natural spline
sp = spline(xi, yi)
```

If you know the start and end slope of the spline, you can add them 
as optional parameters

```fortran
! Create a spline with flat ends
sp = spline(xi, yi, 0.0d0, 0.0d0)
```

To extract a value for the spline use the elemental type bound function `value(x)`

```fortran
y = sp%value(x)
```

Here `x` and `y` are scalar values. But since it is an elemental function, it can work with arrays directly. 
For example

```fortran
xi= [(x0+i*h, i=0, n-1)]
yi = sp%value(xi)
```

where `xi` and `yi` are arrays

Additionaly, you can interpolate the slope and the second derivative using the
`slope(x)` and `slope2(x)` functions. For example

```fortran
y = sp%value(x)
yp = sp%slope(x)
ypp = sp%slope2(x)
```

As a note, if the evaluation point `x` is outside the defined range for the
spline, then a linear iterpolation is done from the ends of the spline.

Finally, if a local extrema point is needed (local minimum or maximum) a bisection
method is used to find the `x` value where the slope is zero using
the `extrema(x_low,x_high)` function

```fortran
x_max = sp%extrema()
y_max = sp%value(x_max)
```

The high and low limits of the extrema function are optional and the
start and end of the spline are used by default.

## License

[MIT](LICENSE) © John Alexiou