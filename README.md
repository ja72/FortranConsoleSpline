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

Additionally, you can interpolate the slope and the second derivative using the
`slope(x)` and `slope2(x)` functions. For example

```fortran
y = sp%value(x)
yp = sp%slope(x)
ypp = sp%slope2(x)
```

As a note, if the evaluation point `x` is outside the defined range for the
spline, then a linear interpolation is done from the ends of the spline.

Finally, if a local extrema point is needed (local minimum or maximum) a bisection
method is used to find the `x` value where the slope is zero using
the `extrema(x_low,x_high)` function

```fortran
x_max = sp%extrema()
y_max = sp%value(x_max)
```

The high and low limits of the extrema function are optional, and the
start and end of the spline are used by default.

## Example

The included example program below defines a spline using 9 points
between 0 and 2 with an interval of 0.25.

The interpolation output has an interval of 0.20.

In the end, the local maximum is found.

### Program


```fortran
! Variables
real(wp), allocatable :: xi(:), yi(:), h, x, y, yp, x0
type(spline) :: sp
integer :: i, n
! compile with /fpconstant
xi = [0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0]
yi = [18.0,18.4921875,18.9375,19.2890625,19.5,19.5234375,19.3125,18.8203125,18.0]

print *, 'Cubic Spline Interpolation Demo'
    
n = 11
h = (xi(size(xi))-xi(1))/(n-1)
x0 = xi(1)
sp = spline(xi, yi)
print *, ""
print '(1x,a6,1x,a18,1x,a18,1x,a18)', "Index", "x", "y", "yp"
do i=0,n-1
    x = x0 + i*h
    y = sp%value(x)
    yp = sp%slope(x)
    print '(1x,i6,1x,g18.11,1x,g18.6,1x,g18.6)', i, x, y, yp
end do
print *, ""
    
x = sp%extrema()
i = sp%indexof(x)
y = sp%value(x)
yp = sp%slope(x)
    
print *, "Local Extrema"
print '(1x,a6,1x,a18,1x,a18,1x,a18)', "Index", "x", "y", "yp"
print '(1x,i6,1x,g18.11,1x,g18.6,1x,g18.6)', i, x, y, yp
```

### Output

The output of the example is as follows:

```text
 Cubic Spline Interpolation Demo

  Index                  x                  y                 yp
      0   0.0000000000            18.0000            2.06799
      1  0.20000000000            18.4009            1.87745
      2  0.40000000000            18.7637            1.73300
      3  0.60000000000            19.0861            1.47687
      4  0.80000000000            19.3398           0.943939
      5   1.0000000000            19.5000           0.209478
      6   1.2000000000            19.5304          -0.461936E-01
      7   1.4000000000            19.4106          -0.938651
      8   1.6000000000            19.1328           -1.85224
      9   1.8000000000            18.6726           -3.07239
     10   2.0000000000            18.0000           -3.50827

 Local Extrema
  Index                  x                  y                 yp
      5   1.1857554913            19.5308           0.738816E-07
```

## License

[MIT](LICENSE) © John Alexiou
