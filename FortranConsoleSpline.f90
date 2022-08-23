
    program FortranConsoleSpline
    use mod_splines
    implicit none

    ! Variables
    real(wp), allocatable :: xi(:), yi(:), h, x, y, yp, x0
    type(spline) :: sp
    integer :: i, n
    ! compile with /fpconstant
    xi = [0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0]
    yi = [18.0,18.4921875,18.9375,19.2890625,19.5,19.5234375,19.3125,18.8203125,18.0]

    print *, 'Cubic Spline Interpolation Demo'
    print *, 'n=', size(yi)
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
    
    
    end program FortranConsoleSpline

