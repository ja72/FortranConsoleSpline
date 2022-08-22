    module mod_splines
    use, intrinsic :: iso_fortran_env
    implicit none

    integer, parameter :: wp = real64
    real(wp), parameter :: big = 1e30_wp, tiny=1/big

    type :: spline
        real(wp), allocatable :: x(:), y(:), y2(:)
    contains
        procedure :: indexof => sp_index_of_x
        procedure :: value => sp_interpolate_value
        procedure :: slope => sp_interpolate_slope
        procedure :: slope2 => sp_interpolate_slope2
        procedure :: extrema => sp_find_local_extrema
    end type

    interface spline
    module procedure :: sp_calculate_from_data
    end interface

    contains

    pure function sp_calculate_from_data(x,y,y1_slope,yn_slope) result(sp)
    ! =====================================================
    ! Input x and y=f(x), n (dimension of x,y), (Ordered)
    ! y1 and yn are the first derivatives of f in the 1st point and the n-th
    ! Output: array y2(n) containing second derivatives of f(x_i)
    ! =====================================================
    
    type(spline) :: sp
    real(wp), intent(in) :: x(:), y(:)
    real(wp) :: y2(size(y))
    real(wp), optional, intent(in) :: y1_slope, yn_slope
    real(wp):: p, qn, sig, un, u(size(y))
    INTEGER:: n, i, j
    
        n = size(y)
        IF (present(y1_slope)) THEN    ! natural spline conditions
            y2(1) = -0.5
            u(1) = (3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-y1_slope)
        ELSE
            y2(1) = 0
            u(1) = 0
        END IF

        DO i = 2, n-1                            ! tridiag. decomposition
            sig = (x(i)-(i-1))/(x(i+1)-x(i-1))
            p = sig*y2(i-1)+2.
            y2(i) = (sig-1.)/p
            u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
        END DO

        IF (present(yn_slope)) THEN   ! natural spline conditions
            qn = 0.5
            un=(3./(x(n)-x(n-1)))*(yn_slope-(y(n)-y(n-1))/(x(n)-x(n-1)))
        ELSE
            qn = 0
            un = 0
        END IF

        y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)

        DO j = n-1, 1, -1          !  backwards substitution tri-diagonale
            y2(j) = y2(j)*y2(j+1)+u(j)
        END DO
        
        sp%x = x
        sp%y = y
        sp%y2 = y2
        
        RETURN
    end function sp_calculate_from_data
    
    elemental function sp_index_of_x(sp,x) result(k_low)
    class(spline), intent(in) :: sp
    real(wp), intent(in) :: x
    integer:: n, k, k_low, k_high
        n = size(sp%y)
        k_low = 1
        k_high = n
        if(x<sp%x(k_low)) then
            return
        elseif (x>sp%x(k_high)) then
            k_low = k_high-1
            return
        end if
        do while(k_high - k_low > 1) 
            k = (k_high + k_low) / 2
            IF (sp%x(k) > x) THEN
                k_high = k
            ELSE
                k_low = k
            END IF
        end do
    end function
    
    elemental function sp_interpolate_value(sp,x) result(y)
    ! =====================================================
    ! Subroutine that does the actual interpolation
    ! Input arrays of x_in and y_in=f(x), spline_res is the result of
    ! the 'spline' subroutine, x is the corresponding value we are looking for
    ! i.e. (time_at_max in hubble), y is the output result
    ! =====================================================
    class(spline), intent(in) :: sp
    real(wp), intent(in) :: x
    real(wp) :: y
    integer:: n, k 
    real(wp):: a, b, c, d, h, t
        n = size(sp%y)
        k= sp%indexof(x)
        h = sp%x(k+1) - sp%x(k)
        IF (h == 0) error STOP "Bad x input"
        t = (x-sp%x(k))/h
        a = 1-t
        b = t
        if( x>=sp%x(k) .and. x<=sp%x(k+1)) then
            ! Cubic inside the interval
            c = (a**3-a)*(h**2)/6
            d = (b**3-b)*(h**2)/6
        else
            ! Linear outside the interval
            c = 0.0_wp
            d = 0.0_wp
        end if
        y = a*sp%y(k)+b*sp%y(k+1)+c*sp%y2(k)+d*sp%y2(k+1)

        RETURN
    end function sp_interpolate_value
    
    elemental function sp_interpolate_slope(sp,x) result(yp)
    ! =====================================================
    ! Subroutine that does the actual interpolation
    ! Input arrays of x_in and y_in=f(x), spline_res is the result of
    ! the 'spline' subroutine, x is the corresponding value we are looking for
    ! i.e. (time_at_max in hubble), yp is the output result slope
    ! =====================================================
    class(spline), intent(in) :: sp
    real(wp), intent(in) :: x
    real(wp) :: yp
    integer:: n, k 
    real(wp):: a, b, c, d, h, t
        n = size(sp%y)
        k= sp%indexof(x)
        
        h = sp%x(k+1) - sp%x(k)
        IF (h == 0) error STOP "Bad x input"
        t = (x-sp%x(k))/h
        a = -1/h
        b = 1/h
        if( x>=sp%x(k) .and. x<=sp%x(k+1)) then
            ! Cubic inside the interval
            c = (1-3*(1-t)**2)*(h/6)
            d = (3*t**2-1)*(h/6)
        else
            ! Linear outside the interval
            c = 0.0_wp
            d = 0.0_wp
        end if
        yp = a*sp%y(k)+b*sp%y(k+1)+c*sp%y2(k)+d*sp%y2(k+1)

        RETURN
    end function sp_interpolate_slope
    
    elemental function sp_interpolate_slope2(sp,x) result(yp2)
    ! =====================================================
    ! Subroutine that does the actual interpolation
    ! Input arrays of x_in and y_in=f(x), spline_res is the result of
    ! the 'spline' subroutine, x is the corresponding value we are looking for
    ! i.e. (time_at_max in hubble), yp is the output result 2nd slope
    ! =====================================================
    class(spline), intent(in) :: sp
    real(wp), intent(in) :: x
    real(wp) :: yp2
    integer:: n, k 
    real(wp):: a, b, c, d, h, t
        n = size(sp%y)
        k= sp%indexof(x)
        
        h = sp%x(k+1) - sp%x(k)
        IF (h == 0) error STOP "Bad x input"
        t = (x-sp%x(k))/h
        a = 0.0_wp
        b = 0.0_wp
        if( x>=sp%x(k) .and. x<=sp%x(k+1)) then
            ! Cubic inside the interval
            c = 1-t
            d = t
        else
            ! Linear outside the interval
            c = 0.0_wp
            d = 0.0_wp
        end if
        yp2 = a*sp%y(k)+b*sp%y(k+1)+c*sp%y2(k)+d*sp%y2(k+1)
        RETURN
    end function sp_interpolate_slope2
    
    pure function sp_find_local_extrema(sp, x_low, x_high) result(x)
    class(spline), intent(in) :: sp
    real(wp) :: x
    real(wp), intent(in), optional :: x_low, x_high
    integer :: n, k1, k2
    real(wp) :: x1, x2, yp1, yp2, h, tol, yp
        n = size(sp%y)
        if(present(x_low)) then
            x1 = x_low
        else
            x1 = sp%x(1)
        end if
        if(present(x_high)) then
            x2 = x_high
        else
            x2 = sp%x(n)
        end if
        h = x2 - x1
        tol = h/(2**23)
        yp1 = sp_interpolate_slope(sp, x1)
        yp2 = sp_interpolate_slope(sp, x2)
        
        if( yp1*yp2 > 0 ) then
            ! no solution
            if( yp1>0 ) then
                x = big
            else
                x = tiny
            end if
        end if
        
        do while (x2-x1>tol)
            x = (x1+x2)/2
            yp = sp_interpolate_slope(sp, x)
            if( yp1*yp > 0) then
                x1 = x
                yp1 = yp
            else
                x2 = x
                yp2 = yp
            end if
        end do
        
    end function

    end module mod_splines