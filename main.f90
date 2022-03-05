module config
    implicit none
    public

    !Integration precision
    integer, parameter :: prec = 1000
    !Threads count
    integer, parameter :: thr = 4
!---------------------------------------------------
    !Draw smoothness
    integer, parameter :: dots = 1000
    !Segment boundaries
    real(8), parameter :: fl = 0.d0, cl = 5.d0
!---------------------------------------------------
end module config

module com
    implicit none
    public

    real(8) :: gr
end module com


program main
    use config
    use omp_lib
    implicit none

    real(8), parameter :: pi = 4*atan(1.d0)

    real(8) :: arg
    integer :: i

    open(unit = 1, file = "output.dat")
        do i = 1, dots+1
            arg = fl + (i-1)*((cl-fl)/(dots))
            write(1, *) arg, phi(arg, 0.d0), phi(arg, pi/6), phi(arg, pi/4), phi(arg, pi/3), phi(arg, pi/2), phi(arg, pi)
        end do
    close(1)

contains
    real(8) function lim0(lr)
        use com
        real(8), intent(in) :: lr
        lim0 = 5.d-1*(lr**2.d0)*exp(-lr)/max(gr, lr)
    end function lim0

    real(8) function lim2(lr)
        use com
        real(8), intent(in) :: lr
        lim2 = 1d-1*(lr**2.d0)*exp(-lr)*(min(gr, lr)**2)/(max(gr, lr)**3)
    end function lim2

    real(8) function inf0(z)
        use com
        real(8), intent(in) :: z
        inf0 = 5.d-1*exp(-1.d0/z)*(z**-4.d0)*min(1.d0/gr, z)
    end function inf0

    real(8) function inf2(z)
        use com
        real(8), intent(in) :: z
        inf2 = 1d-1*(z**-4.d0)*exp(-1.d0/z)*(min(1.d0/gr, z)**3)/(max(1.d0/gr, z)**2)
    end function inf2

    real(8) function phi(rh, theta)
        use com

        real(8), intent(in) :: rh, theta

        abstract interface
            real(8) function func(x)
            real(8), intent(in) :: x
            end function func
        end interface

        procedure(func), pointer :: ptr1, ptr2
        gr = rh
        phi = 0

        ptr1 => lim0
        ptr2 => inf0
        phi = phi + integ(ptr1, 0.d0, 1.d0, prec) + integ(ptr2, 1.d-4, 1.d0, prec)

        ptr1 => lim2
        ptr2 => inf2
        phi = phi + (3*cos(theta)**2-1)*(integ(ptr1, 0.d0, 1.d0, prec) + integ(ptr2, 1.d-4, 1.d0, prec))

    end function phi

    real(8) function integ(fptr, bt, tp, cnt)
        !10 point Gauss-Legendre algorithm

        abstract interface
            real(8) function func(x)
            real(8), intent(in) :: x
            end function func
        end interface

        procedure(func), pointer, intent(in) :: fptr

        real(8), intent(in) :: bt, tp
        integer, intent(in) :: cnt

        real(8) :: x, half, mid, step
        real(8) :: tmp, er, next
        real(8), dimension(10) :: a, w
        integer(kind=8) :: i, j

        data w /0.066671344308688d0, 0.149451349150581d0, 0.219086362515982d0, &
                0.269266719309996d0, 0.295524224714753d0, 0.295524224714753d0, &
                0.269266719309996d0, 0.219086362515982d0, 0.149451349150581d0, &
                0.066671344308688d0/
        data a /0.97390652851717174d0, 0.86506336668898454d0, 0.67940956829902444d0, &
                0.43339539412924721d0, 0.14887433898163122d0, -0.14887433898163119d0, &
                -0.43339539412924721d0, -0.67940956829902444d0, -0.86506336668898454d0, &
                -0.97390652851717174d0/

        step = (tp-bt)/cnt
        half = step/2
        integ=0
        !$omp parallel default(none) private(x, mid, next, er, tmp) &
        !$omp shared(half, cnt, bt, step, a, w, fptr) num_threads(thr) reduction(+: integ)
            er = 0
            !$omp do
            do i = 1, cnt
                x = bt + (i-1)*step
                mid = x + half
                do j = 1, 10
                    next = half*w(j)*fptr(a(j)*half+mid) - er
                    tmp = integ + next
                    er = (tmp - integ) - next
                    integ = tmp
                end do
            end do
            !$omp end do

        !$omp end parallel

    end function integ

end program main
