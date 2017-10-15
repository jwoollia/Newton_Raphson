program newton_raphson
implicit none
integer, parameter :: dp=kind(0.0d0)
real(kind=dp) :: x,fx,fdash,tol,xl,fl,xu,fu,dx,dxold,c
character(len=1) :: aflag
integer :: i
!
tol=1.d-6
print *, " ... tolerance is set to ",tol
print *, " ... do you wish to change tolerance? (y/n) ... "
read(*,*) aflag
if(aflag=="y".or.aflag=="Y") then
    print *, " ... enter value ... "
    read(*,*) tol
end if
values: do
    i=0
    dx=huge(x)
    print '(a)', " ... enter target value ... "
    read(*,*) c
    print '(a)', " ... enter initial estimate of root ... "
    read(*,*) x
    solve: do while (dx>tol)
        i=i+1
        xl=x-2.d0*tol
        fl=func(xl,c)
        xu=x+2.d0*tol
        fu=func(xu,c)
        fdash=(fu-fl)/(xu-xl)
        fx=func(x,c)
        dx=fx/fdash
        x=x-dx
        dx=abs(dx)
        if(i>3.and.dx>dxold) then ! not converging
            print '(a)', " ... algorithm is diverging ... "
            print '(/a)', " ... do you wish another target value? (y/n) ... "
            read(*,*) aflag
            if(aflag=="n".or.aflag=="N") exit values
            cycle values
        end if
        dxold=dx
        print '(f12.6)', x
    end do solve
    print '(a,f12.6)', " ... convergence criterion reached with root ... ",x
    print '(a,i4)', " ... convergence achieved after cycle ",i
    print '(/a)', " ... do you wish another target value? (y/n) ... "
    read(*,*) aflag
    if(aflag=="n".or.aflag=="N") exit values
end do values
!
contains
!
real(kind=dp) function func(u,b)
real(kind=dp) :: u,b
func=u*u-b
return
end function func
!
end program newton_raphson