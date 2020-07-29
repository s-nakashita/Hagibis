program stat_test2
    USE stat
    implicit none
    integer,parameter :: n=100,p=3
    real :: x(p,n),y(n),r(n)
    real :: a,b(p)
    real :: cor(p),R2,thr(p)
    real :: y_(n)
    integer :: i

    call random_seed()
    call random_number(r)
    do i = 1, n
        x(1,i) = (real(i)-50.0)**3/(10.0**3)
        x(2,i) = (real(i)-50.0)*(real(i)-50.0)/10.0/10.0
        x(3,i) = (real(i)-50.0)/10.0
    enddo
    do i = 1, n
        y(i) = 2*x(1,i) + 5*x(2,i) + 3*x(3,i) + (2*r(i) - 1)
    enddo 

    open(10,file="test.txt",status="new")
    do i = 1, n
        write(10,'(2f9.3)') x(3,i), y(i)
    end do
    close(10)
    
    call reg_mul(p,n,x,y,.False.,a,b,cor,R2,thr)
    print *, "a=", a, "b(1)=", b(1), "b(2)=", b(2), "b(3)=", b(3)
    print *, "cor=", cor, "thr=", thr
    print *, "R2=", R2

    do i = 1, n
        y_(i) = a + b(1)*x(1,i) + b(2)*x(2,i) + b(3)*x(3,i)
    end do

    open(10,file="test-r.txt",status="new")
    do i = 1, n
        write(10,'(2f9.3)') x(3,i), y_(i)
    end do
    close(10)

end program