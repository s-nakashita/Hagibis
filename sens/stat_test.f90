program stat_test
    USE stat
    implicit none
    integer,parameter :: n=100
    real :: x(n),y(n),r(n)
    real :: a,b,cor,thr
    real :: y_(n)
    integer :: i

    call random_seed()
    call random_number(r)
    do i = 1, n
        x(i) = i
    enddo
    do i = 1, n
        y(i) = 2.0 + 0.5*x(i) + (2*r(i) - 1)
    enddo 

    open(10,file="test.txt",status="new")
    do i = 1, n
        write(10,'(2f7.3)') x(i), y(i)
    end do
    close(10)
    
    call reg(n,x,y,.False.,a,b,cor,thr)
    print *, "a=", a, "b=", b
    print *, "cor=", cor, "thr=", thr

    do i = 1, n
        y_(i) = a + b*x(i)
    end do

    open(10,file="test-r.txt",status="new")
    do i = 1, n
        write(10,'(2f7.3)') x(i), y_(i)
    end do
    close(10)

end program