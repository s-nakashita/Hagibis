program det_test
    use stat
    implicit none
    integer,parameter :: m=3,n=3
    double precision :: A(m,n),Ai(m-1,n-1)
    double precision :: detA
    integer :: i,j,ii,jj
    
    A(1,:) = (/-1.0d0,2.0d0,3.0d0/)
    A(2,:) = (/4.0d0,5.0d0,6.0d0/)
    A(3,:) = (/2.0d0,4.0d0,6.0d0/)

    call det(m,n,A,detA)
    print*, detA
    ii = 1
    do i = 1, m
        if(i==2) cycle
        jj = 1
        do j = 1, n
            if(j==2) cycle
            print *, ii, jj
            print *, i,j
            Ai(ii,jj) = A(i,j)
            jj = jj + 1
        end do
        ii = ii + 1
    end do
    print*,Ai
end program