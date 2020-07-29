module stat
!    
!   modules for statistical analysis
!
    implicit none

contains
!
! determinant (use lapack routine)
!
subroutine det(m,n,A,detA)
    implicit none
    integer,intent(in) :: m, n
    double precision,intent(in) :: A(m,n)
    double precision,intent(out) :: detA

    double precision,allocatable :: dummy(:,:) ! destroyed
    integer,allocatable :: ipiv(:)
    integer :: info,i

    allocate(ipiv(min(m,n)))
    allocate(dummy(m,n))
    dummy = A
! step1. A = (P)LU
    call dgetrf(m, n, dummy, m, ipiv, info)
! step2. det(A) = det(P)*det(U) (det(L)=1)
    if(info.eq.0)then
        detA = 1.0d0
        do i = 1, min(m,n)
            if(ipiv(i).eq.i)then
                detA = detA * dummy(i,i)
            else
                detA = detA * (-dummy(i,i))
            endif
        enddo
    elseif(info>0)then
        detA = 0.0d0
    else
        print *, "Error in ",info,"th value"
    endif
    return
end subroutine det
!
! calcurate single linear regression yi=a+b*xi
!
subroutine reg(n,xin,yin,lnorm,a,b,cor,thr)
    implicit none
    integer,intent(in) :: n
    real,intent(in) :: xin(n),yin(n)
    logical,intent(in) :: lnorm
    real,intent(out) :: a,b,cor,thr

    real :: x(n),y(n)
    real :: xm,ym,xym,x2m
    real :: vx,vy,vxy
    real :: y_(n)
    real :: se
    integer :: i

    x = xin
    y = yin
    if(lnorm) then
        ! Normalize
        xm = 0.0
        ym = 0.0
        do i = 1, n
            xm = xm + x(i)
            ym = ym + y(i)
        end do
        xm = xm/n
        ym = ym/n
    
        vx = 0.0
        vy = 0.0
        do i = 1, n
            vx = vx + (x(i)-xm)**2
            vy = vy + (y(i)-ym)**2
        end do
        vx = vx/n
        vy = vy/n

        do i = 1, n
            x(i) = (x(i)-xm)/sqrt(vx)
            y(i) = (y(i)-ym)/sqrt(vy)
        end do
    endif

    ! Regression
    xm = 0.0
    ym = 0.0
    xym = 0.0
    x2m = 0.0
    do i = 1, n
        xm = xm + x(i)
        ym = ym + y(i)
        xym = xym + x(i)*y(i)
        x2m = x2m + x(i)*x(i)
    end do
    xm = xm/n
    ym = ym/n
    
    b = (xym - n*xm*ym)/(x2m - n*xm*xm)
    a = ym - b*xm

    ! t-test
    do i = 1, n
        y_(i) = a + b*x(i)
    end do
    se = 0.0
    do i = 1, n
        se = se + (y(i)-y_(i))**2
    end do
    se = sqrt(se/(n-2))

    thr = 0.0
    do i = 1, n
        thr = thr + (x(i)-xm)**2
    end do
    thr = b*sqrt(thr)/se

    vx = 0.0
    vy = 0.0
    vxy = 0.0
    do i = 1, n
        vx = vx + (x(i)-xm)**2
        vy = vy + (y(i)-ym)**2
        vxy = vxy + (x(i)-xm)*(y(i)-ym)
    end do
    cor = vxy/sqrt(vx*vy)

    return
end subroutine reg
!
! calcurate multiple linear regression yi=a+b1*xi1+b2*xi2+...
!
subroutine reg_mul(p,n,xin,yin,lnorm,a,b,cor,R2,thr)
    implicit none
    integer,intent(in) :: p,n
    real,intent(in) :: xin(p,n),yin(n)
    logical,intent(in) :: lnorm
    real,intent(out) :: a,b(p),cor(p),R2,thr(p) 
    !cor(p) is partial correlation, R2 is coefficient of determinant

    real :: x(p,n),y(n)
    real :: xm(p),ym
    real :: vx(p),vy,vxy(p),vxx
    double precision :: R(p+1,p+1),Ri(p,p)
    double precision :: detR,vdetR
    double precision :: detRi(p+1),detR1i(p) !cofactor
    integer :: i,j,k,ii,jj

    x = xin
    y = yin
    if(lnorm) then
! Normalize
        xm = 0.0
        ym = 0.0
        do i = 1, n
            do j = 1, p
                xm(j) = xm(j) + x(j,i)
            end do
            ym = ym + y(i)
        end do
        xm = xm/n
        ym = ym/n
    
        vx = 0.0
        vy = 0.0
        do i = 1, n
            do j = 1, p
                vx(j) = vx(j) + (x(j,i)-xm(j))**2
            end do
            vy = vy + (y(i)-ym)**2
        end do
        vx = vx/n
        vy = vy/n

        do i = 1, n
            do j = 1, p
                x(j,i) = (x(j,i)-xm(j))/sqrt(vx(j))
            end do
            y(i) = (y(i)-ym)/sqrt(vy)
        end do
    endif

! Regression
    xm = 0.0
    ym = 0.0
    do i = 1, n
        do j = 1, p
            xm(j) = xm(j) + x(j,i)
        end do
        ym = ym + y(i)
    end do
    xm = xm/n
    ym = ym/n
    !print*, "xm",xm,"ym",ym

    vx = 0.0
    vy = 0.0
    vxy = 0.0
    do i = 1, n
        do j = 1, p
            vx(j) = vx(j) + (x(j,i)-xm(j))**2
            vxy(j) = vxy(j) + (x(j,i)-xm(j))*(y(i)-ym)
        end do
        vy = vy + (y(i)-ym)**2
    end do
    vx = vx/n
    vy = vy/n
    vxy = vxy/n
    !print*, "vx",vx,"vy",vy,"vxy",vxy

!correlation matrix
    R = 1.0d0
    do j = 1, p
        R(1,j+1) = vxy(j)/(sqrt(vx(j))*sqrt(vy))
    end do

    do i = 1, p-1
        do j = i+1, p
            vxx = 0.0d0
            do k = 1, n
                vxx = vxx + (x(i,k)-xm(i))*(x(j,k)-xm(j))
            enddo
            vxx = vxx/n
            R(i+1,j+1) = vxx/(sqrt(vx(i))*sqrt(vx(j)))
        end do
    end do
    !print *, "R(1)",R(1,:)
    !print *, "R(2)",R(2,:)
    !print *, "R(3)",R(3,:)

    do i = 2, p+1
        do j = 1, i
            R(i,j) = R(j,i)
        end do
    end do
    !print *, "R(1)",R(1,:)
    !print *, "R(2)",R(2,:)
    !print *, "R(3)",R(3,:)

    call det(p+1,p+1,R,detR)
    !print *, "detR",detR
!co-factor
    Ri = R(2:,2:)
    !print *, "Ri(1)",Ri(1,:)
    !print *, "Ri(2)",Ri(2,:)
    call det(p,p,Ri,detRi(1))

    do k = 1, p
        ii = 1
        do i = 1, p+1
            if(i==k+1) cycle
            jj = 1
            do j = 1, p+1
                if(j==k+1) cycle
                Ri(ii,jj) = R(i,j)
                jj = jj + 1
            end do
            ii = ii + 1
        end do
        !print *, "Ri(1)",Ri(1,:)
        !print *, "Ri(2)",Ri(2,:)
        call det(p,p,Ri,detRi(k+1))
    end do
    !print *, "detRi", detRi
    do k = 1, p
        ii = 1
        do i = 2, p+1
            jj = 1
            do j = 1, p+1
                if(j==k+1) cycle
                Ri(ii,jj) = R(i,j)
                jj = jj + 1
            end do
            ii = ii + 1
        end do
        !print *, "Ri(1)",Ri(1,:)
        !print *, "Ri(2)",Ri(2,:)
        call det(p,p,Ri,detR1i(k))
        if(mod(k,2)==1)then
            detR1i(k) = -detR1i(k)
        endif
    end do
    !print *, "detR1i", detR1i

!verificate
    vdetR = R(1,1)*detRi(1)
    do k = 1, p
        vdetR = vdetR + R(1,k+1)*detR1i(k)
    end do
    !print *, "verify",detR-vdetR
! partial regression coefficient
    a = 0.0
    do k = 1, p
        b(k) = -sqrt(vy/vx(k))*detR1i(k)/detRi(1)
        a = a + b(k)*xm(k)
    end do
    a = ym - a

! partial correlation
    do k = 1, p
        cor(k) = -detR1i(k)/sqrt(detRi(1)*detRi(k+1))
    end do

! determinant coefficient
    R2 = 1.0 - detR/detRi(1)

! t-test
    do i = 1, p
        thr(i) = abs(cor(i))*sqrt((n-p-1)/(1-cor(i)**2))
    end do
    
    return
end subroutine reg_mul
end module stat