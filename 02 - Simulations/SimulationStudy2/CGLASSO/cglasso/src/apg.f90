subroutine grad_b(q,p,n,m,f,A,b,x,Tht,g)
implicit none
integer :: q,p,n
double precision :: m(q),f(q),A(n,q),b(n,p),x(q,p),Tht(p,p),g(q,p)
!internal variables
integer :: i
double precision :: MU(n,p),R(n,p),tmp1(q,p),tmp2(q,p)
call DGEMM('N', 'N', n, p, q, 1.d0, A, n, x, q, 0.d0, MU, n)
R = MU - b
call DGEMM('T', 'N', q, p, n, 1.d0, A, n, R, n, 0.d0, tmp1, q)
call DGEMM('N', 'N', q, p, p, 1.d0, tmp1, q, Tht, p, 0.d0, tmp2, q)
do i = 1, p
    g(:,i) = f * tmp2(:,i) / m
end do
end subroutine grad_b

subroutine prox_lasso(p,q,k,x,lam)
implicit none
integer :: p,q,k
double precision :: x(q,p),lam(q,p),temp(q/k,p/k),l(q/k,p/k)
!internal variables
integer :: i,q1,p1
q1 = q / k
p1 = p / k
! loop over all elements
do i = 1, k
    temp = x( ((i-1)*q1 + 1):(i*q1), ((i-1)*p1 + 1):(i*p1) )
    l = lam( ((i-1)*q1 + 1):(i*q1), ((i-1)*p1 + 1):(i*p1) )
    call softmatrix_b(temp, l, p1, q1)
    x( ((i-1)*q1 + 1):(i*q1), ((i-1)*p1 + 1):(i*p1) ) = temp
end do
end subroutine prox_lasso

subroutine prox_grouplasso(p,q,k,x,lam)
implicit none
integer :: p,q,k
double precision :: x(q,p),lam(q,p)
!internal variables
integer :: i,q1,p1
double precision :: nrmB(q/k,p/k),temp(q/k,p/k),l(q/k,p/k)
q1 = q / k
p1 = p / k
nrmB = 0.d0
! loop over all elements
do i = 1, k
    temp = x( ((i-1)*q1 + 1):(i*q1), ((i-1)*p1 + 1):(i*p1) )
    nrmB = nrmB + temp**2
end do
nrmB = sqrt(nrmB)
do i = 1, k
    temp = x( ((i-1)*q1 + 1):(i*q1), ((i-1)*p1 + 1):(i*p1) )
    l = lam( ((i-1)*q1 + 1):(i*q1), ((i-1)*p1 + 1):(i*p1) )
    call softmatrix_b_joint(temp, l, nrmB, p1, q1)
    x( ((i-1)*q1 + 1):(i*q1), ((i-1)*p1 + 1):(i*p1) ) = temp
end do
end subroutine prox_grouplasso

subroutine prox_sparsegrouplasso(p,q,k,x,lam,alpha)
implicit none
integer :: p,q,k
double precision :: x(q,p),lam(q,p),alpha
!internal variables
integer :: i,q1,p1
double precision :: nrmB(q/k,p/k),temp(q/k,p/k),l(q/k,p/k)
q1 = q / k
p1 = p / k
nrmB = 0.d0
! loop over all elements
do i = 1, k
    temp = x( ((i-1)*q1 + 1):(i*q1), ((i-1)*p1 + 1):(i*p1) )
    l = lam( ((i-1)*q1 + 1):(i*q1), ((i-1)*p1 + 1):(i*p1) )
    call softmatrix_b(temp, alpha * l, p1, q1)
    nrmB = nrmB + temp**2
    x( ((i-1)*q1 + 1):(i*q1), ((i-1)*p1 + 1):(i*p1) ) = temp
end do
nrmB = sqrt(nrmB)
do i = 1, k
    temp = x( ((i-1)*q1 + 1):(i*q1), ((i-1)*p1 + 1):(i*p1) )
    l = lam( ((i-1)*q1 + 1):(i*q1), ((i-1)*p1 + 1):(i*p1) )
    call softmatrix_b_joint(temp, (1.d0 - alpha) * l, nrmB, p1, q1)
    x( ((i-1)*q1 + 1):(i*q1), ((i-1)*p1 + 1):(i*p1) ) = temp
end do
end subroutine prox_sparsegrouplasso

subroutine apg(p,q,n,k,nk,fk,A,b,Tht,weights,lambda,alpha,xm,ym,x,beta,maxit,thr,trace,df,nit,conv)
implicit none
integer :: p,q,n,k,maxit,trace,df((p/k)+1,k),nit,conv
double precision :: nk(q),fk(q),A(n,q),b(n,p),Tht(p,p),weights(q,p),lambda,alpha,xm(q/k,k), &
    ym(p/k,k),x(q,p),beta(0:(q/k),p/k,k),thr
!internal variables
integer :: i,q1,p1,i2,i3
double precision :: y(q,p),g(q,p),nrmg,theta,t,xhat(q,p),ghat(q,p),xold(q,p),yold(q,p),pen(q,p), &
    tmp1,tmp2,err1,gold(q,p),that,tmp3(q/k,p/k),tmp4(p/k)

!USE_RESTART = 1 ! use adaptive restart scheme
!ALPHA = 1.0d1 ! step-size growth factor
!BETA = 0.5d0 ! step-size shrinkage factor
!USE_GRA = 0 ! if true uses UN-accelerated proximal gradient descent (typically slower)
!STEP_SIZE = 0 ! starting step-size estimate, if not set then apg makes initial guess
!FIXED_STEP_SIZE = 0

q1 = q / k
p1 = p / k

y = x
call grad_b(q, p, n, nk, fk, A, b, y, Tht, g)
nrmg = sqrt(sum(g**2))

theta = 1.d0

! Initial step size
! Barzilai-Borwein step-size initialization:
t = 1.d0 / nrmg
xhat = x - t * g
call grad_b(q, p, n, nk, fk, A, b, xhat, Tht, ghat)
t = abs(sum((x - xhat)*(g - ghat)) / sum((g - ghat)**2))

! Main loop
do i = 1, maxit
    call rchkusr()

    xold = x
    yold = y

    ! The proximal gradient step (update x)
    x = y - t * g
    pen = t * lambda * weights
    call prox_sparsegrouplasso(p, q, k, x, pen, alpha)

    ! The error for the stopping criterion
    tmp1 = sqrt(sum((y - x)**2))
    tmp2 = sqrt(sum(x**2))
    err1 = tmp1 / max(1.d0, tmp2)
    if (trace.eq.2) then
        call apg_trace(i, err1, t)
    end if
    if (err1.lt.thr) then
        exit
    end if

    ! Update theta for acceleration
    theta = 2.d0 / (1.d0 + sqrt(1.d0 + 4.d0 / (theta**2)))

    ! Update y
    if (sum((y - x) * (x - xold)).gt.0.d0) then
      x = xold
      y = x
      theta = 1.d0
    else
      y = x + (1.d0 - theta) * (x - xold)
    end if

    ! New gradient
    gold = g
    call grad_b(q, p, n, nk, fk, A, b, y, Tht, g)

    ! Update stepsize by TFOCS-style backtracking
    that = 0.5d0 * sum((y - yold)**2) / abs(sum((y - yold) * (gold - g)))
    tmp1 = max(0.5d0 * t, that)
    t = min(1.01d0 * t, tmp1)
end do
nit = i
if(i.eq.maxit) then
    conv = 1
end if
do i = 1, k
    tmp3 = x( ((i-1)*q1 + 1):(i*q1), ((i-1)*p1 + 1):(i*p1) )
    beta(1:q1,:,i) = tmp3

    tmp4 = ym(:,i)
    call DGEMV('T', q1, p1, -1.d0, tmp3, q1, xm(:,i), 1, 1.d0, tmp4, 1)
    beta(0,:,i) = tmp4
    do i2 = 1, p1
        df(i2, i) = 0
        do i3 = 1, q1
            if(beta(i3, i2, i).ne.0.d0) df(i2, i) = df(i2, i) + 1
        end do
    end do
    df(p1 + 1, i) = sum(df(1:p1, i))
end do
end subroutine apg
