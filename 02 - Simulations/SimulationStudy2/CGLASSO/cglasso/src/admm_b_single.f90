!p,q,xtx_n,xtr_n,B_n,Tht_n,wB,lambda(h),maxit_bcd,thr_bcd,conv,subrout,nnit,dfB(:,h,k),trace
subroutine admm_b_single(p,q,qp,xtx_n,xtr_n,B,Tht,W,lmb,maxit,thr,rho,trace,conv,subrout,nit,df)
implicit none
double precision :: xtx_n(q,q),xtr_n(q,p),B(q,p),Tht(p,p),W(q,p),lmb,thr,rho
integer :: p,q,qp,maxit,trace,conv,subrout,nit(2),df(p+1)
!internal variables
double precision :: Gamma(q,p),U(q,p),xtr_n_theta(q,p),theta_xtx_n(qp,qp),Iqp(qp,qp)
double precision :: theta_xtx_n_inv(qp,qp),Gamma2(q,p),tempY(qp),s2,r2,eps1,eps2!,Db
integer :: i,criterion,info,j

Gamma = 0.d0
U = 0.d0
info = 0

Iqp = 0.d0
do i = 1, qp
    Iqp(i,i) = 1.d0
end do

theta_xtx_n = 0.d0
call DGEMM('N', 'N', q, p, p, 1.d0, xtr_n, q, Tht, p, 0.d0, xtr_n_theta, q)
call kron(p, q, xtx_n, Tht, theta_xtx_n)
theta_xtx_n_inv = theta_xtx_n + rho * Iqp
call dpotrf('U', qp, theta_xtx_n_inv, qp, info)
! or
!call inv(qp, theta_xtx_n + rho * Iqp, theta_xtx_n_inv, info)

criterion = 1
i = 0
do while ((criterion.eq.1).and.(i.lt.maxit))
    call rchkusr()

    ! update values
    i = i + 1
    Gamma2 = Gamma

    ! update 'Beta'
    call mattovec(p, q, xtr_n_theta + rho * (Gamma - U), tempY)
    call dpotrs('U', qp, 1, theta_xtx_n_inv, qp, tempY, qp, info)
    ! or one of the following
    !call DGEMV('N', qp, qp, 1.d0, theta_xtx_n_inv, qp, tempY, 1, 0.d0, tempX, 1)
    !call solve(theta_xtx_n + rho * Iqp, tempY, qp, info)
    
    call vectomat(p, q, tempY, B)

    ! update 'Gamma'
    Gamma = B + U
    call softmatrix_b(Gamma, lmb * W / rho, p, q)

    ! update 'U'
    U = U + B - Gamma

    ! stopping criterion
    !Db = sum(abs(Gamma2 - Gamma)) / sum(abs(Gamma))
    s2 = sqrt(sum((rho*(Gamma2 - Gamma))**2))
    r2 = sqrt(sum((B - Gamma)**2))
    eps1 = thr * max(sqrt(sum(B**2)), sqrt(sum(Gamma**2)))
    eps2 = thr * sqrt(sum(U**2))
    criterion = 0
    !if(Db.gt.thr) then
    !    criterion = 1
    !end if
    if ((r2.ge.eps1).or.(s2.ge.eps2)) then
        criterion = 1
    end if

    !if(trace.eq.9) call admm_trace_1(i, Db, rho)
    !if((trace.eq.2).and.(criterion.eq.0)) call admm_trace_1(i, Db, rho)
    if(trace.eq.9) call admm_trace_2(i, r2, eps1, s2, eps2, rho)
    if((trace.eq.2).and.(criterion.eq.0)) call admm_trace_2(i, r2, eps1, s2, eps2, rho)
end do
!call DGEMV('T', q, p, 1.d0, B2, q, xm, 1, 0.d0, ym2, 1)
B = Gamma
nit(1) = 1
nit(2) = i
do i = 1, p
    df(i) = 0
    do j = 1, q
        if(B(j, i).ne.0.d0) df(i) = df(i) + 1
    end do
end do
df(p + 1) = sum(df(1:p))
if(i.eq.maxit) then
    conv = 1
    subrout = 1
    return
end if
end subroutine admm_b_single

subroutine kron(p, q, XtX, theta, thetaXtX)
implicit none
double precision :: XtX(q,q), theta(p,p), thetaXtX(p*q,p*q)
integer :: p, q
integer :: i, j
do i = 1, p
    do j = i, p
        if(theta(i,j).ne.0) then
            !thetaXtX(row + 1:q, col + 1:q) = theta(i,j) * XtX
            thetaXtX(((i-1)*q+1):(i*q), ((j-1)*q+1):(j*q)) = theta(i,j) * XtX
            if(j.ne.i) thetaXtX(((j-1)*q+1):(j*q), ((i-1)*q+1):(i*q)) = thetaXtX(((i-1)*q+1):(i*q), &
                ((j-1)*q+1):(j*q))
        end if
    end do
end do
end subroutine kron

subroutine mattovec(p, q, X, Y)
implicit none
double precision :: X(q,p), Y(q*p)
integer :: p, q
integer :: i, row
row = 0
do i = 1, p
    Y(row + 1:q) = X(:,i)
    row = row + q
end do
end subroutine mattovec

subroutine vectomat(p, q, Y, X)
implicit none
double precision :: X(q,p), Y(q*p)
integer :: p, q
integer :: i, row
row = 0
do i = 1, p
    X(:,i) = Y(row + 1:q)
    row = row + q
end do
end subroutine vectomat

subroutine softmatrix_b(S, Tau, p, q)
implicit none
double precision :: S(q,p), Tau(q,p)
integer :: p, q
!internal variables
integer :: i, j
! loop over all elements
do i = 1, q
    do j = 1, p
        ! soft threshold each element
        call softc(S(i, j), Tau(i, j))
    end do
end do
end subroutine softmatrix_b

subroutine solve(A, b, n, info)
integer :: info, n
double precision :: A(n, n), b(n, 1)
call DPOSV('U', n, 1, A, n, b, n, info)
if (info.ne.0) return
end subroutine solve
