!p,q,xtx_n,xtr_n,B_n,Tht_n,wB,lambda(h),maxit_bcd,thr_bcd,conv,subrout,nnit,dfB(:,h,k),trace
subroutine admm_multilasso(p,q,qp,xtx_n,xtr_n,B,Tht,W,lmb,maxit,thr,conv,subrout,nit,df,trace,rho)
implicit none
double precision :: xtx_n(q,q),xtr_n(q,p),B(0:(q-1),p),Tht(p,p),W(q,p),lmb,thr
integer :: p,q,qp,maxit,conv,subrout,nit(2),df(p+1),trace
!internal variables
double precision :: tau_inc,tau_dec,tol_abs,tol_rel,rho
double precision :: Bold(q,p),x(q,p),u(q,p),xtr_n_theta(q,p),theta_xtx_n(qp,qp)
double precision :: Iqp(qp,qp),theta_xtx_n_inv(qp,qp),tempY(qp)!,tempX(qp)
double precision :: s2,r2,mu,eps1,eps2,rhonew,W2(q,p),B2(q,p),Db
integer :: i,criterion,update,info,j
tau_inc = 2.d0
tau_dec = 2.d0
tol_abs = thr
tol_rel = thr
eps1 = thr
eps2 = thr
mu = 10.d0
!rho = qp !2.d0
B2(1,:) = B(0,:)
B2(2:q,:) = B(1:(q-1),:)
x = B2
u = 0.d0
theta_xtx_n = 0.d0
call DGEMM('N', 'N', q, p, p, 1.d0, xtr_n, q, Tht, p, 0.d0, xtr_n_theta, q)
call kron(p, q, xtx_n, Tht, theta_xtx_n)
Iqp = 0.d0
do i = 1, qp
    Iqp(i,i) = 1.d0
end do
do i = 1, p
    W2(:,i) = W(:,i) * Tht(i,i)
end do
criterion = 1
update = 1
info = 0
i = 0
do while ((criterion.eq.1).and.(i.lt.maxit))
    call rchkusr()

    ! update values
    i = i + 1
    Bold = B2

    ! update 'x'
    !if(update.eq.1) call inv(qp, theta_xtx_n + rho * Iqp, theta_xtx_n_inv, info)
    
    if(update.eq.1) then
        theta_xtx_n_inv = theta_xtx_n + rho * Iqp
        call dpotrf('U', qp, theta_xtx_n_inv, qp, info)
    end if
    call mattovec(p, q, xtr_n_theta + rho * B2 - u, tempY)
    !call DGEMV('N', qp, qp, 1.d0, theta_xtx_n_inv, qp, tempY, 1, 0.d0, tempX, 1)
    !call solve(theta_xtx_n + rho * Iqp, tempY, qp, info)
    call dpotrs('U', qp, 1, theta_xtx_n_inv, qp, tempY, qp, info)
    !call vectomat(p, q, tempX, x)
    call vectomat(p, q, tempY, x)

    ! update 'B'
    B2 = u + rho * x
    call softmatrixc2(B2, lmb * W2, p, q)
    B2 = B2 / rho

    ! update 'u'
    u = u + rho * (x - B2)

    ! calculate new rho
    !Db = maxval(sqrt(sum((rho * (B2 - Bold))**2, dim = 2)))
    Db = sum(abs(B2 - Bold)) / sum(abs(Bold))
    s2 = sqrt(sum((rho * (B2 - Bold))**2))
    r2 = sqrt(sum((x - B2)**2))
    rhonew = rho
    if (r2.gt.(mu*s2)) then
        rhonew = rho * tau_inc
    end if
    if (s2.gt.(mu*r2)) then
        rhonew = rho / tau_dec
    end if

    ! stopping criterion
    ! ADMM criterion
    !eps1 = qp * tol_abs + tol_rel * max(sqrt(sum(x**2)), sqrt(sum(B2**2)))
    !eps2 = qp * tol_abs + tol_rel * sqrt(sum(u**2))
    criterion = 0
    !if ((r2.ge.eps1).or.(s2.ge.eps2)) then
    !    criterion = 1
    !end if
    if(Db.gt.thr) then
        criterion = 1
    end if

    if((trace.eq.2).and.(criterion.eq.0)) call admm_trace_1(i, r2, s2, rho)

    update = 0
    if(rhonew.ne.rho) then
        update = 1
        rho = rhonew
    end if
end do
!call DGEMV('T', q, p, 1.d0, B2, q, xm, 1, 0.d0, ym2, 1)
B(0,:) = B2(1,:)
B(1:(q-1),:) = B2(2:q,:)
nit(1) = 1
nit(2) = i
do i = 1, p
    df(i) = 0
    do j = 1, (q-1)
        if(B(j, i).ne.0.d0) df(i) = df(i) + 1
    end do
end do
df(p + 1) = sum(df(1:p))
if(i.eq.maxit) then
    conv = 1
    subrout = 1
    return
end if
end subroutine admm_multilasso

subroutine softmatrixc2(S, Tau, p, q)
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
end subroutine softmatrixc2
