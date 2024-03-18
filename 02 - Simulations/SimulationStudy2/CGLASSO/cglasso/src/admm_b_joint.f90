subroutine admm_b_joint(p,q,qp,K,fk,xtx_n,xtr_n,B,Tht,W,lmb,alpha,maxit,thr,rho,trace,conv,&
    subrout,nit,df)
implicit none
integer :: p,q,qp,K,maxit,trace,conv,subrout,nit(2),df(p+1,K)
double precision :: fk(K),xtx_n(q,q,K),xtr_n(q,p,K),B(q,p,K),Tht(p,p,K),W(q,p,K),lmb,alpha,thr,rho
!internal variables
double precision :: Gamma(q,p,K),U(q,p,K),xtr_n_theta(q,p,K),theta_xtx_n(qp,qp,K),Iqp(qp,qp)
double precision :: theta_xtx_n_inv(qp,qp,K),Gamma2(q,p,K),normsoftB(q,p),tempY(qp),s2,r2,eps1,eps2!,Db
integer :: i,k2,criterion,info,j

Gamma = 0.d0
U = 0.d0
info = 0

Iqp = 0.d0
do i = 1, qp
    Iqp(i,i) = 1.d0
end do

theta_xtx_n = 0.d0
do k2 = 1, K
    call DGEMM('N', 'N', q, p, p, 1.d0, xtr_n(:,:,k2), q, Tht(:,:,k2), p, 0.d0, xtr_n_theta(:,:,k2), q)
    xtr_n_theta(:,:,k2) = fk(k2) * xtr_n_theta(:,:,k2)
    call kron(p, q, xtx_n(:,:,k2), Tht(:,:,k2), theta_xtx_n(:,:,k2))
    theta_xtx_n(:,:,k2) = fk(k2) * theta_xtx_n(:,:,k2)
    theta_xtx_n_inv(:,:,k2) = theta_xtx_n(:,:,k2) + rho * Iqp
    call dpotrf('U', qp, theta_xtx_n_inv(:,:,k2), qp, info)
    ! or
    !call inv(qp, theta_xtx_n + rho * Iqp, theta_xtx_n_inv, info)
end do

criterion = 1
i = 0
do while ((criterion.eq.1).and.(i.lt.maxit))
    call rchkusr()

    ! update values
    i = i + 1
    Gamma2 = Gamma

    normsoftB = 0.d0
    do k2 = 1, K
        ! update 'Beta'
        call mattovec(p, q, xtr_n_theta(:,:,k2) + rho * (Gamma(:,:,k2) - U(:,:,k2)), tempY)
        call dpotrs('U', qp, 1, theta_xtx_n_inv(:,:,k2), qp, tempY, qp, info)
        ! or one of the following
        !call DGEMV('N', qp, qp, 1.d0, theta_xtx_n_inv, qp, tempY, 1, 0.d0, tempX, 1)
        !call solve(theta_xtx_n + rho * Iqp, tempY, qp, info)
    
        call vectomat(p, q, tempY, B(:,:,k2))

        ! update 'Gamma'
        Gamma(:,:,k2) = B(:,:,k2) + U(:,:,k2)
        call softmatrix_b(Gamma(:,:,k2), alpha * lmb * W(:,:,k2) / rho, p, q)
        normsoftB = normsoftB + Gamma(:,:,k2)**2
    end do
    normsoftB = sqrt(normsoftB)

    !Db = 0.d0
    s2 = 0.d0
    r2 = 0.d0
    eps1 = 0.d0
    eps2 = 0.d0
    do k2 = 1, K
        call softmatrix_b_joint(Gamma(:,:,k2), (1 - alpha) * lmb * W(:,:,k2) / rho, normsoftB, p, q)
        ! update 'U'
        U(:,:,k2) = U(:,:,k2) + B(:,:,k2) - Gamma(:,:,k2)
        ! stopping criterion
        ! Db = Db + fk(k2) * sum(abs(Gamma2(:,:,k2) - Gamma(:,:,k2))) / sum(abs(Gamma(:,:,k2)))
        s2 = s2 + fk(k2) * sqrt(sum((rho*(Gamma2(:,:,k2) - Gamma(:,:,k2)))**2))
        r2 = r2 + fk(k2) * sqrt(sum((B(:,:,k2) - Gamma(:,:,k2))**2))
        eps1 = eps1 + fk(k2) * (thr * max(sqrt(sum(B(:,:,k2)**2)), sqrt(sum(Gamma(:,:,k2)**2))))
        eps2 = eps2 + fk(k2) * (thr * sqrt(sum(U(:,:,k2)**2)))
    end do

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
do k2 = 1, K
    do i = 1, p
        df(i, k2) = 0
        do j = 1, q
            if(B(j, i, k2).ne.0.d0) df(i, k2) = df(i, k2) + 1
        end do
    end do
    df(p + 1, k2) = sum(df(1:p, k2))
end do
if(i.eq.maxit) then
    conv = 1
    subrout = 1
    return
end if
end subroutine admm_b_joint

subroutine softmatrix_b_joint(S, Tau, norm, p, q)
implicit none
double precision :: S(q,p), Tau(q,p), norm(q,p)
integer :: p, q
!internal variables
integer :: i, j
! loop over all elements
do i = 1, q
    do j = 1, p
        ! soft threshold each element
        call soft_group(S(i, j), Tau(i, j), norm(i, j))
    end do
end do
end subroutine softmatrix_b_joint

subroutine soft_group(s, tau, norm)
implicit none
double precision :: s, tau, norm
!internal variables
double precision :: d
! soft-thresholding
d = 0.d0
if(norm.ne.0.d0) then
    if((1.d0 - tau/norm).gt.0) then
        s = (1.d0 - tau/norm) * s
    else
        s = d
    end if
end if
end subroutine soft_group
