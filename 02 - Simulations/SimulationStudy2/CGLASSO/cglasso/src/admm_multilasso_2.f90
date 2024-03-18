subroutine admm_multilasso_2(p,q,K,qp,ym,xm,xtx_n,xtr_n,B,Tht,W,lmb1,lmb2,maxit,thr,rho,conv,subrout,nit,df,trace)
implicit none
double precision :: xm(q,K),ym(p,K),xtx_n(q,q,K),xtr_n(q,p,K)
double precision :: B(q+1,p,K),Tht(p,p,K),W(q,p,K),lmb1,lmb2,thr,rho
integer :: p,q,K,qp,maxit,conv,subrout,nit(2),df(p+1,k),trace
!internal variables
double precision :: Bold(q+1,p,K),x(q,p,K),u(q,p,K)
double precision :: xtr_n_theta(q,p,K),theta_xtx_n(qp,qp,K)
double precision :: Iqp(qp,qp),theta_xtx_n_inv(qp,qp,K),tempY(qp)
double precision :: W2(q,p,K),Db,normsoftB(q,p)
integer :: i,j,k2,criterion,info
x = B(2:(q+1),:,:)
u = 0.d0
theta_xtx_n = 0.d0
Iqp = 0.d0
do i = 1, qp
    Iqp(i,i) = 1.d0
end do

info = 0
do k2 = 1, K
    call DGEMM('N','N',q,p,p,1.d0,xtr_n(:,:,k2),q,Tht(:,:,k2),p,0.d0,xtr_n_theta(:,:,k2),q)
    call kron(p,q,xtx_n(:,:,k2),Tht(:,:,k2),theta_xtx_n(:,:,k2))
    do i = 1, p
        W2(:,i,k2) = W(:,i,k2) * Tht(i,i,k2)
    end do
    theta_xtx_n_inv(:,:,k2) = theta_xtx_n(:,:,k2) + rho * Iqp
    call dpotrf('U', qp, theta_xtx_n_inv(:,:,k2), qp, info)
    !call inv(qp, theta_xtx_n(:,:,k2) + rho * Iqp, theta_xtx_n_inv(:,:,k2), info)
end do

criterion = 1
i = 0
do while ((criterion.eq.1).and.(i.lt.maxit))
    call rchkusr()

    ! update values
    i = i + 1
    Bold = B
    normsoftB = 0.d0
    Db = 0.d0

    do k2 = 1, K
        ! update 'x'
        call mattovec(p, q, xtr_n_theta(:,:,k2) + rho * (B(2:(q+1),:,k2) - u(:,:,k2)), tempY)
        call dpotrs('U', qp, 1, theta_xtx_n_inv(:,:,k2), qp, tempY, qp, info)
        !call DGEMV('N', qp, qp, 1.d0, theta_xtx_n_inv(:,:,k2), qp, tempY, 1, 0.d0, tempY, 1)
        call vectomat(p, q, tempY, x(:,:,k2))

        ! update 'B'
        B(2:(q+1),:,k2) = u(:,:,k2) + x(:,:,k2)
        call softmatrixc2(B(2:(q+1),:,k2), lmb1 / rho * W2(:,:,k2), p, q)

        ! update 'b0'
        call DGEMV('T', q, p, 1.d0, B(2:(q+1),:,k2), q, xm(:,k2), 1, 0.d0, B(1,:,k2), 1)
        B(1,:,k2) = ym(:,k2) - B(1,:,k2)

        ! update 'normsoftB'
        normsoftB = normsoftB + B(2:(q+1),:,k2)*B(2:(q+1),:,k2)
    end do

    normsoftB = sqrt(normsoftB)

    do k2 = 1, K
        call softmatrixc3(B(2:(q+1),:,k2), normsoftB, lmb2 / rho, p, q)
        
        ! update 'u'
        u(:,:,k2) = u(:,:,k2) + x(:,:,k2) - B(2:(q+1),:,k2)

        Db = Db + sum(abs(B(:,:,k2) - Bold(:,:,k2))) / sum(abs(Bold(:,:,k2)))
    end do
    
    criterion = 0
    if(Db.gt.thr) then
        criterion = 1
    end if

    if((trace.eq.1).and.(criterion.eq.0)) call admm_trace_1(i, Db, thr, rho)
    if((trace.eq.2)) call admm_trace_1(i, Db, thr, rho)
end do
!call DGEMV('T', q, p, 1.d0, B2, q, xm, 1, 0.d0, ym2, 1)
nit(1) = 1
nit(2) = i
do k2 = 1, K
    do i = 1, p
        df(i,k2) = 0
        do j = 2, q
            if(B(j,i,k2).ne.0.d0) df(i,k2) = df(i,k2) + 1
        end do
    end do
    df(p+1,k2) = sum(df(1:p,k2))
end do

if(i.eq.maxit) then
    conv = 1
    subrout = 1
    return
end if
end subroutine admm_multilasso_2

subroutine softmatrixc3(A, normA, lmb2, p, q)
implicit none
double precision :: A(q,p),normA(q,p),lmb2
integer :: p,q
!internal variables
double precision :: temp
integer :: i,j
! loop over all elements
do i = 1, q
    do j = 1, p
        if(normA(i,j).eq.0) then
            A(i,j) = 0
        else
            temp = 1.d0 - lmb2 / normA(i,j)
            if(temp.le.0) then
                A(i,j) = 0
            else
                A(i,j) = A(i,j) * temp
            end if
        end if
    end do
end do
end subroutine softmatrixc3
