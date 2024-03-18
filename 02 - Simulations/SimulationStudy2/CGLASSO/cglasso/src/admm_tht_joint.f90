!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Gianluca Sottile
! e-mail: gianluca.sottile@unipa.it
! webpage: bit.ly/gianlucasottile
!
! version: 1.0.0
! Data: Feb 09, 2021
!
! Description
! 'admm' Penalized precision matrix estimation using the ADMM algorithm
subroutine admm_tht_joint(p,K,fk,S,Theta,wTht,lam,alpha,maxit,thr,rho,trace,conv,nit,df)
implicit none
integer :: p,K,maxit,trace,conv,nit,df(K)
double precision :: fk(K),S(p,p,K),Theta(p,p,K),wTht(p,p,K),lam,alpha,thr,rho
!internal variables
double precision :: Z(p,p,K),U(p,p,K),Z2(p,p,K),normsoftTht(p,p),s2,r2,eps1,eps2!,Db
integer :: info,criterion,iter,k2,i,j

Z = 0.d0
U = 0.d0

info = 0

criterion = 1
iter = 0

! loop until convergence
do while ((criterion.eq.1).and.(iter.lt.maxit))
    call rchkusr()
    
    ! update values
    iter = iter + 1
    Z2 = Z
    
    normsoftTht = 0.d0
    do k2 = 1, K
        ! ridge equation (1)
        ! gather eigen values (spectral decomposition)
        call ridgec_joint(fk(k2), S(:,:,k2) - rho * Z(:,:,k2) / fk(k2) + U(:,:,k2) / fk(k2), rho, p, Theta(:,:,k2))
    
        ! penalty equation (2)
        ! soft-thresholding
        Z(:,:,k2) = Theta(:,:,k2) + U(:,:,k2) / rho
        call softmatrix_tht(Z(:,:,k2), alpha * lam * wTht(:,:,k2) / rho, p)
        normsoftTht = normsoftTht + Z(:,:,k2)**2
    end do
    normsoftTht = sqrt(normsoftTht)
    
    !Db = 0.d0
    s2 = 0.d0
    r2 = 0.d0
    eps1 = 0.d0
    eps2 = 0.d0
    do k2 = 1, K
        call softmatrix_tht_joint(Z(:,:,k2), (1 - alpha) * lam * wTht(:,:,k2) / rho, normsoftTht, p)
        ! update U (3)
        U(:,:,k2) = U(:,:,k2) + rho * (Theta(:,:,k2) - Z(:,:,k2))
    
        ! stopping criterion
        ! ADMM criterion
        ! calculate new rho
        !Db = Db + fk(k2) * sum(abs(Z2(:,:,k2) - Z(:,:,k2))) / sum(abs(Z(:,:,k2)))
        s2 = s2 + fk(k2) * sqrt(sum((rho*(Z2(:,:,k2) - Z(:,:,k2)))**2))
        r2 = r2 + fk(k2) * sqrt(sum((Theta(:,:,k2) - Z(:,:,k2))**2))
        eps1 = eps1 + fk(k2) * (p * thr + thr * max(sqrt(sum(Theta(:,:,k2)**2)), sqrt(sum(Z(:,:,k2)**2))))
        eps2 = eps2 + fk(k2) * (p * thr + thr * sqrt(sum(U(:,:,k2)**2)))
    end do

    criterion = 0
    if ((r2.ge.eps1).or.(s2.ge.eps2)) then
        criterion = 1
    end if
    !if(Db.gt.thr) then
    !    criterion = 1
    !end if
    
    !if(trace.eq.9) call admm_trace_1(iter, Db, rho)
    !if((trace.eq.2).and.(criterion.eq.0)) call admm_trace_1(iter, Db, rho)
    if(trace.eq.9) call admm_trace_2(iter, r2, eps1, s2, eps2, rho)
    if((trace.eq.2).and.(criterion.eq.0)) call admm_trace_2(iter, r2, eps1, s2, eps2, rho)

    if (r2.gt.(10.d0 * s2)) then
        rho = rho * 2.d0
    end if
    if (s2.gt.(10.d0 * r2)) then
        rho = rho / 2.d0
    end if
end do

Theta = Z
nit = iter
df = 0
do k2 = 1, K
    do i = 1, p - 1
        do j = i + 1, p
            if(Theta(i, j, k2).ne.0.d0) df(k2) = df(k2) + 1
        end do
    end do
end do
if(iter.eq.maxit) then
    conv = 1
    return
end if
!call inv(p, Z, Sgm, info)
end subroutine admm_tht_joint

subroutine softmatrix_tht_joint(S, Tau, norm, p)
implicit none
double precision :: S(p,p), Tau(p,p), norm(p,p)
integer :: p
!internal variables
integer :: i, j
! loop over all elements
do i = 1, p
    do j = i, p
        ! soft threshold each element
        call soft_group(S(i, j), Tau(i, j), norm(i, j))
        S(j, i) = S(i, j)
    end do
end do
end subroutine softmatrix_tht_joint

subroutine ridgec_joint(fk, S, lam, p, Omega)
implicit none
double precision :: fk, S(p,p), lam, Omega(p,p)
integer :: p
!internal variables
double precision :: w(p), A(p,p)
double precision, dimension(:), allocatable :: work
integer :: lwork, info, i
A = S
lwork = -1
allocate(work(p), stat = info)
call DSYEV('V', 'U', p, A, p, w, work, lwork, info)
lwork = int(work(1))
deallocate(work, stat = info)
allocate(work(lwork), stat = info)
call DSYEV('V', 'U', p, A, p, w, work, lwork, info)
deallocate(work, stat = info)
w = sqrt(fk * (-w + sqrt((w**2) + 4 * lam / fk)) / (2 * lam))
do i = 1, p
    A(:, i) = A(:, i) * w(i)
end do
!Omega = matmul(Omega, transpose(A))
call DGEMM('N', 'T', p, p, p, 1.d0, A, p, A, p, 0.d0, Omega, p)
end subroutine ridgec_joint
