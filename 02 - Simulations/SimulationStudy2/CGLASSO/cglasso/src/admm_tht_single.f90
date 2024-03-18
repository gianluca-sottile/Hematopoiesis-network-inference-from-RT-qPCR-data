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
subroutine admm_tht_single(p,S,Theta,wTht,lam,maxit,thr,rho,trace,conv,nit,df)
implicit none
double precision :: S(p,p),Theta(p,p),wTht(p,p),lam,thr,rho
integer :: p,maxit,trace,conv,nit,df
!internal variables
double precision :: Z(p,p),U(p,p),Z2(p,p),s2,r2,eps1,eps2!,Db
integer :: info,criterion,iter,i,j

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
    
    ! ridge equation (1)
    ! gather eigen values (spectral decomposition)
    call ridgec(S - rho * Z + U , rho, p, Theta)
    
    ! penalty equation (2)
    ! soft-thresholding
    Z = Theta + U / rho
    call softmatrix_tht(Z, lam * wTht / rho, p)
    
    ! update U (3)
    U = U + rho * (Theta - Z)
    
    ! calculate new rho
    !Db = sum(abs(Z2 - Z)) / sum(abs(Z))
    s2 = sqrt(sum((rho*(Z2 - Z))**2))
    r2 = sqrt(sum((Theta - Z)**2))
    
    ! stopping criterion
    ! ADMM criterion
    eps1 = p * thr + thr * max(sqrt(sum(Theta**2)), sqrt(sum(Z**2)))
    eps2 = p * thr + thr * sqrt(sum(U**2))
    criterion = 0
    if ((r2.ge.eps1).or.(s2.ge.eps2)) then
        criterion = 1
    end if
    !if (Db.gt.thr) then
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
do i = 1, p - 1
    do j = i + 1, p
        if(Theta(i, j).ne.0.d0) df = df + 1
    end do
end do
if(iter.eq.maxit) then
    conv = 1
    return
end if
!call inv(p, Z, Sgm, info)
end subroutine admm_tht_single

subroutine softmatrix_tht(S, Tau, p)
implicit none
double precision :: S(p,p), Tau(p,p)
integer :: p
!internal variables
integer :: i, j
! loop over all elements
do i = 1, p
    do j = i, p
        ! soft threshold each element
        call softc(S(i, j), Tau(i, j))
        S(j, i) = S(i, j)
    end do
end do
end subroutine softmatrix_tht

subroutine softc(s, tau)
implicit none
double precision :: s, tau
!internal variables
double precision :: d
! soft-thresholding
d = 0.d0
if((s.gt.0).and.(s.gt.tau)) then
    s = s - tau
else if((s.lt.0).and.(-s.gt.tau)) then
    s = s + tau
else
    s = d
end if
end subroutine softc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Gianluca Sottile
! e-mail: gianluca.sottile@unipa.it
! webpage: bit.ly/gianlucasottile
!
! version: 1.0.0
! Data: Feb 09, 2021
!
! Description
! 'ridgec' gather eigen values of S (spectral decomposition) and compute omega hat for lambda (zero gradient equation)
!
! Input
! S = estimated covariance matrix (p x p)
! lam = tuning parameter (scalar)
!
! Output
! Omega = p x p matrix
subroutine ridgec(S, lam, p, Omega)
implicit none
double precision :: S(p,p), lam, Omega(p,p)
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
w = sqrt((-w + sqrt((w**2) + 4 * lam)) / (2 * lam))
do i = 1, p
    A(:, i) = A(:, i) * w(i)
end do
!Omega = matmul(Omega, transpose(A))
call DGEMM('N', 'T', p, p, p, 1.d0, A, p, A, p, 0.d0, Omega, p)
end subroutine ridgec
