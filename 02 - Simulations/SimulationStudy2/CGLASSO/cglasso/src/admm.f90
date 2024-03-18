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
subroutine admm(p, S, initOmega, initZ, initY, lam, thr, maxit, Sgm, iter, trace)
implicit none
double precision :: S(p,p), initOmega(p,p), initZ(p,p), initY(p,p), lam(p,p)
double precision :: thr, Sgm(p,p)
integer :: p, maxit, iter, trace
!internal variables
double precision :: s2, r2, Z2(p,p), Z(p,p), Y(p,p), alpha, rho, tau_inc, tau_dec, mu
double precision :: Omega(p,p), Tau(p,p), Taum(p,p), Db
integer :: criterion, info!, i, j, count

!thr = 0.d0
!do i = 1, p
!    do j = i + 1, p
!        thr = thr + abs(S(i, j))
!    end do
!end do
!thr = tol_abs * thr / (p**2)

alpha = 1.0d0
mu = 10.0d0
tau_inc = 2.0d0
tau_dec = 2.0d0
rho = 1.0d0

criterion = 1
iter = 0
s2 = 0.d0
r2 = 0.d0
Z2 = initZ
Z = initZ
Y = initY
Omega = initOmega
info = 0

! save values
Tau = lam * alpha
Taum = lam - Tau

if(trace.eq.2) call admm_trace_1_2()
! loop until convergence
do while ((criterion.eq.1).and.(iter.lt.maxit))
    call rchkusr()
    
    ! update values
    iter = iter + 1
    Z = Z2
    
    ! ridge equation (1)
    ! gather eigen values (spectral decomposition)
    call ridgec(S + Y - rho * Z, rho, p, Omega)
    
    ! penalty equation (2)
    ! soft-thresholding
    Z2 = Y + rho * Omega
    call softmatrixc(Z2, Tau, p)
    Z2 = Z2 / (Taum + rho)
    
    ! update Y (3)
    Y = Y + rho * (Omega - Z2)
    
    ! calculate new rho
    !Db = maxval(sqrt(sum((rho*(Z2 - Z))**2, dim = 2)))
    Db = sum(abs(Z2 - Z)) / sum(abs(Z))
    s2 = sqrt(sum((rho*(Z2 - Z))**2))
    r2 = sqrt(sum((Omega - Z2)**2))
    if (r2.gt.(mu*s2)) then
        rho = rho * tau_inc
    end if
    if (s2.gt.(mu*r2)) then
        rho = rho / tau_dec
    end if
    
    ! stopping criterion
    ! ADMM criterion
    !eps1 = p * tol_abs + tol_rel * max(sqrt(sum(Omega**2)), sqrt(sum(Z2**2)))
    !eps2 = p * tol_abs + tol_rel * sqrt(sum(Y**2))
    criterion = 0
    !if ((r2.ge.eps1).or.(s2.ge.eps2)) then
    !    criterion = 1
    !end if
    if (Db.gt.thr) then
        criterion = 1
    end if

    if(trace.eq.9) call admm_trace_1(iter, Db, rho)
    if((trace.eq.2).and.(criterion.eq.0)) call admm_trace_1(iter, Db, rho)
    !if(trace.eq.2) call admm_trace_1(iter, Db, rho)
end do

initOmega = Omega
initZ = Z2
initY = Y
call inv(p, Z2, Sgm, info)
end subroutine admm

subroutine softmatrixc(S, Tau, p)
implicit none
double precision :: S(p,p), Tau(p,p)
integer :: p
!internal variables
integer :: i, j
! loop over all elements
do i = 1, p
    do j = 1, p
        ! soft threshold each element
        call softc(S(i, j), Tau(i, j))
    end do
end do
end subroutine softmatrixc

