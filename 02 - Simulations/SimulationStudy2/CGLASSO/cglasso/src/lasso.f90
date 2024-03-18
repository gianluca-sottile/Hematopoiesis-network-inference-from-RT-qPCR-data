!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Luigi Augugliaro
! e-mail: luigi.augugliaro@unipa.it
! webpage: http://dssm.unipa.it/augugliaro
!
! version: 1.0.0
! Data: July 12, 2018
!
! Description
! 'lasso' solves a lasso problem using a ccd algorithm
!
! Input
! pendiag = integer used to specify if the diagonal elements of the                 (integer)
!           concentration matrix (Tht) are penalized
!           'pendiag = 0' means that the diagonal elements are not penalized
!           'pendiag = 1' means that the diagonal elements are penalized
! p = dimension of the matrix Sgm11                                                 (integer)
!
! Output
! conv = integer used to encode the convergence of the algorithm:                   (integer)
!           '0' convergence is met
!           '1' maximum number of steps is reached
subroutine lasso(p, Sgm11, Sgm12, rhom, niter, thr, b, lnit, conv)
implicit none
integer :: p, niter, lnit, conv
double precision :: Sgm11(p,p), Sgm12(p), rhom(p), thr, b(p), temp(p)
!internal variables
double precision, parameter :: fac = 0.2
integer :: i,m,card
double precision :: bm,dbm,y,dbm_max
card = 0
do m = 1, p
    if(abs(b(m)).gt.0.d0) card = card + 1
end do
if(card.le.int(fac * p)) then
    call DGEMV('N', p, p, 1.d0, Sgm11, p, b, 1, 0.d0, temp, 1)
    Sgm12 = Sgm12 - temp !matmul(Sgm11, b)
else
    do m = 1, p
        if(abs(b(m)).gt.0.d0) Sgm12 = Sgm12 - Sgm11(:, m) * b(m)
    end do
end if
do i = 1, niter
    call rchkusr()

    lnit = i
    dbm_max = 0.d0
    do m = 1, p
        bm = b(m)
        b(m) = 0.d0
        y = Sgm12(m) + Sgm11(m,m) * bm
        if(abs(y).gt.rhom(m)) b(m) = sign(abs(y) - rhom(m), y) / Sgm11(m,m)
        dbm = b(m) - bm
        dbm_max = max(dbm_max, abs(dbm))
! code_conv    
!        dbm_max = dbm_max + (dbm)**2
        Sgm12 = Sgm12 - dbm * Sgm11(:,m)
    end do
! code_conv    
!    dbm_max = sqrt(dbm_max / p)
    if(dbm_max.lt.thr) exit
end do
if(i.eq.niter) conv = 1
end subroutine lasso
