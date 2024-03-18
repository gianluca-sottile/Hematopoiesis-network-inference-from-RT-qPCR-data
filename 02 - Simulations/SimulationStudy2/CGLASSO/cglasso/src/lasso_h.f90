!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Luigi Augugliaro              (R code)
! e-mail: luigi.augugliaro@unipa.it
!
! Author: Gianluca Sottile              (Fortran code)
! e-mail: gianluca.sottile@unipa.it
!
! version: 1.0.0
! Data: December 5, 2018
!
! Description
! The subroutine 'lasso_h' solves a modified lasso problem for 
! a fixed value of the tuning parameter, i.e.
! 
!       b^hat = argmin_b 1 / (2 * n) ||y + vh - Xb||^2_2 + lambda ||b||_1
!
!                               INPUT
!
! q = number of predictors
! ym = mean of the responce variable
! xm = q-dimensional vector of means of the predictors
! xtx = q x q dimensional matrix X'X / n
! xtr = q x q dimensional matrix X'R / n
! v = q-dimensional vector of vixed values
! w = q-dimensional vector of poritive weights
! lmb = value of the tuning parameter
! maxit = maximum number of iteration steps
! thr = convergence threshould
!
!                               OUTPUT
!
! bn = lasso estimate of the regression coefficients. In input bn is used 
!       as starting value
! conv = integer used to encode the convergence of the algorithm:
!       '0' convergence is met
!       '1' maximum number of iterations has been exceeded
! nit = number of iteration steps
subroutine lasso_h(q,ym,xm,xtx,xtr,v,w,lmb,bn,maxit,thr,conv,nit)
implicit none
integer :: q,maxit,conv,nit
double precision :: ym,xm(q),xtx(q,q),xtr(q),v(q),w(q),lmb,bn(0:q),thr
integer :: j
double precision :: bo(0:q),lmbj,bj,dbj,dbmax
bo = bn
conv = 0
nit = 0
do
    call rchkusr()

    nit = nit + 1
    if(nit.gt.maxit) then 
        conv = 1
        exit
    end if
    dbmax = 0.d0
    bn(0) = ym
    do j = 1, q
        lmbj = lmb / xtx(j,j)
        bj = bo(j) + (xtr(j) + v(j)) / xtx(j,j)
        if(abs(bj).lt.(lmbj * w(j))) then
            bn(j) = 0.d0
        else
            bn(j) = bj - lmbj * w(j) * sign(1.d0, bj)
        end if
        dbj = bn(j) - bo(j)
        dbmax = max(dbmax, abs(dbj))
! code_conv     
!        dbmax = dbmax + (dbj)**2
        xtr = xtr - xtx(j,:) * dbj
        bn(0) = bn(0) - xm(j) * bn(j)
    end do
    dbj = bn(0) - bo(0)
    dbmax = max(dbmax, abs(dbj))
! code_conv
!    dbmax = dbmax + (dbj)**2
!    dbmax = sqrt(dbmax  / (q + 1))
    xtr = xtr - xm * dbj
!    if(sum((bn - bo)**2).lt.thr) then
!    if(sum(abs(bn - bo)).lt.thr) then
    if(dbmax.lt.thr) then
        exit
    else
        bo = bn
    end if
end do
end subroutine lasso_h
