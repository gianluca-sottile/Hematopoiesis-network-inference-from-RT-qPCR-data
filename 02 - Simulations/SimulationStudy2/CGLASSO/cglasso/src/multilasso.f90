!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Luigi Augugliaro              (R code)
! e-mail: luigi.augugliaro@unipa.it
!
! Author: Gianluca Sottile              (Fortran code)
! e-mail: gianluca.sottile@unipa.it
!
! version: 1.0.0
! Data: December 12, 2018
!
! Description
! The subroutine 'multilasso_h' solves a multivariate lasso problem for 
! a fixed value of the tuning parameter, i.e.
! 
!       B^hat = argmin_B 1 / (2 * n) tr{(Y - XB)Tht(Y - XB)'} + lambda ||B||_1
!
! Arguments:
!
!                                                  INPUT
!
! p = number of responce variables
! q = number of predictor
! ym = p dimensional vector of the p responce variables
! xm = q dimensional vector of the q predictor
! XtX = q x q dimensional matrix X'X / n
! Tht = p x p precision matrix
! W = q x p dimensional matrix of weights
! lmb = value of the tuning parameter
! maxit = maximum number of iterations
! thr = convergence threshold
! trace = integer for printing out information on video
!
!                                                  OUPUT
!
! xtr_n = q x p dimensional matrix X'R / n
! B = (q + 1) x dimensional matrix of the regression coefficients
! conv = integer used to encode the convergence of the algorithm
!           '0' convergence is met
!           '1' maximum number of iterations has been exceeded
! subrout = integer used to encode the subroutine
! nit = 2 dimensional vector used to store the number of steps until
!           convergence is met:
!               nit(1) = number of ciclies
!               nit(2) = total number of steps
! df = (p + 1) dimensional vector of non-zero regression coefficients
!
subroutine multilasso(p,q,ym,xm,xtx_n,xtr_n,B,Tht,W,lmb,maxit,thr,conv,subrout,nit,df,trace)
implicit none
integer :: p,q,maxit,conv,subrout,nit(2),df(p+1),trace
double precision :: ym(p),xm(q),xtx_n(q,q),xtr_n(q,p),B(0:q,p),Tht(p,p),W(q,p),lmb,thr
! internal variables
integer :: i,j,h,k,i2
double precision :: Db,bo(0:q),bn(0:q),lmbh,vh(q),Db_o
if(trace.eq.2) call trace_cglasso_v2_2_4_1()
Db = 0.d0
bo = 0.d0
bn = 0.d0
lmbh = 0.d0
conv = 0
nit = 0
Db_o = 99999999999.0d0
do
    call rchkusr()

    i = 0
    nit(1) = nit(1) + 1
    Db = 0.d0
    do j = 1, p
        bo = B(:, j)
        bn = bo
!        lmbh = lmb / Tht(j,j)
        lmbh = lmb
        vh = 0.d0
        do k = 1, q
            do h = 1, p
                if(j.ne.h) then
                    if(Tht(h,j).ne.0.d0) vh(k) = vh(k) + xtr_n(k,h) * Tht(h,j)
                end if
            end do
        end do
        vh = vh / Tht(j,j)
        call lasso_h(q,ym(j),xm,xtx_n,xtr_n(:,j),vh,W(:,j),lmbh,bn,maxit,thr/q,conv,i2)
        if(conv.eq.1) then
            subrout = 2
            exit
        end if
        i = i + i2
        if(i.ge.maxit) then
            conv = 1
            exit
        end if
        B(:, j) = bn
        Db = max(Db, sum(abs(bn - bo)) / (q + 1))
! code_conv     
!        Db = Db + sum((bn - bo)**2)
    end do
! code_conv     
!    Db = sqrt(Db / (p * (q + 1)))
    if(conv.eq.1) then
        subrout = 1
        exit
    end if
    if(trace.eq.2) call trace_cglasso_v2_2_4_2(nit(1), i, Db)
    if((Db.le.thr).or.(Db_o.lt.Db)) then
        if(trace.eq.2) call trace_cglasso_v2_2_4_3(thr)
        exit
    end if
    Db_o = Db
end do
nit(2) = i
do h = 1, p
    df(h) = 0
    do k = 1, q
        if(B(k, h).ne.0.d0) df(h) = df(h) + 1
    end do
end do
df(p + 1) = sum(df(1:p))
end subroutine multilasso
