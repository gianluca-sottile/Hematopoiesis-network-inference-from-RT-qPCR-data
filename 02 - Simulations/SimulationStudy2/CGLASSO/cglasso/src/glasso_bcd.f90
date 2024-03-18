!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Luigi Augugliaro
! e-mail: luigi.augugliaro@unipa.it
! webpage: http://dssm.unipa.it/augugliaro
!
! version: 1.0.0
! Data: July 12, 2018
!
! Description
! 'glasso_bcd' implements the block coordinate discent (bcd) algorithm proposed in Friedman et al. (2008)
! to fit a l1-penalized guassian graphical model.
!
! Input
! p = dimension of the matrix S                                                     (integer)
! S = p x p dimensional variance/covariance matrix                                  (double precision)
! rho = p x p dimenaion matrix of the values of the tuning parameter                (double precision)
! maxit = maximum number of steps of the bcd algorithm                              (integer)
! thr = threshold value used to declare the convergence of the bcd algorithm        (double precision)
!
! Output
! Sgm = estiamted l1-penalized variance/covariance matrix                           (double precision)
! Tht = estiamted l1-penalized concentration matrix                                 (double precision)
! nit = number of bcd steps                                                         (integer)
! conv = integer used to encode the convergence of the algorithm:                   (integer)
!           '0' convergence is met
!           '1' maximum number of steps is reached
subroutine glasso_bcd(p,S,rho,maxit,thr,Sgm,Tht,nit,conv,trace)
implicit none
integer :: p,maxit,nit,conv,trace
double precision :: S(p,p),rho(p,p),thr,Sgm(p,p),Tht(p,p)
!internal variables
integer :: m,i,j,p_1,idx(p-1),lnit
double precision :: shr,b(p-1),bmat(p-1,p),Sgm11(p-1,p-1),Sgm12(p-1),Sgm12o(p-1),rhom(p-1),dSgm12
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! shr is the threshold used to declare the convergence of the bcd algorithm, i.e.
! shr = thr * sum_{i < j} |S_{ij}| / (p * (p - 1) / 2)
p_1 = p - 1
shr = 0.d0
do i = 1, p
    do j = i + 1, p
        shr = shr + abs(S(i, j))
    end do
end do
shr = thr * 2.d0 * shr / (p * p_1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! inizializing bmat matrix     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
bmat = 0.d0
do m = 1, p
    j = 0
    do i = 1, p
        if(i.ne.m) then
            j = j + 1
            idx(j) = i
        end if
    end do
    bmat(:, m) = -Tht(idx, m) / Tht(m, m)
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! starting bcd algorithm       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
nit = 0
if(trace.eq.2) call glasso_trace_2_3_1()
do
    call rchkusr()

    nit = nit + 1
    dSgm12 = 0.d0
    do m = 1, p
        j = 0
        do i = 1, p
            if(i.ne.m) then
                j = j + 1
                idx(j) = i
            end if
        end do
        b = bmat(:, m)
        Sgm12o = Sgm(idx, m)
        Sgm12 = S(idx, m)
        Sgm11 = Sgm(idx, idx)
        rhom = rho(idx, m)
        lnit = 0
        call lasso(p_1, Sgm11, Sgm12, rhom, maxit, thr / p_1, b, lnit, conv)
        if(conv.eq.1) return
        bmat(:, m) = b
        Sgm(idx, m) = S(idx, m) - Sgm12
        Sgm(m, idx) = Sgm(idx, m)
        dSgm12 = max(dSgm12, sum(abs(Sgm(idx, m) - Sgm12o)) / p_1)
! code_conv     
!        dSgm12 = dSgm12 + sum((Sgm(idx, m) - Sgm12o)**2)
    end do
! code_conv     
!    dSgm12 = sqrt(dSgm12) / (p * (p - 1))
    if(trace.eq.2) call glasso_trace_2_3_2(nit, lnit, dSgm12)
    if(dSgm12.lt.shr) then
        if(trace.eq.2) call glasso_trace_2_3_3(shr)
        exit
    end if
    if(nit.eq.maxit) then
        conv = 1
        return
    end if
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Convert the bmat matrix to Tht        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do m = 1, p
    j = 0
    do i = 1, p
        if(i.ne.m) then
            j = j + 1
            idx(j) = i
        end if
    end do
    Sgm12 = Sgm(idx, m)
    b = -bmat(:, m)
    Tht(m, m) = 1.d0 / (Sgm(m, m) + dot_product(Sgm12, b))
    Tht(idx, m) = Tht(m, m) * b
end do
end subroutine glasso_bcd




























