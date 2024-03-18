!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Luigi Augugliaro
! e-mail: luigi.augugliaro@unipa.it
!
! version: 1.0.0
! Data: 2020-06-23
!
! conv      [OUT] :: error code
! subrout   [OUT] :: integer encoding the subroutines
!           '1' e_step_v2
!           '2' multilasso
!           '3' glassosub
!           '4' cglasso_v2
subroutine cglasso_v4(n,q,X,p,Y,Id,nP,InfoP,lo,up,wB,pendiag,wTht,ym,yv,nlambda,lambdaratio,lambda,nrho,rhoratio,&
rho,maxit_em,thr_em,maxit_bcd,thr_bcd,Yipt,B,mu,R,S,Sgm,Tht,Adj_yy,dfB,dfTht,ncomp,Ck,pk,Adj_xy,nit,conv,subrout,trace)
implicit none
integer :: n,q,p,Id(n,p),nP,InfoP(nP,0:(p+3)),pendiag,nlambda,nrho,maxit_em,maxit_bcd,Adj_yy(p,p,nlambda,nrho),&
dfB(p+1,nlambda,nrho),dfTht(nlambda,nrho),ncomp(nlambda,nrho),Ck(p,nlambda,nrho),pk(p,nlambda,nrho),&
Adj_xy(q,p,nlambda,nrho),nit(2,nlambda,nrho),conv,subrout,trace
double precision :: X(n,q),Y(n,p),lo(p),up(p),wB(q,p),wTht(p,p),ym(p),yv(p),lambdaratio,lambda(nlambda),rhoratio,&
rho(nrho),thr_em,thr_bcd,Yipt(n,p,nlambda,nrho),B(0:q,p,nlambda,nrho),mu(n,p,nlambda,nrho),R(n,p,nlambda,nrho),&
S(p,p,nlambda,nrho),Sgm(p,p,nlambda,nrho),Tht(p,p,nlambda,nrho)
!!!!!!!!!!!!!!!!!!!!!!
! internal variables !
!!!!!!!!!!!!!!!!!!!!!!
integer :: i,j,h,k,ii,nnit(2),ncomp_n,Ck_n(p),pk_n(p),count,info
double precision :: X2(n,q),xtx_n(q,q),xtr_n(q,p),xm(q),mintp,maxtp,dtp,S_n(p,p),Yipt_n(n,p),R_n(n,p),&
T1o(p),T2o(p,p),T1(p),T2(p,p),B_o(0:q,p),B_n(0:q,p),mu_n(n,p),Sgm_o(p,p),Sgm_n(p,p),Tht_n(p,p),dB,dSgm,&
Yipt_lo(p),Yipt_up(p),z,tmean,Tht_o(p,p),tau!,O_n(p,p),Y_n(p,p)
double precision, external :: rdnorm, rpnorm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! computing starting values !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
nnit = 0
X2 = X
do h = 1, q
    xm(h) = sum(X(:, h)) / n
    X2(:, h) = X2(:, h) - xm(h)
end do
call DGEMM('T', 'N', q, q, n, 1.d0, X2, n, X2, n, 0.d0, xtx_n, q)
!call DGEMM('N', 'N', q, q, n, 1.d0, transpose(X), q, X, n, 0.d0, xtx_n, q)
xtx_n = xtx_n / n
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! computing the observed statistics T1o and T2o !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
T1o = 0.d0
do j = 1, p
    T1o(j) = sum(Y(:, j), mask = Id(:, j).eq.0)
end do
T2o = 0.d0
do i = 1, p
    do j = i, p
        if(any((Id(:, i).eq.0).and.(Id(:, j).eq.0))) then
            T2o(i, j) = sum(Y(:, i) * Y(:, j), mask = (Id(:, i).eq.0).and.(Id(:, j).eq.0))
            T2o(j, i) = T2o(i, j)
        end if
    end do
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! computing Yipt_n,Yipt_lo and Yipt_up !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Yipt_lo = ym - 3.d0 * sqrt(yv)
Yipt_up = ym + 3.d0 * sqrt(yv)
Yipt_n = Y
do j = 1, p
    where (Id(:, j).eq.9) Yipt_n(:, j) = ym(j)
    if (any(Id(:, j).eq.-1)) then
        z = (lo(j) - ym(j)) / sqrt(yv(j))
        tmean = ym(j) - sqrt(yv(j)) * rdnorm(z, 0.0d0, 1.0d0, 0) / rpnorm(z, 0.0d0, 1.0d0, 1, 0)
        if(tmean.le.Yipt_lo(j)) tmean = Yipt_lo(j)
        where (Id(:, j).eq.-1) Yipt_n(:, j) = tmean
    end if
    if (any(Id(:, j).eq.1)) then
        z = (up(j) - ym(j)) / sqrt(yv(j))
        tmean = ym(j) + sqrt(yv(j)) * rdnorm(z, 0.0d0, 1.0d0, 0) / rpnorm(z, 0.0d0, 1.0d0, 0, 0)
        if(tmean.ge.Yipt_up(j)) tmean = Yipt_up(j)
        where (Id(:, j).eq.1) Yipt_n(:, j) = tmean
    end if
    ym(j) = sum(Yipt_n(:, j)) / n
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!
! computing R_n and xtr_n !
!!!!!!!!!!!!!!!!!!!!!!!!!!!
do k = 1, p
    R_n(:, k) = Yipt_n(:, k) - ym(k)
end do
!call DGEMM('N', 'N', q, p, n, 1.d0, transpose(X), q, R_n, n, 0.d0, xtr_n, q)
call DGEMM('T', 'N', q, p, n, 1.d0, X, n, R_n, n, 0.d0, xtr_n, q)
xtr_n = xtr_n / n
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! computing S_n and rho-values !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
maxtp = 0.d0
do k = 1, p
    S_n(k, k) = yv(k)
    do h = k + 1, p
        S_n(h, k) = dot_product(Yipt_n(:, h), Yipt_n(:, k)) / n - ym(h) * ym(k)
        S_n(k, h) = S_n(h, k)
        if(wTht(h, k).gt.0.0d0) then
            if(wTht(h, k).eq.1.0d0) then
                maxtp = max(maxtp, abs(S_n(h, k)))
            else
                maxtp = max(maxtp, abs(S_n(h, k) / wTht(h, k)))
            end if
        end if
    end do
end do
if(rho(1).eq.0.d0) then
    maxtp = maxtp * 1.001d0
    mintp = rhoratio * maxtp
    rho(1) = maxtp
    dtp = (rho(1) - mintp) / (nrho - 1)
    do k = 2, nrho
        rho(k) = rho(k - 1) - dtp
    end do
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Inizializing B_n, mu_n, Sgm_n and Tht_n !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
B_n = 0.d0
B_n(0, :) = ym
mu_n = 0.d0
Sgm_n = 0.d0
Tht_n = 0.d0
do k = 1, p
    mu_n(:, k) = ym(k)
    Sgm_n(k, k) = S_n(k, k)
    if(pendiag.eq.1) then
        if(wTht(k, k).eq.1.0d0) then
            Sgm_n(k, k) = Sgm_n(k, k) + maxtp
        else
            Sgm_n(k, k) = Sgm_n(k, k) + maxtp * wTht(k, k)
        end if
    end if
    Tht_n(k, k) = 1.d0 / Sgm_n(k, k)
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!
! computing lambda-values !
!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(lambda(1).eq.0.d0) then
    maxtp = 0.d0
    do k = 1, p
        do h = 1, q
            if(wB(h, k).gt.0.0d0) then
                if(wB(h, k).eq.1.0d0) then
                    maxtp = max(maxtp, abs(xtr_n(h, k)))
                else
                    maxtp = max(maxtp, abs(xtr_n(h, k)) / wB(h, k))
                end if
            end if
        end do
    end do
    maxtp = maxtp * 1.001d0
    mintp = lambdaratio * maxtp
    lambda(1) = maxtp
    dtp = (lambda(1) - mintp) / (nlambda - 1)
    do k = 2, nlambda
        lambda(k) = lambda(k - 1) - dtp
    end do
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Starting optimization                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do h = 1, nlambda
    do k = 1, nrho
        if(trace.eq.2) call trace_cglasso_v2_2_1(k + nrho * (h - 1), rho(k), lambda(h))
        nnit = 0
        !!!!!!!!!!!!!!!!!!!!!!!!!
        ! starting EM algorithm !
        !!!!!!!!!!!!!!!!!!!!!!!!!
        do ii = 1, maxit_em
            call rchkusr()
        
            B_o = B_n
            Sgm_o = Sgm_n
            Tht_o = Tht_n
            if(trace.eq.2) call trace_cglasso_v2_2_2(ii)
            !!!!!!!!!!!!!!!!!!!
            ! Starting E-step !
            !!!!!!!!!!!!!!!!!!!
            call e_step_v2(n,p,Y,lo,up,nP,InfoP,T1o,T2o,mu_n,Sgm_n,Tht_n,Yipt_lo,Yipt_up,Yipt_n,T1,T2,conv)
            if(conv.ne.0) then
                subrout = 1
                return
            end if
            !!!!!!!!!!!!!!!!!!!
            ! Starting M-step !
            !!!!!!!!!!!!!!!!!!!
            if(trace.eq.2) call trace_cglasso_v2_2_3()
            ym = T1 / n
            do j = 1, p
            !R_n = Yipt_n - mu_n
                R_n(:,j) = Yipt_n(:,j) - ym(j)
            end do
            !call DGEMM('N', 'N', q, p, n, 1.d0, transpose(X) / n, q, R_n, n, 0.d0, xtr_n, q)
            call DGEMM('T', 'N', q, p, n, 1.d0, X2, n, R_n, n, 0.d0, xtr_n, q)
            xtr_n = xtr_n / n
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! fitting multilasso model !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if(trace.eq.2) call trace_cglasso_v2_2_4(lambda(h))
            tau = 1.d0
            call admm_b_single(p,q,p*q,xtx_n,xtr_n,B_n(1:q,:),Tht_n,wB,lambda(h),maxit_bcd,thr_bcd,tau,&
                trace,conv,subrout,nnit,dfB(:,h,k))
            if(conv.ne.0) then
                subrout = 2
                return
            end if
            B_n(0,:) = ym
            call DGEMV('T', q, p, -1.d0, B_n(1:q,:), q, xm, 1, 1.d0, B_n(0,:), 1)
            do j = 1, p
                mu_n(:, j) = B_n(0, j)
                do i = 1, q
                    if(abs(B_n(i, j)).gt.0.d0) mu_n(:, j) = mu_n(:, j) + X(:, i) * B_n(i, j)
                end do
            end do
            R_n = Yipt_n - mu_n
            !!!!!!!!!!!!!!!!
            ! updating S_n !
            !!!!!!!!!!!!!!!!
            !call DGEMM('N', 'N', p, p, n, 1.d0, transpose(Yipt_n), p, mu_n, n, 0.d0, S_n, p)
            call DGEMM('T', 'N', p, p, n, 1.d0, Yipt_n, n, mu_n, n, 0.d0, S_n, p)
            S_n = T2 - S_n - transpose(S_n)
            do j = 1, p
                S_n(j, j) = S_n(j, j) + dot_product(mu_n(:, j), mu_n(:, j))
                do i = j + 1, p
                    S_n(i, j) = S_n(i, j) + dot_product(mu_n(:, i), mu_n(:, j))
                    S_n(j, i) = S_n(i, j)
                end do
            end do
            S_n = S_n / n
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! solving the glasso problem !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if(trace.eq.2) call trace_cglasso_v2_2_5(rho(k))
            !call glassosub(p,S_n,pendiag,rho_gl,maxit_bcd,thr_bcd,Sgm_n,Tht_n,ncomp_n,Ck_n,pk_n,nnit(1),conv,trace)
            call glassosub_v2(p,S_n,pendiag,rho(k),wTht,maxit_bcd,thr_bcd,Tht_n,ncomp_n,Ck_n,pk_n,nnit(1),conv,trace)
            if(conv.ne.0) then
                subrout = 3
                return
            end if
            if(trace.eq.2) call trace_cglasso_v2_2_6()
            nit(1, h, k) = ii
            nit(2, h, k) = nit(2, h, k) + sum(nnit)
            call inv(p, Tht_n, Sgm_n, info)
            
            dB = sqrt(sum((B_o - B_n)**2.0d0)) / (dfB(p + 1, h, k) + p)
            dSgm = 0.d0
            count = 0
            do j = 1, p
                do i = 1, j
                    dSgm = dSgm + (Tht_o(i, j) - Tht_n(i, j))**2.0d0
                    if(abs(Tht_n(i, j)).gt.0.d0) count = count + 1
                end do
            end do
            dSgm = sqrt(dSgm) / count
            
            if(trace.eq.2) call trace_cglasso_v2_2_7(thr_em, dB, dSgm, max(dB, dSgm))
            if((max(dB, dSgm)).le.thr_em) exit
        end do
        if(ii.ge.maxit_em) then
            conv = 1
            subrout = 4
            return
        end if
        if(trace.eq.1) call trace_cglasso_v2_1(k + nrho * (h - 1), rho(k), lambda(h), nit(1, h, k), nit(2, h, k))
        Yipt(:, :, h, k) = Yipt_n
        B(:, :, h, k) = B_n
        mu(:, :, h, k) = mu_n
        R(:, :, h, k) = R_n
        S(:, :, h, k) = S_n
        Sgm(:, :, h, k) = Sgm_n
        Tht(:, :, h, k) = Tht_n
        ncomp(h, k) = ncomp_n
        Ck(:, h, k) = Ck_n
        pk(:, h, k) = pk_n
        do j = 1, p
            do i = j + 1, p
                if (abs(Tht_n(i, j)).gt.0.d0) then
                    dfTht(h, k) = dfTht(h, k) + 1
                    Adj_yy(i, j, h, k) = 1
                    Adj_yy(j, i, h, k) = 1
                end if
            end do
            do i = 1, q
                if (abs(B_n(i, j)).gt.0.d0) Adj_xy(i, j, h, k) = 1
            end do
        end do
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! setting starting values for the next multilasso and glasso problem !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    B_n = B(:, :, h, 1)
    Tht_n = Tht(:, :, h, 1)
    mu_n = mu(:, :, h, 1)
    ncomp_n = ncomp(h, 1)
    Ck_n = Ck(:, h, 1)
    pk_n = pk(:, h, 1)
end do
end subroutine cglasso_v4
