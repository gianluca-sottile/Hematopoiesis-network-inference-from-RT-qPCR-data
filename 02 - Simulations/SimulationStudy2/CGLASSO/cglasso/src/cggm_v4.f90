!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Luigi Augugliaro
! e-mail: luigi.augugliaro@unipa.it
!
! version: 1.0.0
! Data: 2020-07-09
!
! conv      [OUT] :: error code
! subrout   [OUT] :: integer encoding the subroutines
!           '1' e_step_v2
!           '2' multilasso
!           '3' glassosub
!           '4' cglasso_v2
subroutine cggm_v4(n,q,X,p,Y,Id,nP,InfoP,lo,up,ym,yv,wB,pendiag,wTht,ntp,lambda,rho,maxit_em,thr_em,maxit_bcd,thr_bcd,&
Yipt,B,mu,R,S,Sgm,Tht,nit,conv,subrout,trace)
implicit none
integer :: n,q,p,Id(n,p),nP,InfoP(nP,0:(p+3)),pendiag,ntp,maxit_em,maxit_bcd,nit(2,ntp),conv,subrout,trace
double precision :: X(n,q),Y(n,p),lo(p),up(p),wB(q,p),wTht(p,p),ym(p),yv(p),rho(ntp),lambda(ntp),thr_em,thr_bcd,&
Yipt(n,p,1,1),B(0:q,p,1,1),mu(n,p,1,1),R(n,p,1,1),S(p,p,1,1),Sgm(p,p,1,1),Tht(p,p,1,1)
!!!!!!!!!!!!!!!!!!!!!!
! internal variables !
!!!!!!!!!!!!!!!!!!!!!!
integer :: i,j,h,ii,nnit(2),dfB(p+1),ncomp_n,Ck_n(p),pk_n(p),count
double precision :: xtx_n(q,q),xtr_n(q,p),xm(q),T1o(p),T2o(p,p),T1(p),T2(p,p),Yipt_n(n,p),Yipt_lo(p),&
Yipt_up(p),S_n(p,p),R_n(n,p),mu_n(n,p),B_o(0:q,p),B_n(0:q,p),Sgm_o(p,p),Sgm_n(p,p),Tht_n(p,p),rho_gl(p,p),dB,dSgm,&
Tht_o(p,p)
double precision, parameter :: big = huge(0.0d0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initializing starting values !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Yipt_lo = ym - 3.d0 * sqrt(yv)
Yipt_up = ym + 3.d0 * sqrt(yv)
B_o = B(:, :, 1, 1)
B_n = B_o
Sgm_o = Sgm(:, :, 1, 1)
Sgm_n = Sgm_o
Yipt_n = Yipt(:, :, 1, 1)
S_n = S(:, :,  1, 1)
R_n = R(:, :, 1, 1)
mu_n = mu(:, :, 1, 1)
Tht_n = Tht(:, :, 1, 1)
xtx_n = 0.0d0
!call DGEMM('N', 'N', q, q, n, 1.d0, transpose(X), q, X, n, 0.0d0, xtx_n, q)
call DGEMM('T', 'N', q, q, n, 1.d0, X, n, X, n, 0.0d0, xtx_n, q)
xtx_n = xtx_n / n
!call DGEMM('N', 'N', q, p, n, 1.d0, transpose(X) / n, q, R_n, n, 0.0d0, xtr_n, q)
call DGEMM('T', 'N', q, p, n, 1.d0, X, n, R_n, n, 0.0d0, xtr_n, q)
xtr_n = xtr_n / n
where (wB.eq.0.0d0) wB = big
!where (wTht.eq.0.0d0) wTht = big
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! computing starting values !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
nnit = 0
do h = 1, q
    xm(h) = sum(X(:, h)) / n
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! computing the observed statistics T1o and T2o !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
T1o = 0.0d0
do j = 1, p
    T1o(j) = sum(Y(:, j), mask = Id(:, j).eq.0)
end do
T2o = 0.0d0
do i = 1, p
    do j = i, p
        if(any((Id(:, i).eq.0).and.(Id(:, j).eq.0))) then
            T2o(i, j) = sum(Y(:, i) * Y(:, j), mask = (Id(:, i).eq.0).and.(Id(:, j).eq.0))
            T2o(j, i) = T2o(i, j)
        end if
    end do
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Starting optimization                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do h = 2, ntp
    if(trace.eq.2) call trace_cglasso_v2_2_1(h, rho(h), lambda(h))
    nnit = 0
!    rho_gl = rho(h) * wTht
    where (wTht.ne.0.0d0) rho_gl = rho(h)
    where (wTht.eq.0.0d0) rho_gl = big

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
            where (wB.eq.big) wB = 0.0d0
!            where (wTht.eq.big) wTht = 0.0d0
            return
        end if
        !!!!!!!!!!!!!!!!!!!
        ! Starting M-step !
        !!!!!!!!!!!!!!!!!!!
        if(trace.eq.2) call trace_cglasso_v2_2_3()
        ym = T1 / n
        R_n = Yipt_n - mu_n
        !call DGEMM('N', 'N', q, p, n, 1.d0, transpose(X) / n, q, R_n, n, 0.0d0, xtr_n, q)
        call DGEMM('T', 'N', q, p, n, 1.d0, X, n, R_n, n, 0.0d0, xtr_n, q)
        xtr_n = xtr_n / n
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! fitting multilasso model !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(trace.eq.2) call trace_cglasso_v2_2_4(lambda(h))
        call multilasso(p,q,ym,xm,xtx_n,xtr_n,B_n,Tht_n,wB,lambda(h),maxit_bcd,thr_bcd,conv,subrout,nnit,dfB,trace)
        if(conv.ne.0) then
            subrout = 2
            where (wB.eq.big) wB = 0.0d0
!            where (wTht.eq.big) wTht = 0.0d0
            return
        end if
        do j = 1, p
            mu_n(:, j) = B_n(0, j)
            do i = 1, q
                if(abs(B_n(i, j)).gt.0.0d0) mu_n(:, j) = mu_n(:, j) + X(:, i) * B_n(i, j)
            end do
        end do
        R_n = Yipt_n - mu_n
        !!!!!!!!!!!!!!!!
        ! updating S_n !
        !!!!!!!!!!!!!!!!
        !call DGEMM('N', 'N', p, p, n, 1.d0, transpose(Yipt_n), p, mu_n, n, 0.0d0, S_n, p)
        call DGEMM('T', 'N', p, p, n, 1.d0, Yipt_n, n, mu_n, n, 0.0d0, S_n, p)
        S_n = T2 - S_n - transpose(S_n)
        do j = 1, p
            S_n(j, j) = S_n(j, j) + dot_product(mu_n(:, j), mu_n(:, j))
            do i = j + 1, p
                S_n(i, j) = S_n(i, j) + dot_product(mu_n(:, i), mu_n(:, j))
                S_n(j, i) = S_n(i, j)
            end do
        end do
        S_n = S_n / n
        !!!!!!!!!!!!!!!!!!!!!!
        ! inizializing Sgm_n !
        !!!!!!!!!!!!!!!!!!!!!!
        Sgm_n = S_n
        do j = 1, p
            do i = j + 1, p
                if(abs(Tht_n(i, j)).gt.0.0d0) then
                    Sgm_n(i, j) = Sgm_n(i, j) + rho_gl(i, j) * sign(1.d0, Tht_n(i, j))
                    Sgm_n(j, i) = Sgm_n(i, j)
                end if
            end do
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! solving the glasso problem !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(trace.eq.2) call trace_cglasso_v2_2_5(rho(h))
        call glassosub_v2(p,S_n,pendiag,rho_gl,maxit_bcd,thr_bcd,Sgm_n,Tht_n,ncomp_n,Ck_n,pk_n,nnit(1),conv,trace)
        if(conv.ne.0) then
            subrout = 3
            where (wB.eq.big) wB = 0.0d0
!            where (wTht.eq.big) wTht = 0.0d0
            return
        end if
        if(trace.eq.2) call trace_cglasso_v2_2_6()
        nit(1, h) = ii
        nit(2, h) = nit(2, h) + sum(nnit)
        !Sgm_n = (Sgm_n + transpose(Sgm_n)) / 2.0d0
        !Tht_n = (Tht_n + transpose(Tht_n)) / 2.0d0
!        dB = maxval(abs(B_o - B_n))
! code_conv
        !dB = sqrt(sum((B_o - B_n)**2.0d0) / (p * (q + 1)))
        dB = sqrt(sum((B_o - B_n)**2.0d0)) / (dfB(p + 1) + p)
        dSgm = 0.0d0
        count = 0
        do j = 1, p
            do i = 1, j
!                dSgm = max(dSgm, abs(Sgm_o(i, j) - Sgm_n(i, j)))
! code_conv
                !dSgm = dSgm + (Sgm_o(i, j) - Sgm_n(i, j))**2.0d0
                !if(abs(Sgm_n(j, i)).gt.0.d0) count = count + 1
                dSgm = dSgm + (Tht_o(i, j) - Tht_n(i, j))**2.0d0
                if(abs(Tht_n(i, j)).gt.0.d0) count = count + 1
            end do
        end do
! code_conv
        !dSgm = sqrt(dSgm / (p * (p + 1) * 0.5d0))
        dSgm = sqrt(dSgm) / count
        !count = count + dfB(p + 1) + p
        if(trace.eq.2) call trace_cglasso_v2_2_7(thr_em, dB, dSgm, max(dB, dSgm))
!        if(max(dB, dSgm).le.thr_em) exit
! code_conv        
        !if((0.5d0 * (dB + dSgm)).le.thr_em) exit
        if(max(dB, dSgm).le.thr_em) exit
    end do
    if(ii.ge.maxit_em) then
        conv = 1
        subrout = 4
        where (wB.eq.big) wB = 0.0d0
!        where (wTht.eq.big) wTht = 0.0d0
        return
    end if
    if(trace.eq.1) call  trace_cglasso_v2_1(h, rho(h), lambda(h), nit(1, h), nit(2, h))
end do
Yipt(:, :, 1, 1) = Yipt_n
B(:, :, 1, 1) = B_n
mu(:, :, 1, 1) = mu_n
R(:, :, 1, 1) = R_n
S(:, :, 1, 1) = S_n
Sgm(:, :, 1, 1) = Sgm_n
Tht(:, :, 1, 1) = Tht_n
where (wB.eq.big) wB = 0.0d0
!where (wTht.eq.big) wTht = 0.0d0
end subroutine cggm_v4
