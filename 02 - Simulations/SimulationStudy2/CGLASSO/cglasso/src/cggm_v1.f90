!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Luigi Augugliaro
! e-mail: luigi.augugliaro@unipa.it
!
! version: 1.0.0
! Data: 2020-07-10
!
! conv      [OUT] :: error code
! subrout   [OUT] :: integer encoding the subroutines
!           '1' e_step_v2
!           '2' glassosub
!           '3' cglasso_v1
subroutine cggm_v1(n,p,Y,Id,nP,InfoP,lo,up,ym,yv,pendiag,wTht,ntp,rho,maxit_em,thr_em,maxit_bcd,thr_bcd,Yipt,B,mu,R,S,&
Sgm,Tht,nit,conv,subrout,trace)
implicit none
integer :: n,p,Id(n,p),nP,InfoP(nP,0:(p+3)),pendiag,ntp,maxit_em,maxit_bcd,nit(2,ntp),conv,subrout,trace
double precision :: Y(n,p),lo(p),up(p),wTht(p,p),ym(p),yv(p),rho(ntp),thr_em,thr_bcd,Yipt(n,p,1,1),B(1,p,1,1),&
mu(n,p,1,1),R(n,p,1,1),S(p,p,1,1),Sgm(p,p,1,1),Tht(p,p,1,1)
!!!!!!!!!!!!!!!!!!!!!!
! internal variables !
!!!!!!!!!!!!!!!!!!!!!!
integer :: i,j,h,ii,nnit,ncomp_n,Ck_n(p),pk_n(p),count
double precision :: T1o(p),T2o(p,p),T1(p),T2(p,p),Yipt_n(n,p),Yipt_lo(p),Yipt_up(p),S_n(p,p),R_n(n,p),mu_n(n,p),B_o(p),&
B_n(p),Sgm_o(p,p),Sgm_n(p,p),Tht_n(p,p),rho_gl(p,p),dB,dSgm
double precision, parameter :: big = huge(0.0d0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initializing starting values !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
nnit = 0
Yipt_lo = ym - 3.d0 * sqrt(yv)
Yipt_up = ym + 3.d0 * sqrt(yv)
B_o = B(1, :, 1, 1)
B_n = B_o
Sgm_o = Sgm(:, :, 1, 1)
Sgm_n = Sgm_o
Yipt_n = Yipt(:, :, 1, 1)
S_n = S(:, :,  1, 1)
R_n = R(:, :, 1, 1)
mu_n = mu(:, :, 1, 1)
Tht_n = Tht(:, :, 1, 1)
!where (wTht.eq.0.0d0) wTht = big
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
    if(trace.eq.2) call trace_cglasso_v1_2_1(h, rho(h))
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
        if(trace.eq.2) call trace_cglasso_v1_2_2(ii)
        !!!!!!!!!!!!!!!!!!!
        ! Starting E-step !
        !!!!!!!!!!!!!!!!!!!
        call e_step_v2(n,p,Y,lo,up,nP,InfoP,T1o,T2o,mu_n,Sgm_n,Tht_n,Yipt_lo,Yipt_up,Yipt_n,T1,T2,conv)
        if(conv.ne.0) then
            subrout = 1
!            where (wTht.eq.big) wTht = 0.0d0
            return
        end if
        !!!!!!!!!!!!!!!!!!!
        ! Starting M-step !
        !!!!!!!!!!!!!!!!!!!
        if(trace.eq.2) call trace_cglasso_v1_2_3()
        B_n = T1 / n
        do j = 1, p
            mu_n(:, j) = B_n(j)
            S_n(:, j) = T2(:, j) / n - B_n * B_n(j)
        end do
        R_n = Yipt_n - mu_n
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
        if(trace.eq.2) call trace_cglasso_v1_2_4()
        call glassosub(p,S_n,pendiag,rho_gl,maxit_bcd,thr_bcd,Sgm_n,Tht_n,ncomp_n,Ck_n,pk_n,nnit,conv,trace)
        if(conv.ne.0) then
            subrout = 3
!            where (wTht.eq.big) wTht = 0.0d0
            return
        end if
        if(trace.eq.2) call trace_cglasso_v1_2_5()
        nit(1, h) = ii
        nit(2, h) = nit(2, h) + nnit
        Sgm_n = (Sgm_n + transpose(Sgm_n)) / 2.0d0
        Tht_n = (Tht_n + transpose(Tht_n)) / 2.0d0
!        dB = maxval(abs(B_o - B_n))
! code_conv
        !dB = sqrt(sum((B_o - B_n)**2.0d0) / p)
        dB = sqrt(sum((B_o - B_n)**2.0d0)) / p
        dSgm = 0.0d0
        count = 0
        do j = 1, p
            do i = 1, j
!                dSgm = max(dSgm, abs(Sgm_o(i, j) - Sgm_n(i, j)))
! code_conv
                dSgm = dSgm + (Sgm_o(i, j) - Sgm_n(i, j))**2.0d0
                if(abs(Sgm_n(j, i)).gt.0.d0) count = count + 1
            end do
        end do
! code_conv
        !dSgm = sqrt(dSgm / (p * (p + 1) * 0.5d0))
        dSgm = sqrt(dSgm) / count
        !count = count + p
        if(trace.eq.2) call trace_cglasso_v1_2_6(thr_em, dB, dSgm, max(dB, dSgm))
!        if(max(dB, dSgm).le.thr_em) exit
! code_conv        
        !if((0.5d0 * (dB + dSgm)).le.thr_em) exit
        if(max(dB, dSgm).le.thr_em) exit
    end do
    if(ii.ge.maxit_em) then
        conv = 1
        subrout = 4
!        where (wTht.eq.big) wTht = 0.0d0
        return
    end if
    if(trace.eq.1) call trace_cglasso_v1_1(h, rho(h), nit(1, h), nit(2, h))
end do
Yipt(:, :, 1, 1) = Yipt_n
B(1, :, 1, 1) = B_n
mu(:, :, 1, 1) = mu_n
R(:, :, 1, 1) = R_n
S(:, :, 1, 1) = S_n
Sgm(:, :, 1, 1) = Sgm_n
Tht(:, :, 1, 1) = Tht_n
!where (wTht.eq.big) wTht = 0.0d0
end subroutine cggm_v1
