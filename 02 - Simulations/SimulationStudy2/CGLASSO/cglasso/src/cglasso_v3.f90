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
!           '2' glassosub
!           '3' cglasso_v1
subroutine cglasso_v3(n,p,Y,Id,nP,InfoP,lo,up,pendiag,wTht,ym,yv,nrho,rhoratio,rho,maxit_em,thr_em,maxit_bcd,thr_bcd,&
Yipt,B,mu,R,S,Sgm,Tht,Adj_yy,dfB,dfTht,ncomp,Ck,pk,nit,conv,subrout,trace)
implicit none
integer :: n,p,Id(n,p),nP,InfoP(nP,0:(p+3)),pendiag,nrho,maxit_em,maxit_bcd,Adj_yy(p,p,1,nrho),dfB(p+1,1,nrho),dfTht(1,nrho),&
ncomp(1,nrho),Ck(p,1,nrho),pk(p,1,nrho),nit(2,1,nrho),conv,subrout,trace
double precision :: Y(n,p),lo(p),up(p),wTht(p,p),ym(p),yv(p),rhoratio,rho(nrho),thr_em,thr_bcd,Yipt(n,p,1,nrho),&
B(1,p,1,nrho),mu(n,p,1,nrho),R(n,p,1,nrho),S(p,p,1,nrho),Sgm(p,p,1,nrho),Tht(p,p,1,nrho)
!!!!!!!!!!!!!!!!!!!!!!
! internal variables !
!!!!!!!!!!!!!!!!!!!!!!
integer :: i,j,h,k,ii,nnit,ncomp_n,Ck_n(p),pk_n(p),count,info
double precision :: T1o(p),T2o(p,p),Yipt_n(n,p),Yipt_lo(p),Yipt_up(p),mintp,maxtp,dtp,S_n(p,p),R_n(n,p),B_o(p),B_n(p),&
T1(p),T2(p,p),mu_n(n,p),Sgm_o(p,p),Sgm_n(p,p),Tht_n(p,p),dB,dSgm,z,tmean,Tht_o(p,p)
double precision, external :: rdnorm, rpnorm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! computing starting values !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
nnit = 0
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! computing Yipt_n, Yipt_lo and Yipt_up !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
do j = 1, p
    R_n(:, j) = Yipt_n(:, j) - ym(j)
end do
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
B_n = ym
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
do k = 1, nrho
    if(trace.eq.2) call trace_cglasso_v1_2_1(k, rho(k))
    nnit = 0
    do ii = 1, maxit_em
        call rchkusr()
    
        B_o = B_n
        Sgm_o = Sgm_n
        Tht_o = Tht_n
        if(trace.eq.2) call trace_cglasso_v1_2_2(ii)
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
        if(trace.eq.2) call trace_cglasso_v1_2_3()
        B_n = T1 / n
        do j = 1, p
            mu_n(:, j) = B_n(j)
            S_n(:, j) = T2(:, j) / n - B_n * B_n(j)
        end do
        R_n = Yipt_n - mu_n
        Sgm_n = S_n
        !do j = 1, p
        !    do i = j + 1, p
        !        if(abs(Tht_n(i, j)).gt.0.d0) then
        !            Sgm_n(i, j) = Sgm_n(i, j) + rho_gl(i, j) * sign(1.d0, Tht_n(i, j))
        !            Sgm_n(j, i) = Sgm_n(i, j)
        !        end if
        !    end do
        !end do
        if(trace.eq.2) call trace_cglasso_v1_2_4()
        call glassosub_v2(p,S_n,pendiag,rho(k),wTht,maxit_bcd,thr_bcd,Tht_n,ncomp_n,Ck_n,pk_n,nnit,conv,trace)
        if(conv.ne.0) then
            subrout = 2
            return
        end if
        if(trace.eq.2) call trace_cglasso_v1_2_5()
        nit(1, 1, k) = ii
        nit(2, 1, k) = nit(2, 1, k) + nnit
        call inv(p, Tht_n, Sgm_n, info)
        !Sgm_n = (Sgm_n + transpose(Sgm_n)) / 2.0d0
        !Tht_n = (Tht_n + transpose(Tht_n)) / 2.0d0
!        dB = maxval(abs(B_o - B_n))
! code_conv
        !dB = sqrt(sum((B_o - B_n)**2.0d0) / p)
        dB = sqrt(sum((B_o - B_n)**2.0d0)) / p
        dSgm = 0.d0
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
        !count = count + p
        if(trace.eq.2) call trace_cglasso_v1_2_6(thr_em, dB, dSgm, max(dB, dSgm))
!        if(max(dB, dSgm).le.thr_em) exit
! code_conv        
        !if((0.5d0 * (dB + dSgm)).le.thr_em) exit
        if((max(dB, dSgm)).le.thr_em) exit
    end do
    if(ii.ge.maxit_em) then
        conv = 1
        subrout = 3
        return
    end if
    if(trace.eq.1) call  trace_cglasso_v1_1(k, rho(k), nit(1, 1, k), nit(2, 1, k))
    Yipt(:, :, 1, k) = Yipt_n
    B(1, :, 1, k) = B_n
    mu(:, :, 1, k) = mu_n
    R(:, :, 1, k) = R_n
    S(:, :, 1, k) = S_n
    Sgm(:, :, 1, k) = Sgm_n
    Tht(:, :, 1, k) = Tht_n
    ncomp(1, k) = ncomp_n
    Ck(:, 1, k) = Ck_n
    pk(:, 1, k) = pk_n
    dfB(p + 1, 1, nrho) = 0
    do j = 1, p
        do i = j + 1, p
            if(abs(Tht_n(i, j)).gt.0.d0) then
                dfTht(1, k) = dfTht(1, k) + 1
                Adj_yy(i, j, 1, k) = 1
                Adj_yy(j, i, 1, k) = 1
            end if
        end do
    end do
end do
end subroutine cglasso_v3
