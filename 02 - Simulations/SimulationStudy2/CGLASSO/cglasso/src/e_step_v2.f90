subroutine e_step_v2(n,p,Y,lo,up,nP,InfoP,T1o,T2o,mu,Sgm,Tht,Yipt_lo,Yipt_up,Yipt,T1,T2,conv)
! conv = integer used to encode errors
!       '-1' error in memory allocation
!        '2' error in computing tmean and tvar
!        '3' Tht_mm inversion failed
implicit none
integer :: n, p,nP,InfoP(nP,0:(p+3)),conv
double precision :: Y(n,p),lo(p),up(p),T1o(p),T2o(p,p),mu(n,p),Sgm(p,p),Tht(p,p),Yipt_lo(p),Yipt_up(p),Yipt(n,p),T1(p),T2(p,p)
! internal variables
double precision, external :: rdnorm, rpnorm
integer :: i,j,ii,nmar,nlc,nrc,nmis,no,idP(nP + 1)
double precision :: z,ratio
double precision, dimension (:), allocatable :: mu_m_given_o,tmean,tvar
double precision, dimension (:,:), allocatable :: B,Sgm_m_given_o

T1 = T1o
T2 = T2o
idP(1:nP) = InfoP(:, 0)
idP(nP + 1) = n + 1
if (idP(1).eq. n + 1) return

do ii = 1, nP
    call rchkusr()

    nmar = InfoP(ii, p + 1)
    nlc = InfoP(ii, p + 2)
    nrc = InfoP(ii, p + 3)
    nmis = nmar + nlc + nrc
    no = p - nmis
    allocate(B(nmis,max(no,1)),Sgm_m_given_o(nmis,nmis),mu_m_given_o(nmis),tmean(nmis),tvar(nmis),stat=conv)
    if(conv.ne.0) then
        conv = -1
        return
    end if
    if (no.ne.0) then
        call inv(nmis, Tht(InfoP(ii, 1:nmis), InfoP(ii, 1:nmis)), Sgm_m_given_o, conv)
        if(conv.ne.0) then
            conv = 3
            return
        end if
!        B = matmul(Sgm_m_given_o, Tht(InfoP(ii, 1:nmis), InfoP(ii, (nmis + 1):p)))
        call DGEMM('N', 'N', nmis, no, nmis, 1.d0, Sgm_m_given_o, nmis, Tht(InfoP(ii, 1:nmis), InfoP(ii, (nmis + 1):p)), &
                nmis, 0.d0, B, nmis)
    else
        Sgm_m_given_o = Sgm(InfoP(ii, 1:nmis), InfoP(ii, 1:nmis))
    end if

    do i = idP(ii), idP(ii + 1) - 1
    
        mu_m_given_o = mu(i, InfoP(ii, 1:nmis))
        if (no.ne.0) then
            do j = 1, no
                mu_m_given_o = mu_m_given_o - B(:, j) * (Y(i, InfoP(ii, nmis + j)) - mu(i, InfoP(ii, nmis + j)))
            end do
        end if

        do j = 1, nmar
            tmean(j) = mu_m_given_o(j)
            if (tmean(j).le.Yipt_lo(InfoP(ii, j))) tmean(j) = Yipt_lo(InfoP(ii, j))
            if (tmean(j).ge.Yipt_up(InfoP(ii, j))) tmean(j) = Yipt_up(InfoP(ii, j))
            tvar(j) = Sgm_m_given_o(j, j)
        end do
        do j = nmar + 1, nmar + nlc
            z = (lo(InfoP(ii, j)) - mu_m_given_o(j)) / sqrt(Sgm_m_given_o(j, j))
            ratio = rdnorm(z, 0.0d0, 1.0d0, 0) / rpnorm(z, 0.0d0, 1.0d0, 1, 0)
            tmean(j) = mu_m_given_o(j) - sqrt(Sgm_m_given_o(j, j)) * ratio
            if(tmean(j).ge.Yipt_lo(InfoP(ii, j))) then
                tvar(j) = Sgm_m_given_o(j, j) * (1.0d0 - z * ratio - ratio**2)
            else
                tmean(j) = Yipt_lo(InfoP(ii, j))
                tvar(j) = Sgm_m_given_o(j, j)
            end if
        end do
        do j = nmar + nlc + 1, nmis
            z = (up(InfoP(ii, j)) - mu_m_given_o(j)) / sqrt(Sgm_m_given_o(j, j))
            ratio = rdnorm(z, 0.0d0, 1.0d0, 0) / rpnorm(z, 0.0d0, 1.0d0, 0, 0)
            tmean(j) = mu_m_given_o(j) + sqrt(Sgm_m_given_o(j, j)) * ratio
            if (tmean(j).le.Yipt_up(InfoP(ii, j))) then
                tvar(j) = Sgm_m_given_o(j, j) * (1.0d0 + z * ratio - ratio**2)
            else
                tmean(j) = Yipt_up(InfoP(ii, j))
                tvar(j) = Sgm_m_given_o(j, j)
            end if
        end do

        Yipt(i, InfoP(ii, 1:nmis)) = tmean
        T1(InfoP(ii, 1:nmis)) = T1(InfoP(ii, 1:nmis)) + tmean
        do j = 1, nmis
            T2(InfoP(ii, 1:nmis), InfoP(ii, j)) = T2(InfoP(ii, 1:nmis), InfoP(ii, j)) + tmean * tmean(j)
            T2(InfoP(ii, j), InfoP(ii, j)) = T2(InfoP(ii, j), InfoP(ii, j)) + tvar(j)
            if(no.ne.0) T2(InfoP(ii, (nmis + 1):p), InfoP(ii, j)) = T2(InfoP(ii, (nmis + 1):p), InfoP(ii, j)) + &
                        Y(i, InfoP(ii, (nmis + 1):p)) * tmean(j)
        end do
        if(no.ne.0) T2(InfoP(ii, 1:nmis), InfoP(ii, (nmis + 1):p)) = transpose(T2(InfoP(ii, (nmis + 1):p), InfoP(ii, 1:nmis)))

    end do
    deallocate(B, Sgm_m_given_o, mu_m_given_o, tmean, tvar, stat = conv)
    if(conv.ne.0) then
        conv = -1
        return
    end if
end do
end subroutine e_step_v2
