subroutine fitmcgm(n,p,Y,lo,up,R,nstp,eps,ym,yv,conv)
implicit none
integer :: n,p,R(n,p),nstp,conv
double precision :: Y(n,p),lo(p),up(p),eps,ym(p),yv(p)
! internal variables
integer :: i,j,k,no,nrc,nlc,nmar
double precision :: ym_o,ym_n,yv_o,yv_n,T1_o,T2_o,T1_n,T2_n,z,ratio,tmean,tvar,dm,dv
double precision, dimension(:), allocatable :: yobs
double precision, parameter :: tol = 0.001d0
double precision, external :: rdnorm, rpnorm
do j = 1, p
    nrc = count(R(:, j).eq.1)
    nlc = count(R(:, j).eq.-1)
    nmar = count(R(:, j).eq.9)
    no = n - nrc - nlc - nmar
    if(no.eq.n) then
        ym(j) = sum(Y(:, j)) / n
        yv(j) = sum(Y(:, j)**2) / n - ym(j)**2
    else
        allocate(yobs(1:no), stat = conv)
        if(conv.ne.0) then
            conv = -1
            return
        end if
        k = 0
        do i = 1, n
            if(R(i, j).eq.0) then
                k = k + 1
                yobs(k) = Y(i, j)
            end if
        end do
        T1_o = sum(yobs)
        T2_o = sum(yobs**2)
        ym(j) = T1_o / no
        yv(j) = T2_o / no - ym(j)**2
        ym_o = ym(j)
        yv_o = yv(j)
        ym_n = 0.d0
        yv_n = 0.d0
        do i = 1, nstp
            T1_n = T1_o
            T2_n = T2_o
            if(nlc.ne.0) then
                z = (lo(j) - ym_o) / sqrt(yv_o)
                if(rdnorm(z, 0.0d0, 1.0d0, 0).gt.(tol * rpnorm(z, 0.0d0, 1.0d0, 0, 0))) then
                    ratio = rdnorm(z, 0.0d0, 1.0d0, 0) / rpnorm(z, 0.0d0, 1.0d0, 1, 0)
                    tmean = ym_o - sqrt(yv_o) * ratio
                    tvar = yv_o * (1.d0 - z * ratio - ratio**2)
                else
                    tmean = ym_o
                    tvar = yv_o
                end if
                T1_n = T1_n + nlc * tmean
                T2_n = T2_n + nlc * (tvar + tmean**2)
            end if
            if(nrc.ne.0) then
                z = (up(j) - ym_o) / sqrt(yv_o)
                if(rdnorm(z, 0.0d0, 1.0d0, 0).gt.(tol * rpnorm(z, 0.0d0, 1.0d0, 1, 0))) then
                    ratio = rdnorm(z, 0.0d0, 1.0d0, 0) / rpnorm(z, 0.0d0, 1.0d0, 0, 0)
                    tmean = ym_o + sqrt(yv_o) * ratio
                    tvar = yv_o * (1.d0 + z * ratio - ratio**2)
                else
                    tmean = ym_o
                    tvar = yv_o
                end if
                T1_n = T1_n + nrc * tmean
                T2_n = T2_n + nrc * (tvar + tmean**2)
            end if
            if(nmar.ne.0) then
                T1_n = T1_n + nmar * ym_o
                T2_n = T2_n + nmar * (yv_o + ym_o**2)
            end if
            ym_n = T1_n / n
            yv_n = T2_n / n - ym_n**2
            dm = abs(ym_n - ym_o)
            dv = abs(yv_n - yv_o)
            if(max(dm, dv).le.eps) exit
            ym_o = ym_n
            yv_o = yv_n
        end do
        if(i.ge.nstp) then
            conv = 1
            return
        end if
        ym(j) = ym_n
        yv(j) = yv_n
        deallocate(yobs, stat = conv)
        if(conv.ne.0) then
            conv = -1
            return
        end if
    end if
end do
end subroutine fitmcgm
