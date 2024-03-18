subroutine predict(newrho, newlambda, nrho, rho, nlambda, lambda, d1, d2, Min, Mout)
implicit none
integer :: nrho, nlambda, d1, d2
double precision :: newrho, newlambda, rho(nrho), lambda(nlambda), Min(d1, d2, nlambda, nrho), Mout(d1, d2)
!!!!!!!!!!!!!!!!!!!!!!
! internal variables !
!!!!!!!!!!!!!!!!!!!!!!
integer :: i, ii, jj
double precision :: dlt
if (newrho.ge.rho(1).and.newlambda.ge.lambda(1)) then
    Mout = Min(:, :, 1, 1)
    return
end if
ii = 0
do i = 1, nrho - 1
    if (rho(i + 1).le.newrho.and.newrho.lt.rho(i)) then
        ii = i
        exit
    end if
end do
jj = 0
do i = 1, nlambda - 1
    if (lambda(i + 1).le.newlambda.and.newlambda.lt.lambda(i)) then
        jj = i
        exit
    end if
end do
Mout = Min(:, :, jj + 1, ii + 1)
if (nrho.gt.1) then
    dlt = 0.0d0
    if (ii.gt.0) dlt = (newrho - rho(ii + 1)) / (rho(ii) - rho(ii + 1))
    if (abs(dlt).gt.0.0d0) Mout = Mout + (Min(:, :, jj + 1, ii) - Min(:, :, jj + 1, ii + 1)) * dlt
end if
if (nlambda.gt.1) then
    dlt = 0.0d0
    if (jj.gt.0) dlt = (newlambda - lambda(jj + 1)) / (lambda(jj) - lambda(jj + 1))
    if (abs(dlt).gt.0.0d0) Mout = Mout + (Min(:, :, jj, ii + 1) - Min(:, :, jj + 1, ii + 1)) * dlt
end if
end subroutine predict
