!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Luigi Augugliaro
! e-mail: luigi.augugliaro@unipa.it
! webpage: http://dssm.unipa.it/augugliaro
!
! version: 2.0.0
! Data: June 29, 2020
!
! DESCRIPTION
! 'setup' returns matrices R and Y. Rows of the Matrix Y are ordered according to the pattern of missing data
!
! Arguments
!
!   n   [IN]        = sample size                                                                                       (integer)
!   p   [IN]        = number of variables                                                                               (integer)
!   Y   [IN/OUT]    = n x p dimensional matrix. In output the rows are ordered according to the pattern of missing data (double)
!   lo  [IN]        = p dimensional vector of left censored values                                                      (double)
!   up  [IN]        = p dimensional vector of right censored values                                                     (double)
!   Yna [IN]        = double coding the Not Availabe data                                                               (double)
!   R   [OUT]       = n x p dimensional matrix coding information about the missing data, that is:                      (integer)
!       'R(i, j) = 0' (i = startmis...n and j = 1...p) means that Y[i, j] NotAvaliable
!       'R(i, j) = 1' (i = startmis...n and j = 1...p) means that Y[i, j] is a left censored data
!       'R(i, j) = 2' (i = startmis...n and j = 1...p) means that Y[i, j] is a right censored data
!       'R(i, j) = 3' (i = startmis...n and j = 1...p) means that Y[i, j] is observed
!   startmis    [OUT]   = the starting row of the censored values                                                       (integer)
!   order   [OUT]   = n-dimensional vector                                                                              (integer)
subroutine setup(n, p, Y, lo, up, Yna, R, startmis, order)
implicit none
integer :: n, p, R(0:n, 0:p), startmis, order(n)
double precision :: Y(n, p), lo(p), up(p), Yna
! internal variables
integer :: i,j,jj, tmp1(p)
double precision :: tmp2(p)
R(:, 1:p) = 3
do j = 1, p
    where (Y(:, j).eq.Yna) R(1:n, j) = 0
    where (Y(:, j).eq.lo(j)) R(1:n, j) = 1
    where (Y(:, j).eq.up(j)) R(1:n, j) = 2
end do
jj = -1
do
    jj = jj + 1
    if(jj.le.n) R(jj, 0) = 1
    if(jj.ge.n) exit
    j = jj
    do i = j + 1, n
        if(all(R(j, 1:p).eq.R(i, 1:p))) then
            jj = jj + 1
            tmp1 = R(jj, 1:p)
            R(jj, 1:p) = R(i, 1:p)
            R(i, 1:p) = tmp1
            tmp1(1) = order(jj)
            order(jj) = order(i)
            order(i) = tmp1(1)
            tmp2 = Y(jj, :)
            Y(jj, :) = Y(i, :)
            Y(i, :) = tmp2
        end if
    end do
end do
startmis = 0
do i = 1, n
if (any(R(i,1:p).ne.3)) then
        startmis = i
        exit
    end if
end do
end subroutine setup
