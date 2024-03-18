!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Luigi Augugliaro
! e-mail: luigi.augugliaro@unipa.it
! webpage: http://dssm.unipa.it/augugliaro
!
! version: 1.0.0
! Data: July 12, 2018
!
! Description
! 'Find_ConnectedComp' is used to identify the connected components; it implements the
! results given in Witten et al. (2011).
!
! Input
! p = dimension of the matrix S                                                     (integer)
! S = p x p dimensional variance/covariance matrix                                  (double precision)
! rho = value of the tuning parameter                                               (double precision)
!
! Output
! k = the number of connected components identified using the results given in      (integer)
!       Witten et al., (2011)
! Ck = a p-dimensional vector of integers used to indentify the connected           (integer)
!       components
! pk = a p-dimensional vector of integers spencifying the number of vertices        (integer)
!       belinging to a given connected component, i.e., 'pk(i)' is the number
!       of vertices of the i-th component
subroutine Find_ConnectedComp(p,S,rho,k,Ck,pk)
integer :: p,k,Ck(p),pk(p)
double precision :: S(p,p),rho(p,p)
!internal variables
integer :: i,j,v,Ck_k(p),pk_k,nc(0:p),ncount,maxncount
logical :: Adj(p,p),connected(p)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computing the adjacency matrix                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Adj = .false.
ncount = 0
do j = 1, p
    Adj(j,j) = .true.
    do i = j + 1, p
        if(abs(S(i,j)).gt.rho(i,j)) then
            ncount = ncount + 1
            Adj(i,j) = .true.
            Adj(j,i) = .true.
        end if
    end do
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 'ncount = 0' => all vertices are isolated              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(ncount.eq.0) then
    k = p
    pk = 1
    Ck = (/ (i, i = 1, p) /)
    return
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 'ncount = maxncount' => all vertices are in the same connected component !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
maxncount = (p - 1) * p / 2
if(ncount.eq.maxncount) then
    k = 1
    pk(k) = p
    pk(2:p) = 0
    Ck = (/ (i, i = 1, p) /)
    return
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! we us the dfs algorithm to find the connected components !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
connected = .false.
k = 0
nc(0) = 0
Ck_k = 0
pk = 0
do v = 1, p
    if(.not.connected(v)) then
        call dfs(v,p,Adj,Ck_k,pk_k)
        k = k + 1
        connected(Ck_k(1:pk_k)) = .true.
        pk(k) = pk_k
        nc(k) = nc(k - 1) + pk_k
        Ck((nc(k - 1) + 1):nc(k)) = Ck_k(1:pk_k)
    end if
    if(nc(k).eq.p) exit
end do
end subroutine Find_ConnectedComp
