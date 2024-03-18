!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Luigi Augugliaro & Gianluca Sottile
! e-mail: luigi.augugliaro@unipa.it, gianluca.sottile@unipa.it
! webpage: http://dssm.unipa.it/augugliaro http://bit.ly/gianlucasottile
!
! version: 1.0.0
! Data: June 19, 2023
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
subroutine graph_adjacency(p,N,fk,S,W,rho,alpha,k,Ck,pk)
integer :: p,N,k,Ck(p),pk(p)
double precision :: fk(N),S(p,p,N),W(p,p,N),rho
!internal variables
integer :: k2,i,j,v,Ck_k(p),pk_k,nc(0:p),ncount,maxncount
double precision :: tempsum(p,p,N),normsoftTht(p,p),rho2(p,p)
logical :: Adj(p,p),connected(p)
rho2 = (1 - alpha) * rho
normsoftTht = 0.d0
do k2 = 1, N
    tempsum(:,:,k2) = fk(k2) * abs(S(:,:,k2)) - alpha * rho * W(:,:,k2)
    do i = 1, p
        do j = i, p
            if(tempsum(i,j,k2).lt.0.d0) then
                tempsum(i,j,k2) = 0.d0
                tempsum(j,i,k2) = 0.d0
            end if
        end do
    end do
    !call softmatrix_tht(tempsum(:,:,k2), alpha * rho * W(:,:,k2), p)
    normsoftTht = normsoftTht + tempsum(:,:,k2)**2
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computing the adjacency matrix                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Adj = .false.
ncount = 0
do j = 1, p
    Adj(j,j) = .true.
    do i = j + 1, p
        if(normsoftTht(i,j).gt.(rho2(i,j)**2)) then
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
end subroutine graph_adjacency
