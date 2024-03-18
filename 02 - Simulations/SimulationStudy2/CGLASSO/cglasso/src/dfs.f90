!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Luigi Augugliaro
! e-mail: luigi.augugliaro@unipa.it
! webpage: http://dssm.unipa.it/augugliaro
!
! version: 1.0.0
! Data: July 12, 2018
!
! Description
! 'dfs' implements the depth-first search algorithm to identify the vertices connected with a given vertex.
!
! Input
! v = index of a given vertex                                                       (integer)
! p = dimension of the adjacency matrix Adj                                         (integer)
! Adj = p x p adjacency matrix                                                      (logical)
!
! Output
! Ck = a p-dimensional vector of integers used to indentify the connected           (integer)
!       components
! pk = a p-dimensional vector of integers pencifying the number of verices          (integer)
!       belinging to a given connected component, i.e., 'pk(i)' is the number
!       of vertices of the i-th component
subroutine dfs(v,p,Adj,Ck,pk)
implicit none
integer :: v,p,Ck(p),pk
logical :: Adj(p,p)
!internal variables
logical :: visited(p),isolated
integer :: stack(p),i,j,k
visited = .false.
visited(v) = .true.
pk = 1
Ck(pk) = v
stack = 0
k = 1
stack(k) = v
do
    isolated = .true.
    i = stack(k)
    do j = 1, p
        if(.not.visited(j).and.Adj(i,j)) then
            visited(j) = .true.
            pk = pk + 1
            Ck(pk) = j
            k = k + 1
            stack(k) = j
            isolated = .false.
            exit
        end if
    end do
    if(isolated) k = k - 1
    if(k.eq.0) exit
end do
end subroutine dfs
