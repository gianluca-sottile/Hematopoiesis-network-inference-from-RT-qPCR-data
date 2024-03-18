!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Luigi Augugliaro & Gianluca Sottile
! e-mail: luigi.augugliaro@unipa.it gianluca.sottile@unipa.it
!
! version: 1.0.0
! Data: June 19, 2023
!
! Description
! 'glassosub' implements the block coordinate discent (bcd) algorithm proposed in Friedman et al. (2008)
! to fit a l1-penalized guassian graphical model. We also implement the results given in Witten et al.(2011)
! to reduce the computational burdern of the algorithm.
!
! Input
! p = dimension of the matrix S                                                     (integer)
! S = p x p dimensional variance/covariance matrix                                  (double precision)
! pendiag = integer used to specify if the diagonal elements of the                 (integer)
!           concentration matrix (Tht) are penalized
!           'pendiag = 0' means that the diagonal elements are not penalized
!           'pendiag = 1' means that the diagonal elements are penalized
! rho = p x p dimensional matrix of the values of the tuning parameter              (double precision)
! maxit = maximum number of steps of the bcd algorithm                              (integer)
! thr = threshold value used to declare the convergence of the bcd algorithm        (double precision)
!
! Output
! Sgm = estiamted l1-penalized variance/covariance matrix                           (double precision)
! Tht = estiamted l1-penalized concentration matrix                                 (double precision)
!
!   NOTE: Sgm and Tht are also used in input as starting values for the bcd algorithm
!
! k = output of the subroutine 'Find_ConnectedComp'; 'k' is the number              (integer)
!       of connected components identified using the results given in
!       Witten et al., (2011)
! Ck = output of the subroutine 'Find_ConnectedComp'; 'Ck' is a p-dimensional       (integer)
!       vector of integers used to indentify the connected components
! pk = output of the subroutine 'Find_ConnectedComp'; 'pk' is a p-dimensional       (integer)
!       vector of integers pencifying the number of verices belonging to a
!       connected component, i.e., 'pk(i)' is the number of vertices of the
!       i-th component
! nit = number of steps of the bcd algorithm                                        (integer)
! conv = integer used to encode the convergence of the algorithm:                   (integer)
!          '-1' error in memory allocation
!           '0' convergence is met
!           '1' maximum number of steps is reached
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine admm_tht_sub(p,N,fk,S,wTht,pendiag,rho,alpha,maxit,thr,Tht,k,Ck,pk,nit,df,conv,trace)
implicit none
integer :: p,N,pendiag,maxit,k,Ck(p),pk(p),nit,df(N),conv,trace
double precision :: fk(N),S(p,p,N),wTht(p,p,N),rho,alpha,thr,Tht(p,p,N)
! internal variables
integer :: i,p_k1,p_k2,k2,m,nnit,j
double precision :: tau
integer, dimension(:), allocatable :: idx1,idx2
double precision, dimension(:,:,:), allocatable :: wTht_k,S_k,Tht_k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Finding the connected components             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!call Find_ConnectedComp(p, S, rho, k, Ck, pk)
call graph_adjacency(p, N, fk, S, wTht, rho, alpha, k, Ck, pk)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Starting optimization                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
nit = 0
conv = 0
do i = 1, k
    call rchkusr()
    
    !if(trace.eq.2 .or. trace.eq.9) call glasso_trace_2_2(i, k)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! extraction of the indices associated with the k-th connected component !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    p_k1 = pk(1)
    p_k2 = p - p_k1
    allocate(idx1(1:p_k1), idx2(1:p_k2), stat = conv)
    if(conv.ne.0) then
        conv = -1
        return
    end if
    idx1 = Ck(1:p_k1)
    if (p_k2.ne.0) then
        idx2 = Ck((p_k1 + 1):p)
    end if
    if(p_k1 .eq. 1) then
        do k2 = 1, N
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! if p_k1 is equal to 1 then we are working with an isolated vertix      !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if(pendiag.eq.1) S(idx1, idx1, k2) = S(idx1, idx1, k2) + alpha * rho * wTht(idx1, idx1, k2)
            Tht(idx1, idx1, k2) = 1.d0 / S(idx1, idx1, k2)
        end do
    else
        !if(trace.eq.9) call glasso_trace_2_2(i, k)
        if(trace.eq.2 .or. trace.eq.9) call glasso_trace_2_2(i, k)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! glasso model is applied to the k-th connected component                !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        allocate(wTht_k(1:p_k1,1:p_k1,1:N), S_k(1:p_k1,1:p_k1,1:N), Tht_k(1:p_k1,1:p_k1,1:N), stat=conv)
        if(conv.ne.0) then
            conv = -1
            return
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! inizializing submatrices passed to glasso                              !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do k2 = 1, N
            wTht_k(:, :, k2) = wTht(idx1, idx1, k2)
            S_k(:, :, k2) = S(idx1, idx1, k2)
            if(pendiag.eq.1) then
                do m = 1, p_k1
                    S_k(m, m, k2) = S_k(m, m, k2) + alpha * rho * wTht_k(m, m, k2)
                end do
            end if
            Tht_k(:, :, k2) = Tht(idx1, idx1, k2)
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! fitting glasso model to the k-th connected component                   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        tau = 2.d0
        call admm_tht_joint(p_k1, N, fk, S_k, Tht_k, wTht_k, rho, alpha, maxit, thr, tau, trace, conv, nnit, df)
        if(conv.ne.0) return
        nit = nit + nnit
        do k2 = 1, N
            Tht(idx1, idx1, k2) = Tht_k(:, :, k2)
        end do
        deallocate(wTht_k, S_k, Tht_k, stat = conv)
        if(conv.ne.0) then
            conv = -1
            return
        end if
    end if
    if (p_k2.ne.0) then
        do k2 = 1, N
            Tht(idx1, idx2, k2) = 0.d0
            Tht(idx2, idx1, k2) = 0.d0
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! shifting Ck and pk vectors                                             !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Ck(1:(p_k2)) = idx2
        Ck((p_k2 + 1):p) = idx1
        pk(1:(k - 1)) = pk(2:k)
        pk(k) = p_k1
    end if
    deallocate(idx1, idx2, stat = conv)
    if(conv.ne.0) then
        conv = -1
        return
    end if
end do
df = 0
do k2 = 1, N
    do i = 1, p - 1
        do j = i + 1, p
            if(Tht(i, j, k2).ne.0.d0) df(k2) = df(k2) + 1
        end do
    end do
end do
end subroutine admm_tht_sub
