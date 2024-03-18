!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Luigi Augugliaro
! e-mail: luigi.augugliaro@unipa.it
!
! version: 1.0.0
! Data: July 12, 2018
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
subroutine glassosub_v2(p,S,pendiag,rho,wTht,maxit,thr,Tht,k,Ck,pk,nit,conv,trace)
implicit none
integer :: p,pendiag,maxit,k,Ck(p),pk(p),nit,conv,trace
double precision :: S(p,p),rho,wTht(p,p),thr,Tht(p,p)
! internal variables
integer :: i,m,p_k1,p_k2,nnit,df
double precision :: tau
integer, dimension(:), allocatable :: idx1, idx2
double precision, dimension(:,:), allocatable :: wTht_k,S_k,Tht_k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Finding the connected components             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call Find_ConnectedComp(p, S, rho * wTht, k, Ck, pk)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Starting optimization                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
nit = 0
conv = 0
do i = 1, k
    call rchkusr()

    !if(trace.eq.2) call glasso_trace_2_2(i, k)
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
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! if p_k1 is equal to 1 then we are working with an isolated vertix      !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(pendiag.eq.1) S(idx1, idx1) = S(idx1, idx1) + rho * wTht(idx1, idx1)
        Tht(idx1, idx1) = 1.d0 / S(idx1, idx1)
        !if(trace.eq.2) then
        !    call admm_trace_1_2()
        !    call admm_trace_1(1, 0.d0, 0.d0)
        !end if
    else
        !if(trace.eq.9) call glasso_trace_2_2(i, k)
        !if(trace.eq.2) call glasso_trace_2_2(i, k)
        if(trace.eq.2 .or. trace.eq.9) call glasso_trace_2_2(i, k)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! glasso model is applied to the k-th connected component                !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        allocate(wTht_k(1:p_k1,1:p_k1), S_k(1:p_k1,1:p_k1), Tht_k(1:p_k1,1:p_k1), stat=conv)
        if(conv.ne.0) then
            conv = -1
            return
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! inizializing submatrices passed to glasso                              !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        wTht_k = wTht(idx1, idx1)
        S_k = S(idx1, idx1)
        if(pendiag.eq.1) then
            do m = 1, p_k1
                S_k(m, m) = S_k(m, m) + rho * wTht_k(m, m)
            end do
        end if
        Tht_k = Tht(idx1, idx1)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! fitting glasso model to the k-th connected component                   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !call glasso_bcd(p_k1, S_k, rho_k, maxit, thr, Sgm_k, Tht_k, nnit, conv, trace)
        tau = 2.d0
        call admm_tht_single(p_k1, S_k, Tht_k, wTht_k, rho, maxit, thr,&
            tau, trace, conv, nnit, df)
        !call admm(p_k1, S_k, Tht_k, Y_k, rho_k, thr, maxit, Sgm_k, nnit, trace)
        if(conv.ne.0) return
        nit = nit + nnit
        Tht(idx1, idx1) = Tht_k
        deallocate(wTht_k, S_k, Tht_k, stat = conv)
        if(conv.ne.0) then
            conv = -1
            return
        end if
    end if
    if (p_k2.ne.0) then
        Tht(idx1, idx2) = 0.d0
        Tht(idx2, idx1) = 0.d0
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
end subroutine glassosub_v2
