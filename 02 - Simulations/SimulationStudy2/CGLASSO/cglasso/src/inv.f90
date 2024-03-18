subroutine inv(n, A, invA, info)
integer :: info, n
double precision :: A(n, n), invA(n, n)
integer :: i, j
invA = A
call DPOTRF('U', n, invA, n, info)
if (info.ne.0) return
call DPOTRI('U', n, invA, n, info)
if (info.ne.0) return
do j = 1, n
    do i = j + 1, n
        invA(i, j) = invA(j, i)
    end do
end do
end subroutine inv
