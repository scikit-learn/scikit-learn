function logsum(N,X) result(Xsum)
    ! calculate log(sum exp(x))
    implicit none
    integer, intent(in) :: N ! length of X
    real(8), intent(in) :: X(N) ! data
    real(8) :: Xsum, maxv
    integer :: i
    maxv = maxval(X)
    Xsum = 0.0d0
    do i = 1, N
        Xsum = Xsum + exp(X(i) - maxv)
    enddo
    Xsum = log(Xsum) + maxv
end function

subroutine forward_F(T,N,lnpi,lnA,lnf,lnalpha)
    implicit none
    integer, intent(in) :: T,N
    real(8), intent(in) :: lnpi(N), lnA(N,N), lnf(T,N)
    real(8), intent(out) :: lnalpha(T,N)
    integer :: tt,i,j
    real(8) :: temp(N)
    real(8), external :: logsum
    
    lnalpha(1,:) = lnpi(:) + lnf(1,:)
    
    do tt = 2, T
        do j = 1, N
            temp(:) = lnalpha(tt-1,:) + lnA(:,j)
            lnalpha(tt,j) = logsum(N,temp) + lnf(tt,j)
        enddo
    enddo
end subroutine

subroutine backward_F(T,N,lnpi,lnA,lnf,lnbeta)
    implicit none
    integer, intent(in) :: T,N
    real(8), intent(in) :: lnpi(N), lnA(N,N), lnf(T,N)
    real(8), intent(out) :: lnbeta(T,N)
    integer :: tt,i,j
    real(8) :: temp(N)
    real(8), external :: logsum
    
    lnbeta(T,:) = 0.0d0
    
    do tt = T-1,1,-1
        do i = 1, N
            temp(:) = lnA(i,:) + lnf(tt+1,:) + lnbeta(tt+1,:)
            lnbeta(tt,i) = logsum(N,temp)
        enddo
    enddo
end subroutine

subroutine compute_lnEta_F(T,N,lnalpha,lnA,lnbeta,lnf,lnP_f,lneta)
    implicit none
    integer, intent(in) :: T,N
    real(8), intent(in) :: lnA(N,N), lnf(T,N), lnalpha(T,N), lnbeta(T,N), lnP_f
    real(8), intent(out) :: lneta(T-1,N,N)
    integer :: i,j,tt

    do tt = 1, T-1
        do i = 1, N
            do j = 1, N
                lneta(tt,i,j) = lnalpha(tt,i) + lnA(i,j) &
                    + lnf(tt+1,j) + lnbeta(tt+1,j) - lnP_f
            enddo
        enddo
    enddo
end subroutine

subroutine viterbi_F(T,N,lnpi,lnA,lnf,z,lnP)
    implicit none
    integer, intent(in) :: T,N
    real(8), intent(in) :: lnpi(N), lnA(N,N), lnf(T,N)
    integer(8), intent(out) :: z(T)
    real(8), intent(out) :: lnP
    integer :: tt,i,j,imax(1)
    real(8) :: lndelta(T,N)
    real(8) :: temp(N)

    ! Initialization
    lndelta(1,:) = lnpi + lnf(1,:)

    ! Induction
    do tt = 2, T
        do j = 1, N
            temp(:) = lndelta(tt-1,:) + lnA(:,j)
            lndelta(tt,j) = maxval(temp(:)) + lnf(tt,j)
        enddo
    enddo

    lnP = maxval(lndelta(T,:))

    ! Traceback
    imax = maxloc(lndelta(T,:))
    z(T) = imax(1)
    
    do tt = T-1, 1, -1
        imax = maxloc(lndelta(tt,:) + lnA(:,z(tt+1)))
        z(tt) = imax(1)
    enddo

    ! make numbering consistent to that of python
    z(:) = z(:) - 1

end subroutine
