module constants
    implicit none
    real(8), parameter :: pi = 3.14159265358979323846d0
    real(8), parameter :: a0 = 0.529e-10d0  ! Bohr radius in meters
end module constants

module functions
    use constants
    implicit none

contains

    ! Function to compute S3^T
    function S3_T(lambda) result(S3)
        real(8), intent(in) :: lambda
        real(8) :: S3
        
        ! Computing the terms in the formula for S3^T
        S3 = -(1.0d0 / (60.0d0 * pi)) * ( &
               5.0d0 - (lambda + 5.0d0 / lambda) * atan(lambda) & 
               - 2.0d0 / lambda * asin(lambda / sqrt(1.0d0 + lambda**2)) + &
               (2.0d0 / (lambda * sqrt(2.0d0 + lambda**2))) * ( &
               pi / 2.0d0 - atan(lambda / sqrt(2.0d0 + lambda**2))) &
              )
    end function S3_T

    ! Function to compute a3_T given n
    function a3_T(n) result(a3)
        real(8), intent(in) :: n
        real(8) :: a3, lambda, S3
        
        ! Compute lambda
        lambda = (3.0d0**(1.0d0/6.0d0)) * (pi**(5.0d0/6.0d0)) * (a0**0.5d0) * (n**(1.0d0/6.0d0))
        
        ! Get S3^T
        S3 = S3_T(lambda)
        
        ! Compute a3^T
        a3 = 2.0d0 * (1.0d0 / 3.0d0)**(1.0d0/3.0d0) * (pi**(-2.0d0/3.0d0)) * &
             (3.0d0 / (4.0d0 * pi * n))**(2.0d0/3.0d0) * S3
    end function a3_T

end module functions

program main
    use functions
    implicit none
    real(8) :: n, result
    integer :: i
    real(8), parameter :: n_min = 1.0d27, n_max = 1.0d31
    integer, parameter :: num_points = 200
    real(8) :: n_values(num_points)

    ! Generate n values in a log scale
    do i = 1, num_points
        n_values(i) = n_min * (n_max / n_min)**(real(i-1) / (num_points-1))
    end do

    ! Loop over all values of n and calculate a3_T
    do i = 1, num_points
        n = n_values(i)
        result = a3_T(n)
        print *, "n = ", n, " a3_T = ", result
    end do

end program main
