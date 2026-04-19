module CG_method_mod
    use iso_fortran_env, only: real64

    implicit none
    integer, parameter :: dp = real64

    !public :: CG_method_1 ! Using the mat A of Ax = b explicitly
    public :: CG_method ! Using the mat A of Ax = b implicitly, i.e., using the function to calculate Ax.
    public :: mat_to_vec
    public :: index_3D_to_1D
    


contains
    subroutine mat_to_vec(A, x) ! reshape the 3D array A to a 1D array x, where A is the discretized laplacian operator in 3D, size is (n_x*n_y*n_z, n_x*n_y*n_z), and x is the vectorized form of A, size is (n_x*n_y*n_z * n_x*n_y*n_z)
        implicit none
        real(dp), intent(in) :: A(:,:,:)
        real(dp), intent(out) :: x(:)
        
        integer :: i, j, k, n_x, n_y, n_z
        n_x = size(A, 1)
        n_y = size(A, 2)
        n_z = size(A, 3)
        if (size(x) /= n_x*n_y*n_z) stop "mat_to_vec: size of x must be equal to n_x*n_y*n_z"
        x = reshape(A, [n_x*n_y*n_z]) 
    end subroutine mat_to_vec

    subroutine index_3D_to_1D( n_x, n_y, n_z, index)
        implicit none
        integer :: i, j, k
        integer, intent(in) :: n_x, n_y, n_z
        integer, intent(out) :: index(n_x,n_y,n_z) ! index(i,j,k) gives the corresponding index in the 1D array for the 3D grid point (i,j,k)

        if (i < 1 .or. i > n_x .or. j < 1 .or. j > n_y .or. k < 1 .or. k > n_z) then
            stop "index_3D_to_1D: i, j, k must be within the bounds of the grid"
        end if
        do k = 1, n_z
            do j = 1, n_y
                do i = 1, n_x
                    index(i,j,k) = (k-1)*n_x*n_y + (j-1)*n_x + i
                end do
            end do
        end do
    end subroutine index_3D_to_1D

    subroutine initialize_CG_method(n_x, n_y, n_z, h_x, a_pot, A)
        implicit none
        integer, intent(in) :: n_x, n_y, n_z
        real(dp), intent(in) :: h_x, a_pot
        real(dp), intent(out) :: A(n_x*n_y*n_z, n_x*n_y*n_z)
        
        integer :: i, j, k, row
        real(dp) :: h2inv
        
        h2inv = 1.0_dp / (12.0_dp * h_x * h_x)
        A = 0.0_dp

        ! 3次元格子でループを回す
        do k = 1, n_z
            do j = 1, n_y
                do i = 1, n_x
                    ! 現在のグリッド点の1次元インデックス（行番号）
                    row = (k-1)*n_x*n_y + (j-1)*n_x + i
                    
                    ! --- 対角成分 ---
                    A(row, row) = -90.0_dp * h2inv - 1.0_dp / a_pot**2
                    
                    ! --- X方向の隣接点 ---
                    if (i+1 <= n_x) A(row + 1, row)     = 16.0_dp * h2inv
                    if (i-1 >= 1)   A(row - 1, row)     = 16.0_dp * h2inv
                    if (i+2 <= n_x) A(row + 2, row)     = -1.0_dp * h2inv
                    if (i-2 >= 1)   A(row - 2, row)     = -1.0_dp * h2inv
                    
                    ! --- Y方向の隣接点 ---
                    if (j+1 <= n_y) A(row + n_x, row)   = 16.0_dp * h2inv
                    if (j-1 >= 1)   A(row - n_x, row)   = 16.0_dp * h2inv
                    if (j+2 <= n_y) A(row + 2*n_x, row) = -1.0_dp * h2inv
                    if (j-2 >= 1)   A(row - 2*n_x, row) = -1.0_dp * h2inv
                    
                    ! --- Z方向の隣接点 ---
                    if (k+1 <= n_z) A(row + n_x*n_y, row)   = 16.0_dp * h2inv
                    if (k-1 >= 1)   A(row - n_x*n_y, row)   = 16.0_dp * h2inv
                    if (k+2 <= n_z) A(row + 2*n_x*n_y, row) = -1.0_dp * h2inv
                    if (k-2 >= 1)   A(row - 2*n_x*n_y, row) = -1.0_dp * h2inv
                end do
            end do
        end do
    end subroutine initialize_CG_method

    





    subroutine CG_method(A,b,x)
        implicit none
        real(dp), intent(in) :: A(:,:)
        real(dp), intent(in) :: b(:)
        real(dp), intent(out) :: x(:)
        real(dp) :: tol ! tolerance for convergence
        integer  :: max_iter
        real(dp), allocatable :: temp(:), r(:), p(:)
        real(dp) :: alpha, beta

        integer :: n
        n = size(b)
        max_iter = 10000
        tol = 1.0e-20_dp
        allocate(temp(n), r(n), p(n))
        if (size(A,1) /= n .or. size(A,2) /= n) stop "CG_method: size of A must be (n,n)"

        ! Initial guess

    end subroutine CG_method

    




end module CG_method_mod