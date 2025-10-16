PROGRAM test_multi_precision_linear_algebra
  USE multi_precision_float_mod, ONLY: mpf, &
      new_mpf_from_string, mpf_value_equal, mpf_abs, new_mpf_from_mpi_exp, mpf_to_string
  USE multi_precision_integer_mod, ONLY: mpi, new_mpi_from_integer, OPERATOR(==)
  USE multi_precision_linear_algebra_mod
  IMPLICIT NONE
  
  WRITE(*, '(/A)') "Running tests for MultiPrecisionLinearAlgebra module..."
  WRITE(*, '(/A)') "===================================================="

  ! CALL test_vector_construction()
  ! CALL test_matrix_construction()
  ! CALL test_vector_arithmetic()
  ! CALL test_matrix_arithmetic()
  ! CALL test_level1_blas()
  ! CALL test_level2_blas()
  ! CALL test_level3_blas()

  CALL test_large_dot_product()

  WRITE(*, '(/A)') "===================================================="
  WRITE(*, '(A/)') "All MultiPrecisionLinearAlgebra tests completed."

CONTAINS

  SUBROUTINE assert(condition, message)
    LOGICAL, INTENT(IN) :: condition
    CHARACTER(LEN=*), INTENT(IN) :: message
    IF (condition) THEN
      WRITE(*, '(A, A)') "  [PASS] ", message
    ELSE
      WRITE(*, '(A, A)') "  [FAIL] ", message
      STOP "Test assertion failed."
    END IF
  END SUBROUTINE assert

  SUBROUTINE test_vector_construction()
    TYPE(mpf), DIMENSION(3) :: mpf_elements
    TYPE(mpf_vector) :: vec
    TYPE(mpi) :: mpi_val_1, mpi_val_2, mpi_val_3

    WRITE(*, '(A)') "Testing block floating-point vector construction..."

    mpf_elements(1) = new_mpf_from_string("12.0") ! 3 * 2^2
    mpf_elements(2) = new_mpf_from_string("1.5")  ! 3 * 2^-1
    mpf_elements(3) = new_mpf_from_string("24.0") ! 3 * 2^3

    vec = new_mpf_vector_from_mpfs(mpf_elements)

    CALL assert(vec%n == 3, "Vector Construction: Size")
    CALL assert(vec%exponent == -1, "Vector Construction: Common exponent")

    mpi_val_1 = new_mpi_from_integer(24_8) ! 12.0 scaled to exp=-1 -> 12.0 * 2^1 = 24
    mpi_val_2 = new_mpi_from_integer(3_8)  ! 1.5 scaled to exp=-1 -> 1.5 * 2^1 = 3
    mpi_val_3 = new_mpi_from_integer(48_8) ! 24.0 scaled to exp=-1 -> 24.0 * 2^1 = 48

    CALL assert(vec%mantissas(1) == mpi_val_1, "Vector Construction: Mantissa 1 scaled correctly")
    CALL assert(vec%mantissas(2) == mpi_val_2,  "Vector Construction: Mantissa 2 scaled correctly")
    CALL assert(vec%mantissas(3) == mpi_val_3, "Vector Construction: Mantissa 3 scaled correctly")

  END SUBROUTINE test_vector_construction

  SUBROUTINE test_matrix_construction()
    TYPE(mpf), DIMENSION(2,2) :: mpf_elements
    TYPE(mpf_matrix) :: mat
    
    WRITE(*, '(A)') "Testing block floating-point matrix construction..."

    mpf_elements(1,1) = new_mpf_from_string("1.0")
    mpf_elements(1,2) = new_mpf_from_string("2.0")
    mpf_elements(2,1) = new_mpf_from_string("0.5") ! 1 * 2^-1
    mpf_elements(2,2) = new_mpf_from_string("4.0")

    mat = new_mpf_matrix_from_mpfs(mpf_elements)
    CALL assert(mat%m == 2, "Matrix Construction: Rows")
    CALL assert(mat%n == 2, "Matrix Construction: Columns")
    CALL assert(mat%exponent == -1, "Matrix Construction: Common exponent")
    ! 1.0 -> exp 0, mant 1. Shifted to exp -1 -> mant 2
    ! 0.5 -> exp -1, mant 1. Shifted to exp -1 -> mant 1
    CALL assert(mat%mantissas(1,1) == new_mpi_from_integer(2_8), "Matrix Construction: Mantissa (1,1) scaled correctly")
    CALL assert(mat%mantissas(2,1) == new_mpi_from_integer(1_8), "Matrix Construction: Mantissa (2,1) scaled correctly")

  END SUBROUTINE test_matrix_construction

  SUBROUTINE test_vector_arithmetic()
    TYPE(mpf_vector) :: vec1, vec2, vec_res
    TYPE(mpf) :: dot_res, expected_res, scalar

    WRITE(*, '(A)') "Testing vector arithmetic..."

    vec1 = new_mpf_vector_from_mpfs([new_mpf_from_string("1.5"), new_mpf_from_string("2.0")])
    vec2 = new_mpf_vector_from_mpfs([new_mpf_from_string("3.0"), new_mpf_from_string("4.25")])

    PRINT *, "--- Inside test_vector_arithmetic ---"
    
    ! Dot Product
    dot_res = vec1 * vec2
    expected_res = new_mpf_from_string("13.0")
    CALL assert(mpf_value_equal(dot_res, expected_res), "Vector Dot Product")

    ! Addition
    vec_res = vec1 + vec2
    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(vec_res%mantissas(1), vec_res%exponent), new_mpf_from_string("4.5")), "Vector Addition: Element 1")
    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(vec_res%mantissas(2), vec_res%exponent), new_mpf_from_string("6.25")), "Vector Addition: Element 2")

    ! Scalar Multiplication
    scalar = new_mpf_from_string("2.0")
    vec_res = vec1 * scalar
    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(vec_res%mantissas(1), vec_res%exponent), new_mpf_from_string("3.0")), "Vector-Scalar Mul: Element 1")
    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(vec_res%mantissas(2), vec_res%exponent), new_mpf_from_string("4.0")), "Vector-Scalar Mul: Element 2")

  END SUBROUTINE test_vector_arithmetic

  SUBROUTINE test_matrix_arithmetic()
      TYPE(mpf_matrix) :: mat1, mat2, mat_res
      TYPE(mpf_vector) :: vec, vec_res

      WRITE(*, '(A)') "Testing matrix arithmetic..."

      mat1 = new_mpf_matrix_from_mpfs(RESHAPE([new_mpf_from_string("1.0"), new_mpf_from_string("3.0"), &
                                               new_mpf_from_string("2.0"), new_mpf_from_string("4.0")], [2,2]))
      vec = new_mpf_vector_from_mpfs([new_mpf_from_string("5.0"), new_mpf_from_string("6.0")])

      ! Matrix-Vector product
      vec_res = mat1 * vec
      ! Expected: [1*5+2*6, 3*5+4*6] = [17, 39]
      CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(vec_res%mantissas(1), vec_res%exponent), new_mpf_from_string("17.0")), "Matrix-Vector Mul: Element 1")
      CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(vec_res%mantissas(2), vec_res%exponent), new_mpf_from_string("39.0")), "Matrix-Vector Mul: Element 2")

      ! Matrix-Matrix product
      mat2 = new_mpf_matrix_from_mpfs(RESHAPE([new_mpf_from_string("1.0"), new_mpf_from_string("0.0"), &
                                               new_mpf_from_string("0.0"), new_mpf_from_string("1.0")], [2,2]))
      mat_res = mat1 * mat2
      CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(mat_res%mantissas(1,1), mat_res%exponent), new_mpf_from_string("1.0")), "Matrix-Matrix Mul: Identity (1,1)")
      CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(mat_res%mantissas(1,2), mat_res%exponent), new_mpf_from_string("2.0")), "Matrix-Matrix Mul: Identity (1,2)")

      ! Matrix Addition
      mat_res = mat1 + mat1
      CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(mat_res%mantissas(1,1), mat_res%exponent), new_mpf_from_string("2.0")), "Matrix Addition")

  END SUBROUTINE test_matrix_arithmetic

  SUBROUTINE test_level1_blas()
    TYPE(mpf_vector) :: vec
    TYPE(mpf) :: asum_res, nrm2_res, expected_res
    INTEGER :: idx

    WRITE(*, '(A)') "Testing Level 1 BLAS routines..."

    vec = new_mpf_vector_from_mpfs([new_mpf_from_string("-3.0"), new_mpf_from_string("4.0")])

    ! ASUM
    asum_res = mpf_vector_asum(vec)
    expected_res = new_mpf_from_string("7.0")
    CALL assert(mpf_value_equal(asum_res, expected_res), "ASUM: |-3| + |4| = 7")
    
    ! IAMAX / IAMIN
    vec = new_mpf_vector_from_mpfs([new_mpf_from_string("1.0"), new_mpf_from_string("-5.0"), &
                                    new_mpf_from_string("2.0")])
    idx = mpf_vector_iamax(vec)
    CALL assert(idx == 2, "IAMAX: Index of max absolute value is 2")
    idx = mpf_vector_iamin(vec)
    CALL assert(idx == 1, "IAMIN: Index of min absolute value is 1")

  END SUBROUTINE test_level1_blas

  SUBROUTINE test_level2_blas()
    TYPE(mpf_matrix) :: mat
    TYPE(mpf_vector) :: vec, vec_res, b_vec

    WRITE(*, '(A)') "Testing Level 2 BLAS routines..."

    ! --- TRMV Test ---
    ! A = [1 2]  x = [3]  A*x = [1*3+2*4] = [11]
    !     [0 4]      [4]        [0*3+4*4]   [16]
    mat = new_mpf_matrix_from_mpfs(RESHAPE([new_mpf_from_string("1.0"), new_mpf_from_string("0.0"), &
                                            new_mpf_from_string("2.0"), new_mpf_from_string("4.0")], [2,2]))
    vec = new_mpf_vector_from_mpfs([new_mpf_from_string("3.0"), new_mpf_from_string("4.0")])
    CALL mpf_trmv(mat, vec, 'U')
    vec_res = new_mpf_vector_from_mpfs([new_mpf_from_string("11.0"), new_mpf_from_string("16.0")])
    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(vec%mantissas(1), vec%exponent), new_mpf_from_string("11.0")), "TRMV (Upper): Element 1")
    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(vec%mantissas(2), vec%exponent), new_mpf_from_string("16.0")), "TRMV (Upper): Element 2")

    ! --- TRSV Test ---
    ! A*x = b.  x = inv(A)*b.  Let's solve A*x = [11, 16]
    ! We should get back x = [3, 4]
    b_vec = new_mpf_vector_from_mpfs([new_mpf_from_string("11.0"), new_mpf_from_string("16.0")])
    CALL mpf_trsv(mat, b_vec, 'U')
    vec_res = new_mpf_vector_from_mpfs([new_mpf_from_string("3.0"), new_mpf_from_string("4.0")])
    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(b_vec%mantissas(1), b_vec%exponent), new_mpf_from_string("3.0")), "TRSV (Upper): Element 1")
    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(b_vec%mantissas(2), b_vec%exponent), new_mpf_from_string("4.0")), "TRSV (Upper): Element 2")

    ! Lower triangular
    mat = new_mpf_matrix_from_mpfs(RESHAPE([new_mpf_from_string("3.0"), new_mpf_from_string("1.0"), &
                                            new_mpf_from_string("0.0"), new_mpf_from_string("2.0")], [2,2]))
    b_vec = new_mpf_vector_from_mpfs([new_mpf_from_string("6.0"), new_mpf_from_string("8.0")])
    ! 3x1 = 6 -> x1=2.  1*x1 + 2*x2 = 8 -> 1*2 + 2*x2 = 8 -> 2*x2=6 -> x2=3
    CALL mpf_trsv(mat, b_vec, 'L')
    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(b_vec%mantissas(1), b_vec%exponent), new_mpf_from_string("2.0")), "TRSV (Lower): Element 1")
    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(b_vec%mantissas(2), b_vec%exponent), new_mpf_from_string("3.0")), "TRSV (Lower): Element 2")

  END SUBROUTINE test_level2_blas

  SUBROUTINE test_level3_blas()
    TYPE(mpf_matrix) :: mat_a, mat_b, mat_res

    WRITE(*, '(A)') "Testing Level 3 BLAS routines..."

    ! --- TRMM Test (Upper) ---
    ! A = [1 2], B = [3 4], A*B = [1*3+2*5, 1*4+2*6] = [13, 16]
    !     [0 4]      [5 6]         [0*3+4*5, 0*4+4*6]   [20, 24]
    mat_a = new_mpf_matrix_from_mpfs(RESHAPE([new_mpf_from_string("1.0"), new_mpf_from_string("0.0"), &
                                              new_mpf_from_string("2.0"), new_mpf_from_string("4.0")], [2,2]))
    mat_b = new_mpf_matrix_from_mpfs(RESHAPE([new_mpf_from_string("3.0"), new_mpf_from_string("5.0"), &
                                              new_mpf_from_string("4.0"), new_mpf_from_string("6.0")], [2,2]))
    CALL mpf_trmm(mat_a, mat_b, 'U')
    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(mat_b%mantissas(1,1), mat_b%exponent), new_mpf_from_string("13.0")), "TRMM (Upper): C(1,1)")
    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(mat_b%mantissas(2,1), mat_b%exponent), new_mpf_from_string("20.0")), "TRMM (Upper): C(2,1)")
    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(mat_b%mantissas(1,2), mat_b%exponent), new_mpf_from_string("16.0")), "TRMM (Upper): C(1,2)")
    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(mat_b%mantissas(2,2), mat_b%exponent), new_mpf_from_string("24.0")), "TRMM (Upper): C(2,2)")

    ! --- TRSM Test (Lower) ---
    ! A = [2 0], X = [1 2], B = A*X = [2  4]
    !     [1 3]      [3 4]           [10 14]
    ! We solve A*X=B and expect to get X back.
    mat_a = new_mpf_matrix_from_mpfs(RESHAPE([new_mpf_from_string("2.0"), new_mpf_from_string("1.0"), &
                                              new_mpf_from_string("0.0"), new_mpf_from_string("3.0")], [2,2]))
    mat_b = new_mpf_matrix_from_mpfs(RESHAPE([new_mpf_from_string("2.0"), new_mpf_from_string("10.0"), &
                                              new_mpf_from_string("4.0"), new_mpf_from_string("14.0")], [2,2]))
    CALL mpf_trsm(mat_a, mat_b, 'L')
    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(mat_b%mantissas(1,1), mat_b%exponent), new_mpf_from_string("1.0")), "TRSM (Lower): X(1,1)")
    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(mat_b%mantissas(2,1), mat_b%exponent), new_mpf_from_string("3.0")), "TRSM (Lower): X(2,1)")
    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(mat_b%mantissas(1,2), mat_b%exponent), new_mpf_from_string("2.0")), "TRSM (Lower): X(1,2)")
    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(mat_b%mantissas(2,2), mat_b%exponent), new_mpf_from_string("4.0")), "TRSM (Lower): X(2,2)")

  END SUBROUTINE test_level3_blas

  SUBROUTINE test_large_dot_product()
    INTEGER, PARAMETER :: n = 10

    TYPE(mpf_vector) :: vec1_mpf, vec2_mpf
    TYPE(mpf) :: dot_res_mpf
    TYPE(mpf), ALLOCATABLE :: mpf_elements1(:), mpf_elements2(:)

    REAL(16), ALLOCATABLE :: vec1_q(:), vec2_q(:)
    REAL(16) :: dot_res_q, rel_error

    REAL :: rand_val
    CHARACTER(LEN=50) :: str_val
    INTEGER :: i

    WRITE(*, '(A, I0, A)') "Testing large dot product (n=", n, ")..."

    ALLOCATE(vec1_q(n), vec2_q(n))
    ALLOCATE(mpf_elements1(n), mpf_elements2(n))

    CALL RANDOM_SEED()

    ! 1. Generate random vectors for both mpf and REAL(16)
    DO i = 1, n
      CALL RANDOM_NUMBER(rand_val)
      vec1_q(i) = (REAL(rand_val, 16) - 0.5_16) * 200.0_16
      WRITE(str_val, '(F0.40)') vec1_q(i)
      mpf_elements1(i) = new_mpf_from_string(TRIM(str_val))
      WRITE(*, '(A,I2,A,F45.35)') "DEBUG: vec1_q(", i, ") = ", vec1_q(i)
      WRITE(*, '(A,I2,A,A)')     "DEBUG: vec1_mpf(", i, ") = ", mpf_to_string(mpf_elements1(i))

      CALL RANDOM_NUMBER(rand_val)
      vec2_q(i) = (REAL(rand_val, 16) - 0.5_16) * 200.0_16
      WRITE(str_val, '(F0.40)') vec2_q(i)
      mpf_elements2(i) = new_mpf_from_string(TRIM(str_val))
      WRITE(*, '(A,I2,A,F45.35)') "DEBUG: vec2_q(", i, ") = ", vec2_q(i)
      WRITE(*, '(A,I2,A,A)')     "DEBUG: vec2_mpf(", i, ") = ", mpf_to_string(mpf_elements2(i))
    END DO

    ! 2. Compute dot products
    vec1_mpf = new_mpf_vector_from_mpfs(mpf_elements1)
    vec2_mpf = new_mpf_vector_from_mpfs(mpf_elements2)
    dot_res_mpf = vec1_mpf * vec2_mpf
    dot_res_q = DOT_PRODUCT(vec1_q, vec2_q)

    ! 3. Compare the results by converting the REAL(16) result to an mpf
    WRITE(str_val, '(F0.40)') dot_res_q

    WRITE(*, '(A, A)') "  MPF Result:    ", mpf_to_string(dot_res_mpf)
    WRITE(*, '(A, F36.15)') "  REAL(16) Result: ", dot_res_q

    CALL assert(mpf_value_equal(dot_res_mpf, new_mpf_from_string(TRIM(str_val))), "Large Dot Product: mpf vs REAL(16)")

  END SUBROUTINE test_large_dot_product

END PROGRAM test_multi_precision_linear_algebra