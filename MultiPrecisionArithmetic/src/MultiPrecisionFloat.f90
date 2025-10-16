MODULE multi_precision_float_mod
  USE multi_precision_integer_mod, ONLY : mpi, COEFFS_LIMIT, MULTI_PRECISION_BASE, &
                                     normalize_mpi, new_mpi_from_coeffs, mpi_shift_bits_left, &
                                     mpi_abs, new_mpi_from_integer, mpi_to_string, &
                                     OPERATOR(==), OPERATOR(<), OPERATOR(+), OPERATOR(-), &
                                     OPERATOR(*), new_mpi_from_string, mpi_is_zero, &
                                     mpi_div_rem, mpi_sign, mpi_shift_bits_right
  IMPLICIT NONE

  TYPE mpf
    TYPE(mpi)              :: mantissa
    INTEGER                :: exponent
  END TYPE mpf

  INTERFACE OPERATOR(==)
    MODULE PROCEDURE mpf_equal
  END INTERFACE

  INTERFACE OPERATOR(<)
    MODULE PROCEDURE mpf_less
  END INTERFACE

  INTERFACE OPERATOR(+)
    MODULE PROCEDURE mpf_add
  END INTERFACE

  INTERFACE OPERATOR(-)
    MODULE PROCEDURE mpf_subtract
    MODULE PROCEDURE mpf_unary_negate
  END INTERFACE

  INTERFACE OPERATOR(*)
    MODULE PROCEDURE mpf_multiply
  END INTERFACE

  INTERFACE OPERATOR(/)
    MODULE PROCEDURE mpf_divide
  END INTERFACE

  PRIVATE
  PUBLIC :: mpf, OPERATOR(==), OPERATOR(<), OPERATOR(+), &
            OPERATOR(-), OPERATOR(*), OPERATOR(/), new_mpf_from_mpi_exp, &
            new_mpf_from_integer, mpf_is_zero, &
            normalize_mpf_float, mpf_to_string, mpf_value_equal, &
            mpf_abs, mpf_scale_up_by_base_power, new_mpf_from_string

CONTAINS

  FUNCTION new_mpf_from_mpi_exp(mantissa_in, exponent_in) RESULT(mpf_out)
    TYPE(mpi), INTENT(IN)              :: mantissa_in
    INTEGER, OPTIONAL, INTENT(IN)        :: exponent_in
    TYPE(mpf)           :: mpf_out

    mpf_out%mantissa = mantissa_in
    IF (PRESENT(exponent_in)) THEN
      mpf_out%exponent = exponent_in
    ELSE
      mpf_out%exponent = 0
    END IF
    CALL normalize_mpf_float(mpf_out)
  END FUNCTION new_mpf_from_mpi_exp
  
SUBROUTINE normalize_mpf_float(mpf_val)
    TYPE(mpf), INTENT(INOUT) :: mpf_val
    TYPE(mpi) :: mantissa_abs
    INTEGER :: num_trailing_zeros, i
    LOGICAL :: was_negative

    was_negative = mpi_sign(mpf_val%mantissa)
    
    CALL normalize_mpi(mpf_val%mantissa)

    IF (mpf_is_zero(mpf_val)) THEN
      mpf_val%exponent = 0
      RETURN
    END IF

    mantissa_abs = mpi_abs(mpf_val%mantissa)

    num_trailing_zeros = 0
    DO i = 1, COEFFS_LIMIT
      IF (mantissa_abs%coeffs(i) == 0_8) THEN
        num_trailing_zeros = num_trailing_zeros + 32
      ELSE
        num_trailing_zeros = num_trailing_zeros + TRAILZ(mantissa_abs%coeffs(i))
        EXIT
      END IF
    END DO

    IF (num_trailing_zeros > 0) THEN
      mantissa_abs = mpi_shift_bits_right(mantissa_abs, num_trailing_zeros)
      CALL normalize_mpi(mantissa_abs)
      mpf_val%exponent = mpf_val%exponent + num_trailing_zeros
    END IF

    IF (was_negative) THEN
      mpf_val%mantissa = -mantissa_abs
    ELSE
      mpf_val%mantissa = mantissa_abs
    END IF

END SUBROUTINE normalize_mpf_float

  ! Creates a HighPrecisionFloat from an Integer.
  FUNCTION new_mpf_from_integer(x_in) RESULT(mpf_out)
    INTEGER(KIND=8), INTENT(IN) :: x_in
    TYPE(mpf)  :: mpf_out
    TYPE(mpi)                 :: hpi_val
    hpi_val = new_mpi_from_integer(x_in)
    mpf_out = new_mpf_from_mpi_exp(hpi_val, 0) ! Exponent is 0 for integers
  END FUNCTION new_mpf_from_integer


  ! Returns the absolute value of a MultiPrecisionFloat.
  FUNCTION mpf_abs(mpf_in) RESULT(mpf_out)
    TYPE(mpf), INTENT(IN) :: mpf_in
    TYPE(mpf)             :: mpf_out
    mpf_out = new_mpf_from_mpi_exp(mpi_abs(mpf_in%mantissa), mpf_in%exponent)
  END FUNCTION mpf_abs


  ! Checks if two MultiPrecisionFloat numbers are equal (a == b)
  FUNCTION mpf_equal(a, b) RESULT(res)
    TYPE(mpf), INTENT(IN) :: a, b
    LOGICAL                                :: res
    res = (a%mantissa == b%mantissa) .AND. (a%exponent == b%exponent)
  END FUNCTION mpf_equal


  ! Compares two HighPrecisionFloat numbers for less than (a < b).
  FUNCTION mpf_less(a, b) RESULT(res)
    TYPE(mpf), INTENT(IN) :: a, b
    LOGICAL                                :: res
    TYPE(mpi)                            :: scaled_a_mantissa, scaled_b_mantissa
    INTEGER                                :: common_exponent, diff_a, diff_b

    IF (mpi_sign(a%mantissa) .NEQV. mpi_sign(b%mantissa)) THEN
      res = mpi_sign(a%mantissa) ! If a is negative and b is not, a < b
      RETURN
    ELSE IF (mpf_is_zero(a)) THEN
      res = .FALSE. ! Both are zero, so not less than
      RETURN
    END IF

    ! To correctly compare values, we must align their exponents first.
    common_exponent = MIN(a%exponent, b%exponent)
    diff_a = a%exponent - common_exponent ! bit shift for a
    diff_b = b%exponent - common_exponent ! bit shift for b

    scaled_a_mantissa = mpi_shift_bits_left(a%mantissa, diff_a)
    scaled_b_mantissa = mpi_shift_bits_left(b%mantissa, diff_b)

    res = (scaled_a_mantissa < scaled_b_mantissa)
  END FUNCTION mpf_less

  ! Compares the numerical value of two HPFs by aligning their exponents.
  ! This is different from the built-in `==` which does a structural comparison.
  LOGICAL FUNCTION mpf_value_equal(a, b)
    TYPE(mpf), INTENT(IN) :: a, b
    INTEGER :: common_exp, diff_a, diff_b
    TYPE(mpi) :: scaled_a, scaled_b

    IF (mpf_is_zero(a) .AND. mpf_is_zero(b)) THEN
        mpf_value_equal = .TRUE.
        RETURN
    END IF

    ! To correctly compare values, we must align their exponents first.
    common_exp = MIN(a%exponent, b%exponent)
    diff_a = a%exponent - common_exp
    diff_b = b%exponent - common_exp

    scaled_a = mpi_shift_bits_left(a%mantissa, diff_a)
    scaled_b = mpi_shift_bits_left(b%mantissa, diff_b)

    mpf_value_equal = (scaled_a == scaled_b)
  END FUNCTION mpf_value_equal

  ! Adds two MultiPrecisionFloat numbers, aligning their exponents.
  FUNCTION mpf_add(a, b) RESULT(mpf_sum)
    TYPE(mpf), INTENT(IN) :: a, b
    TYPE(mpf)             :: mpf_sum
    TYPE(mpi)                            :: scaled_a_mantissa, scaled_b_mantissa, result_mantissa
    INTEGER                                :: common_exponent
    INTEGER                                :: diff_a, diff_b

    IF (mpf_is_zero(a)) THEN
      mpf_sum = b
      RETURN
    END IF
    IF (mpf_is_zero(b)) THEN
      mpf_sum = a
      RETURN
    END IF

    common_exponent = MIN(a%exponent, b%exponent)
    diff_a = a%exponent - common_exponent ! bit shift for a
    diff_b = b%exponent - common_exponent ! bit shift for b

    ! Scale mantissas by shifting bits to align exponents.
    scaled_a_mantissa = mpi_shift_bits_left(a%mantissa, diff_a)
    scaled_b_mantissa = mpi_shift_bits_left(b%mantissa, diff_b)

    ! Perform HighPrecisionInt addition.
    result_mantissa = scaled_a_mantissa + scaled_b_mantissa
    
    mpf_sum = new_mpf_from_mpi_exp(result_mantissa, common_exponent)
  END FUNCTION mpf_add


  ! Scales a MultiPrecisionFloat by 2^power.
  FUNCTION mpf_scale_up_by_base_power(mpf_in, power) RESULT(mpf_out)
    TYPE(mpf), INTENT(IN) :: mpf_in
    INTEGER, INTENT(IN)                  :: power
    TYPE(mpf)             :: mpf_out

    ! Scaling a float by a power of 2 is equivalent to
    ! adding to its exponent. The result is then normalized by the constructor.
    mpf_out = new_mpf_from_mpi_exp(mpf_in%mantissa, mpf_in%exponent + power * 32)
  END FUNCTION mpf_scale_up_by_base_power


  ! Unary negation operator for MultiPrecisionFloat.
  FUNCTION mpf_unary_negate(mpf_in) RESULT(mpf_out)
    TYPE(mpf), INTENT(IN) :: mpf_in
    TYPE(mpf)             :: mpf_out

    mpf_out = new_mpf_from_mpi_exp(-mpf_in%mantissa, mpf_in%exponent)
  END FUNCTION mpf_unary_negate


  ! Subtraction operator for MultiPrecisionFloat, implemented as a + (-b).
  FUNCTION mpf_subtract(a, b) RESULT(mpf_diff)
    TYPE(mpf), INTENT(IN) :: a, b
    TYPE(mpf)             :: mpf_diff
    mpf_diff = mpf_add(a, mpf_unary_negate(b))
  END FUNCTION mpf_subtract


  ! Multiplies two MultiPrecisionFloat numbers.
  FUNCTION mpf_multiply(a, b) RESULT(mpf_prod)
    TYPE(mpf), INTENT(IN) :: a, b
    TYPE(mpf)             :: mpf_prod
    TYPE(mpi)                            :: result_mantissa
    INTEGER                                :: result_exponent
    TYPE(mpi)                            :: zero_mpi
    IF (mpf_is_zero(a) .OR. mpf_is_zero(b)) THEN
      zero_mpi = new_mpi_from_integer(0_8)
      mpf_prod = new_mpf_from_mpi_exp(zero_mpi, 0)
      RETURN
    END IF

    result_mantissa = a%mantissa * b%mantissa
    result_exponent = a%exponent + b%exponent
    
    mpf_prod = new_mpf_from_mpi_exp(result_mantissa, result_exponent)
  END FUNCTION mpf_multiply

  FUNCTION mpf_divide(a, b) RESULT(mpf_quotient)
    TYPE(mpf), INTENT(IN) :: a, b
    TYPE(mpf)             :: mpf_quotient
    TYPE(mpi)                            :: scaled_a_mantissa, result_mantissa, remainder_ignored
    INTEGER                                :: result_exponent, scale_bits
    INTEGER                                :: i, bit_pos
    INTEGER(KIND=8)                        :: temp_val, MASK32

    MASK32 = INT(Z'FFFFFFFF', KIND=8)

    IF (mpf_is_zero(b)) THEN
      STOP "Error: Division by zero MultiPrecisionFloat"
    END IF
    IF (mpf_is_zero(a)) THEN
      mpf_quotient = new_mpf_from_integer(0_8)
      RETURN
    END IF

    temp_val = 0
    bit_pos = -1
    DO i = COEFFS_LIMIT, 1, -1
      temp_val = IAND(ABS(a%mantissa%coeffs(i)), MASK32)
      IF (temp_val /= 0_8) THEN
        ! Find position of MSB within this 32-bit word
        ! Use a simple loop since LEADZ is problematic
        DO bit_pos = 31, 0, -1
          IF (BTEST(temp_val, bit_pos)) EXIT
        END DO
        bit_pos = (i - 1) * 32 + bit_pos
        EXIT
      END IF
    END DO

    ! Calculate how much to scale up
    ! Target: shift so MSB is near the top
    scale_bits = MAX(0, (COEFFS_LIMIT * 32 - 8) - bit_pos)

    ! Scale the numerator
    scaled_a_mantissa = mpi_shift_bits_left(a%mantissa, scale_bits)

    ! Perform division
    CALL mpi_div_rem(scaled_a_mantissa, b%mantissa, result_mantissa, remainder_ignored)

    ! Compute result exponent
    result_exponent = a%exponent - b%exponent - scale_bits

    mpf_quotient = new_mpf_from_mpi_exp(result_mantissa, result_exponent)

  END FUNCTION mpf_divide

  !> @brief Checks if a high-precision float is zero.
  ELEMENTAL FUNCTION mpf_is_zero(val) RESULT(is_zero)
    TYPE(mpf), INTENT(IN) :: val
    LOGICAL :: is_zero
    is_zero = mpi_is_zero(val%mantissa)
  END FUNCTION mpf_is_zero

FUNCTION new_mpf_from_string(s_in) RESULT(mpf_out)
  CHARACTER(LEN=*), INTENT(IN) :: s_in
  TYPE(mpf)   :: mpf_out
  CHARACTER(LEN=:), ALLOCATABLE :: cmd
  CHARACTER(LEN=2048) :: buffer
  CHARACTER(LEN=512) :: temp_file
  CHARACTER(LEN=10) :: coeffs_limit_str
  INTEGER :: stat, unit_num, io_stat
  REAL :: r
  INTEGER(KIND=8) :: c(4)
  INTEGER :: e, sign_flag

  ! Use standard-compliant random number generation for temporary file names
  CALL RANDOM_NUMBER(r)
  WRITE(temp_file, '(A,I0,A)') 'mpf_temp_', INT(r*1000000), '.txt'
  WRITE(coeffs_limit_str, '(I0)') COEFFS_LIMIT
  cmd = './mpf_from_string.sh "' // TRIM(s_in) // '" "' // TRIM(temp_file) // '" ' // TRIM(coeffs_limit_str)
  
  CALL EXECUTE_COMMAND_LINE(TRIM(cmd), WAIT=.TRUE., EXITSTAT=stat)
  
  IF (stat /= 0) THEN
    PRINT *, "ERROR: Julia command failed with exit status:", stat
    mpf_out%mantissa = new_mpi_from_coeffs([INT(0,8), INT(0,8), &
                                            INT(0,8), INT(0,8)])
    mpf_out%exponent = 0
    RETURN
  END IF
  
  ! Read the output from file
  OPEN(NEWUNIT=unit_num, FILE=TRIM(temp_file), STATUS='OLD', ACTION='READ', IOSTAT=io_stat)
  
  IF (io_stat /= 0) THEN
    PRINT *, "ERROR: Could not open temporary file:", TRIM(temp_file)
    mpf_out%mantissa = new_mpi_from_coeffs([INT(0,8), INT(0,8), INT(0,8), INT(0,8)])
    mpf_out%exponent = 0
    RETURN
  END IF
  
  READ(unit_num, '(A)', IOSTAT=io_stat) buffer
  CLOSE(unit_num, STATUS='DELETE')
  
  IF (io_stat /= 0) THEN
    PRINT *, "ERROR: Could not read from temporary file"
    mpf_out%mantissa = new_mpi_from_coeffs([INT(0,8), INT(0,8), &
                                            INT(0,8), INT(0,8)])
    mpf_out%exponent = 0
    RETURN
  END IF
  
  ! Parse comma-separated values: c0,c1,c2,c3,exponent,sign_flag
  READ(buffer, *, IOSTAT=io_stat) c(1), c(2), c(3), c(4), e, sign_flag
  
  IF (io_stat /= 0) THEN
    PRINT *, "ERROR: Could not parse output:", TRIM(buffer)
    mpf_out%mantissa = new_mpi_from_coeffs([INT(0,8), INT(0,8), &
                                            INT(0,8), INT(0,8)])
    mpf_out%exponent = 0
    RETURN
  END IF
  
  mpf_out%mantissa = new_mpi_from_coeffs(c)
  mpf_out%exponent = e
  
  IF (sign_flag == 1) THEN
    mpf_out%mantissa = -mpf_out%mantissa 
  END IF

  CALL normalize_mpf_float(mpf_out)

END FUNCTION new_mpf_from_string

FUNCTION mpf_to_string(mpf_val) RESULT(str)
  TYPE(mpf), INTENT(IN) :: mpf_val
  CHARACTER(LEN=:), ALLOCATABLE :: str, cmd
  CHARACTER(LEN=2048) :: buffer
  CHARACTER(LEN=512) :: temp_file
  CHARACTER(LEN=32) :: coeff_str(4), exp_str, sign_str
  INTEGER :: stat, unit_num, io_stat, i
  REAL :: r
  INTEGER(KIND=8) :: c(4)
  INTEGER :: sign_flag
  TYPE(mpf) :: mpf_work
  
  mpf_work = mpf_val
  sign_flag = 0
  
  IF (mpi_sign(mpf_val%mantissa)) THEN
    sign_flag = 1
    mpf_work%mantissa = mpi_abs(mpf_val%mantissa)
  END IF
  
  IF (mpf_is_zero(mpf_work)) THEN
      str = "0.0"
      RETURN
  END IF

  c = mpf_work%mantissa%coeffs(1:4)
  
  DO i = 1, 4
    WRITE(coeff_str(i), '(I0)') c(i)
  END DO
  WRITE(exp_str, '(I0)') mpf_work%exponent
  WRITE(sign_str, '(I0)') sign_flag
  
  CALL RANDOM_NUMBER(r)
  WRITE(temp_file, '(A,I0,A)') 'mpf_temp_', INT(r*1000000), '.txt'

  cmd = './mpf_to_string.sh ' // &
      TRIM(coeff_str(1)) // ' ' // TRIM(coeff_str(2)) // ' ' // &
      TRIM(coeff_str(3)) // ' ' // TRIM(coeff_str(4)) // ' ' // &
      TRIM(exp_str)      // ' ' // TRIM(sign_str)     // ' "' // TRIM(temp_file) // '"'
  
  CALL EXECUTE_COMMAND_LINE(TRIM(cmd), WAIT=.TRUE., EXITSTAT=stat)  
  
  IF (stat /= 0) THEN
    str = "ERROR: Julia conversion failed"
    RETURN
  END IF
  
  OPEN(NEWUNIT=unit_num, FILE=TRIM(temp_file), STATUS='OLD', ACTION='READ', IOSTAT=io_stat)
  
  IF (io_stat /= 0) THEN
    str = "ERROR: Could not open output file"
    RETURN
  END IF
  
  READ(unit_num, '(A)', IOSTAT=io_stat) buffer
  CLOSE(unit_num, STATUS='DELETE') 
  
  IF (io_stat /= 0) THEN
    str = "ERROR: Could not read output"
    RETURN
  END IF
  
  str = TRIM(buffer)
  
END FUNCTION mpf_to_string

END MODULE multi_precision_float_mod