MODULE high_precision_float_mod
  USE high_precision_integer_mod, ONLY : high_precision_int, new_hpi_from_coeffs, &
                                     normalize_hpi, hpi_abs, new_hpi_from_integer, &
                                     hpi_to_string, hpi_equal, hpi_less,  &
                                     hpi_add, hpi_subtract, hpi_multiply, &
                                     hpi_unary_negate, HIGH_PRECISION_BASE, new_hpi_from_string, &
                                     hpi_is_zero, hpi_power, hpi_divide, hpi_div_rem, &
                                     hpi_scale_up_by_base_power
  IMPLICIT NONE

  ! Precision for division operations. This determines how many extra digits (in base B)
  ! are calculated. A larger value increases precision at the cost of performance.
  INTEGER, PARAMETER :: DEFAULT_DIVISION_PRECISION = 50
  ! Number of decimal places to generate when converting a float to a string.
  INTEGER, PARAMETER :: TO_STRING_DECIMAL_PLACES = 50

  ! The main derived type for a high-precision float.
  TYPE high_precision_float
    TYPE(high_precision_int) :: mantissa
    INTEGER                  :: exponent
  END TYPE high_precision_float

  ! Interface for operator overloading
  INTERFACE OPERATOR(==)
    MODULE PROCEDURE hpf_equal
  END INTERFACE

  INTERFACE OPERATOR(<)
    MODULE PROCEDURE hpf_less
  END INTERFACE

  INTERFACE OPERATOR(+)
    MODULE PROCEDURE hpf_add
  END INTERFACE

  INTERFACE OPERATOR(-)
    MODULE PROCEDURE hpf_subtract
    MODULE PROCEDURE hpf_unary_negate
  END INTERFACE

  INTERFACE OPERATOR(*)
    MODULE PROCEDURE hpf_multiply
  END INTERFACE

  INTERFACE OPERATOR(/)
    MODULE PROCEDURE hpf_divide
  END INTERFACE

  PRIVATE
  PUBLIC :: high_precision_float, OPERATOR(==), OPERATOR(<), OPERATOR(+), &
            OPERATOR(-), OPERATOR(*), OPERATOR(/), &
            new_hpf_from_hpi_exp, new_hpf_from_hpf, &
            normalize_hpf_float, new_hpf_from_integer, new_hpf_from_real64, DEFAULT_DIVISION_PRECISION, &
            hpf_abs, &
            print_hpf_float, hpf_scale_up_by_base_power, new_hpf_from_string, hpf_to_string

CONTAINS

  ! Inner constructor
  FUNCTION new_hpf_from_hpi_exp(mantissa_in, exponent_in) RESULT(hpf_out)
    TYPE(high_precision_int), INTENT(IN) :: mantissa_in
    INTEGER, OPTIONAL, INTENT(IN)        :: exponent_in
    TYPE(high_precision_float)           :: hpf_out

    hpf_out%mantissa = mantissa_in
    IF (PRESENT(exponent_in)) THEN
      hpf_out%exponent = exponent_in
    ELSE
      hpf_out%exponent = 0
    END IF
    CALL normalize_hpf_float(hpf_out)
  END FUNCTION new_hpf_from_hpi_exp

  ! Copy constructor
  FUNCTION new_hpf_from_hpf(hpf_in) RESULT(hpf_out)
    TYPE(high_precision_float), INTENT(IN) :: hpf_in
    TYPE(high_precision_float)             :: hpf_out
    hpf_out%mantissa = hpf_in%mantissa
    hpf_out%exponent = hpf_in%exponent
  END FUNCTION new_hpf_from_hpf

  SUBROUTINE normalize_hpf_float(hpf)
    TYPE(high_precision_float), INTENT(INOUT) :: hpf

    ! 1. Normalize Mantissa: Ensures the HighPrecisionInt mantissa is canonical.
    CALL normalize_hpi(hpf%mantissa)

    ! If mantissa is zero, the exponent must also be zero for canonical form
    IF (hpf%mantissa%sign == 0_1) THEN
      hpf%exponent = 0
      RETURN
    END IF

    ! Adjust exponent by removing least significant zero coefficients from mantissa.
    ! This ensures hpf.mantissa.coeffs(1) is non-zero unless the mantissa is [0].
    DO WHILE (SIZE(hpf%mantissa%coeffs) > 1 .AND. hpf%mantissa%coeffs(1) == 0_8)
      ! Julia's popfirst! removes the first element and shifts others.
      ! Fortran: reallocate or shift elements.
      hpf%mantissa%coeffs = hpf%mantissa%coeffs(2:SIZE(hpf%mantissa%coeffs))
      hpf%exponent = hpf%exponent + 1 ! Compensate for removing B^0 term, effectively multiplying by B
    END DO

  END SUBROUTINE normalize_hpf_float

  ! Creates a HighPrecisionFloat from an Integer.
  FUNCTION new_hpf_from_integer(x_in) RESULT(hpf_out)
    INTEGER(KIND=8), INTENT(IN) :: x_in
    TYPE(high_precision_float)  :: hpf_out
    TYPE(high_precision_int)    :: hpi_val

    hpi_val = new_hpi_from_integer(x_in)
    hpf_out = new_hpf_from_hpi_exp(hpi_val, 0) ! Exponent is 0 for integers
  END FUNCTION new_hpf_from_integer

  ! Returns the absolute value of a HighPrecisionFloat.
  FUNCTION hpf_abs(hpf_in) RESULT(hpf_out)
    TYPE(high_precision_float), INTENT(IN) :: hpf_in
    TYPE(high_precision_float)             :: hpf_out
    hpf_out = new_hpf_from_hpi_exp(hpi_abs(hpf_in%mantissa), hpf_in%exponent)
  END FUNCTION hpf_abs

  ! Checks if two HighPrecisionFloat numbers are equal (a == b)
  FUNCTION hpf_equal(a, b) RESULT(res)
    TYPE(high_precision_float), INTENT(IN) :: a, b
    LOGICAL                                :: res
    res = hpi_equal(a%mantissa, b%mantissa) .AND. (a%exponent == b%exponent)
  END FUNCTION hpf_equal

  ! Compares two HighPrecisionFloat numbers for less than (a < b).
  FUNCTION hpf_less(a, b) RESULT(res)
    TYPE(high_precision_float), INTENT(IN) :: a, b
    LOGICAL                                :: res
    TYPE(high_precision_int)               :: scaled_a_mantissa, scaled_b_mantissa
    INTEGER                                :: common_exponent, diff_a, diff_b

    IF (a%mantissa%sign /= b%mantissa%sign) THEN
      res = (a%mantissa%sign < b%mantissa%sign)
      RETURN
    ELSE IF (a%mantissa%sign == 0_1) THEN
      res = .FALSE. ! Both are zero, so not less than
      RETURN
    END IF

    ! To correctly compare values, we must align their exponents first.
    common_exponent = MIN(a%exponent, b%exponent)
    diff_a = a%exponent - common_exponent
    diff_b = b%exponent - common_exponent

    scaled_a_mantissa = hpi_scale_up_by_base_power(a%mantissa, diff_a)
    scaled_b_mantissa = hpi_scale_up_by_base_power(b%mantissa, diff_b)

    res = hpi_less(scaled_a_mantissa, scaled_b_mantissa)
  END FUNCTION hpf_less

  ! Adds two HighPrecisionFloat numbers, aligning their exponents.
  FUNCTION hpf_add(a, b) RESULT(hpf_sum)
    TYPE(high_precision_float), INTENT(IN) :: a, b
    TYPE(high_precision_float)             :: hpf_sum
    TYPE(high_precision_int)               :: scaled_a_mantissa, scaled_b_mantissa, result_mantissa
    INTEGER                                :: common_exponent
    INTEGER                                :: diff_a, diff_b

    IF (a%mantissa%sign == 0_1) THEN
      hpf_sum = new_hpf_from_hpf(b)
      RETURN
    END IF
    IF (b%mantissa%sign == 0_1) THEN
      hpf_sum = new_hpf_from_hpf(a)
      RETURN
    END IF

    common_exponent = MIN(a%exponent, b%exponent)
    diff_a = a%exponent - common_exponent
    diff_b = b%exponent - common_exponent

    ! Scale mantissas using `hpf_scale_up_by_base_power` to align exponents.
    scaled_a_mantissa = hpi_scale_up_by_base_power(a%mantissa, diff_a)
    scaled_b_mantissa = hpi_scale_up_by_base_power(b%mantissa, diff_b)

    ! Perform HighPrecisionInt addition.
    result_mantissa = hpi_add(scaled_a_mantissa, scaled_b_mantissa)
    
    hpf_sum = new_hpf_from_hpi_exp(result_mantissa, common_exponent)
  END FUNCTION hpf_add

  ! Scales a HighPrecisionFloat by HIGH_PRECISION_BASE^power.
  FUNCTION hpf_scale_up_by_base_power(hpf_in, power) RESULT(hpf_out)
    TYPE(high_precision_float), INTENT(IN) :: hpf_in
    INTEGER, INTENT(IN)                  :: power
    TYPE(high_precision_float)             :: hpf_out

    ! Scaling a float by a power of the internal base is equivalent to
    ! adding to its exponent. The result is then normalized by the constructor.
    hpf_out = new_hpf_from_hpi_exp(hpf_in%mantissa, hpf_in%exponent + power)
  END FUNCTION hpf_scale_up_by_base_power

  ! Unary negation operator for HighPrecisionFloat.
  FUNCTION hpf_unary_negate(hpf_in) RESULT(hpf_out)
    TYPE(high_precision_float), INTENT(IN) :: hpf_in
    TYPE(high_precision_float)             :: hpf_out

    hpf_out = new_hpf_from_hpi_exp(hpi_unary_negate(hpf_in%mantissa), hpf_in%exponent)
  END FUNCTION hpf_unary_negate

  ! Subtraction operator for HighPrecisionFloat, implemented as a + (-b).
  FUNCTION hpf_subtract(a, b) RESULT(hpf_diff)
    TYPE(high_precision_float), INTENT(IN) :: a, b
    TYPE(high_precision_float)             :: hpf_diff
    hpf_diff = hpf_add(a, hpf_unary_negate(b))
  END FUNCTION hpf_subtract

  ! Multiplies two HighPrecisionFloat numbers.
  FUNCTION hpf_multiply(a, b) RESULT(hpf_prod)
    TYPE(high_precision_float), INTENT(IN) :: a, b
    TYPE(high_precision_float)             :: hpf_prod
    TYPE(high_precision_int)               :: result_mantissa
    INTEGER                                :: result_exponent

    IF (a%mantissa%sign == 0_1 .OR. b%mantissa%sign == 0_1) THEN
      hpf_prod = new_hpf_from_hpi_exp(new_hpi_from_integer(0_8), 0)
      RETURN
    END IF

    result_mantissa = hpi_multiply(a%mantissa, b%mantissa)
    result_exponent = a%exponent + b%exponent
    
    hpf_prod = new_hpf_from_hpi_exp(result_mantissa, result_exponent)
  END FUNCTION hpf_multiply

  ! Divides two HighPrecisionFloat numbers using high-precision integer division.
  FUNCTION hpf_divide(a, b) RESULT(hpf_quotient)
    TYPE(high_precision_float), INTENT(IN) :: a, b
    TYPE(high_precision_float)             :: hpf_quotient
    TYPE(high_precision_int)               :: scaled_a_mantissa, result_mantissa
    INTEGER                                :: result_exponent

    IF (hpi_is_zero(b%mantissa)) THEN
      STOP "Error: Division by zero HighPrecisionFloat"
    END IF
    IF (hpi_is_zero(a%mantissa)) THEN
      hpf_quotient = new_hpf_from_integer(0_8)
      RETURN
    END IF

    scaled_a_mantissa = hpi_scale_up_by_base_power(a%mantissa, DEFAULT_DIVISION_PRECISION)
    result_mantissa = hpi_divide(scaled_a_mantissa, b%mantissa)
    result_exponent = a%exponent - b%exponent - DEFAULT_DIVISION_PRECISION
    hpf_quotient = new_hpf_from_hpi_exp(result_mantissa, result_exponent)
  END FUNCTION hpf_divide

  ! Displays a HighPrecisionFloat by converting it to BigFloat for decimal
  ! representation, and showing its internal mantissa and exponent.
  SUBROUTINE print_hpf_float(io, hpf)
    CLASS(*), INTENT(IN)               :: io ! Generic IO unit
    TYPE(high_precision_float), INTENT(IN) :: hpf

    CHARACTER(LEN=:), ALLOCATABLE      :: value_str
    CHARACTER(LEN=256)                 :: mantissa_coeffs_str
    CHARACTER(LEN=20)                  :: coeff_val_str
    INTEGER                            :: i
    
    SELECT TYPE (io)
    TYPE IS (INTEGER)
        IF (hpi_is_zero(hpf%mantissa)) THEN
          WRITE(io, '(A)') "HighPrecisionFloat(0.0, mantissa=HighPrecisionInt(0, coeffs=[0]), exponent=0)"
        ELSE
          value_str = hpf_to_string(hpf)
    
          ! Build the coefficient string representation
          mantissa_coeffs_str = "["
          DO i = 1, SIZE(hpf%mantissa%coeffs)
            WRITE(coeff_val_str, '(I0)') hpf%mantissa%coeffs(i); mantissa_coeffs_str = TRIM(mantissa_coeffs_str) // TRIM(ADJUSTL(coeff_val_str))
            IF (i < SIZE(hpf%mantissa%coeffs)) THEN
              mantissa_coeffs_str = TRIM(mantissa_coeffs_str) // ", "
            END IF
          END DO
          mantissa_coeffs_str = TRIM(mantissa_coeffs_str) // "]"
    
          WRITE(io, '(A, A, A, A, A, A, A, I0, A)') &
            "HighPrecisionFloat(", TRIM(value_str), &
            ", mantissa=HighPrecisionInt(", hpi_to_string(hpf%mantissa), &
            ", coeffs=", TRIM(mantissa_coeffs_str), &
            "), exponent=", hpf%exponent, ")"
        END IF
    CLASS DEFAULT
        ! Handle cases where io is not an integer. For this library's purpose,
        ! we can assume it's an error.
        STOP "Error in print_hpf_float: I/O unit must be an INTEGER."
    END SELECT
  END SUBROUTINE print_hpf_float

  ! Converts a HighPrecisionFloat to a decimal string representation.
  FUNCTION hpf_to_string(hpf) RESULT(str_out)
    TYPE(high_precision_float), INTENT(IN) :: hpf
    CHARACTER(LEN=:), ALLOCATABLE          :: str_out
    TYPE(high_precision_float)             :: abs_hpf
    TYPE(high_precision_int)               :: num, den, q, r, ten
    INTEGER                                :: i ! Loop index for decimal places

    IF (hpi_is_zero(hpf%mantissa)) THEN
      str_out = "0.0"
      RETURN
    END IF

    ! Work with the absolute value first, then apply sign at the end.
    abs_hpf = hpf_abs(hpf)

    ! Unified logic: always treat the number as a fraction num/den
    ! The mantissa is the numerator. The denominator is HIGH_PRECISION_BASE raised
    ! to the power of the negative exponent.
    IF (abs_hpf%exponent <= 0) THEN
      ! If exponent is 0 or negative, the number is mantissa / B^(-exponent)
      num = abs_hpf%mantissa
      den = hpi_power(new_hpi_from_integer(HIGH_PRECISION_BASE), -abs_hpf%exponent)
    ELSE
      ! If exponent is positive, the number is mantissa * B^(exponent) / 1
      num = hpi_multiply(abs_hpf%mantissa, hpi_power(new_hpi_from_integer(HIGH_PRECISION_BASE), abs_hpf%exponent))
      den = new_hpi_from_integer(1_8)
    END IF

    ! Convert the fraction num/den to a decimal string
    CALL hpi_div_rem(num, den, q, r)
    str_out = hpi_to_string(q)

    IF (hpi_is_zero(r) .AND. TO_STRING_DECIMAL_PLACES > 0) THEN
      ! No remainder, but decimal places requested. Add ".0" for consistency.
      str_out = TRIM(str_out) // ".0"
    ELSE IF (.NOT. hpi_is_zero(r)) THEN ! There is a remainder, so we must generate the fractional part.
      str_out = TRIM(str_out) // "."
      ten = new_hpi_from_integer(10_8)
      DO i = 1, TO_STRING_DECIMAL_PLACES
        IF (hpi_is_zero(r)) EXIT
        CALL hpi_div_rem(hpi_multiply(r, ten), den, q, r)
        str_out = TRIM(str_out) // hpi_to_string(q) ! Append the digit
      END DO
    END IF

    IF (hpf%mantissa%sign < 0) THEN
      str_out = "-" // TRIM(str_out)
    END IF
  END FUNCTION hpf_to_string

  ! Creates a HighPrecisionFloat from a REAL(KIND=8) by converting to a string first.
  FUNCTION new_hpf_from_real64(x_in) RESULT(hpf_out)
    REAL(KIND=8), INTENT(IN)   :: x_in
    TYPE(high_precision_float) :: hpf_out
    CHARACTER(LEN=100)         :: s_val

    WRITE(s_val, '(ES34.25E3)') x_in
    hpf_out = new_hpf_from_string(TRIM(s_val))
  END FUNCTION new_hpf_from_real64

  ! Creates a HighPrecisionFloat from a decimal string (e.g., "-123.45e-6").
  FUNCTION new_hpf_from_string(s_in) RESULT(hpf_out)
    CHARACTER(LEN=*), INTENT(IN) :: s_in
    TYPE(high_precision_float)   :: hpf_out
    CHARACTER(LEN=LEN(s_in))     :: temp_s, num_part, exp_part, int_str
    INTEGER                      :: e_pos, dot_pos, num_frac_digits, exp_val, den_power, i
    TYPE(high_precision_int)     :: num_hpi, den_hpi
    LOGICAL                      :: is_negative = .FALSE.


    is_negative = .FALSE. ! Explicit assignment for robust initialization
    temp_s = TRIM(ADJUSTL(s_in))
    IF (LEN_TRIM(temp_s) == 0) THEN
      hpf_out = new_hpf_from_integer(0_8); RETURN
    END IF

    IF (temp_s(1:1) == '-') THEN
      is_negative = .TRUE.
      temp_s = temp_s(2:)
    ELSE IF (temp_s(1:1) == '+') THEN
      temp_s = temp_s(2:)
    END IF

    e_pos = INDEX(temp_s, 'e', .TRUE.) ! Case-insensitive search for 'e'
    IF (e_pos == 0) e_pos = INDEX(temp_s, 'E', .TRUE.)

    IF (e_pos > 0) THEN
      num_part = temp_s(:e_pos-1)
      exp_part = temp_s(e_pos+1:)
      READ(exp_part, *, IOSTAT=i) exp_val
      IF (i /= 0) THEN
         STOP "Invalid exponent in string"
      END IF
    ELSE
      num_part = temp_s
      exp_val = 0
    END IF

    dot_pos = INDEX(num_part, '.')
    IF (dot_pos > 0) THEN
      num_frac_digits = LEN_TRIM(num_part) - dot_pos
      int_str = num_part(:dot_pos-1) // num_part(dot_pos+1:)
    ELSE
      num_frac_digits = 0
      int_str = num_part
    END IF

    IF (is_negative) THEN
      int_str = "-" // TRIM(int_str)
    END IF

    num_hpi = new_hpi_from_string(int_str)
    IF (hpi_is_zero(num_hpi)) THEN
      hpf_out = new_hpf_from_integer(0_8); RETURN
    END IF

    den_power = num_frac_digits - exp_val

    IF (den_power == 0) THEN
      hpf_out = new_hpf_from_hpi_exp(num_hpi)
    ELSE IF (den_power > 0) THEN
      den_hpi = hpi_power(new_hpi_from_integer(10_8), den_power)
      hpf_out = hpf_divide(new_hpf_from_hpi_exp(num_hpi), new_hpf_from_hpi_exp(den_hpi))
    ELSE ! den_power < 0
      num_hpi = hpi_multiply(num_hpi, hpi_power(new_hpi_from_integer(10_8), -den_power))
      hpf_out = new_hpf_from_hpi_exp(num_hpi)
    END IF

  END FUNCTION new_hpf_from_string

END MODULE high_precision_float_mod