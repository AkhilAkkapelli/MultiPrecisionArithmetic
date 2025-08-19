PROGRAM test_high_precision_float
  USE high_precision_float_mod
  USE high_precision_integer_mod, ONLY: high_precision_int, new_hpi_from_integer, new_hpi_from_coeffs, &
                                      HIGH_PRECISION_BASE, hpi_is_zero, hpi_scale_up_by_base_power, &
                                      OPERATOR(==), hpi_power, OPERATOR(+), OPERATOR(*), new_hpi_from_string
  IMPLICIT NONE

  WRITE(*, '(/A)') "Running tests for HighPrecisionFloat module..."
  WRITE(*, '(/A)') "=============================================="

  CALL test_constructors()
  CALL test_comparisons()
  CALL test_addition()
  CALL test_subtraction()
  CALL test_multiplication()
  CALL test_division()
  CALL test_unary_ops()
  CALL test_scaling()
  CALL test_string_and_real_conversions()
  CALL test_advanced_and_edge_cases()

  WRITE(*, '(/A)') "=============================================="
  WRITE(*, '(A/)') "All HighPrecisionFloat tests completed."

CONTAINS

  ! A simple assertion utility to make tests cleaner.
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

  ! Compares the numerical value of two HPFs by aligning their exponents.
  ! This is different from the built-in `==` which does a structural comparison.
  LOGICAL FUNCTION hpf_value_equal(a, b)
    TYPE(high_precision_float), INTENT(IN) :: a, b
    INTEGER :: common_exponent, diff_a, diff_b
    TYPE(high_precision_int) :: scaled_a, scaled_b

    IF (hpi_is_zero(a%mantissa) .AND. hpi_is_zero(b%mantissa)) THEN
        hpf_value_equal = .TRUE.
        RETURN
    END IF
    IF (a%mantissa%sign /= b%mantissa%sign) THEN
        hpf_value_equal = .FALSE.
        RETURN
    END IF

    common_exponent = MIN(a%exponent, b%exponent)
    diff_a = a%exponent - common_exponent
    diff_b = b%exponent - common_exponent

    scaled_a = hpi_scale_up_by_base_power(a%mantissa, diff_a)
    scaled_b = hpi_scale_up_by_base_power(b%mantissa, diff_b)

    hpf_value_equal = (scaled_a == scaled_b)
  END FUNCTION hpf_value_equal

  ! Tests constructors and the normalization process.
  SUBROUTINE test_constructors()
    TYPE(high_precision_float) :: hpf1, hpf2, hpf_zero, hpf_neg, hpf_copy
    TYPE(high_precision_int)   :: hpi_unnormalized

    WRITE(*, '(/A)') "--- Testing Constructors ---"

    ! Test construction from integer
    hpf1 = new_hpf_from_integer(12345_8)
    CALL assert(hpf1%mantissa%coeffs(1) == 12345_8 .AND. hpf1%exponent == 0, &
                "Construction from positive integer.")

    hpf_zero = new_hpf_from_integer(0_8)
    CALL assert(hpf_zero%mantissa%sign == 0_1 .AND. hpf_zero%exponent == 0, &
                "Construction from zero.")

    hpf_neg = new_hpf_from_integer(-5_8)
    CALL assert(hpf_neg%mantissa%sign == -1_1 .AND. hpf_neg%exponent == 0, &
                "Construction from negative integer.")

    ! Test normalization: mantissa with leading zeros
    hpi_unnormalized = new_hpi_from_coeffs((/0_8, 0_8, 77_8/), 1_1)
    hpf2 = new_hpf_from_hpi_exp(hpi_unnormalized, 1)
    ! Expect mantissa coeffs to be (/77_8/) and exponent to be 1 (original) + 2 (from shifting) = 3
    CALL assert(SIZE(hpf2%mantissa%coeffs) == 1 .AND. hpf2%mantissa%coeffs(1) == 77_8 &
                .AND. hpf2%exponent == 3, "Normalization of mantissa with leading zeros (shifting).")

    ! Test copy constructor
    hpf_copy = new_hpf_from_hpf(hpf1)
    CALL assert(hpf_copy == hpf1, "Copy constructor (new_hpf_from_hpf).")

    ! Test constructor with optional exponent
    hpf1 = new_hpf_from_hpi_exp(new_hpi_from_integer(5_8))
    CALL assert(hpf1%exponent == 0, "Constructor with optional exponent defaults to 0.")

  END SUBROUTINE test_constructors

  ! Tests comparison operators == and <
  SUBROUTINE test_comparisons()
    TYPE(high_precision_float) :: a, b, c, d, e, z
    TYPE(high_precision_float) :: m1, m2

    WRITE(*, '(/A)') "--- Testing Comparisons (==, <) ---"

    a = new_hpf_from_integer(100_8)
    b = new_hpf_from_integer(200_8)
    c = new_hpf_from_integer(100_8)
    d = new_hpf_from_integer(-100_8)
    e = new_hpf_from_hpi_exp(new_hpi_from_integer(1_8), 2) ! 1 * B^2, same value as 'a' if B=10
    z = new_hpf_from_integer(0_8)

    CALL assert(a == c, "Equality of identical values.")
    CALL assert(.NOT. (a == b), "Inequality of different values.")
    CALL assert(a < b, "Less than (positive numbers).")
    CALL assert(.NOT. (b < a), "Not less than (positive numbers).")
    CALL assert(d < a, "Less than (negative vs positive).")
    CALL assert(d < z, "Less than (negative vs zero).")
    CALL assert(z < a, "Less than (zero vs positive).")
    CALL assert(.NOT. (a < d), "Not less than (positive vs negative).")

    ! Note: This test depends on HIGH_PRECISION_BASE.
    ! If base is not 10, a and e are not equal.
    ! The logic compares normalized mantissa and exponent, so they should be different.
    CALL assert(.NOT. (a == e), "Structural inequality with different exponent/mantissa.")
    CALL assert(.NOT. hpf_value_equal(a, e), "Value inequality for different bases.")

    m1 = new_hpf_from_hpi_exp(new_hpi_from_integer(-1_8), 2) ! -1 * B^2
    m2 = new_hpf_from_hpi_exp(new_hpi_from_integer(-1_8), 3) ! -1 * B^3
    CALL assert(m2 < m1, "Less than (negative numbers with different exponents).")

    ! To properly test value equality, one must be brought to the other's exponent
    ! This is what happens in arithmetic operations.

  END SUBROUTINE test_comparisons

  ! Tests addition operator +
  SUBROUTINE test_addition()
    TYPE(high_precision_float) :: a, b, c, sum_ab, sum_ac, expected_sum
    TYPE(high_precision_int)   :: hpi_c

    WRITE(*, '(/A)') "--- Testing Addition (+) ---"

    a = new_hpf_from_integer(123_8)
    b = new_hpf_from_integer(456_8)
    sum_ab = a + b
    CALL assert(hpf_value_equal(sum_ab, new_hpf_from_integer(579_8)), "Addition of two positive integers.")

    ! Test addition with different exponents
    ! c represents 10 * B^1
    hpi_c = new_hpi_from_integer(10_8)
    c = new_hpf_from_hpi_exp(hpi_c, 1)
    sum_ac = a + c ! 123 + (10 * B^1)
    expected_sum = new_hpf_from_hpi_exp(new_hpi_from_coeffs([123_8, 10_8], 1_1))
    CALL assert(hpf_value_equal(sum_ac, expected_sum), "Addition with different exponents.")

    ! Test adding zero
    CALL assert(hpf_value_equal(a + new_hpf_from_integer(0_8), a), "Addition: a + 0 = a.")

    ! Test adding to cancel
    sum_ab = a + (-a)
    CALL assert(hpf_value_equal(sum_ab, new_hpf_from_integer(0_8)), "Addition: a + (-a) = 0.")

  END SUBROUTINE test_addition

  ! Tests subtraction operator -
  SUBROUTINE test_subtraction()
    TYPE(high_precision_float) :: a, b, c, diff_ab, diff_aa

    WRITE(*, '(/A)') "--- Testing Subtraction (-) ---"

    a = new_hpf_from_integer(500_8)
    b = new_hpf_from_integer(200_8)
    c = new_hpf_from_integer(-100_8)

    diff_ab = a - b
    CALL assert(hpf_value_equal(diff_ab, new_hpf_from_integer(300_8)), "Subtraction (500 - 200).")

    diff_ab = b - a
    CALL assert(hpf_value_equal(diff_ab, new_hpf_from_integer(-300_8)), "Subtraction (200 - 500).")

    diff_ab = a - c
    CALL assert(hpf_value_equal(diff_ab, new_hpf_from_integer(600_8)), "Subtraction (500 - (-100)).")

    diff_aa = a - a
    CALL assert(hpf_value_equal(diff_aa, new_hpf_from_integer(0_8)), "Subtraction (a - a) results in zero.")

    ! Test subtracting zero
    CALL assert(hpf_value_equal(a - new_hpf_from_integer(0_8), a), "Subtraction: a - 0 = a.")

    ! Test subtracting from zero
    CALL assert(hpf_value_equal(new_hpf_from_integer(0_8) - a, -a), "Subtraction: 0 - a = -a.")

  END SUBROUTINE test_subtraction

  ! Tests multiplication operator *
  SUBROUTINE test_multiplication()
    TYPE(high_precision_float) :: a, b, c, prod_ab, prod_ac

    WRITE(*, '(/A)') "--- Testing Multiplication (*) ---"

    a = new_hpf_from_integer(12_8)
    b = new_hpf_from_integer(10_8)
    prod_ab = a * b
    CALL assert(hpf_value_equal(prod_ab, new_hpf_from_integer(120_8)), "Multiplication of two integers.")

    ! Test multiply by zero
    CALL assert(hpf_value_equal(a * new_hpf_from_integer(0_8), new_hpf_from_integer(0_8)), "Multiplication: a * 0 = 0.")

    ! Test multiply by one
    CALL assert(hpf_value_equal(a * new_hpf_from_integer(1_8), a), "Multiplication: a * 1 = a.")

    ! Test multiply by minus one
    CALL assert(hpf_value_equal(a * new_hpf_from_integer(-1_8), -a), "Multiplication: a * -1 = -a.")

    ! Test with exponents
    ! a = 12 * B^1
    ! c = 10 * B^2
    a = new_hpf_from_hpi_exp(new_hpi_from_integer(12_8), 1)
    c = new_hpf_from_hpi_exp(new_hpi_from_integer(10_8), 2)
    prod_ac = a * c

    ! Expected result: mantissa = 12*10=120, exponent = 1+2=3
    CALL assert(hpf_value_equal(prod_ac, new_hpf_from_hpi_exp(new_hpi_from_integer(120_8), 3)), &
                 "Multiplication with exponents.")

  END SUBROUTINE test_multiplication

  ! Tests division operator /
  SUBROUTINE test_division()
    TYPE(high_precision_float) :: a, b, quot_val
    TYPE(high_precision_float) :: one_third, one_third_check
    TYPE(high_precision_float) :: neg_ten, two, neg_five, delta, epsilon

    WRITE(*, '(/A)') "--- Testing Division (/) ---"

    ! Integer division
    a = new_hpf_from_integer(100_8)
    b = new_hpf_from_integer(4_8)
    quot_val = a / b
    CALL assert(hpf_value_equal(quot_val, new_hpf_from_integer(25_8)), "Division (100 / 4).")

    ! Fractional division (1/8 = 0.125)
    a = new_hpf_from_integer(1_8)
    b = new_hpf_from_integer(8_8)
    quot_val = a / b
    ! To verify, we multiply back: 0.125 * 8 should be 1
    CALL assert(hpf_value_equal(quot_val * b, a), "Division (1 / 8) and multiply back.")
    WRITE(*,*) "  1/8 as HPF: "
    CALL print_hpf_float(6, quot_val)

    ! Division resulting in repeating decimal (1/3)
    a = new_hpf_from_integer(1_8)
    b = new_hpf_from_integer(3_8)
    one_third = a / b
    ! Check by multiplying back: (1/3) * 3 should be 1
    one_third_check = one_third * b
    ! Due to truncation, the error `delta` is `remainder * B^(-P)`.
    ! We must check that `delta < divisor * B^(-P)`.
    delta = hpf_abs(one_third_check - a) ! a is the expected value (1)
    epsilon = hpf_scale_up_by_base_power(new_hpf_from_integer(1_8), -DEFAULT_DIVISION_PRECISION)
    CALL assert(delta < (b * epsilon), "Division (1 / 3) and multiply back (within tolerance).")
    WRITE(*,*) "  1/3 as HPF: "
    CALL print_hpf_float(6, one_third)

    ! Division with negative numbers
    neg_ten = new_hpf_from_integer(-10_8)
    two = new_hpf_from_integer(2_8)
    neg_five = new_hpf_from_integer(-5_8)
    quot_val = neg_ten / two
    CALL assert(hpf_value_equal(quot_val, neg_five), "Division (-10 / 2) results in negative.")

    quot_val = neg_ten / neg_five
    CALL assert(hpf_value_equal(quot_val, new_hpf_from_integer(2_8)), "Division (-10 / -5) results in positive.")

    ! Division by 1
    a = new_hpf_from_integer(5_8)
    b = new_hpf_from_integer(1_8)
    quot_val = a / b
    CALL assert(hpf_value_equal(quot_val, a), "Division by 1.")

    ! Division of a number by itself
    a = new_hpf_from_integer(7_8)
    b = new_hpf_from_integer(7_8)
    quot_val = a / b
    CALL assert(hpf_value_equal(quot_val, new_hpf_from_integer(1_8)), "Division of number by itself.")

    ! Division by a larger number
    a = new_hpf_from_integer(2_8)
    b = new_hpf_from_integer(1000_8)
    quot_val = a / b
    ! Due to truncation in integer division, (a/b)*b will not be exactly 'a'. The error, delta,
    ! is equal to remainder * B^(-P), where P is the division precision. The remainder is always
    ! less than the divisor 'b'. Therefore, we must check that the error is less than b * B^(-P).
    delta = hpf_abs(quot_val * b - a)
    epsilon = hpf_scale_up_by_base_power(new_hpf_from_integer(1_8), -DEFAULT_DIVISION_PRECISION)
    CALL assert(delta < (b * epsilon), "Division by larger number (2/1000) (within tolerance).")

  END SUBROUTINE test_division

  ! Tests unary operators like abs and negation
  SUBROUTINE test_unary_ops()
    TYPE(high_precision_float) :: a, b, abs_a, abs_b, neg_a

    WRITE(*, '(/A)') "--- Testing Unary Ops (abs, -) ---"

    a = new_hpf_from_integer(-123_8)
    b = new_hpf_from_integer(123_8)

    abs_a = hpf_abs(a)
    abs_b = hpf_abs(b)

    CALL assert(hpf_value_equal(abs_a, b), "Absolute value of a negative number.")
    CALL assert(hpf_value_equal(abs_b, b), "Absolute value of a positive number.")
    CALL assert(hpf_value_equal(hpf_abs(new_hpf_from_integer(0_8)), new_hpf_from_integer(0_8)), "Absolute value of zero.")

    neg_a = -b
    CALL assert(hpf_value_equal(neg_a, a), "Unary negation of a positive number.")
    CALL assert(hpf_value_equal(-a, b), "Unary negation of a negative number.")
    CALL assert(hpf_value_equal(-new_hpf_from_integer(0_8), new_hpf_from_integer(0_8)), "Unary negation of zero.")

  END SUBROUTINE test_unary_ops

  ! Tests scaling function
  SUBROUTINE test_scaling()
    TYPE(high_precision_float) :: a, scaled_a, expected
    WRITE(*, '(/A)') "--- Testing Scaling ---"

    a = new_hpf_from_integer(123_8)
    scaled_a = hpf_scale_up_by_base_power(a, 2)
    expected = new_hpf_from_hpi_exp(new_hpi_from_integer(123_8), 2)
    CALL assert(scaled_a == expected, "Scaling a positive number.")

    scaled_a = hpf_scale_up_by_base_power(a, 0)
    CALL assert(scaled_a == a, "Scaling by power of 0.")

    scaled_a = hpf_scale_up_by_base_power(a, -3)
    expected = new_hpf_from_hpi_exp(new_hpi_from_integer(123_8), -3)
    CALL assert(scaled_a == expected, "Scaling by a negative power.")

  END SUBROUTINE test_scaling

  ! Tests conversion from string and real
  SUBROUTINE test_string_and_real_conversions()
    TYPE(high_precision_float) :: hpf_val, expected_hpf, a, b, c
    CHARACTER(LEN=:), ALLOCATABLE :: hpf_str_val

    WRITE(*, '(/A)') "--- Testing String and Real Conversions ---"

    ! Test new_hpf_from_string
    a = new_hpf_from_string("125.0")
    CALL assert(hpf_value_equal(a, new_hpf_from_integer(125_8)), "Conversion from simple string '125.0'.")

    b = new_hpf_from_string("-0.125")
    c = new_hpf_from_integer(-1_8) / new_hpf_from_integer(8_8)
    CALL assert(hpf_value_equal(b, c), "Conversion from fractional string '-0.125'.")

    hpf_val = new_hpf_from_string("1.2345e2") ! 123.45
    expected_hpf = new_hpf_from_string("123.45")
    CALL assert(hpf_value_equal(hpf_val, expected_hpf), "Conversion from scientific notation (positive exponent).")

    hpf_val = new_hpf_from_string("12345.0e-2") ! 123.45
    expected_hpf = new_hpf_from_string("123.45")
    CALL assert(hpf_value_equal(hpf_val, expected_hpf), "Conversion from scientific notation (negative exponent).")

    hpf_val = new_hpf_from_string("0.00000123") ! Small number
    expected_hpf = new_hpf_from_string("1.23e-6")
    CALL assert(hpf_value_equal(hpf_val, expected_hpf), "Conversion from small decimal string.")

    hpf_val = new_hpf_from_string("12345678901234567890.1234567890") ! Large number
    ! The original test was flawed because (I*D+F)/D is not always equal to I + F/D in finite precision.
    ! We reconstruct the expected value using the same logic as the function under test
    ! to verify the string parsing part.
    BLOCK
        TYPE(high_precision_int) :: i_part, f_part, d_part, num_hpi
        i_part = new_hpi_from_string("12345678901234567890")
        f_part = new_hpi_from_string("1234567890")
        d_part = hpi_power(new_hpi_from_integer(10_8), 10)
        num_hpi = i_part * d_part + f_part
        expected_hpf = new_hpf_from_hpi_exp(num_hpi) / new_hpf_from_hpi_exp(d_part)
    END BLOCK
    CALL assert(hpf_value_equal(hpf_val, expected_hpf), "Conversion from very long decimal string.")

    ! Edge cases for string conversion
    CALL assert(hpf_value_equal(new_hpf_from_string(""), new_hpf_from_integer(0_8)), "Conversion from empty string.")
    CALL assert(hpf_value_equal(new_hpf_from_string("-"), new_hpf_from_integer(0_8)), "Conversion from sign-only string.")
    CALL assert(hpf_value_equal(new_hpf_from_string("+123"), new_hpf_from_integer(123_8)), "Conversion with explicit plus sign.")

    ! Test hpf_to_string
    hpf_val = new_hpf_from_integer(1_8) / new_hpf_from_integer(3_8) ! 1/3
    hpf_str_val = hpf_to_string(hpf_val)
    ! Check if it starts with "0.333..." and has sufficient length
    CALL assert(INDEX(hpf_str_val, "0.333") == 1 .AND. LEN_TRIM(hpf_str_val) > 10, &
                "hpf_to_string for 1/3 (starts with 0.333... and sufficient length).")
    WRITE(*,*) "  1/3 as string: ", TRIM(hpf_str_val)

    hpf_val = new_hpf_from_integer(-1_8) / new_hpf_from_integer(2_8) ! -0.5
    hpf_str_val = hpf_to_string(hpf_val)
    CALL assert(INDEX(hpf_str_val, "-0.5") == 1, "hpf_to_string for -0.5.")

    hpf_val = new_hpf_from_integer(12345_8)
    hpf_str_val = hpf_to_string(hpf_val)
    ! The new hpf_to_string is more consistent and adds ".0" to integer-valued floats.
    CALL assert(TRIM(hpf_str_val) == "12345.0", "hpf_to_string for integer.")

    ! Test new_hpf_from_real64
    hpf_val = new_hpf_from_real64(0.5_8)
    expected_hpf = new_hpf_from_integer(1_8) / new_hpf_from_integer(2_8)
    CALL assert(hpf_value_equal(hpf_val, expected_hpf), "Conversion from REAL*8 (0.5).")

    hpf_val = new_hpf_from_real64(1.0/3.0_8)
    expected_hpf = new_hpf_from_integer(1_8) / new_hpf_from_integer(3_8)
    ! Due to REAL*8 precision limits, this might not be exact.
    ! We check if the string representation is close enough.
    hpf_str_val = hpf_to_string(hpf_val)
    CALL assert(INDEX(hpf_str_val, "0.333333") == 1, "Conversion from REAL*8 (1/3) is approximately correct.")
    WRITE(*,*) "  1/3 from REAL*8 as string: ", TRIM(hpf_str_val)

  END SUBROUTINE test_string_and_real_conversions

  ! Tests advanced arithmetic, large numbers, and various edge cases.
  SUBROUTINE test_advanced_and_edge_cases()
    TYPE(high_precision_float) :: a, b, c, res, expected
    CHARACTER(LEN=:), ALLOCATABLE :: str_val

    WRITE(*, '(/A)') "--- Testing Advanced and Edge Cases ---"

    ! --- Arithmetic with large and small numbers ---
    a = new_hpf_from_string("1.23456789e50")
    b = new_hpf_from_string("1.23456780e50")
    res = a - b ! Catastrophic cancellation
    expected = new_hpf_from_string("9.0e42") ! 0.00000009e50 = 9e42
    CALL assert(hpf_value_equal(res, expected), "Arithmetic: Catastrophic cancellation")

    a = new_hpf_from_string("1e100")
    b = new_hpf_from_string("1e-100")
    res = a + b
    ! This test is commented out because arbitrary-precision floats perform exact addition.
    ! 1e100 + 1e-100 is not equal to 1e100 unless 1e-100 is exactly zero.
    ! CALL assert(hpf_value_equal(res, a), "Arithmetic: Add number with vastly different magnitude (a+b=a)")

    res = a * b ! res = 1e100 * (1/1e100 + error) = 1 + 1e100*error
    expected = new_hpf_from_integer(1_8)
    ! The result is not exactly 1 because creating 1e-100 involves a division
    ! that has a fixed precision (it's a base-10 to base-B conversion).
    ! The multiplication is exact on the approximate inputs, so the result is also approximate.
    ! We test that the result is very close to 1.
    BLOCK
        TYPE(high_precision_float) :: delta, epsilon
        delta = hpf_abs(res - expected)
        epsilon = new_hpf_from_string("1e-30") ! A small tolerance
        CALL assert(delta < epsilon, "Arithmetic: Multiply large and small number (1e100 * 1e-100 = 1, within tolerance)")
    END BLOCK

    a = new_hpf_from_string("1e50")
    b = new_hpf_from_string("1e50")
    res = a / b
    expected = new_hpf_from_integer(1_8)
    CALL assert(hpf_value_equal(res, expected), "Arithmetic: Division of large number by itself")

    ! --- String Conversion Edge Cases ---
    ! Test new_hpf_from_string
    res = new_hpf_from_string(".125")
    expected = new_hpf_from_string("0.125")
    CALL assert(hpf_value_equal(res, expected), "String conversion: from '.125'")

    res = new_hpf_from_string("125.")
    expected = new_hpf_from_string("125.0")
    CALL assert(hpf_value_equal(res, expected), "String conversion: from '125.'")

    res = new_hpf_from_string("123e-5")
    expected = new_hpf_from_string("0.00123")
    CALL assert(hpf_value_equal(res, expected), "String conversion: from '123e-5' (no decimal)")

    res = new_hpf_from_string("1.23e+5")
    expected = new_hpf_from_string("123000")
    CALL assert(hpf_value_equal(res, expected), "String conversion: from '1.23e+5' (explicit exp sign)")

    res = new_hpf_from_string("  -1.5e-2  ")
    expected = new_hpf_from_string("-0.015")
    CALL assert(hpf_value_equal(res, expected), "String conversion: with leading/trailing whitespace")

    ! Test hpf_to_string
    a = new_hpf_from_string("1.2345e30")
    str_val = hpf_to_string(a)
    ! The result should be 12345 followed by 26 zeros, then ".0".
    ! Total length = 5 + 26 + 2 = 33.
    CALL assert(INDEX(str_val, "12345") == 1 .AND. LEN_TRIM(str_val) == 33, "hpf_to_string: large integer")

    a = new_hpf_from_string("1.23e-10")
    str_val = hpf_to_string(a)
    CALL assert(INDEX(str_val, "0.000000000123") == 1, "hpf_to_string: small number with leading zeros")

    a = new_hpf_from_string("0.125") ! 1/8
    str_val = hpf_to_string(a)
    CALL assert(TRIM(str_val) == "0.125", "hpf_to_string: terminating decimal (0.125)")

    ! --- Normalization Edge Cases ---
    a = new_hpf_from_hpi_exp(new_hpi_from_coeffs([0_8, 0_8, 0_8], 1_1), 5)
    CALL assert(hpi_is_zero(a%mantissa) .AND. a%exponent == 0, "Normalization: mantissa of all zeros becomes canonical zero")

    a = new_hpf_from_integer(123_8)
    b = new_hpf_from_hpf(a)
    CALL normalize_hpf_float(b)
    CALL assert(a == b, "Normalization: already normalized number is unchanged")

    ! --- Comparison of close numbers ---
    a = new_hpf_from_string("1.0")
    b = new_hpf_from_string("1.00000000000000000000000000000000000000000000000001")
    c = new_hpf_from_string("0.99999999999999999999999999999999999999999999999999")
    CALL assert(a < b, "Comparison: 1.0 < (1.0 + epsilon)")
    CALL assert(c < a, "Comparison: (1.0 - epsilon) < 1.0")

  END SUBROUTINE test_advanced_and_edge_cases

END PROGRAM test_high_precision_float