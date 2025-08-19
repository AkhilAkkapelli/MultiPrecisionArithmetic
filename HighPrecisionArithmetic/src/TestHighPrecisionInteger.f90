PROGRAM test_high_precision_integer
  USE high_precision_integer_mod
  IMPLICIT NONE


  CALL test_constructors()
  CALL test_conversions()
  CALL test_addition()
  CALL test_subtraction()
  CALL test_multiplication()
  CALL test_comparisons()
  CALL test_from_string()
  CALL test_is_zero()
  CALL test_power()
  CALL test_division_and_rem()
  CALL test_scaling()
  CALL test_multiplication_overflow()

CONTAINS

  ! =============================================================================
  ! Test Helper Subroutines
  ! =============================================================================


!! @return     The factorial of `n` as INTEGER(KIND=8).
  SUBROUTINE check(description, condition)
    CHARACTER(LEN=*), INTENT(IN) :: description
    LOGICAL, INTENT(IN)          :: condition
    IF (condition) THEN
      WRITE(*, '(A, A)') "[PASS] ", description
    ELSE
      WRITE(*, '(A, A)') "[FAIL] ", description
    END IF
  END SUBROUTINE check

  SUBROUTINE check_string(description, actual, expected)
    CHARACTER(LEN=*), INTENT(IN) :: description, actual, expected
    CALL check(description, TRIM(actual) == TRIM(expected))
    IF (TRIM(actual) /= TRIM(expected)) THEN
        WRITE(*, '(2A)') "       Expected: ", TRIM(expected)
        WRITE(*, '(2A)') "       Actual:   ", TRIM(actual)
    END IF
  END SUBROUTINE check_string

  SUBROUTINE check_hpi(description, hpi, expected_coeffs, expected_sign)
    CHARACTER(LEN=*), INTENT(IN)             :: description
    TYPE(high_precision_int(*)), INTENT(IN)  :: hpi
    INTEGER(KIND=8), DIMENSION(:), INTENT(IN) :: expected_coeffs
    INTEGER(KIND=1), INTENT(IN)              :: expected_sign
    LOGICAL                                  :: is_ok
    INTEGER                                  :: i

    is_ok = (hpi%sign == expected_sign) .AND. (hpi%ncoeffs == SIZE(expected_coeffs))
    IF (is_ok .AND. hpi%ncoeffs > 0) THEN
        is_ok = ALL(hpi%coeffs(1:hpi%ncoeffs) == expected_coeffs)
    END IF

    CALL check(description, is_ok)
    IF (.NOT. is_ok) THEN
        CALL print_hpi(hpi)
        WRITE(*, '(A, I0)') "       Expected sign: ", expected_sign
        WRITE(*, '(A)', ADVANCE='NO') "       Expected coeffs: ["
        IF (SIZE(expected_coeffs) > 0) THEN
            WRITE(*, '(I0, A)', ADVANCE='NO') (expected_coeffs(i), ", ", i=1, SIZE(expected_coeffs)-1)
            WRITE(*, '(I0)', ADVANCE='NO') expected_coeffs(SIZE(expected_coeffs))
        END IF
        WRITE(*, '(A)') "]"
    END IF
  END SUBROUTINE check_hpi

  ! =============================================================================
  ! Test Suites
  ! =============================================================================

  SUBROUTINE test_constructors()
    TYPE(high_precision_int(:)), ALLOCATABLE :: hpi
    PRINT *
    PRINT *, "--- Testing Constructors ---"

    ! From integer
    hpi = new_hpi_from_integer(0_8)
    CALL check_hpi("From integer: 0", hpi, [0_8], 0_1)

    hpi = new_hpi_from_integer(12345_8)
    CALL check_hpi("From integer: small positive", hpi, [12345_8], 1_1)

    hpi = new_hpi_from_integer(-54321_8)
    CALL check_hpi("From integer: small negative", hpi, [54321_8], -1_1)

    hpi = new_hpi_from_integer(HIGH_PRECISION_BASE)
    CALL check_hpi("From integer: base", hpi, [0_8, 1_8], 1_1)

    hpi = new_hpi_from_integer(HIGH_PRECISION_BASE + 123_8)
    CALL check_hpi("From integer: base + 123", hpi, [123_8, 1_8], 1_1)

    ! From coefficients
    hpi = new_hpi_from_coeffs([1_8, 2_8, 3_8], 1_1)
    CALL check_hpi("From coeffs: positive", hpi, [1_8, 2_8, 3_8], 1_1)

    hpi = new_hpi_from_coeffs([5_8, 0_8, 0_8], -1_1)
    CALL check_hpi("From coeffs: with trailing zeros (should normalize)", hpi, [5_8], -1_1)

    hpi = new_hpi_from_coeffs([0_8], 1_1)
    CALL check_hpi("From coeffs: zero", hpi, [0_8], 0_1)
  END SUBROUTINE test_constructors

  SUBROUTINE test_conversions()
    TYPE(high_precision_int(:)), ALLOCATABLE :: hpi
    PRINT *
    PRINT *, "--- Testing Conversions ---"

    ! To String
    hpi = new_hpi_from_integer(1234567890123456789_8)
    CALL check_string("To string: large positive", hpi_to_string(hpi), "1234567890123456789")

    hpi = new_hpi_from_integer(-987654321_8)
    CALL check_string("To string: negative", hpi_to_string(hpi), "-987654321")

    hpi = new_hpi_from_integer(0_8)
    CALL check_string("To string: zero", hpi_to_string(hpi), "0")

    hpi = new_hpi_from_coeffs([0_8, 0_8, 1_8], 1_1) ! 2^64
    CALL check_string("To string: 2^64", hpi_to_string(hpi), "18446744073709551616")

    ! To Integer
    hpi = new_hpi_from_integer(123_8)
    CALL check("To integer: small positive", hpi_to_integer(hpi) == 123_8)

    hpi = new_hpi_from_integer(-456_8)
    CALL check("To integer: small negative", hpi_to_integer(hpi) == -456_8)

    hpi = new_hpi_from_integer(HUGE(0_8))
    CALL check("To integer: max int8", hpi_to_integer(hpi) == HUGE(0_8))

    hpi = new_hpi_from_integer(-HUGE(0_8)-1_8)
    CALL check("To integer: min int8", hpi_to_integer(hpi) == -HUGE(0_8)-1_8)

    hpi = new_hpi_from_coeffs([0_8, 0_8, 1_8], 1_1) ! 2^64, should overflow
    CALL check("To integer: overflow positive", hpi_to_integer(hpi) == HUGE(0_8))

    hpi = new_hpi_from_coeffs([1_8, 0_8, 1_8], -1_1) ! -(2^64+1), should overflow
    CALL check("To integer: overflow negative", hpi_to_integer(hpi) == -HUGE(0_8)-1_8)
  END SUBROUTINE test_conversions

  SUBROUTINE test_addition()
    TYPE(high_precision_int(:)), ALLOCATABLE :: a, b, res
    PRINT *
    PRINT *, "--- Testing Addition (+) ---"

    a = new_hpi_from_integer(100_8)
    b = new_hpi_from_integer(23_8)
    res = a + b
    CALL check_hpi("pos + pos (no carry)", res, [123_8], 1_1)

    a = new_hpi_from_integer(HIGH_PRECISION_BASE - 10_8)
    b = new_hpi_from_integer(20_8)
    res = a + b
    CALL check_hpi("pos + pos (with carry)", res, [10_8, 1_8], 1_1)

    a = new_hpi_from_integer(-100_8)
    b = new_hpi_from_integer(-23_8)
    res = a + b
    CALL check_hpi("neg + neg", res, [123_8], -1_1)

    a = new_hpi_from_integer(100_8)
    b = new_hpi_from_integer(-23_8)
    res = a + b
    CALL check_hpi("pos + neg, |pos|>|neg|", res, [77_8], 1_1)

    a = new_hpi_from_integer(23_8)
    b = new_hpi_from_integer(-100_8)
    res = a + b
    CALL check_hpi("pos + neg, |pos|<|neg|", res, [77_8], -1_1)

    a = new_hpi_from_integer(100_8)
    b = new_hpi_from_integer(-100_8)
    res = a + b
    CALL check_hpi("a + (-a) = 0", res, [0_8], 0_1)

    a = new_hpi_from_integer(12345_8)
    b = new_hpi_from_integer(0_8)
    res = a + b
    CALL check_hpi("a + 0 = a", res, [12345_8], 1_1)

    a = new_hpi_from_coeffs([1_8, 1_8], 1_1) ! 2^32+1
    b = new_hpi_from_coeffs([HIGH_PRECISION_BASE - 1_8], 1_1) ! 2^32-1
    res = a + b ! (2^32+1) + (2^32-1) = 2*2^32 = 2^33
    CALL check_hpi("multi-coeff add", res, [0_8, 2_8], 1_1)

    ! Large number addition (> 2^64)
    a = new_hpi_from_coeffs([5_8, 0_8, 1_8], 1_1) ! 2^64 + 5
    b = new_hpi_from_coeffs([10_8, 0_8, 1_8], 1_1) ! 2^64 + 10
    res = a + b ! 2*2^64 + 15 = 2^65 + 15
    CALL check_hpi("large number add (>2^64)", res, [15_8, 0_8, 2_8], 1_1)
  END SUBROUTINE test_addition

  SUBROUTINE test_subtraction()
    TYPE(high_precision_int(:)), ALLOCATABLE :: a, b, res
    PRINT *
    PRINT *, "--- Testing Subtraction (-) ---"

    a = new_hpi_from_integer(100_8)
    b = new_hpi_from_integer(23_8)
    res = a - b
    CALL check_hpi("pos - pos, a > b", res, [77_8], 1_1)

    a = new_hpi_from_integer(23_8)
    b = new_hpi_from_integer(100_8)
    res = a - b
    CALL check_hpi("pos - pos, a < b", res, [77_8], -1_1)

    a = new_hpi_from_integer(100_8)
    b = new_hpi_from_integer(100_8)
    res = a - b
    CALL check_hpi("a - a = 0", res, [0_8], 0_1)

    a = new_hpi_from_integer(-50_8)
    b = new_hpi_from_integer(20_8)
    res = a - b
    CALL check_hpi("neg - pos", res, [70_8], -1_1)

    a = new_hpi_from_integer(50_8)
    b = new_hpi_from_integer(-20_8)
    res = a - b
    CALL check_hpi("pos - neg", res, [70_8], 1_1)

    a = new_hpi_from_coeffs([0_8, 1_8], 1_1) ! 2^32
    b = new_hpi_from_integer(1_8)
    res = a - b
    CALL check_hpi("multi-coeff subtract (borrow)", res, [HIGH_PRECISION_BASE - 1_8], 1_1)

    ! Large number subtraction (> 2^64)
    a = new_hpi_from_coeffs([10_8, 0_8, 2_8], 1_1) ! 2^65 + 10
    b = new_hpi_from_coeffs([5_8, 0_8, 1_8], 1_1)  ! 2^64 + 5
    res = a - b ! (2^65-2^64) + (10-5) = 2^64 + 5
    CALL check_hpi("large number subtract (>2^64)", res, [5_8, 0_8, 1_8], 1_1)
  END SUBROUTINE test_subtraction

  SUBROUTINE test_multiplication()
    TYPE(high_precision_int(:)), ALLOCATABLE :: a, b, res
    PRINT *
    PRINT *, "--- Testing Multiplication (*) ---"

    a = new_hpi_from_integer(123_8)
    b = new_hpi_from_integer(0_8)
    res = a * b
    CALL check_hpi("a * 0 = 0", res, [0_8], 0_1)

    a = new_hpi_from_integer(123_8)
    b = new_hpi_from_integer(1_8)
    res = a * b
    CALL check_hpi("a * 1 = a", res, [123_8], 1_1)

    a = new_hpi_from_integer(123_8)
    b = new_hpi_from_integer(-1_8)
    res = a * b
    CALL check_hpi("a * -1 = -a", res, [123_8], -1_1)

    a = new_hpi_from_integer(1000_8)
    b = new_hpi_from_integer(2000_8)
    res = a * b
    CALL check_hpi("pos * pos", res, [2000000_8], 1_1)

    a = new_hpi_from_integer(1000_8)
    b = new_hpi_from_integer(-2000_8)
    res = a * b
    CALL check_hpi("pos * neg", res, [2000000_8], -1_1)

    a = new_hpi_from_integer(-1000_8)
    b = new_hpi_from_integer(-2000_8)
    res = a * b
    CALL check_hpi("neg * neg", res, [2000000_8], 1_1)

    a = new_hpi_from_integer(HIGH_PRECISION_BASE - 1_8)
    b = new_hpi_from_integer(2_8)
    res = a * b
    CALL check_hpi("mult creates carry", res, [HIGH_PRECISION_BASE - 2_8, 1_8], 1_1)

    a = new_hpi_from_coeffs([1_8, 1_8], 1_1) ! 2^32+1
    b = new_hpi_from_coeffs([1_8, 1_8], 1_1) ! (2^32+1)^2 = 2^64 + 2*2^32 + 1
    res = a * b
    CALL check_hpi("multi-coeff * multi-coeff", res, [1_8, 2_8, 1_8], 1_1)

    ! Large number multiplication (> 2^64)
    a = new_hpi_from_coeffs([0_8, 0_8, 1_8], 1_1) ! 2^64
    b = new_hpi_from_coeffs([0_8, 1_8], 1_1)     ! 2^32
    res = a * b ! 2^64 * 2^32 = 2^96
    CALL check_hpi("large number mult (2^64 * 2^32)", res, [0_8, 0_8, 0_8, 1_8], 1_1)
  END SUBROUTINE test_multiplication

  ! This test specifically targets the integer overflow bug in the naive hpi_multiply.
  SUBROUTINE test_multiplication_overflow()
    TYPE(high_precision_int(:)), ALLOCATABLE :: a, b, res, expected
    PRINT *
    PRINT *, "--- Testing Multiplication Overflow Case ---"

    ! Test (2^32 - 1) * (2^32 - 1) = 2^64 - 2*2^32 + 1
    ! This will overflow a standard INTEGER(KIND=8) during intermediate calculation.
    a = new_hpi_from_coeffs([MASK32], 1_1)
    b = new_hpi_from_coeffs([MASK32], 1_1)
    res = a * b
    expected = new_hpi_from_coeffs([1_8, HIGH_PRECISION_BASE - 2_8], 1_1)
    CALL check("multiplication (2^32-1)^2", res == expected)

    ! Test (2^64 - 1) * (2^32 - 1) = 2^96 - 2^64 - 2^32 + 1
    ! This tests overflow with multi-coefficient numbers.
    a = new_hpi_from_coeffs([MASK32, MASK32], 1_1)
    b = new_hpi_from_coeffs([MASK32], 1_1)
    res = a * b
    expected = new_hpi_from_coeffs([1_8, MASK32, MASK32 - 1_8], 1_1)
    CALL check("multiplication (2^64-1)*(2^32-1)", res == expected)
  END SUBROUTINE test_multiplication_overflow

  SUBROUTINE test_comparisons()
    TYPE(high_precision_int(:)), ALLOCATABLE :: p5, p10, p100, m5, m10, m100, z
    PRINT *
    PRINT *, "--- Testing Comparisons (==, <) ---"

    p5 = new_hpi_from_integer(5_8)
    p10 = new_hpi_from_integer(10_8)
    m5 = new_hpi_from_integer(-5_8)
    z = new_hpi_from_integer(0_8)

    ! Equality
    CALL check("a == a", p5 == p5)
    CALL check("a == b (same value)", p5 == new_hpi_from_integer(5_8))
    CALL check("a /= b", .NOT. (p5 == p10))
    CALL check("a /= -a", .NOT. (p5 == m5))
    CALL check("0 == 0", z == new_hpi_from_integer(0_8))

    ! Less Than
    CALL check("5 < 10", p5 < p10)
    CALL check("NOT 10 < 5", .NOT. (p10 < p5))
    CALL check("NOT 5 < 5", .NOT. (p5 < p5))
    CALL check("-5 < 5", m5 < p5)
    CALL check("NOT 5 < -5", .NOT. (p5 < m5))
    CALL check("-5 < 0", m5 < z)
    CALL check("0 < 5", z < p5)

    ! Test with different lengths
    p10 = new_hpi_from_integer(10_8)
    p100 = new_hpi_from_integer(HIGH_PRECISION_BASE)
    CALL check("len(a) < len(b) => a < b", p10 < p100)
    CALL check("NOT len(b) < len(a)", .NOT. (p100 < p10))

    m10 = -p10
    m100 = -p100
    CALL check("neg, len(a) < len(b) => a > b", .NOT. (m10 < m100))

  END SUBROUTINE test_comparisons

SUBROUTINE test_from_string() ! <--- Defined here
    TYPE(high_precision_int(:)), ALLOCATABLE :: hpi
    PRINT *
    PRINT *, "--- Testing new_hpi_from_string ---"

    hpi = new_hpi_from_string("0")
    CALL check_hpi("From string: '0'", hpi, [0_8], 0_1)

    hpi = new_hpi_from_string("-0")
    CALL check_hpi("From string: '-0'", hpi, [0_8], 0_1)

    hpi = new_hpi_from_string("000")
    CALL check_hpi("From string: '000'", hpi, [0_8], 0_1)

    hpi = new_hpi_from_string("34")
    CALL check_hpi("From string: '34'", hpi, [34_8], 1_1)

    hpi = new_hpi_from_string("12345")
    CALL check_hpi("From string: small positive", hpi, [12345_8], 1_1)

    hpi = new_hpi_from_string("-54321")
    CALL check_hpi("From string: small negative", hpi, [54321_8], -1_1)

    hpi = new_hpi_from_string("+98765")
    CALL check_hpi("From string: with plus sign", hpi, [98765_8], 1_1)

    ! Test with a number that requires multiple coefficients (2^32 + 1)
    CALL check_string("From string: 2^32+1 to string", hpi_to_string(new_hpi_from_string("4294967297")), "4294967297")
    hpi = new_hpi_from_string("4294967297") ! 2^32 + 1
    CALL check_hpi("From string: 2^32+1", hpi, [1_8, 1_8], 1_1)

    ! Test with a very large number (10^30)
    ! 10^30 = 1 * (10^9)^3.
    ! This will be 1 followed by 30 zeros.
    CALL check_string("From string: 10^30 to string", hpi_to_string(new_hpi_from_string("1" // REPEAT("0", 30))), "1" // REPEAT("0", 30))
    hpi = new_hpi_from_string("1" // REPEAT("0", 30)) ! log2(10^30) is approx 99.6, so it needs up to coeff for 2^96
    CALL check("From string: 10^30 (check size)", hpi%ncoeffs == 4)
  END SUBROUTINE test_from_string

  SUBROUTINE test_is_zero()
    TYPE(high_precision_int(:)), ALLOCATABLE :: hpi
    PRINT *
    PRINT *, "--- Testing hpi_is_zero ---"

    hpi = new_hpi_from_integer(0_8)
    CALL check("is_zero: 0", hpi_is_zero(hpi))

    hpi = new_hpi_from_coeffs([0_8, 0_8, 0_8], 1_1) ! Will be normalized to zero
    CALL check("is_zero: from zero coeffs", hpi_is_zero(hpi))

    hpi = new_hpi_from_integer(1_8)
    CALL check("is_zero: non-zero positive", .NOT. hpi_is_zero(hpi))

    hpi = new_hpi_from_integer(-1_8)
    CALL check("is_zero: non-zero negative", .NOT. hpi_is_zero(hpi))
  END SUBROUTINE test_is_zero

  SUBROUTINE test_power()
    TYPE(high_precision_int(:)), ALLOCATABLE :: base, res
    PRINT *
    PRINT *, "--- Testing Power (hpi_power) ---"

    ! base^0 = 1
    base = new_hpi_from_integer(123_8)
    res = hpi_power(base, 0)
    CALL check_hpi("base^0 = 1", res, [1_8], 1_1)

    ! 1^exp = 1
    base = new_hpi_from_integer(1_8)
    res = hpi_power(base, 100)
    CALL check_hpi("1^exp = 1", res, [1_8], 1_1)

    ! 0^exp = 0
    base = new_hpi_from_integer(0_8)
    res = hpi_power(base, 100)
    CALL check_hpi("0^exp = 0", res, [0_8], 0_1)

    ! Simple power
    base = new_hpi_from_integer(2_8)
    res = hpi_power(base, 10)
    CALL check_hpi("2^10 = 1024", res, [1024_8], 1_1)

    ! Negative base, odd power
    base = new_hpi_from_integer(-2_8)
    res = hpi_power(base, 3)
    CALL check_hpi("(-2)^3 = -8", res, [8_8], -1_1)

    ! Negative base, even power
    base = new_hpi_from_integer(-2_8)
    res = hpi_power(base, 4)
    CALL check_hpi("(-2)^4 = 16", res, [16_8], 1_1)

    ! Power resulting in multi-coeff
    base = new_hpi_from_integer(2_8**16) ! 65536
    res = hpi_power(base, 2) ! (2^16)^2 = 2^32
    CALL check_hpi("(2^16)^2 = 2^32", res, [0_8, 1_8], 1_1)

    ! Large power
    base = new_hpi_from_integer(2_8)
    res = hpi_power(base, 100) ! 2^100 = 16 * (2^32)^3
    CALL check_hpi("2^100 (large power)", res, [0_8, 0_8, 0_8, 16_8], 1_1)

    ! -1 to odd power
    base = new_hpi_from_integer(-1_8)
    res = hpi_power(base, 99)
    CALL check_hpi("(-1)^99 = -1", res, [1_8], -1_1)

    ! -1 to even power
    base = new_hpi_from_integer(-1_8)
    res = hpi_power(base, 100)
    CALL check_hpi("(-1)^100 = 1", res, [1_8], 1_1)
  END SUBROUTINE test_power

  SUBROUTINE test_division_and_rem()
    TYPE(high_precision_int(:)), ALLOCATABLE :: a, b, q, r, check_val
    PRINT *
    PRINT *, "--- Testing Division and Remainder (hpi_divide, hpi_div_rem) ---"

    ! --- hpi_div_rem ---
    a = new_hpi_from_integer(100_8)
    b = new_hpi_from_integer(9_8)
    CALL hpi_div_rem(a, b, q, r)
    CALL check_hpi("div_rem: 100/9, quotient", q, [11_8], 1_1)
    CALL check_hpi("div_rem: 100/9, remainder", r, [1_8], 1_1)
    check_val = (q * b) + r
    CALL check("div_rem: (q*b)+r == a (pos/pos)", check_val == a)

    a = new_hpi_from_integer(-100_8)
    b = new_hpi_from_integer(9_8)
    CALL hpi_div_rem(a, b, q, r)
    CALL check_hpi("div_rem: -100/9, quotient", q, [11_8], -1_1)
    CALL check_hpi("div_rem: -100/9, remainder", r, [1_8], -1_1)
    check_val = (q * b) + r
    CALL check("div_rem: (q*b)+r == a (neg/pos)", check_val == a)

    a = new_hpi_from_integer(100_8)
    b = new_hpi_from_integer(-9_8)
    CALL hpi_div_rem(a, b, q, r)
    CALL check_hpi("div_rem: 100/-9, quotient", q, [11_8], -1_1)
    CALL check_hpi("div_rem: 100/-9, remainder", r, [1_8], 1_1)
    check_val = (q * b) + r
    CALL check("div_rem: (q*b)+r == a (pos/neg)", check_val == a)

    a = new_hpi_from_integer(-100_8)
    b = new_hpi_from_integer(-9_8)
    CALL hpi_div_rem(a, b, q, r)
    CALL check_hpi("div_rem: -100/-9, quotient", q, [11_8], 1_1)
    CALL check_hpi("div_rem: -100/-9, remainder", r, [1_8], -1_1)
    check_val = (q * b) + r
    CALL check("div_rem: (q*b)+r == a (neg/neg)", check_val == a)

    a = new_hpi_from_integer(5_8)
    b = new_hpi_from_integer(10_8)
    CALL hpi_div_rem(a, b, q, r)
    CALL check_hpi("div_rem: 5/10, quotient", q, [0_8], 0_1)
    CALL check_hpi("div_rem: 5/10, remainder", r, [5_8], 1_1)
    check_val = (q * b) + r
    CALL check("div_rem: (q*b)+r == a (num < den)", check_val == a)

    ! Large number division
    a = hpi_power(new_hpi_from_integer(2_8), 100) ! 2^100
    b = new_hpi_from_integer(3_8)
    CALL hpi_div_rem(a, b, q, r)
    check_val = (q * b) + r
    CALL check("div_rem: (q*b)+r == a (large number)", check_val == a)
    CALL check_hpi("div_rem: 2^100 / 3, remainder is 1", r, [1_8], 1_1)

    ! Large / Large with exact result: (2^64-1)/(2^32-1) = 2^32+1
    a = hpi_power(new_hpi_from_integer(2_8), 64) - new_hpi_from_integer(1_8)
    b = hpi_power(new_hpi_from_integer(2_8), 32) - new_hpi_from_integer(1_8)
    CALL hpi_div_rem(a, b, q, r)
    CALL check("div_rem: exact large division remainder is 0", hpi_is_zero(r))
    CALL check("div_rem: exact large division quotient check", q == (hpi_power(new_hpi_from_integer(2_8), 32) + new_hpi_from_integer(1_8)))

    ! Large / Large with remainder: (2^100+1)/(2^50+1) -> q=2^50-1, r=2
    a = hpi_power(new_hpi_from_integer(2_8), 100) + new_hpi_from_integer(1_8)
    b = hpi_power(new_hpi_from_integer(2_8), 50) + new_hpi_from_integer(1_8)
    CALL hpi_div_rem(a, b, q, r)
    check_val = (q * b) + r
    CALL check("div_rem: (q*b)+r == a (large/large)", check_val == a)
    CALL check("div_rem: large/large quotient check", q == (hpi_power(new_hpi_from_integer(2_8), 50) - new_hpi_from_integer(1_8)))
    CALL check_hpi("div_rem: large/large remainder check", r, [2_8], 1_1)

    ! Very large / large with exact result: (2^128-1)/(2^64+1) = 2^64-1
    a = hpi_power(new_hpi_from_integer(2_8), 128) - new_hpi_from_integer(1_8)
    b = hpi_power(new_hpi_from_integer(2_8), 64) + new_hpi_from_integer(1_8)
    CALL hpi_div_rem(a, b, q, r)
    CALL check("div_rem: very large exact division remainder is 0", hpi_is_zero(r))
    CALL check("div_rem: very large exact division quotient check", q == (hpi_power(new_hpi_from_integer(2_8), 64) - new_hpi_from_integer(1_8)))

    ! --- hpi_divide (/) ---
    a = new_hpi_from_integer(100_8); b = new_hpi_from_integer(10_8); q = a / b
    CALL check_hpi("divide: 100/10", q, [10_8], 1_1)

    a = new_hpi_from_integer(0_8); b = new_hpi_from_integer(10_8); q = a / b
    CALL check_hpi("divide: 0/10", q, [0_8], 0_1)

    a = new_hpi_from_integer(123_8); b = new_hpi_from_integer(1_8); q = a / b
    CALL check_hpi("divide: a/1", q, [123_8], 1_1)

    a = new_hpi_from_integer(123_8); b = new_hpi_from_integer(123_8); q = a / b
    CALL check_hpi("divide: a/a", q, [1_8], 1_1)

    a = new_hpi_from_coeffs([0_8, 0_8, 1_8], 1_1) ! 2^64
    b = new_hpi_from_integer(2_8)
    q = a / b ! 2^63
    CALL check_hpi("divide: 2^64 / 2 = 2^63", q, [0_8, ISHFT(HIGH_PRECISION_BASE, -1)], 1_1)
  END SUBROUTINE test_division_and_rem

  SUBROUTINE test_scaling()
    TYPE(high_precision_int(:)), ALLOCATABLE :: hpi, res
    PRINT *
    PRINT *, "--- Testing Scaling (hpi_scale_up_by_base_power) ---"

    hpi = new_hpi_from_integer(123_8)
    res = hpi_scale_up_by_base_power(hpi, 0)
    CALL check_hpi("scale by 0", res, [123_8], 1_1)

    res = hpi_scale_up_by_base_power(hpi, 1)
    CALL check_hpi("scale by 1", res, [0_8, 123_8], 1_1)

    res = hpi_scale_up_by_base_power(hpi, 3)
    CALL check_hpi("scale by 3", res, [0_8, 0_8, 0_8, 123_8], 1_1)

    hpi = new_hpi_from_integer(0_8)
    res = hpi_scale_up_by_base_power(hpi, 5)
    CALL check_hpi("scale zero", res, [0_8], 0_1)

    ! Scale a multi-coefficient number
    hpi = new_hpi_from_coeffs([123_8, 456_8], 1_1)
    res = hpi_scale_up_by_base_power(hpi, 2)
    CALL check_hpi("scale multi-coeff number", res, [0_8, 0_8, 123_8, 456_8], 1_1)

    ! Scale a negative number
    hpi = new_hpi_from_integer(-789_8)
    res = hpi_scale_up_by_base_power(hpi, 1)
    CALL check_hpi("scale negative number", res, [0_8, 789_8], -1_1)

    ! Scale a large number (> 2^64)
    hpi = new_hpi_from_coeffs([123_8, 0_8, 1_8], 1_1) ! 2^64 + 123
    res = hpi_scale_up_by_base_power(hpi, 2)
    CALL check_hpi("scale large number (>2^64)", res, [0_8, 0_8, 123_8, 0_8, 1_8], 1_1)
  END SUBROUTINE test_scaling

END PROGRAM test_high_precision_integer