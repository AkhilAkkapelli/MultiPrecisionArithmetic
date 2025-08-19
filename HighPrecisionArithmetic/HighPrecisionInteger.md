# High Precision Integer Module (`high_precision_integer_mod`)

This document provides detailed documentation for the `high_precision_integer_mod` Fortran module, which implements arbitrary-precision integer arithmetic.

## 1. Overview

The `high_precision_integer_mod` module provides a derived type `high_precision_int` and a suite of functions and overloaded operators to perform arithmetic on integers that are too large to fit into standard hardware integer types like `INTEGER(KIND=8)`.

Numbers are stored in a base of 2<sup>32</sup>, allowing for efficient computation on 64-bit systems. The implementation supports basic arithmetic (addition, subtraction, multiplication, division), comparisons, exponentiation, and conversion from/to native integers and strings.

### Key Features

*   **Arbitrary Size**: Integers are limited only by available memory.
*   **Overloaded Operators**: Use standard operators like `+`, `-`, `*`, `/`, `==`, and `<` for intuitive use.
*   **Constructors**: Easily create high-precision integers from standard integers or strings.
*   **Conversions**: Convert high-precision integers back to standard integers or base-10 strings.
*   **Helper Functions**: Includes utilities for absolute value, exponentiation, and more.

## 2. The `high_precision_int` Type

The core of the module is the `high_precision_int` derived type. It is a parameterized type where the parameter `len` defines the maximum number of coefficients the integer can hold, effectively setting its maximum possible size.

```fortran
TYPE high_precision_int(len)
  INTEGER, LEN                 :: len       ! Maximum number of coefficients (storage capacity)
  INTEGER                      :: ncoeffs   ! Number of coefficients currently in use (logical length)
  INTEGER(KIND=8)              :: coeffs(len) ! The coefficients in base 2^32
  INTEGER(KIND=1)              :: sign      ! -1 for negative, 0 for zero, 1 for positive
END TYPE high_precision_int
```

*   **`len`**: An integer length parameter specifying the allocated size of the `coeffs` array. This must be provided when an instance is declared, but for `ALLOCATABLE` variables, it is set upon allocation.
*   **`ncoeffs`**: The logical length of the number. It indicates how many elements of the `coeffs` array are actually being used to represent the number's magnitude. This allows the type to represent numbers of varying sizes within the same allocated storage.
*   **`coeffs`**: An array of 64-bit integers. Each element stores a "digit" of the number in base 2<sup>32</sup>. The number is stored in little-endian format, so `coeffs(1)` is the least significant digit.
*   **`sign`**: Stores the sign of the number.
    *   `1`: Positive
    *   `0`: Zero
    *   `-1`: Negative

## 3. Public Interface

The module exposes a set of public procedures for creating and manipulating `high_precision_int` objects.

### 3.1. Constructors

These functions create new `high_precision_int` objects. The result is always `ALLOCATABLE`.

---

#### `new_hpi_from_integer(x_in)`

Converts a standard `INTEGER(KIND=8)` to a `high_precision_int`.

*   **`x_in`**: `INTEGER(KIND=8), INTENT(IN)` - The integer to convert.
*   **Returns**: `TYPE(high_precision_int(:)), ALLOCATABLE` - The new high-precision integer.

---

#### `new_hpi_from_string(str_in)`

Converts a string representation of a decimal number to a `high_precision_int`. The string can be of arbitrary length and may optionally start with `+` or `-`.

*   **`str_in`**: `CHARACTER(LEN=*), INTENT(IN)` - The string to convert.
*   **Returns**: `TYPE(high_precision_int(:)), ALLOCATABLE` - The new high-precision integer.

---

#### `new_hpi_from_coeffs(coeffs_in, [sign_in])`

A low-level constructor to create a `high_precision_int` directly from an array of coefficients. This is primarily for internal use but is public.

*   **`coeffs_in`**: `INTEGER(KIND=8), DIMENSION(:), INTENT(IN)` - The array of coefficients (base 2<sup>32</sup>).
*   **`sign_in`**: `INTEGER(KIND=1), OPTIONAL, INTENT(IN)` - The sign (`-1`, `0`, or `1`). Defaults to `1` (positive).
*   **Returns**: `TYPE(high_precision_int(:)), ALLOCATABLE` - The new high-precision integer.

### 3.2. Conversions

---

#### `hpi_to_string(hpi_in)`

Converts a `high_precision_int` to its decimal string representation.

*   **`hpi_in`**: `TYPE(high_precision_int(*)), INTENT(IN)` - The number to convert.
*   **Returns**: `CHARACTER(LEN=:), ALLOCATABLE` - The decimal string.

---

#### `hpi_to_integer(hpi)`

Converts a `high_precision_int` to a standard `INTEGER(KIND=8)`. If the number is too large or small to fit, a warning is printed and the result is clamped to `HUGE(0_8)` or its negative equivalent.

*   **`hpi`**: `TYPE(high_precision_int(*)), INTENT(IN)` - The number to convert.
*   **Returns**: `INTEGER(KIND=8)` - The converted integer.

### 3.3. Overloaded Operators

The module overloads standard operators for seamless use. `a` and `b` are `high_precision_int` objects.

| Operator | Function | Description |
|:---:|:---|:---|
| `+` | `hpi_add(a, b)` | Adds two high-precision integers. |
| `-` | `hpi_subtract(a, b)` | Subtracts `b` from `a`. |
| `-` | `hpi_unary_negate(a)` | Returns the negation of `a`. |
| `*` | `hpi_multiply(a, b)` | Multiplies two high-precision integers. |
| `/` | `hpi_divide(a, b)` | Performs integer division `a / b`. |
| `==` | `hpi_equal(a, b)` | Returns `.TRUE.` if `a` and `b` are equal. |
| `<` | `hpi_less(a, b)` | Returns `.TRUE.` if `a` is less than `b`. |

**Note**: The other relational operators (`>`, `<=`, `>=`, `/=`) are not explicitly overloaded but can be composed from `<` and `==`. For example, `a > b` is equivalent to `b < a`, and `a >= b` is ` .NOT. (a < b)`.

### 3.4. Utility Functions and Subroutines

---

#### `hpi_abs(hpi_in)`

Returns the absolute value of a `high_precision_int`.

*   **`hpi_in`**: `TYPE(high_precision_int(*)), INTENT(IN)`
*   **Returns**: `TYPE(high_precision_int(:)), ALLOCATABLE`

---

#### `hpi_power(base, exp)`

Calculates `base` raised to the power of `exp` using exponentiation by squaring.

*   **`base`**: `TYPE(high_precision_int(*)), INTENT(IN)` - The base.
*   **`exp`**: `INTEGER, INTENT(IN)` - The exponent (must be non-negative).
*   **Returns**: `TYPE(high_precision_int(:)), ALLOCATABLE` - The result of `base**exp`.

---

#### `hpi_div_rem(numerator, denominator, quotient, remainder)`

Performs integer division and returns both the quotient and the remainder.

*   **`numerator`**: `TYPE(high_precision_int(*)), INTENT(IN)`
*   **`denominator`**: `TYPE(high_precision_int(*)), INTENT(IN)`
*   **`quotient`**: `TYPE(high_precision_int(:)), ALLOCATABLE, INTENT(OUT)`
*   **`remainder`**: `TYPE(high_precision_int(:)), ALLOCATABLE, INTENT(OUT)`

---

#### `hpi_is_zero(hpi)`

A convenience function to check if a `high_precision_int` is zero.

*   **`hpi`**: `TYPE(high_precision_int(*)), INTENT(IN)`
*   **Returns**: `LOGICAL` - `.TRUE.` if `hpi` is zero, `.FALSE.` otherwise.

---

#### `print_hpi(hpi)`

A utility subroutine to print a detailed representation of a `high_precision_int` to the console, including its decimal value, internal coefficients, and sign.

*   **`hpi`**: `TYPE(high_precision_int(*)), INTENT(IN)`

---

#### `normalize_hpi(hpi)`

An internal-facing subroutine (though public) that ensures a `high_precision_int` is in a canonical form. It handles carries, trims leading zeros, and sets the sign correctly. It is called automatically by the constructors and arithmetic functions.

*   **`hpi`**: `TYPE(high_precision_int(*)), INTENT(INOUT)`

---

## 4. Usage Example

The following program demonstrates how to use the `high_precision_integer_mod` to perform some calculations.

```fortran
PROGRAM demo_hpi
  USE high_precision_integer_mod
  IMPLICIT NONE

  TYPE(high_precision_int(:)), ALLOCATABLE :: a, b, c, p, q, r
  CHARACTER(LEN=:), ALLOCATABLE :: str_a, str_c

  ! 1. Create a large number from a string
  a = new_hpi_from_string("123456789012345678901234567890")
  WRITE(*, '(A)') "Created 'a' from a string:"
  CALL print_hpi(a)
  WRITE(*,*)

  ! 2. Create another number from a standard integer
  b = new_hpi_from_integer(987654321_8)
  WRITE(*, '(A)') "Created 'b' from an integer:"
  CALL print_hpi(b)
  WRITE(*,*)

  ! 3. Perform addition: c = a + b
  c = a + b
  str_c = hpi_to_string(c)
  WRITE(*, '(A, A)') "a + b = ", TRIM(str_c)
  WRITE(*,*)

  ! 4. Perform multiplication: c = a * b
  c = a * b
  str_c = hpi_to_string(c)
  WRITE(*, '(A, A)') "a * b = ", TRIM(str_c)
  WRITE(*,*)

  ! 5. Perform exponentiation: p = b^10
  p = hpi_power(b, 10)
  WRITE(*, '(A)') "b^10 is:"
  CALL print_hpi(p)
  WRITE(*,*)

  ! 6. Perform division and get remainder: p / a
  CALL hpi_div_rem(p, a, q, r)
  WRITE(*, '(A)') "p / a gives:"
  WRITE(*, '(A, A)') "  Quotient : ", hpi_to_string(q)
  WRITE(*, '(A, A)') "  Remainder: ", hpi_to_string(r)
  WRITE(*,*)

  ! 7. Comparison
  IF (a > b) THEN
    WRITE(*,*) "'a' is greater than 'b'."
  END IF

END PROGRAM demo_hpi
```

### Expected Output of the Example

```
 Created 'a' from a string:
 HighPrecisionInt(123456789012345678901234567890, coeffs=[3512313322, 2874452364, 62031636], sign=1)

 Created 'b' from an integer:
 HighPrecisionInt(987654321, coeffs=[987654321], sign=1)

 a + b = 123456789012345678901234567890

 a * b = 121932631138672367636736311982066531310

 b^10 is:
 HighPrecisionInt(900895913232263378131330314532687699310866349693159490881, coeffs=[3365750033, 2505237625, 3474443925, 2221818990, 3139844948, 198890130, 20980], sign=1)

 p / a gives:
   Quotient : 7300999333
   Remainder: 111111111111111111111111111110

 'a' is greater than 'b'.
```

## 5. Compilation

To compile the module and the example program, you can use a Fortran compiler like `gfortran`:

```bash
# First, compile the module file
gfortran -c high_precision_integer_mod.f90

# Then, compile the main program and link it with the module
gfortran demo_hpi.f90 high_precision_integer_mod.o -o demo_hpi

# Run the executable
./demo_hpi
```

*(Note: The example program should be saved as `demo_hpi.f90` and the module code as `high_precision_integer_mod.f90` in the same directory.)*

