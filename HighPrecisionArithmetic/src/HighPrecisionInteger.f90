MODULE high_precision_integer_mod
  IMPLICIT NONE


  TYPE high_precision_int(len)
    INTEGER, LEN    :: len         ! L : total capacity
    INTEGER         :: ncoeffs = 0 ! l : significant length
    INTEGER(KIND=8) :: coeffs(len) ! c_1, c_2, ..., c_L: base-B digits
    INTEGER(KIND=1) :: sign        ! s: sign
  END TYPE high_precision_int


  INTERFACE OPERATOR(+)
    MODULE PROCEDURE hpi_add
  END INTERFACE

  INTERFACE OPERATOR(-)
    MODULE PROCEDURE hpi_unary_negate
    MODULE PROCEDURE hpi_subtract
  END INTERFACE 

  INTERFACE OPERATOR(*)
    MODULE PROCEDURE hpi_multiply
  END INTERFACE 

  INTERFACE OPERATOR(/)
    MODULE PROCEDURE hpi_divide
  END INTERFACE

  INTERFACE OPERATOR(==)
    MODULE PROCEDURE hpi_equal
  END INTERFACE

  INTERFACE OPERATOR(<)
    MODULE PROCEDURE hpi_less
  END INTERFACE OPERATOR(<)

  PUBLIC :: high_precision_int, new_hpi_from_coeffs, new_hpi_from_integer, new_hpi_from_string, &
            normalize_hpi, hpi_abs, hpi_to_integer, hpi_to_string, print_hpi, hpi_is_zero, &
            hpi_power, hpi_divide, hpi_div_rem, hpi_scale_up_by_base_power

CONTAINS

  SUBROUTINE normalize_hpi(hpi)
    TYPE(high_precision_int(*)), INTENT(INOUT) :: hpi

    INTEGER(KIND=8), PARAMETER :: MASK32 = INT(Z'FFFFFFFF', KIND=8) 

    INTEGER(KIND=8) :: carry             ! Overflow from previous digit
    INTEGER(KIND=8) :: current_coeff_val ! Intermediate value for digit + carry
    INTEGER         :: i                 ! Index into coefficient array 
    INTEGER         :: last_nonzero_idx  ! Tracks final non-zero digit index

    IF (hpi%ncoeffs == 0) THEN
      hpi%ncoeffs   = 1
      hpi%coeffs(1) = 0_8
      hpi%sign      = 0_1
      RETURN
    END IF

    carry = 0_8
    i = 1
    last_nonzero_idx = 1

    DO WHILE (i <= hpi%ncoeffs .OR. carry /= 0_8)
      IF (i > hpi%len) THEN
        STOP "FATAL ERROR in normalize_hpi: Overflow. Integer exceeds allocated length capacity."
      END IF

      IF (i <= hpi%ncoeffs) THEN
        current_coeff_val = hpi%coeffs(i) + carry
      ELSE
        current_coeff_val = carry
        hpi%ncoeffs = i
      END IF

      hpi%coeffs(i) = IAND(current_coeff_val, MASK32)
      carry = ISHFT(current_coeff_val, -32) 

      IF (hpi%coeffs(i) /= 0_8) THEN
        last_nonzero_idx = i
      END IF

      i = i + 1
    END DO

    hpi%ncoeffs = last_nonzero_idx

    IF (hpi%ncoeffs == 1 .AND. hpi%coeffs(1) == 0_8) THEN
      hpi%sign = 0_1
    ELSE IF (hpi%sign == 0_1) THEN
      hpi%sign = 1_1
    END IF
  END SUBROUTINE normalize_hpi

  ! FUNCTION new_hpi_from_coeffs(coeffs_in, sign_in, len_in) RESULT(hpi_out)
  !   INTEGER(KIND=8), DIMENSION(:), INTENT(IN) :: coeffs_in
  !   INTEGER(KIND=1), OPTIONAL, INTENT(IN)     :: sign_in
  !   INTEGER, OPTIONAL, INTENT(IN)            :: len_in

  !   TYPE(high_precision_int(:)), ALLOCATABLE  :: hpi_out

  !   INTEGER(KIND=8), PARAMETER :: HIGH_PRECISION_BASE = 2_8**32

  !   INTEGER(KIND=1)            :: actual_sign
  !   INTEGER                    :: min_required_len, final_len
  !   INTEGER                    :: n_in

  !   n_in = SIZE(coeffs_in)

  !   ! 1. Calculate the *minimum* required length to store the input coefficients after normalization.
  !   IF (n_in == 0) THEN
  !     min_required_len = 1
  !   ELSE
  !     min_required_len = n_in
  !     IF (coeffs_in(n_in) >= HIGH_PRECISION_BASE) THEN
  !       min_required_len = n_in + 1
  !     END IF
  !   END IF

  !   ! 2. Determine the final allocation length.
  !   IF (PRESENT(len_in)) THEN
  !     final_len = len_in
  !     IF (final_len < min_required_len) THEN
  !       WRITE(*, '(A, I0, A, I0, A)') "FATAL ERROR in new_hpi_from_coeffs: Provided 'len' (", final_len, &
  !                                     ") is too small. At least ", min_required_len, " is required for the input coefficients."
  !       STOP
  !     END IF
  !   ELSE
  !     final_len = min_required_len
  !   END IF

  !   ALLOCATE(high_precision_int(len=final_len) :: hpi_out)
  !   hpi_out%coeffs = 0_8

  !   ! 3. Set sign.
  !   IF (PRESENT(sign_in)) THEN
  !     actual_sign = sign_in
  !   ELSE
  !     actual_sign = 1_1
  !   END IF
  !   hpi_out%sign = actual_sign

  !   ! 4. Copy coefficients and set initial logical length.
  !   IF (n_in > 0) THEN
  !     hpi_out%ncoeffs = n_in
  !     hpi_out%coeffs(1:n_in) = coeffs_in
  !   ELSE
  !     hpi_out%ncoeffs = 1
  !   END IF

  !   ! 5. Normalize the number.
  !   CALL normalize_hpi(hpi_out)

  ! END FUNCTION new_hpi_from_coeffs

  FUNCTION new_hpi_from_integer(x_in) RESULT(hpi_out)
    INTEGER(KIND=8), INTENT(IN)              :: x_in

    TYPE(high_precision_int(:)), ALLOCATABLE :: hpi_out

    INTEGER          :: min_required_len
    INTEGER(KIND=1)  :: final_sign
    INTEGER(KIND=8)  :: mag_x

    INTEGER(KIND=8), PARAMETER :: HIGH_PRECISION_BASE = 2_8**32
    INTEGER(KIND=8), PARAMETER :: MASK32 = INT(Z'FFFFFFFF', KIND=8)
    INTEGER(KIND=8), PARAMETER :: MOST_NEGATIVE_I8 = -HUGE(0_8) - 1_8

    IF (x_in == 0_8) THEN
      min_required_len = 1
    ELSE IF (x_in == MOST_NEGATIVE_I8) THEN
      min_required_len = 2
    ELSE IF (ABS(x_in) < HIGH_PRECISION_BASE) THEN
      min_required_len = 1
    ELSE
      min_required_len = 2
    END IF

    ALLOCATE(high_precision_int(len=min_required_len) :: hpi_out)
    hpi_out%coeffs = 0_8 
    hpi_out%ncoeffs = min_required_len

    IF (x_in == 0_8) THEN
      hpi_out%sign = 0_1
    ELSE IF (x_in == MOST_NEGATIVE_I8) THEN
      hpi_out%sign = -1_1
      hpi_out%coeffs(1) = 0_8
      hpi_out%coeffs(2) = ISHFT(HIGH_PRECISION_BASE, -1)
    ELSE
      hpi_out%sign = SIGN(1_1, x_in)
      mag_x = ABS(x_in)
      hpi_out%coeffs(1) = IAND(mag_x, MASK32)
      IF (min_required_len == 2) THEN
        hpi_out%coeffs(2) = ISHFT(mag_x, -32)
      END IF
    END IF

  END FUNCTION new_hpi_from_integer

  FUNCTION hpi_to_integer(hpi) RESULT(x)
    TYPE(high_precision_int(*)), INTENT(IN) :: hpi

    INTEGER(KIND=8)                         :: x

    INTEGER(KIND=8) :: abs_val_as_int8
    INTEGER(KIND=8) :: c0, c1

    INTEGER(KIND=8), PARAMETER :: MAX_POS_INT8 = 2_8**31 - 1_8
    INTEGER(KIND=8), PARAMETER :: MAX_NEG_INT8 = 2_8**31

    IF (hpi%sign == 0_1) THEN
      x = 0_8
      RETURN
    END IF

    SELECT CASE (hpi%ncoeffs)
    CASE (1)
      abs_val_as_int8 = hpi%coeffs(1)
    CASE (2)
      c0 = hpi%coeffs(1)
      c1 = hpi%coeffs(2)

      IF (hpi%sign == 1_1 .AND. c1 > MAX_POS_INT8) THEN
        STOP "FATAL ERROR in hpi_to_integer: Positive value is too large to fit in an INTEGER(KIND=8)."
      ELSE IF (hpi%sign == -1_1 .AND. (c1 > ISHFT(MAX_NEG_INT8, 1) .OR. (c1 == MAX_NEG_INT8 .AND. c0 > 0_8))) THEN
        STOP "FATAL ERROR in hpi_to_integer: Negative value is too small to fit in an INTEGER(KIND=8)."
      END IF

      abs_val_as_int8 = c0 + ISHFT(c1, 32)
    CASE DEFAULT
      STOP "FATAL ERROR in hpi_to_integer: Value has more than 2 coefficients and cannot fit in an INTEGER(KIND=8)."
    END SELECT

    x = abs_val_as_int8 * hpi%sign

  END FUNCTION hpi_to_integer

  FUNCTION hpi_to_string(hpi_in) RESULT(str_out)
    TYPE(high_precision_int(*)), INTENT(IN) :: hpi_in
    CHARACTER(LEN=:), ALLOCATABLE          :: str_out

    ! --- Local variables ---
    TYPE(high_precision_int(:)), ALLOCATABLE :: temp_hpi
    ! Process number in chunks of 9 digits for efficiency. 10**9 fits in a 32-bit integer.
    INTEGER, PARAMETER                  :: NUM_DECIMAL_DIGITS = 9
    INTEGER(KIND=8), PARAMETER          :: DECIMAL_BASE = 10_8**NUM_DECIMAL_DIGITS
    CHARACTER(LEN=NUM_DECIMAL_DIGITS)   :: chunk_str
    
    ! Variables for inlined division
    INTEGER(KIND=8)                     :: current_val, remainder
    INTEGER                             :: i, new_len
    
    ! Array to store chunks of digits
    INTEGER(KIND=8), ALLOCATABLE        :: digit_chunks(:)
    INTEGER                             :: chunk_count, max_chunks
    INTEGER                             :: final_str_len, current_pos, start_pos
    CHARACTER(LEN=:), ALLOCATABLE       :: first_chunk_trimmed
    INTEGER                             :: current_hpi_ncoeffs
    INTEGER(KIND=8), ALLOCATABLE :: temp_storage(:)

    IF (hpi_in%sign == 0_1) THEN
      str_out = "0"
      RETURN
    END IF

    temp_hpi = hpi_abs(hpi_in) ! Work with the absolute value
    current_hpi_ncoeffs = temp_hpi%ncoeffs

    ! --- Step 1: Decompose into base 10**9 chunks ---
    ! Pre-allocate digit_chunks to a generous estimated size.
    max_chunks = INT(REAL(current_hpi_ncoeffs) * 9.633 / REAL(NUM_DECIMAL_DIGITS)) + 2
    ALLOCATE(digit_chunks(max_chunks))
    chunk_count = 0

    DO WHILE (temp_hpi%sign /= 0_1)
      chunk_count = chunk_count + 1
      ! If our estimation was too low (rare), grow the chunks array.
      IF (chunk_count > max_chunks) THEN
        max_chunks = chunk_count + 10
        ALLOCATE(temp_storage(max_chunks))
        temp_storage(1:chunk_count-1) = digit_chunks(1:chunk_count-1)
        CALL move_alloc(from=temp_storage, to=digit_chunks)
      END IF

      remainder = 0_8

      ! Perform division in-place on temp_hpi%coeffs using its logical length.
      DO i = current_hpi_ncoeffs, 1, -1
        current_val = temp_hpi%coeffs(i) + ISHFT(remainder, 32)
        temp_hpi%coeffs(i) = current_val / DECIMAL_BASE ! Overwrite with quotient
        remainder = MODULO(current_val, DECIMAL_BASE)
      END DO

      digit_chunks(chunk_count) = remainder

      ! Find the new logical length by trimming leading zeros from the quotient.
      new_len = current_hpi_ncoeffs
      DO WHILE (new_len > 1 .AND. temp_hpi%coeffs(new_len) == 0_8)
        new_len = new_len - 1
      END DO

      IF (new_len == 1 .AND. temp_hpi%coeffs(1) == 0_8) THEN
        temp_hpi%sign = 0_1
      ELSE
        ! Instead of reallocating, just update the logical length for the next iteration.
        current_hpi_ncoeffs = new_len
        temp_hpi%ncoeffs = new_len
      END IF
    END DO
    DEALLOCATE(temp_hpi)

    ! --- Step 2: Assemble the final string from the chunks ---
    ! Calculate final string length to allocate only once.
    WRITE(chunk_str, '(I0)') digit_chunks(chunk_count)
    first_chunk_trimmed = TRIM(ADJUSTL(chunk_str))
    final_str_len = LEN(first_chunk_trimmed)
    IF (chunk_count > 1) THEN
      final_str_len = final_str_len + (chunk_count - 1) * NUM_DECIMAL_DIGITS
    END IF
    start_pos = 1
    IF (hpi_in%sign == -1_1) THEN
      final_str_len = final_str_len + 1
      start_pos = 2
    END IF
    ALLOCATE(CHARACTER(LEN=final_str_len) :: str_out)

    ! Fill the pre-allocated string from right to left.
    current_pos = final_str_len
    DO i = 1, chunk_count - 1
      WRITE(chunk_str, '(I9.9)') digit_chunks(i)
      str_out(current_pos - NUM_DECIMAL_DIGITS + 1 : current_pos) = chunk_str
      current_pos = current_pos - NUM_DECIMAL_DIGITS
    END DO
    str_out(start_pos:current_pos) = first_chunk_trimmed
    DEALLOCATE(digit_chunks)

    ! Place the sign directly.
    IF (hpi_in%sign == -1_1) THEN
      str_out(1:1) = "-"
    END IF

  END FUNCTION hpi_to_string

  
  FUNCTION new_hpi_from_string(str_in) RESULT(hpi_out)
    CHARACTER(LEN=*), INTENT(IN) :: str_in

    TYPE(high_precision_int(:)), ALLOCATABLE :: hpi_out
    
    CHARACTER(LEN=:), ALLOCATABLE :: num_str_trimmed
    INTEGER(KIND=1)               :: input_sign
    INTEGER                       :: i, j, len_str, chunk_start, chunk_end
    INTEGER(KIND=8)               :: chunk_val
    INTEGER, PARAMETER            :: NUM_DECIMAL_DIGITS = 9
    INTEGER(KIND=8), PARAMETER    :: DECIMAL_CHUNK_BASE = 10_8**NUM_DECIMAL_DIGITS
    INTEGER                       :: numeric_len, first_chunk_len
    INTEGER, PARAMETER            :: ICHAR_ZERO = ICHAR('0') ! Pre-calculate ICHAR('0')
    TYPE(high_precision_int(:)), ALLOCATABLE :: hpi_chunk_val_hpi


    num_str_trimmed = TRIM(ADJUSTL(str_in))
    len_str = LEN(num_str_trimmed)
    IF (len_str == 0) THEN
        hpi_out = new_hpi_from_integer(0_8)
        RETURN
    END IF

    input_sign = 1_1
    i = 1
    IF (num_str_trimmed(1:1) == '-') THEN
        input_sign = -1_1
        i = 2
    ELSE IF (num_str_trimmed(1:1) == '+') THEN
        i = 2
    END IF

    ! Trim leading zeros from the numeric part (e.g., "00123" -> "123")
    DO WHILE (i <= len_str .AND. num_str_trimmed(i:i) == '0')
        i = i + 1
    END DO
    IF (i > len_str) THEN ! String was "0", "-0", "000", etc.
        hpi_out = new_hpi_from_integer(0_8)
        RETURN
    END IF

    ! Correctly calculate the length of the first chunk.
    ! It can be shorter than NUM_DECIMAL_DIGITS.
    numeric_len = len_str - (i - 1)
    first_chunk_len = MODULO(numeric_len, NUM_DECIMAL_DIGITS)
    IF (first_chunk_len == 0 .AND. numeric_len > 0) THEN
        first_chunk_len = NUM_DECIMAL_DIGITS
    END IF

    ! Initialize result with the first chunk
    chunk_start = i
    chunk_end = i + first_chunk_len - 1
    ! Manual string-to-integer conversion for the chunk. This is faster than using READ.
    chunk_val = 0_8 ! Initialize chunk_val
    DO j = chunk_start, chunk_end
        IF (num_str_trimmed(j:j) < '0' .OR. num_str_trimmed(j:j) > '9') STOP "ERROR in new_hpi_from_string: Invalid character in numeric string."
        chunk_val = chunk_val * 10_8 + (ICHAR(num_str_trimmed(j:j)) - ICHAR_ZERO)
    END DO
    hpi_out = new_hpi_from_integer(chunk_val)
    i = chunk_end + 1

    ! Allocate hpi_chunk_val_hpi once outside the loop for reuse.
    ! A single 32-bit coefficient is enough for 10^9, so len=1 is sufficient.
    ALLOCATE(high_precision_int(len=1) :: hpi_chunk_val_hpi)
    hpi_chunk_val_hpi%coeffs = 0_8
    hpi_chunk_val_hpi%ncoeffs = 1
    hpi_chunk_val_hpi%sign = 1_1 ! Chunk values are always positive

    ! Process subsequent chunks, which will all be NUM_DECIMAL_DIGITS long.
    DO WHILE (i <= len_str) ! Loop for subsequent chunks
        ! Multiply current hpi_out by DECIMAL_CHUNK_BASE (this shifts the existing number).
        hpi_out = hpi_out * hpi_decimal_chunk_base
        
        chunk_start = i
        chunk_end = i + NUM_DECIMAL_DIGITS - 1
        ! Manual string-to-integer conversion for the chunk. This is faster than using READ.
        chunk_val = 0_8 ! Initialize chunk_val
        DO j = chunk_start, chunk_end
            IF (num_str_trimmed(j:j) < '0' .OR. num_str_trimmed(j:j) > '9') STOP "ERROR in new_hpi_from_string: Invalid character in numeric string."
            chunk_val = chunk_val * 10_8 + (ICHAR(num_str_trimmed(j:j)) - ICHAR_ZERO)
        END DO

        ! Add the new chunk_val
        hpi_out = hpi_out + new_hpi_from_integer(chunk_val)

        i = chunk_end + 1
    END DO

    ! Apply final sign
    IF (input_sign == -1_1) THEN
        hpi_out = -hpi_out
    END IF
    DEALLOCATE(hpi_chunk_val_hpi) ! Deallocate the reusable HPI chunk

 END FUNCTION new_hpi_from_string

  ! Returns the absolute value of a HighPrecisionInt.
  FUNCTION hpi_abs(hpi_in) RESULT(hpi_out)
    TYPE(high_precision_int(*)), INTENT(IN) :: hpi_in
    TYPE(high_precision_int(:)), ALLOCATABLE :: hpi_out

    hpi_out = new_hpi_from_coeffs(hpi_in%coeffs(1:hpi_in%ncoeffs), ABS(hpi_in%sign))
  END FUNCTION hpi_abs

  ! Subtracts the magnitudes of two numbers represented by coefficient vectors.
  ! Computes ||a| - |b|| and returns the result coefficients and a flag indicating if |b| > |a|.
  SUBROUTINE abs_subtract(a_coeffs, b_coeffs, result_coeffs, is_negative_diff)
    INTEGER(KIND=8), DIMENSION(:), INTENT(IN)  :: a_coeffs, b_coeffs
    INTEGER(KIND=8), ALLOCATABLE, INTENT(OUT)  :: result_coeffs(:) ! Resulting coefficients
    LOGICAL, INTENT(OUT)                       :: is_negative_diff ! True if |b| > |a|

    INTEGER                                    :: len_a, len_b
    LOGICAL                                    :: a_mag_is_larger  ! True if |a| >= |b|
    INTEGER                                    :: i

    len_a = SIZE(a_coeffs)
    len_b = SIZE(b_coeffs)

    ! 1. Determine Larger Magnitude.
    a_mag_is_larger = .TRUE.
    IF (len_a < len_b) THEN
      a_mag_is_larger = .FALSE.
    ELSE IF (len_a == len_b) THEN
      DO i = len_a, 1, -1 ! Compare from most significant digit downwards
        IF (a_coeffs(i) < b_coeffs(i)) THEN
          a_mag_is_larger = .FALSE.
          EXIT
        ELSE IF (a_coeffs(i) > b_coeffs(i)) THEN
          EXIT
        END IF
      END DO
    END IF

    ! 2. Perform subtraction by calling the internal helper with the correct operands.
    IF (a_mag_is_larger) THEN
      CALL subtract_core(a_coeffs, b_coeffs, result_coeffs)
    ELSE
      CALL subtract_core(b_coeffs, a_coeffs, result_coeffs)
    END IF

    ! 3. Finalize: Set sign indicator and handle the zero result case.
    IF (SIZE(result_coeffs) == 1 .AND. result_coeffs(1) == 0_8) THEN
      is_negative_diff = .FALSE.
    ELSE
      is_negative_diff = .NOT. a_mag_is_larger
    END IF

CONTAINS

    ! Internal subroutine to subtract magnitudes, where |op1| >= |op2|.
    ! This routine performs the actual subtraction and trims the result.
    SUBROUTINE subtract_core(op1, op2, result)
      INTEGER(KIND=8), DIMENSION(:), INTENT(IN)  :: op1, op2
      INTEGER(KIND=8), ALLOCATABLE, INTENT(OUT) :: result(:)

      INTEGER(KIND=8), PARAMETER :: HIGH_PRECISION_BASE = 2_8**32

      INTEGER       :: len_op1, len_op2, i, last_idx
      INTEGER(KIND=8) :: borrow

      len_op1 = SIZE(op1)
      len_op2 = SIZE(op2)
      ALLOCATE(result(len_op1))
      borrow = 0_8

      ! Subtract common length portion.
      DO i = 1, len_op2
        IF (op1(i) < op2(i) + borrow) THEN
          result(i) = HIGH_PRECISION_BASE + op1(i) - op2(i) - borrow
          borrow = 1_8
        ELSE
          result(i) = op1(i) - op2(i) - borrow
          borrow = 0_8
        END IF
      END DO

      ! Propagate borrow through remaining digits of the larger number.
      DO i = len_op2 + 1, len_op1
        IF (op1(i) < borrow) THEN
          result(i) = HIGH_PRECISION_BASE + op1(i) - borrow
          borrow = 1_8
        ELSE
          result(i) = op1(i) - borrow
          borrow = 0_8
        END IF
      END DO

      ! Trim leading zeros.
      last_idx = len_op1
      DO WHILE (last_idx > 1 .AND. result(last_idx) == 0_8)
        last_idx = last_idx - 1
      END DO

      ! Reallocate to the correct size if trimming occurred.
      IF (last_idx < len_op1) THEN
        result = result(1:last_idx)
      END IF
    END SUBROUTINE subtract_core

  END SUBROUTINE abs_subtract

  ! --- Operator Implementations ---

  ! Checks whether two HighPrecisionInt numbers are equal (a == b).
  FUNCTION hpi_equal(a, b) RESULT(res)
    TYPE(high_precision_int(*)), INTENT(IN) :: a, b
    LOGICAL                                 :: res ! Result of the comparison
    INTEGER                                 :: i   ! Loop index

    IF (a%sign /= b%sign) THEN
      res = .FALSE.
      RETURN
    END IF
    IF (a%sign == 0_1) THEN ! Both are zero (b%sign must also be 0), so they are equal.
      res = .TRUE.
      RETURN
    END IF
    IF (a%ncoeffs /= b%ncoeffs) THEN ! Different logical lengths means different magnitudes.
      res = .FALSE.
      RETURN
    END IF

    res = .TRUE.
    DO i = 1, a%ncoeffs
      IF (a%coeffs(i) /= b%coeffs(i)) THEN ! First differing coefficient determines inequality.
        res = .FALSE.
        RETURN
      END IF
    END DO
  END FUNCTION hpi_equal

  ! Checks if a HighPrecisionInt is zero.
  FUNCTION hpi_is_zero(hpi) RESULT(is_zero)
    TYPE(high_precision_int(*)), INTENT(IN) :: hpi
    LOGICAL                                 :: is_zero

    is_zero = (hpi%sign == 0_1)
  END FUNCTION hpi_is_zero

  ! Compares two HighPrecisionInt numbers for less than (a < b).
  FUNCTION hpi_less(a, b) RESULT(res)
    TYPE(high_precision_int(*)), INTENT(IN) :: a, b
    LOGICAL                                 :: res  ! Result of the comparison

    LOGICAL                              :: flip ! Flag to flip comparison logic for negative numbers
    INTEGER                              :: na, nb ! Number of coefficients
    INTEGER                              :: i      ! Loop index
    INTEGER(KIND=8)                      :: c1, c2 ! Current coefficients being compared

    ! 1. Sign Comparison: Different signs immediately determine the result.
    IF (a%sign /= b%sign) THEN
      res = (a%sign < b%sign)
      RETURN
    ELSE IF (a%sign == 0_1) THEN ! Both are zero, they are equal, not less than.
      res = .FALSE.
      RETURN
    END IF

    ! 'flip' is true if numbers are negative, indicating magnitude comparison should be reversed.
    flip = (a%sign == -1_1)

    na = a%ncoeffs
    nb = b%ncoeffs
    ! 2. Magnitude Comparison (different lengths).
    IF (na /= nb) THEN
      ! IEOR (bitwise XOR) handles the flip logic:
      ! (na < nb) is the base comparison. If 'flip' is true (negative numbers), invert the result.
      res = IEOR(flip, (na < nb))
      RETURN
    END IF

    ! 2. Magnitude Comparison (equal lengths).
    ! Compare coefficients from most significant (rightmost in coeffs array) downwards.
    DO i = na, 1, -1
      c1 = a%coeffs(i)
      c2 = b%coeffs(i)
      IF (c1 /= c2) THEN
        res = IEOR(flip, (c1 < c2))
        RETURN
      END IF
    END DO

    res = .FALSE. ! If all coefficients are identical, a=b, so a is not less than b.
  END FUNCTION hpi_less

  ! Unary negation operator for HighPrecisionInt.
  FUNCTION hpi_unary_negate(hpi_in) RESULT(hpi_out)
    TYPE(high_precision_int(*)), INTENT(IN) :: hpi_in
    TYPE(high_precision_int(:)), ALLOCATABLE :: hpi_out

    IF (hpi_in%sign == 0_1) THEN
      hpi_out = new_hpi_from_integer(0_8) ! Negating zero is zero
    ELSE
      hpi_out = new_hpi_from_coeffs(hpi_in%coeffs(1:hpi_in%ncoeffs), -hpi_in%sign) ! Flip the sign
    END IF
  END FUNCTION hpi_unary_negate

  ! Adds two HighPrecisionInt numbers.
  ! Handles same-sign addition and different-sign subtraction efficiently.
  FUNCTION hpi_add(a, b) RESULT(hpi_sum)
    TYPE(high_precision_int(*)), INTENT(IN) :: a, b
    TYPE(high_precision_int(:)), ALLOCATABLE :: hpi_sum

    INTEGER(KIND=8), PARAMETER :: MASK32 = INT(Z'FFFFFFFF', KIND=8)

    INTEGER(KIND=8)                       :: current_sum, carry       ! For sums and carries
    INTEGER                               :: max_ncoeffs, i           ! Loop variables
    INTEGER(KIND=8), ALLOCATABLE          :: result_coeffs_raw(:)     ! Raw coefficients before normalization
    INTEGER(KIND=8)                       :: val_a, val_b             ! Current coefficient values from a and b
    LOGICAL                               :: is_negative_diff         ! Flag from abs_subtract
    INTEGER(KIND=1)                       :: final_sign               ! Final sign of the result
    INTEGER                               :: na, nb                   ! Number of coefficients

    ! 1. Zero Check: If either operand is zero, return the other operand.
    IF (a%sign == 0_1) THEN
      hpi_sum = b ! Allocatable assignment handles copying
      RETURN
    END IF
    IF (b%sign == 0_1) THEN
      hpi_sum = a ! Allocatable assignment handles copying
      RETURN
    END IF

    ! 2. Same Sign Addition (add magnitudes).
    IF (a%sign == b%sign) THEN
      na = a%ncoeffs
      nb = b%ncoeffs
      max_ncoeffs = MAX(na, nb)
      ALLOCATE(result_coeffs_raw(max_ncoeffs + 1)) ! Allocate space for potential final carry
      result_coeffs_raw = 0_8 ! Initialize with zeros
      
      carry = 0_8

      DO i = 1, max_ncoeffs
        IF (i <= na) THEN
            val_a = a%coeffs(i)
        ELSE
            val_a = 0_8
        END IF

        IF (i <= nb) THEN
            val_b = b%coeffs(i)
        ELSE
            val_b = 0_8
        END IF

        current_sum = val_a + val_b + carry

        result_coeffs_raw(i) = IAND(current_sum, MASK32) ! Extract lower 32 bits
        carry = ISHFT(current_sum, -32)                  ! Extract upper 32 bits (right shift by 32)
      END DO

      ! Handle any final carry.
      IF (carry > 0_8) THEN
        result_coeffs_raw(max_ncoeffs + 1) = carry
      ELSE
        ! If no final carry, resize array to max_ncoeffs to avoid unnecessary trailing zero.
        result_coeffs_raw = result_coeffs_raw(1:max_ncoeffs)
      END IF
      
      hpi_sum = new_hpi_from_coeffs(result_coeffs_raw, a%sign)
      DEALLOCATE(result_coeffs_raw)

    ! 3. Different Sign Subtraction (subtract magnitudes).
    ELSE IF (a%sign == 1_1 .AND. b%sign == -1_1) THEN ! a is positive, b is negative (a + (-|b|))
      CALL abs_subtract(a%coeffs(1:a%ncoeffs), b%coeffs(1:b%ncoeffs), result_coeffs_raw, is_negative_diff)
      IF (is_negative_diff) THEN
          final_sign = -1_1
      ELSE
          final_sign =  1_1
      END IF

      hpi_sum = new_hpi_from_coeffs(result_coeffs_raw, final_sign)
      DEALLOCATE(result_coeffs_raw)

    ELSE ! a is negative, b is positive ((-|a|) + b)
      CALL abs_subtract(b%coeffs(1:b%ncoeffs), a%coeffs(1:a%ncoeffs), result_coeffs_raw, is_negative_diff) ! Subtract |a| from |b|
      IF (is_negative_diff) THEN
          final_sign = -1_1
      ELSE
          final_sign =  1_1
      END IF
      hpi_sum = new_hpi_from_coeffs(result_coeffs_raw, final_sign)
      DEALLOCATE(result_coeffs_raw)
    END IF

  END FUNCTION hpi_add

  ! Subtraction operator for HighPrecisionInt. Implemented as a + (-b).
  FUNCTION hpi_subtract(a, b) RESULT(hpi_diff)
    TYPE(high_precision_int(*)), INTENT(IN) :: a, b
    TYPE(high_precision_int(:)), ALLOCATABLE :: hpi_diff

    ! Subtraction is implemented by adding the negation of the second operand.
    hpi_diff = a + (-b)
  END FUNCTION hpi_subtract

  ! Multiplies two HighPrecisionInt numbers using long multiplication.
  FUNCTION hpi_multiply(a, b) RESULT(hpi_prod)
    TYPE(high_precision_int(*)), INTENT(IN) :: a, b
    TYPE(high_precision_int(:)), ALLOCATABLE :: hpi_prod

    INTEGER(KIND=8), PARAMETER :: MASK32 = INT(Z'FFFFFFFF', KIND=8)

    INTEGER                              :: len_a, len_b, result_len ! Lengths of coefficient arrays and result
    INTEGER(KIND=8), ALLOCATABLE         :: result(:)                ! Result coefficients
    INTEGER                              :: i, j, k                  ! Loop indices
    INTEGER(KIND=8)                      :: ai, prod, lo, hi, carry  ! Intermediate calculation variables

    ! 1. Handle Zero: If either operand is zero, the product is zero.
    IF (hpi_is_zero(a) .OR. hpi_is_zero(b)) THEN
      hpi_prod = new_hpi_from_integer(0_8)
      RETURN
    END IF

    len_a = a%ncoeffs
    len_b = b%ncoeffs
    result_len = len_a + len_b ! Maximum possible length of the product
    ALLOCATE(result(result_len))
    result = 0_8 ! Initialize result coefficients to zero

    ! 2. Multiply and Accumulate with Carry.
    ! Nested loops for multiplying each coefficient of 'a' by each coefficient of 'b'.
    DO i = 1, len_a
      ai = a%coeffs(i)
      carry = 0_8
      DO j = 1, len_b
        k = i + j - 1 ! Current index in the result array
        
        ! Calculate partial product: a_i * b_j + existing_result[k] + carry_from_prev_digit
        prod = ai * b%coeffs(j) + result(k) + carry
        
        lo = IAND(prod, MASK32)     ! Lower 32 bits of the product (current digit)
        hi = ISHFT(prod, -32)       ! Upper 32 bits (carry to the next position)

        result(k) = lo
        carry = hi
      END DO
      ! Propagate any final carry from the current row (ai * b) to higher positions.
      k = i + len_b
      DO WHILE (carry /= 0_8)
        ! This loop extends the result array if the product leads to more digits than initially allocated.
        ! In general, result_len = len_a + len_b should be sufficient, but defensive check.
        IF (k > SIZE(result)) THEN
            ! If result array needs to be dynamically grown, reallocate.
            result = RESHAPE(result, (/k/), ORDER=(/1/))
            result(k) = 0_8 ! Initialize new element
        END IF
        
        prod = result(k) + carry
        result(k) = IAND(prod, MASK32)
        carry = ISHFT(prod, -32)
        k = k + 1
      END DO
    END DO

    ! 3. Finalize: Trim leading zeros from the result.
    DO WHILE (SIZE(result) > 1 .AND. result(SIZE(result)) == 0_8)
      result = result(1:SIZE(result) - 1)
    END DO

    hpi_prod = new_hpi_from_coeffs(result, a%sign * b%sign) ! Set final sign
    DEALLOCATE(result)

  END FUNCTION hpi_multiply

  ! Displays a HighPrecisionInt to the console.
  ! It converts the HighPrecisionInt to an INTEGER(KIND=8) for decimal representation
  ! (if it fits) and also shows its internal coefficient representation.
  SUBROUTINE print_hpi(hpi)
    TYPE(high_precision_int(*)), INTENT(IN) :: hpi
    CHARACTER(LEN=256)                   :: coeffs_str        ! String for coefficients display
    INTEGER                              :: i                 ! Loop index
    CHARACTER(LEN=20)                    :: coeff_val_str     ! Temporary string for individual coefficient
    CHARACTER(LEN=:), ALLOCATABLE        :: decimal_str       ! Full decimal representation

    ! Build the coefficient string representation
    coeffs_str = "["
    DO i = 1, hpi%ncoeffs
      WRITE(coeff_val_str, '(I0)') hpi%coeffs(i) ! Convert coefficient to string
      coeffs_str = TRIM(coeffs_str) // TRIM(ADJUSTL(coeff_val_str))
      IF (i < hpi%ncoeffs) THEN
        coeffs_str = TRIM(coeffs_str) // ", "
      END IF
    END DO
    coeffs_str = TRIM(coeffs_str) // "]"

    ! Get the full decimal representation
    decimal_str = hpi_to_string(hpi)
    
    ! Print the final formatted string including the sign
    WRITE(*, '(A, A, A, A, A, I0, A)') "HighPrecisionInt(", TRIM(decimal_str), ", coeffs=", TRIM(coeffs_str), ", sign=", hpi%sign, ")"
  END SUBROUTINE print_hpi

  ! Scales a HighPrecisionInt by HIGH_PRECISION_BASE^power (for power >= 0).
  ! Prepends 'power' zeros to the 'coeffs' vector. Returns a new HighPrecisionInt.
  FUNCTION hpi_scale_up_by_base_power(hpi_in, power) RESULT(hpi_out)
    TYPE(high_precision_int(*)), INTENT(IN) :: hpi_in
    INTEGER, INTENT(IN)                  :: power
    TYPE(high_precision_int(:)), ALLOCATABLE :: hpi_out
    INTEGER(KIND=8), ALLOCATABLE         :: new_coeffs(:)
    INTEGER                              :: i, old_ncoeffs

    IF (power < 0) THEN
      STOP "Error: Scaling power must be non-negative for hpi_scale_up_by_base_power."
    END IF
    IF (power == 0) THEN
      hpi_out = hpi_in ! Allocatable assignment creates a copy
      RETURN
    END IF
    IF (hpi_in%sign == 0_1) THEN
      hpi_out = new_hpi_from_integer(0_8)
      RETURN
    END IF

    old_ncoeffs = hpi_in%ncoeffs
    ALLOCATE(new_coeffs(old_ncoeffs + power))
    new_coeffs = 0_8 ! Initialize with zeros

    ! Copy original coefficients to the higher indices
    DO i = 1, old_ncoeffs
      new_coeffs(i + power) = hpi_in%coeffs(i + 1 - 1)
    END DO
    hpi_out = new_hpi_from_coeffs(new_coeffs, hpi_in%sign)
    DEALLOCATE(new_coeffs)
  END FUNCTION hpi_scale_up_by_base_power

  ! Calculates base^exp for non-negative integer exponents.
  FUNCTION hpi_power(base, exp) RESULT(res)
    TYPE(high_precision_int(*)), INTENT(IN) :: base
    INTEGER, INTENT(IN)                  :: exp
    TYPE(high_precision_int(:)), ALLOCATABLE :: res
    INTEGER                              :: n
    TYPE(high_precision_int(:)), ALLOCATABLE :: p

    IF (exp < 0) STOP "hpi_power: exponent must be non-negative."
    IF (exp == 0) THEN
      res = new_hpi_from_integer(1_8)
      RETURN
    END IF
    IF (hpi_is_zero(base)) THEN
      res = new_hpi_from_integer(0_8)
      RETURN
    END IF
    IF (hpi_equal(base, new_hpi_from_integer(1_8))) THEN
      res = new_hpi_from_integer(1_8)
      RETURN
    END IF

    ! Exponentiation by squaring
    res = new_hpi_from_integer(1_8)
    p = base
    n = exp
    DO WHILE (n > 0)
      IF (MOD(n, 2) == 1) THEN
        res = res * p
      END IF
      p = p * p
      n = n / 2
    END DO
  END FUNCTION hpi_power

  ! Performs integer division, providing both quotient and remainder.
  SUBROUTINE hpi_div_rem(numerator, denominator, quotient, remainder)
    TYPE(high_precision_int(*)), INTENT(IN)  :: numerator, denominator
    TYPE(high_precision_int(:)), ALLOCATABLE, INTENT(OUT) :: quotient, remainder

    TYPE(high_precision_int(:)), ALLOCATABLE :: num_abs, den_abs, q_abs, r_abs
    INTEGER(KIND=1)                          :: q_sign, r_sign

    IF (hpi_is_zero(denominator)) STOP "hpi_div_rem: Division by zero."

    IF (hpi_is_zero(numerator)) THEN
      quotient = new_hpi_from_integer(0_8)
      remainder = new_hpi_from_integer(0_8)
      RETURN
    END IF

    num_abs = hpi_abs(numerator)
    den_abs = hpi_abs(denominator)

    IF (num_abs < den_abs) THEN
      quotient = new_hpi_from_integer(0_8)
      remainder = numerator ! Allocatable assignment
      RETURN
    END IF

    CALL div_rem_magnitude(num_abs, den_abs, q_abs, r_abs)

    q_sign = numerator%sign * denominator%sign
    r_sign = numerator%sign ! Remainder has the sign of the numerator

    quotient = new_hpi_from_coeffs(q_abs%coeffs(1:q_abs%ncoeffs), q_sign)
    remainder = new_hpi_from_coeffs(r_abs%coeffs(1:r_abs%ncoeffs), r_sign)
    DEALLOCATE(num_abs, den_abs, q_abs, r_abs)

  CONTAINS
    ! Internal subroutine for long division of positive integers.
    SUBROUTINE div_rem_magnitude(num, den, q, r)
      TYPE(high_precision_int(*)), INTENT(IN)  :: num, den ! Assumes num > 0, den > 0
      TYPE(high_precision_int(:)), ALLOCATABLE, INTENT(OUT) :: q, r      ! q = num / den, r = num % den

      INTEGER(KIND=8), PARAMETER :: HIGH_PRECISION_BASE = 2_8**32

      TYPE(high_precision_int(:)), ALLOCATABLE :: current_rem, shifted_den, test_prod
      INTEGER(KIND=8), ALLOCATABLE :: q_coeffs(:)
      INTEGER :: len_diff, i
      INTEGER(KIND=8) :: q_digit, low, high, mid

      len_diff = num%ncoeffs - den%ncoeffs
      ALLOCATE(q_coeffs(len_diff + 1)); q_coeffs = 0_8
      current_rem = num

      DO i = len_diff, 0, -1
        shifted_den = hpi_scale_up_by_base_power(den, i)

        ! Performance Optimization: The original repeated subtraction is replaced with a binary search.
        ! This is dramatically faster when the quotient digit is large.
        IF (current_rem < shifted_den) THEN
            q_digit = 0_8
        ELSE
            ! Binary search for the quotient digit `q_digit` such that
            ! q_digit * shifted_den <= current_rem < (q_digit + 1) * shifted_den
            low = 1_8 ! We know q_digit is at least 1 from the check above.
            high = HIGH_PRECISION_BASE - 1_8
            q_digit = 1_8 ! Smallest possible non-zero guess.

            DO WHILE (low <= high)
                mid = low + ISHFT(high - low, -1) ! mid = (low+high)/2 to avoid overflow
                test_prod = shifted_den * new_hpi_from_integer(mid)
                IF (current_rem < test_prod) THEN
                    high = mid - 1
                ELSE
                    q_digit = mid ! This is a better potential answer, try for a larger one.
                    low = mid + 1
                END IF
            END DO
        END IF

        q_coeffs(i + 1) = q_digit

        ! Update remainder: current_rem = current_rem - q_digit * shifted_den
        IF (q_digit > 0_8) THEN
            current_rem = current_rem - (shifted_den * new_hpi_from_integer(q_digit))
        END IF
      END DO

      q = new_hpi_from_coeffs(q_coeffs, 1_1)
      r = current_rem
      DEALLOCATE(q_coeffs)
    END SUBROUTINE div_rem_magnitude
  END SUBROUTINE hpi_div_rem

  ! Performs integer division, returning only the quotient.
  FUNCTION hpi_divide(a, b) RESULT(q)
    TYPE(high_precision_int(*)), INTENT(IN) :: a, b
    TYPE(high_precision_int(:)), ALLOCATABLE :: q
    TYPE(high_precision_int(:)), ALLOCATABLE :: r ! Dummy remainder

    CALL hpi_div_rem(a, b, q, r)
  END FUNCTION hpi_divide

END MODULE high_precision_integer_mod