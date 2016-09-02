    MODULE mod_space

      USE felib_globals, ONLY : wp

      IMPLICIT NONE

      INTRINSIC associated, present, size

      INTEGER, PRIVATE :: allocate_count
      INTEGER, PRIVATE :: deallocate_count
      INTEGER, PRIVATE :: i

      LOGICAL debug

      INTEGER, PARAMETER :: maximum_number_complex_arrays = 20
      LOGICAL, SAVE , DIMENSION (maximum_number_complex_arrays) :: &
        which_complex_array_used
      DATA (which_complex_array_used(i),i=1,maximum_number_complex_arrays)/ &
        maximum_number_complex_arrays* .FALSE./

      INTEGER, PARAMETER :: maximum_number_complex_vectors = 20
      LOGICAL, SAVE , DIMENSION (maximum_number_complex_vectors) :: &
        which_complex_vector_used
      DATA (which_complex_vector_used(i),i=1,maximum_number_complex_vectors)/ &
        maximum_number_complex_vectors* .FALSE./

      INTEGER, PARAMETER :: maximum_number_real_arrays = 20
      LOGICAL, SAVE , DIMENSION (maximum_number_real_arrays) :: &
        which_real_array_used
      DATA (which_real_array_used(i),i=1,maximum_number_real_arrays)/ &
        maximum_number_real_arrays* .FALSE./

      INTEGER, PARAMETER :: maximum_number_real_vectors = 20
      LOGICAL, SAVE , DIMENSION (maximum_number_real_vectors) :: &
        which_real_vector_used
      DATA (which_real_vector_used(i),i=1,maximum_number_real_vectors)/ &
        maximum_number_real_vectors* .FALSE./

      INTEGER, PARAMETER :: maximum_number_integer_arrays = 20
      LOGICAL, SAVE , DIMENSION (maximum_number_integer_arrays) :: &
        which_integer_array_used
      DATA (which_integer_array_used(i),i=1,maximum_number_integer_arrays)/ &
        maximum_number_integer_arrays* .FALSE./

      INTEGER, PARAMETER :: maximum_number_integer_vectors = 20
      LOGICAL, SAVE , DIMENSION (maximum_number_integer_vectors) :: &
        which_integer_vector_used
      DATA (which_integer_vector_used(i),i=1,maximum_number_integer_vectors)/ &
        maximum_number_integer_vectors* .FALSE./

      TYPE complex_array
        COMPLEX (wp), POINTER, DIMENSION (:,:) :: array
      END TYPE complex_array

      TYPE complex_vector
        COMPLEX (wp), POINTER, DIMENSION (:) :: vector
      END TYPE complex_vector

      TYPE real_array
        REAL (wp), POINTER, DIMENSION (:,:) :: array
      END TYPE real_array

      TYPE real_vector
        REAL (wp), POINTER, DIMENSION (:) :: vector
      END TYPE real_vector

      TYPE integer_array
        INTEGER, POINTER, DIMENSION (:,:) :: array
      END TYPE integer_array

      TYPE integer_vector
        INTEGER, POINTER, DIMENSION (:) :: vector
      END TYPE integer_vector

      TYPE (complex_array), DIMENSION (maximum_number_complex_arrays) :: ca
      TYPE (complex_vector), DIMENSION (maximum_number_complex_vectors) :: cv
      TYPE (real_array), DIMENSION (maximum_number_real_arrays) :: ra
      TYPE (real_vector), DIMENSION (maximum_number_real_vectors) :: rv
      TYPE (integer_array), DIMENSION (maximum_number_integer_arrays) :: ia
      TYPE (integer_vector), DIMENSION (maximum_number_integer_vectors) :: iv


      INTERFACE create
        MODULE PROCEDURE create_complex_array
        MODULE PROCEDURE create_complex_vector
        MODULE PROCEDURE create_real_array
        MODULE PROCEDURE create_real_vector
        MODULE PROCEDURE create_integer_array
        MODULE PROCEDURE create_integer_vector
      END INTERFACE

      INTERFACE destroy
        MODULE PROCEDURE destroy_complex_array
        MODULE PROCEDURE destroy_complex_vector
        MODULE PROCEDURE destroy_real_array
        MODULE PROCEDURE destroy_real_vector
        MODULE PROCEDURE destroy_integer_array
        MODULE PROCEDURE destroy_integer_vector
      END INTERFACE

    CONTAINS

      SUBROUTINE destroy_complex_array(array)
        IMPLICIT NONE
        COMPLEX (wp), POINTER, DIMENSION (:,:) :: array
        INTEGER :: next_array, error

        DO next_array = 1, maximum_number_complex_arrays
          IF (associated(ca(next_array)%array,array)) THEN
            which_complex_array_used(next_array) = .FALSE.
            DEALLOCATE (ca(next_array)%array,STAT=error)
            NULLIFY (array)
            IF (error==0) THEN
              deallocate_count = deallocate_count + 1
            END IF
          END IF
        END DO

      END SUBROUTINE destroy_complex_array

      SUBROUTINE destroy_complex_vector(vector)
        IMPLICIT NONE
        COMPLEX (wp), POINTER, DIMENSION (:) :: vector

        INTEGER :: next_vector, error

        DO next_vector = 1, maximum_number_complex_vectors
          IF (associated(cv(next_vector)%vector,vector)) THEN
            which_complex_vector_used(next_vector) = .FALSE.
            DEALLOCATE (cv(next_vector)%vector,STAT=error)
            NULLIFY (vector)
            IF (error==0) THEN
              deallocate_count = deallocate_count + 1
            END IF
          END IF
        END DO

      END SUBROUTINE destroy_complex_vector

      SUBROUTINE destroy_real_array(array)
        IMPLICIT NONE
        REAL (wp), POINTER, DIMENSION (:,:) :: array
        INTEGER :: next_array, error

        DO next_array = 1, maximum_number_real_arrays
          IF (associated(ra(next_array)%array,array)) THEN
            which_real_array_used(next_array) = .FALSE.
            DEALLOCATE (ra(next_array)%array,STAT=error)
            NULLIFY (array)
            IF (error==0) THEN
              deallocate_count = deallocate_count + 1
            END IF
          END IF
        END DO

      END SUBROUTINE destroy_real_array

      SUBROUTINE destroy_real_vector(vector)
        IMPLICIT NONE
        REAL (wp), POINTER, DIMENSION (:) :: vector

        INTEGER :: next_vector, error

        DO next_vector = 1, maximum_number_real_vectors
          IF (associated(rv(next_vector)%vector,vector)) THEN
            which_real_vector_used(next_vector) = .FALSE.
            DEALLOCATE (rv(next_vector)%vector,STAT=error)
            NULLIFY (vector)
            IF (error==0) THEN
              deallocate_count = deallocate_count + 1
            END IF
          END IF
        END DO

      END SUBROUTINE destroy_real_vector

      SUBROUTINE destroy_integer_array(array)
        IMPLICIT NONE
        INTEGER, POINTER, DIMENSION (:,:) :: array
        INTEGER :: next_array, error

        DO next_array = 1, maximum_number_integer_arrays
          IF (associated(ia(next_array)%array,array)) THEN
            which_integer_array_used(next_array) = .FALSE.
            DEALLOCATE (ia(next_array)%array,STAT=error)
            NULLIFY (array)
            IF (error==0) THEN
              deallocate_count = deallocate_count + 1
            END IF
          END IF
        END DO

      END SUBROUTINE destroy_integer_array

      SUBROUTINE destroy_integer_vector(vector)
        IMPLICIT NONE
        INTEGER, POINTER, DIMENSION (:) :: vector
        INTEGER :: next_vector, error

        DO next_vector = 1, maximum_number_integer_vectors
          IF (associated(iv(next_vector)%vector,vector)) THEN
            which_integer_vector_used(next_vector) = .FALSE.
            DEALLOCATE (iv(next_vector)%vector,STAT=error)
            NULLIFY (vector)
            IF (error==0) THEN
              deallocate_count = deallocate_count + 1
            END IF
          END IF
        END DO

      END SUBROUTINE destroy_integer_vector

      SUBROUTINE destroy_all

        IMPLICIT NONE

        INTEGER :: next_array, next_vector, error

        DO next_array = 1, maximum_number_real_arrays
          IF (which_real_array_used(next_array)) THEN
            which_real_array_used(next_array) = .FALSE.
            IF (associated(ra(next_array)%array)) THEN
              DEALLOCATE (ra(next_array)%array,STAT=error)
              IF (error==0) THEN
                deallocate_count = deallocate_count + 1
              END IF
            END IF
          END IF
        END DO

        DO next_vector = 1, maximum_number_real_vectors
          IF (which_real_vector_used(next_vector)) THEN
            which_real_vector_used(next_vector) = .FALSE.
            IF (associated(rv(next_vector)%vector)) THEN
              DEALLOCATE (rv(next_vector)%vector,STAT=error)
              IF (error==0) THEN
                deallocate_count = deallocate_count + 1
              END IF
            END IF
          END IF
        END DO

        DO next_array = 1, maximum_number_integer_arrays
          IF (which_integer_array_used(next_array)) THEN
            which_integer_array_used(next_array) = .FALSE.
            IF (associated(ia(next_array)%array)) THEN
              DEALLOCATE (ia(next_array)%array,STAT=error)
              IF (error==0) THEN
                deallocate_count = deallocate_count + 1
              END IF
            END IF
          END IF
        END DO

        DO next_vector = 1, maximum_number_integer_vectors
          IF (which_integer_vector_used(next_vector)) THEN
            which_integer_vector_used(next_vector) = .FALSE.
            IF (associated(iv(next_vector)%vector)) THEN
              DEALLOCATE (iv(next_vector)%vector,STAT=error)
              IF (error==0) THEN
                deallocate_count = deallocate_count + 1
              END IF
            END IF
          END IF
        END DO

      END SUBROUTINE destroy_all

      SUBROUTINE check(array)

        IMPLICIT NONE

        REAL (wp), POINTER, OPTIONAL, DIMENSION (:) :: array
        INTEGER :: i, count

        PRINT *, ' '
        PRINT *, 'Check Space'
        PRINT *, ' '
        IF (present(array)) THEN
        ELSE
          PRINT *, 'DEALLOCATED_COUNT=', deallocate_count
          PRINT *, 'ALLOCATED_COUNT=', allocate_count

          PRINT *, 'Real Vectors'
          count = 0
          DO i = 1, maximum_number_real_vectors
            WRITE (*,'(i3,a1,l2)',advance='no') i, ':', &
              which_real_vector_used(i)
            IF (which_real_vector_used(i)) count = count + 1
          END DO
          PRINT *, '1'
          PRINT *, count, ' real vectors set'

          count = 0
          PRINT *, ' '
          PRINT *, 'Real Arrays'
          DO i = 1, maximum_number_real_arrays
            WRITE (*,'(i3,a1,l2)',advance='no') i, ':', &
              which_real_array_used(i)
            IF (which_real_array_used(i)) count = count + 1
          END DO
          PRINT *, ' '
          PRINT *, count, ' real arrays set'
          PRINT *, ' '

          PRINT *, 'Integer Vectors'
          count = 0
          DO i = 1, maximum_number_integer_vectors
            WRITE (*,'(i3,a1,l2)',advance='no') i, ':', &
              which_integer_vector_used(i)
            IF (which_integer_vector_used(i)) count = count + 1
          END DO
          PRINT *, ' '
          PRINT *, count, ' integer vectors set'

          count = 0
          PRINT *, ' '
          PRINT *, 'Integer Arrays'
          DO i = 1, maximum_number_integer_arrays
            WRITE (*,'(i3,a1,l2)',advance='no') i, ':', &
              which_integer_array_used(i)
            IF (which_integer_array_used(i)) count = count + 1
          END DO
          PRINT *, ' '
          PRINT *, count, ' integer arrays set'
          PRINT *, ' '

        END IF

      END SUBROUTINE check

      RECURSIVE SUBROUTINE create_complex_array(array,size1,size2)

        IMPLICIT NONE

        INTEGER :: size1, size2
        COMPLEX (wp), POINTER, DIMENSION (:,:) :: array

        INTEGER :: next_array, error

        debug = .FALSE.

! Check that array has not been allocated

        IF (associated(array)) THEN
          IF (debug) PRINT *, 'Array has association'
          IF (size(array)/=size1*size2) THEN
            IF (debug) PRINT *, 'Array incorrect size'
            CALL destroy(array)
            CALL create(array,size1,size2)
          ELSE
          END IF
        ELSE

! Scan through available arrays and select first free

          IF (debug) PRINT *, 'Array not associated'
          next_array = 1
          DO
            IF ( .NOT. which_complex_array_used(next_array)) EXIT
            next_array = next_array + 1
            IF (next_array>maximum_number_complex_arrays) THEN
              RETURN
            END IF
          END DO
          IF (debug) PRINT *, 'Allocating array no. ', next_array
          ALLOCATE (ca(next_array)%array(size1,size2),STAT=error)
          array => ca(next_array) %array
          which_complex_array_used(next_array) = .TRUE.
          IF (error==0) THEN
            allocate_count = allocate_count + 1
          END IF
        END IF

      END SUBROUTINE create_complex_array

      RECURSIVE SUBROUTINE create_complex_vector(vector,size1)

        IMPLICIT NONE

        INTEGER :: size1
        COMPLEX (wp), POINTER, DIMENSION (:) :: vector

        INTEGER :: next_vector, error

        debug = .FALSE.

        IF (associated(vector)) THEN
          IF (debug) PRINT *, 'Array has association'
          IF (size(vector)/=size1) THEN
            CALL destroy(vector)
            CALL create(vector,size1)
          ELSE
          END IF
        ELSE
          IF (debug) PRINT *, 'Array not associated'
          next_vector = 1
          DO
            IF ( .NOT. which_complex_vector_used(next_vector)) EXIT
            next_vector = next_vector + 1
            IF (next_vector>maximum_number_complex_vectors) THEN
              RETURN
            END IF
          END DO
          IF (debug) PRINT *, 'Allocating vector no. ', next_vector
          ALLOCATE (cv(next_vector)%vector(size1),STAT=error)
          vector => cv(next_vector) %vector
          which_complex_vector_used(next_vector) = .TRUE.
          IF (error==0) THEN
            allocate_count = allocate_count + 1
          END IF
        END IF

      END SUBROUTINE create_complex_vector

      RECURSIVE SUBROUTINE create_real_array(array,size1,size2)

        IMPLICIT NONE

        INTEGER :: size1, size2
        REAL (wp), POINTER, DIMENSION (:,:) :: array

        INTEGER :: next_array, error

        debug = .FALSE.

! Check that array has not been allocated

        IF (associated(array)) THEN
          IF (debug) PRINT *, 'Array has association'
          IF (size(array)/=size1*size2) THEN
            IF (debug) PRINT *, 'Array incorrect size'
            CALL destroy(array)
            CALL create(array,size1,size2)
          ELSE
          END IF
        ELSE

! Scan through available arrays and select first free

          IF (debug) PRINT *, 'Array not associated!'
          next_array = 1
          DO
            IF ( .NOT. which_real_array_used(next_array)) EXIT
            next_array = next_array + 1
            IF (next_array>maximum_number_real_arrays) THEN
              IF (debug) PRINT *, 'Maximum number of real arrays exceeded (', &
                next_array, ')'
              RETURN
            END IF
          END DO
          IF (debug) PRINT *, 'Allocating array no. ', next_array
          ALLOCATE (ra(next_array)%array(size1,size2),STAT=error)
          array => ra(next_array) %array
          which_real_array_used(next_array) = .TRUE.
          IF (error==0) THEN
            allocate_count = allocate_count + 1
          END IF
        END IF

      END SUBROUTINE create_real_array

      RECURSIVE SUBROUTINE create_real_vector(vector,size1)

        IMPLICIT NONE

        INTEGER :: size1
        REAL (wp), POINTER, DIMENSION (:) :: vector

        INTEGER :: next_vector, error

        debug = .FALSE.

        IF (associated(vector)) THEN
          IF (debug) PRINT *, 'Vector has association'
          IF (size(vector)/=size1) THEN
            CALL destroy(vector)
            CALL create(vector,size1)
          ELSE
          END IF
        ELSE
          IF (debug) PRINT *, 'Vector not associated'
          next_vector = 1
          DO
            IF ( .NOT. which_real_vector_used(next_vector)) EXIT
            next_vector = next_vector + 1
            IF (next_vector>maximum_number_real_vectors) THEN
              RETURN
            END IF
          END DO
          IF (debug) PRINT *, 'Allocating vector no. ', next_vector
          ALLOCATE (rv(next_vector)%vector(size1),STAT=error)
          vector => rv(next_vector) %vector
          which_real_vector_used(next_vector) = .TRUE.
          IF (error==0) THEN
            allocate_count = allocate_count + 1
          END IF
        END IF

      END SUBROUTINE create_real_vector

      RECURSIVE SUBROUTINE create_integer_array(array,size1,size2)

        IMPLICIT NONE

        INTEGER :: size1, size2
        INTEGER, POINTER, DIMENSION (:,:) :: array
        INTEGER :: next_array, error

! Check DEBUG option

        debug = .FALSE.

! Check that array has not been allocated

        IF (associated(array)) THEN
          IF (debug) PRINT *, 'Array associated'
          IF (size(array)/=size1*size2) THEN
            IF (debug) PRINT *, 'Array has wrong size'
            CALL destroy(array)
            CALL create(array,size1,size2)
          ELSE
          END IF
        ELSE

! Scan through available arrays and select first free

          IF (debug) PRINT *, 'Array not associated'
          next_array = 1
          DO
            IF ( .NOT. which_integer_array_used(next_array)) EXIT
            next_array = next_array + 1
            IF (next_array>maximum_number_integer_arrays) THEN
              RETURN
            END IF
          END DO
          IF (debug) PRINT *, 'Allocating array no. ', next_array
          ALLOCATE (ia(next_array)%array(size1,size2),STAT=error)
          array => ia(next_array) %array
          which_integer_array_used(next_array) = .TRUE.
          IF (error==0) THEN
            allocate_count = allocate_count + 1
          END IF
        END IF

      END SUBROUTINE create_integer_array

      RECURSIVE SUBROUTINE create_integer_vector(vector,size1)

        IMPLICIT NONE

        INTEGER :: size1
        INTEGER, POINTER, DIMENSION (:) :: vector
        INTEGER :: next_vector, error


        IF (associated(vector)) THEN
          IF (size(vector)/=size1) THEN
            CALL destroy(vector)
            CALL create(vector,size1)
          ELSE
          END IF
        ELSE
          next_vector = 1
          DO
            IF ( .NOT. which_integer_vector_used(next_vector)) EXIT
            next_vector = next_vector + 1
            IF (next_vector>maximum_number_integer_vectors) THEN
              RETURN
            END IF
          END DO
          ALLOCATE (iv(next_vector)%vector(size1),STAT=error)
          vector => iv(next_vector) %vector
          which_integer_vector_used(next_vector) = .TRUE.
          IF (error==0) THEN
            allocate_count = allocate_count + 1
          END IF
        END IF

      END SUBROUTINE create_integer_vector

    END MODULE mod_space

