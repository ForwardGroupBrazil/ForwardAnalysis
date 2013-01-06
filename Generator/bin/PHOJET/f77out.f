C FORTRAN implementation of file I/O
C
      SUBROUTINE F77OUT(STRING)
      CHARACTER*(*)     STRING
      INTEGER         IUNIT
      COMMON /FILEIO/ IUNIT
 12   FORMAT(A)
      WRITE(IUNIT,12) STRING
      END

      SUBROUTINE F77OPN(IU,NAME)
      INTEGER           IU
      CHARACTER*(*)        NAME
      INTEGER         IUNIT
      COMMON /FILEIO/ IUNIT
      CALL F77CLS()
      IF(IU.NE.6)OPEN(UNIT=IU,FILE=NAME,STATUS='UNKNOWN')
      IUNIT = IU
      END

      SUBROUTINE F77CLS()
      INTEGER         IUNIT
      COMMON /FILEIO/ IUNIT
      IF(IUNIT.GT.0)THEN
         CLOSE(IUNIT)
      ENDIF
      END

      SUBROUTINE F77RWD()
      INTEGER         IUNIT
      COMMON /FILEIO/ IUNIT
      IF(IUNIT.GT.0)THEN
         REWIND(IUNIT)
      ENDIF
      END