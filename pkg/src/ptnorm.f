      SUBROUTINE MYUNORM(PROB, Z, LENGTH)
      DOUBLE PRECISION PROB(*), Z(*)
      DOUBLE PRECISION UNORM
      INTEGER LENGTH
      DO I = 1, LENGTH, 1
         PROB(I) = UNORM(Z(I))
      END DO
      RETURN
      END

      SUBROUTINE MYBNORM(PROB, H1, H2, CORREL, LENGTH)
      DOUBLE PRECISION PROB(*), H1(*), H2(*), CORREL
      DOUBLE PRECISION BNORM
      INTEGER LENGTH
      DO I = 1, LENGTH, 1
         PROB(I) = BNORM(H1(I), H2(I), CORREL)
      END DO
      RETURN
      END


      SUBROUTINE MYTNORM(PROB, H1, H2, H3, CORREL, LENGTH)
      DOUBLE PRECISION PROB(*), H1(*), H2(*), H3(*), CORREL(3)
      DOUBLE PRECISION TNORM
      INTEGER LENGTH
      DO I = 1, LENGTH, 1
         PROB(I) = TNORM(H1(I), H2(I), H3(I), CORREL)
      END DO
      RETURN
      END
      
      
      DOUBLE PRECISION FUNCTION UNORM(X)
*     
*     Normal distribution probabilities accurate to 1d-15.
*     Reference: J.L. Schonfelder, Math Comp 32(1978), pp 1232-1240. 
*
      DOUBLE PRECISION X, XABS, EXPONENTIAL, BUILD, CUMNORM
      XABS = ABS(X)
      IF (ABSX .GT. 37) THEN
         CUMNORM = 0
      ELSE
         EXPONENTIAL = EXP(- XABS ** 2 / 2)
         IF (XABS .LT. 7.07106781186547) THEN
            BUILD = 3.52624965998911E-2 * XABS + 0.700383064443688
            BUILD = BUILD * XABS + 6.37396220353165
            BUILD = BUILD * XABS + 33.912866078383
            BUILD = BUILD * XABS + 112.079291497871
            BUILD = BUILD * XABS + 221.213596169931
            BUILD = BUILD * XABS + 220.206867912376
            CUMNORM = EXPONENTIAL * BUILD
            BUILD = 8.83883476483184E-02 * XABS + 1.75566716318264
            BUILD = BUILD * XABS + 16.064177579207
            BUILD = BUILD * XABS + 86.7807322029461
            BUILD = BUILD * XABS + 296.564248779674
            BUILD = BUILD * XABS + 637.333633378831
            BUILD = BUILD * XABS + 793.826512519948
            BUILD = BUILD * XABS + 440.413735824752
            CUMNORM = CUMNORM / BUILD
         ELSE
            BUILD = XABS + 0.65
            BUILD = XABS + 4 / BUILD
            BUILD = XABS + 3 / BUILD
            BUILD = XABS + 2 / BUILD
            BUILD = XABS + 1 / BUILD
            CUMNORM = EXPONENTIAL / BUILD / 2.506628274631
         END IF
      END IF
      IF (X .GT. O) THEN
         CUMNORM = 1 - CUMNORM
      END IF
      UNORM = CUMNORM
      END
      
      DOUBLE PRECISION FUNCTION BNORM(H1, H2, R)
      DOUBLE PRECISION R, X(5), W(5)
      DOUBLE PRECISION H1, H2, H3, H12, R1, R2, RR
      DOUBLE PRECISION H8, R3, H7, H5, H6, AA, AB, LH
      DOUBLE PRECISION BCUM, UNORM
      INTEGER I
      DATA X/
     &     0.04691008,
     &     0.23076534,
     &     0.5,
     &     0.76923466,
     &     0.95308992/
      DATA W/
     &     0.018854042,
     &     0.038088059,
     &     0.0452707394,
     &     0.038088059,
     &     0.018854042/

      LH = 0.
      H12 = (H1 ** 2  + H2 ** 2) / 2
      IF (ABS(R) .GE. 0.7) THEN
         R2 = 1 - R ** 2
         R3 = SQRT(R2)
         IF (R .LT. 0) H2 = - H2
         H3 = H1 * H2
         H7 = EXP(- H3 / 2)
         IF (ABS(R) .LT. 1) THEN
            H6 = ABS(H1 - H2)
            H5 = H6 ** 2 / 2
            H6 = H6 / R3
            AA = 0.5 - H3 / 8
            AB = 3 - 2 * AA * H5
            LH = 0.13298076 * H6 * AB * (1 - UNORM(H6)) -
     &           EXP(- H5 / R2) * (AB + AA * R2) * 0.053051647
            DO I = 1, 5
               R1 = R3 * X(I)
               RR = R1 ** 2
               R2 = SQRT(1 - RR)
               IF (H7 .EQ. 0) THEN
                  H8 = 0
               ELSE
                  H8 = EXP(- H3 / (1 + R2)) / R2 / H7
               END IF
               LH = LH - W(I) * EXP(- H5 / RR) * (H8 - 1 - AA * RR)
            END DO
         END IF
         BCUM = LH * R3 * H7 + UNORM(MIN(H1, H2))
         IF (R .LT. 0) THEN
            BCUM = UNORM(H1) - BCUM
         END IF
      ELSE
         H3 = H1 * H2
         IF (R .NE. 0) THEN
            DO I = 1, 5
               R1 = R * X(I)
               R2 = 1 - R1 ** 2
               LH = LH + W(I) * EXP((R1 * H3 - H12) / R2) / SQRT(R2)
            END DO
         END IF
         BCUM = UNORM(H1) * UNORM(H2) + R * LH
      END IF
      BNORM = BCUM
      END      

      DOUBLE PRECISION FUNCTION TNORM(H1, H2, H3, R)
      DOUBLE PRECISION X(5), W(5), R(3)
      DOUBLE PRECISION H1, H2, H3, R12, R13, R23
      DOUBLE PRECISION DEL, FAC
      DOUBLE PRECISION HP1, HP2, Z12, Z13, XI, XIS
      DOUBLE PRECISION UNORM, BNORM
      INTEGER I
      DATA X/
     &     0.04691008,
     &     0.23076534,
     &     0.5,
     &     0.76923466,
     &     0.95308992/
      DATA W/
     &     0.018854042,
     &     0.038088059,
     &     0.0452707394,
     &     0.038088059,
     &     0.018854042/

      R12 = R(1)
      R13 = R(2)
      R23 = R(3)
      TNORM= 0.0

      DO I = 1, 5
         XI = X(I)
         XIS = XI ** 2
         Z12 = EXP(- 0.5 * (H1 ** 2 + H2 ** 2 - 
     &        2 * XI * R12 * H1 * H2) / 
     &        (1. - XIS * R12 ** 2)) / SQRT(1. - XIS * R12 ** 2)
         Z13 = EXP(- 0.5 * (H1 ** 2 + H3 ** 2 - 
     &        2 * XI * R13 * H1 * H3) / 
     &        (1. - XIS * R13 ** 2)) / SQRT(1. - XIS * R13 ** 2)
         DEL = 1. - XIS * R12 ** 2 - XIS * R13 ** 2 - 
     &        R23 ** 2 + 2 * XIS * R12 * R13 * R23
         FAC = SQRT(DEL)
         HP1 = (H3 * (1. - XIS * R12 ** 2) - 
     &        H1 * (XI * R13 - XI * R12 * R23)  - 
     &        H2 * (R23 - XIS * R12 * R13)) / FAC / 
     &        SQRT(1. - XIS * R12 ** 2)
         HP2 = (H2 * (1. - XIS * R13 ** 2) -
     &        H1 * (XI * R12 - XI * R13 * R23) - 
     &        H3 * (R23 - XIS * R12 * R13)) / FAC /
     &        SQRT(1. - XIS * R13 ** 2)
         TNORM = TNORM + W(I) * Z12 * UNORM(HP1) * R12
         TNORM = TNORM + W(I) * Z13 * UNORM(HP2) * R13
      END DO
      TNORM = TNORM + UNORM(H1) * BNORM(H2, H3, R23)
      RETURN
      END
