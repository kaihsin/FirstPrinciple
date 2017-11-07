!
PROGRAM ONE_D_EULER
!
! Euler algorithm for the motion of a particle subject to an external
! force f(x) = -x. An equal time step DT is used. The position, 
! velocity and energy of the particle are written out.
!
  IMPLICIT NONE
  INTEGER, PARAMETER :: N=80001,IN= 80
  INTEGER :: I
  REAL*8 :: PI,DT
  REAL*8, DIMENSION (N):: T,V,X,E
!
! Assign constants, initial position, and initial velocity
!
  PI   = 4.D0*DATAN(1.D0)
  DT   = 0.0025D0*PI
  X(1) = 0.D0
  T(1) = 0.D0
  V(1) = 1.D0
!
! Recursion for position and velocity at later time
!
  DO I = 1, N-1
    T(I+1) = DT*I
    X(I+1) = X(I)+V(I)*DT
    V(I+1) = V(I)-X(I)*DT
  END DO
  DO I = 1, N
    E(I) = 0.5D0*(V(I)**2+X(I)**2)
  END DO
!
! Write the position, velocity and energy every IN steps
!
  WRITE (6,"(4F18.8)") (T(I)/PI,X(I),V(I),E(I),I=1,N,IN)
END PROGRAM ONE_D_EULER 
