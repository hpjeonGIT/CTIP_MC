!
!############################################################################
SUBROUTINE POSITION_UPDATE(Nset, q, h)
USE DATAFMT
!
TYPE(AMNT):: Nset
TYPE(PTCL):: q(Nset%Npt_all)
REAL(DP)  :: h
!
INTEGER(GI):: i, id, Nall
REAL(DP):: xm, x(3), box(3), F_local, Fmax
!
!
Nall = Nset%Npt_all
Fmax = 1.D-8
DO i=1, Nall
   F_local = q(i)%ff(1)**2 + q(i)%ff(2)**2 + q(i)%ff(3)**2
   F_local = DSQRT(F_local)
   Fmax = MAX(Fmax, F_local)
END DO
!
box(:) = Nset%box(:)
DO i=1, Nset%Npt_all
   q(i)%xx(:) = q(i)%xx(:) + q(i)%ff(:)*h/Fmax
   q(i)%fr(:) = Nset%R(1,:)*q(i)%xx(1) + Nset%R(2,:)*q(i)%xx(2) + &
        &  Nset%R(3,:)*q(i)%xx(3)
   q(i)%fr(:) = MOD(q(i)%fr(:) + 1.0D0, 1.0D0)
   q(i)%xx(:) =  Nset%L(1,:)*q(i)%fr(1) + &
        & Nset%L(2,:)*q(i)%fr(2) + Nset%L(3,:)*q(i)%fr(3)
END DO
!
PRINT *, "# Fmax = ", Fmax
!
RETURN
END SUBROUTINE POSITION_UPDATE
!
!############################################################################
SUBROUTINE VOLUME_CONTROL(Nset, sys, q)
USE DATAFMT
!
TYPE(AMNT):: Nset
TYPE(STAT):: sys
TYPE(PTCL):: q(Nset%Npt_all)
!
REAL(DP) :: P_short, P_thermo, dx(3), df(3)
INTEGER(GI):: i, j, k, Nall
!
Nall = Nset%Npt_all
P_short = 0.0D0
DO i = 1, Nall - 1
   DO j = 1+1, Nall
      dx(:) = q(i)%xx(:) - q(j)%xx(:)
      df(:) = q(i)%ff_short(:) - q(j)%ff_short(:)
      DO k = 1, 3
         P_short = P_short + dx(k)*df(k)
      END DO
   END DO
END DO
P_short = - P_short /3.D0 /Nset%V

P_thermo  = P_short + sys%Ees/3.D0/Nset%V

sys%press = P_thermo
!PRINT *, P_thermo, " = pressure"
RETURN
END SUBROUTINE VOLUME_CONTROL
