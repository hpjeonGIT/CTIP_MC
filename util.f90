!
!#############################################################################
SUBROUTINE PRINT_STAT(Nset, sys, Nloop)
USE DATAFMT
IMPLICIT NONE
TYPE(AMNT):: Nset
TYPE(STAT):: sys
INTEGER(GI):: Nloop
!
REAL(DP):: E_total, E_es, E_eam, Temperature
E_es = sys%Eself + sys%Edir + sys%Erec + sys%Ectip
sys%Ees = E_es
E_eam = sys%Epair + sys%Egroup
E_total = E_es + E_eam
sys%E_old = sys%Eall
sys%Eall = E_total
WRITE(15,100) Nloop, E_total, E_es, E_eam, sys%Epair, sys%Egroup, sys%press
100 FORMAT (I4, 6(1X, ES11.4))
RETURN
END SUBROUTINE PRINT_STAT
!
!#############################################################################
SUBROUTINE PRINT_SNAPSHOT(Nset, Nloop, q)
USE DATAFMT
IMPLICIT NONE
TYPE(AMNT):: Nset
TYPE(PTCL):: q(Nset%Npt_all)
INTEGER(GI):: Nloop
!
INTEGER(GI):: openstatus, nframe, i, id
CHARACTER(LEN=256):: dummy, SNAPFILE
!
nframe = Nloop/Nset%freq_snap
WRITE(SNAPFILE,110) nframe
110 FORMAT("image", I3.3, ".xyz")
!
OPEN(UNIT=21, FILE=SNAPFILE, ACTION="WRITE", &
     & FORM = "FORMATTED", IOSTAT = openstatus)
IF (openstatus >0) THEN
   PRINT '("=== cannot open ", A, "===")', SNAPFILE ; STOP
END IF

WRITE(21,*) Nset%Npt_all
WRITE(21,*) "Frame = ", Nframe, " energy =  0.0 "

DO i=1, Nset%Npt_all
   id = q(i)%id
   SELECT CASE (id)
      CASE (1)
         dummy = "Al"
      CASE (2)
         dummy = "Ni"
      CASE (3)
         dummy = "Fe"
      CASE(4)
         dummy = "O "
      CASE(5)
         dummy = "H "
      CASE DEFAULT
         CLOSE(20)
         PRINT *, id
         STOP "=== Particle id error ==="
   END SELECT
   WRITE(21,300) dummy, q(i)%xx(:), q(i)%q
END DO
200 FORMAT(A3, 3(ES13.6, 1X), 6(ES11.4, 1X), ES13.6)
300 FORMAT(A3, 4(ES11.4, 1X))
CLOSE(21)
RETURN
END SUBROUTINE PRINT_SNAPSHOT
