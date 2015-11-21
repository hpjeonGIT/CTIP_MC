PROGRAM SIMOX !****************************************************************
!
!
![1] Zhou et al, PRB v.69, pp.035402, (2004)
![2] Hasnaoui et al, Suf. Sci. v.579, pp.47-57 (2005)
!
!*** Note
! The code is using the feature of F95/2003. Some of F90 features may not work.
! Required external library
! - FFT for particle mesh Ewald
USE DATAFMT
IMPLICIT NONE
!
TYPE(AMNT) :: Nset
TYPE(TCKT) :: tag
TYPE(METL) :: eam(Nmetal)
TYPE(POTN) :: param(Natom)
TYPE(OXID) :: oxygn
TYPE(STAT) :: sys
TYPE(CORR) :: pair
REAL(DP)   :: h
INTEGER(GI):: Nloop
TYPE(PTCL), POINTER:: q(:)
REAL(DP),   POINTER:: psi(:,:,:)
REAL(SP)   :: secnds, time0, time1, time2
time0 = 0
time1 = secnds(time0)
!
CALL READ_INPUT(Nset, tag, sys, q, psi, pair, h)
CALL POTENTIAL_SETUP(eam, param, oxygn)
CALL PME_SELF(Nset, sys, param, q, psi)
!
Nloop = 1
h = 0.1D0
sys%E_old = 0.0D0
sys%Eall = 0.0D0
OPEN(UNIT=15, FILE="temporal_status.dat")
WRITE(15,*) "# Loop (ea) - E_total - E_Ewald -  E_eam -  ", &
     & "Epair   -  Egroup (eV)  - Pressure (eV/A^3)"
!
DO WHILE (Nloop  <=  Nset%freq_max)
   IF (MOD(Nloop,10) == 1) CALL CTIP_SOLVER(Nset, param, sys, q, psi)
   CALL FORCE(Nset, eam, oxygn, sys, param, pair, q)
   CALL PME(Nset, sys, q, psi)
   CALL POSITION_UPDATE(Nset, q, h)
   IF (sys%Eall <  sys%E_old) THEN 
      h = 1.2D0*h
      PRINT *, "h is increased", h
   ELSE
      h = h * 0.2D0
      PRINT *, "h is shortened", h
   END IF
   !
   IF (MOD(Nloop, Nset%freq_stat) == 0) CALL PRINT_STAT(Nset, sys, Nloop)
   IF (MOD(Nloop, Nset%freq_snap) == 0) CALL PRINT_SNAPSHOT(Nset, Nloop, q)
   CALL VOLUME_CONTROL(Nset, sys, q)
   Nloop = Nloop + 1
END DO
!
CALL PRINT_SNAPSHOT(Nset, Nloop, q)
!CALL ELASTIC_CONSTANT
CALL EVACUATE_MEMORY(q, psi)
CLOSE(15)
time2 = secnds(time1)
PRINT '("# Wall time =", ES11.4, " seconds with", I8, " loops")', time2, Nloop-1
STOP
!
CONTAINS
!
!###############################################################################
SUBROUTINE READ_INPUT(Nset, tag, sys, q, psi, pair, h)
USE DATAFMT
IMPLICIT NONE
TYPE(AMNT):: Nset
TYPE(TCKT):: tag
TYPE(STAT):: sys
TYPE(CORR):: pair
REAL(DP)  :: h
TYPE(PTCL), POINTER:: q(:)
REAL(DP),   POINTER:: psi(:,:,:)
!
INTEGER(GI):: openstatus, istatus, i, id
CHARACTER(256):: dummy, FILENAME
!
! Initialize tag variables
tag%regular_run = .FALSE.
tag%new_o2      = .FALSE.
tag%new_charge  = .FALSE.
tag%snapshot    = .FALSE.
tag%thermostat  = .FALSE.
tag%restart     = .FALSE.
tag%status      = .FALSE.
tag%print_corr  = .FALSE.
!
! Read simulation parameters **************************************************
OPEN(UNIT=10, FILE='config.prm', STATUS = "OLD", ACTION="READ", &
     & FORM = "FORMATTED", IOSTAT = openstatus)
IF (openstatus >0) STOP "#### cannot  open config.prm ###"
READ(10,*) dummy
READ(10,*) dummy
READ(10,*) dummy
CALL ANY_2_UPPER(dummy)
IF (dummy == 'REGULAR') THEN
   tag%regular_run = .TRUE.
ELSE IF (dummy == 'OXIDATION') THEN
   tag%regular_run = .FALSE.
ELSE
   STOP "=== Simulation type is unknown ==="
END IF
READ(10,*) dummy
READ(10,*) sys%time_max, sys%time_snap, sys%time_rdf, sys%time_stat,h
Nset%freq_max  = NINT(sys%time_max /h)
Nset%freq_snap = NINT(sys%time_snap/h)
Nset%freq_rdf  = NINT(sys%time_rdf /h)
Nset%freq_stat = NINT(sys%time_stat /h)
h = h/ TFM
READ(10,*) dummy
READ(10,*) sys%box(:)
Nset%box(:) = sys%box(:)
Nset%V = Nset%box(1)*Nset%box(2)*Nset%box(3)
READ(10,*) dummy
READ(10,*) Nset%fft(:)
ALLOCATE(psi(Nset%fft(1),Nset%fft(2),Nset%fft(3)), STAT=istatus)
READ(10,*) dummy
READ(10,*) dummy, sys%T_given, Nset%alpha
sys%T_given = sys%T_given*kB ! K -> eV
CALL ANY_2_UPPER(dummy)
IF (dummy == 'YES') THEN
   tag%thermostat = .TRUE.
ELSE IF (dummy == 'NO') THEN
   tag%thermostat = .FALSE.
ELSE
   STOP "=== Temperature control is unknown ==="
END IF
READ(10,*) dummy
READ(10,*) dummy
CALL ANY_2_UPPER(dummy)
IF (dummy == 'YES') THEN
   tag%restart = .TRUE.
ELSE IF (dummy == 'NO') THEN
   tag%restart = .FALSE.
ELSE
   STOP "=== Restart option is unknown ==="
END IF
READ(10,*) dummy
READ(10,*) sys%Rcut
sys%Rcut2 = sys%Rcut**2
sys%a = 3.8D0/sys%Rcut
READ(10,*) dummy
READ(10,*) Nset%freq_new
READ(10,*) dummy
READ(10,*) Nset%z_new
READ(10,*) dummy
READ(10,*) Nset%q_crit
CLOSE(10)
!
! Read particle data **********************************************************
FILENAME = 'input.xyz'
OPEN(UNIT=11, FILE=FILENAME, STATUS = "OLD", ACTION="READ", &
     & FORM = "FORMATTED", IOSTAT = openstatus)
IF (openstatus >0) THEN
   PRINT *, '("=== cannot open ", A, "===")', FILENAME; STOP
END IF
!
! Atomic data configuration
! Al=1 Ni=2 Cu=3 O=4 H=5
Nset%Npt(:) = 0
READ(11,*) Nset%Npt_all
ALLOCATE(q(Nset%Npt_all), STAT=istatus)
IF(istatus /=0) STOP "=== Particle allocation error ==="
READ(11,*) dummy
DO i=1, Nset%Npt_all
   READ(11,*) dummy, q(i)%xx(:), q(i)%q
   SELECT CASE (dummy)
   CASE ("Al")
      id = 1
   CASE ("Ni")
      id = 2
   CASE ("Cu")
      id = 3
   CASE ("O")
      id = 4
   CASE ("H")
      id = 5
   CASE DEFAULT
      STOP "=== Particle type error ==="
   END SELECT
   Nset%Npt(id) = Nset%Npt(id) + 1
   q(i)%id = id
END DO
CLOSE(11)
id = 0; DO i=1,Natom; id = id + Nset%Npt(i); END DO
IF (id /= Nset%Npt_all) THEN
   DEALLOCATE(q, STAT = istatus)
   STOP "=== Number of particles mismatch ==="
END IF
!
! Read lattice vectors *********************************************************
FILENAME = 'latticevector.dat'
OPEN(UNIT=33, FILE=FILENAME, STATUS = "OLD", ACTION="READ", &
     & FORM = "FORMATTED", IOSTAT = openstatus)
IF (openstatus >0) THEN
   PRINT *, '("=== cannot open ", A, "===")', FILENAME; STOP
END IF
READ(33,*) dummy
READ(33,*) Nset%L(1,:)
READ(33,*) Nset%L(2,:)
READ(33,*) Nset%L(3,:)
READ(33,*) dummy
READ(33,*) Nset%R(1,:)
READ(33,*) Nset%R(2,:)
READ(33,*) Nset%R(3,:)
CLOSE(33)
!
! Coordinate translation into fractional coordinate
DO i=1, Nset%Npt_all
   q(i)%fr(:) = Nset%R(1,:)*q(i)%xx(1) + Nset%R(2,:)*q(i)%xx(2) + &
        &       Nset%R(3,:)*q(i)%xx(3)
   q(i)%fr(:) = MOD(q(i)%fr(:) + 1.0D0, 1.0D0)
END DO
!
! Correlation data setup ******************************************************
pair%dr = sys%Rcut / DBLE(Nrdf)
pair%V  = Nset%V
pair%rdf(:,:,:) = 0
pair%Npt_old    = Nset%Npt_all
!
RETURN 
END SUBROUTINE READ_INPUT
!
!##############################################################################
SUBROUTINE EVACUATE_MEMORY(q, psi)
USE DATAFMT
IMPLICIT NONE
!
TYPE(PTCL), POINTER:: q(:)
REAL(DP),   POINTER:: psi(:,:,:)
INTEGER(GI):: istatus
DEALLOCATE(q, STAT = istatus)
IF(istatus /=0) STOP "=== Particle deallocation error ==="
DEALLOCATE(psi, STAT = istatus)
IF(istatus /=0) STOP "=== PME coefficient deallocation error ==="
!
RETURN
END SUBROUTINE EVACUATE_MEMORY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END PROGRAM SIMOX
!
!http://coding.derkeiler.com/Archive/Fortran/comp.lang.fortran/2005-03/0761.html
!###############################################################################
SUBROUTINE ANY_2_UPPER(txt_string)
USE DATAFMT, ONLY: GI
IMPLICIT NONE
CHARACTER(LEN=*):: txt_string
INTEGER(GI):: i, nlen, id
nlen = LEN(txt_string)
DO i=1, nlen
   id = ichar(txt_string(i:i))
   IF (id >= 97 .AND. id < 122) txt_string(i:i) = CHAR(id-32)
END DO
RETURN
END SUBROUTINE ANY_2_UPPER
!
!###############################################################################
SUBROUTINE TAG_SET(Nset, tag, Nloop)
USE DATAFMT, ONLY: GI, DP, AMNT, TCKT, STAT
IMPLICIT NONE
!
TYPE(AMNT) :: Nset
TYPE(TCKT) :: tag
INTEGER(GI):: Nloop
!
IF (Nset%freq_new < 0) THEN
   tag%new_charge = .FALSE.
ELSE
   IF (MOD(Nloop, Nset%freq_new) == 1) THEN
      tag%new_charge = .TRUE.
   ELSE 
      tag%new_charge = .FALSE.
   END IF
END IF
IF (MOD(Nloop, Nset%freq_snap) == 0) THEN
   tag%snapshot = .TRUE.
ELSE
   tag%snapshot = .FALSE.
END IF
IF (MOD(Nloop, Nset%freq_rdf) == 0) THEN
   tag%print_corr = .TRUE.
ELSE
   tag%print_corr = .FALSE.
END IF
IF (MOD(Nloop, Nset%freq_stat) == 0) THEN
   tag%status = .TRUE.
ELSE
   tag%status = .FALSE.
END IF
RETURN
END SUBROUTINE TAG_SET
