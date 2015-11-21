!
!##############################################################################
SUBROUTINE CTIP_SOLVER(Nset, param, sys, q, psi)
USE DATAFMT
IMPLICIT NONE
!
TYPE(AMNT):: Nset
TYPE(POTN):: param(Natom)
TYPE(STAT):: sys
TYPE(PTCL):: q(Nset%Npt_all)
REAL(DP)  :: psi(Nset%fft(1),Nset%fft(2),Nset%fft(3))
!
REAL(DP):: XJZ(Nset%Npt_all,3), q_i(Nset%Npt_all), g_sum, g_sum_old, d_sum
REAL(DP):: d(Nset%Npt_all), p(Nset%Npt_all), g(Nset%Npt_all), Xall, q_sum
REAL(DP):: afb_i, xi_i, EXP_rxi_i, EXP_rxi_i_r, a, ar, box(3)
REAL(DP):: afb_j, xi_j, EXP_rxi_j, EXP_rxi_j_r, dx(3), rx(3), r2, r
REAL(DP):: xir, fifj, kc_const, alpha, V, Eself, c_pme, Esum, dq
REAL(DP):: asym_mat(Nset%Npt_all, Nset%Npt_all), phi(Nset%Npt_all)
REAL(DP):: sym_mat(Nset%Npt_all,Nset%Npt_all), qmx(Nset%Npt_all, 2)
REAL(DP):: u(3), z(3), AX(4,3), Qlocal(Nset%Npt_all,4,4,4), xfft(3), Error
INTEGER(GI):: Nall, i, j, k, l, id_i, id_j, grid(Nset%Npt_all,4,4,4,3)
INTEGER(GI):: fft(3), nn(3)
REAL(DP), PARAMETER:: E_crit = 1.D-3
!
! Preparation of PME data array
Nall = Nset%Npt_all; Xall = DBLE(Nall)
V = Nset%V; a = sys%a; kc_const = EPS*2.D0*a/SQRTPI
box(:) = Nset%box(:); fft(:) = Nset%fft(:); xfft(:) = DBLE(fft(:))
DO i=1, Nall
   id_i = q(i)%id
   q_i(i) = q(i)%q
   XJZ(i,1) = param(id_i)%chi
   XJZ(i,2) = param(id_i)%J - kc_const
   XJZ(i,3) = param(id_i)%Z
   qmx(i,1) = param(id_i)%qmin; qmx(i,2) = param(id_i)%qmax
   !
   !u(:) = xfft(:)*(q(i)%xx(:) + box(:)*0.5D0)/box(:)
   u(:) = xfft(:)*q(i)%fr(:)

   z(:) = DINT(u(:)) - u(:) + 1.D0
   AX(1,:)  =       z(:)**3                                      /6.D0
   z(:) = z(:) + 1.D0
   AX(2,:) = (-3.D0*z(:)**3 + 12.D0*z(:)**2 - 12.D0*z(:) +  4.D0)/6.D0
   z(:) = z(:) + 1.D0
   AX(3,:) = ( 3.D0*z(:)**3 - 24.D0*z(:)**2 + 60.D0*z(:) - 44.D0)/6.D0
   z(:) = z(:) + 1.D0
   AX(4,:) = (     -z(:)**3 + 12.D0*z(:)**2 - 48.D0*z(:) + 64.D0)/6.D0
   DO j=1,4; nn(1) = j
      DO k=1,4; nn(2) = k
         DO l=1,4; nn(3) = l
            grid(i,j,k,l,:) = MOD(INT(u(:)) + nn(:) - 2 + fft(:), fft(:)) + 1
            Qlocal(i,j,k,l) = AX(j,1)*AX(k,2)*AX(l,3)
         END DO
      END DO
   END DO
END DO
c_pme = EPS * xfft(1)*xfft(2)*xfft(3)/V/PI
!
! Preparation of Coulomb integration data
! asym_mat(i,j) = [j|fi] - [fi|fj]
! sym_mat(i,j) = [fi|fj] - 1/rij + erfc(arij)/rij
asym_mat(:,:) = 0.0D0; sym_mat(:,:) = 0.0D0
DO i=1, Nall-1
   id_i  = q(i)%id
   xi_i  = param(id_i)%xi
   DO j=i+1, Nall
      id_j  = q(j)%id
      xi_j  = param(id_j)%xi


      dx(:) = q(i)%fr(:) - q(j)%fr(:)
      dx(:) = dx(:) - DNINT(dx(:))
      rx(:) = Nset%L(1,:)*dx(1) + Nset%L(2,:)*dx(2) + Nset%L(3,:)*dx(3)
      r2 = rx(1)**2 + rx(2)**2 + rx(3)**2

      r = DSQRT(r2)
      ar = a*r
      EXP_rxi_j = DEXP(-2.D0*xi_j*r); EXP_rxi_j_r = EXP_rxi_j/r
      afb_j = -xi_j*EXP_rxi_j - EXP_rxi_j_r
      IF( id_i == id_j) THEN
         EXP_rxi_i = EXP_rxi_j; EXP_rxi_i_r = EXP_rxi_j_r
         afb_i = afb_j
         xir = xi_i*r
         fifj = - EXP_rxi_j*(1.D0 + xir*(param(id_i)%fafb(id_i,2) + &
              & xir*(param(id_i)%fafb(id_i,3) + &
              & xir*param(id_i)%fafb(id_i,4))))/r
      ELSE
         EXP_rxi_i = DEXP(-2.D0*xi_i*r); EXP_rxi_i_r = EXP_rxi_i/r
         afb_i = -xi_i*EXP_rxi_i - EXP_rxi_i_r
         fifj = - EXP_rxi_i* &
              & (param(id_i)%fafb(id_j,1) + param(id_i)%fafb(id_j,3)/r) &
              & - EXP_rxi_j* &
              & (param(id_i)%fafb(id_j,2) + param(id_i)%fafb(id_j,4)/r)
      END IF
      ! They are short-range and cut-off radius might be implemented
      asym_mat(i,j) = EPS*(afb_i - fifj)
      asym_mat(j,i) = EPS*(afb_j - fifj)
      sym_mat (i,j) = EPS*(fifj + (1.D0 - DERF(ar))/r)
      sym_mat (j,i) = sym_mat(i,j)
   END DO
END DO
!
! Iterative conjugate gradient solver
Error = 1.0; p(:) = 0.0D0; g_sum = 1.0D0; k = 0
DO WHILE (Error > E_crit)
   CALL PHI_MESH(Nall, fft, psi, Qlocal, q_i, grid, phi); phi(:) = phi(:)*c_pme
   CALL dUdq(Nall, XJZ, asym_mat, sym_mat, q_i, qmx, phi, g)
   g_sum_old = g_sum; g_sum = 0.0D0
   DO i=1, Nall
      g_sum = g_sum + g(i)
   END DO
   d(:) = g(:) + p(:)*g_sum/g_sum_old
   d_sum = 0.0D0
   DO i=1, Nall
      d_sum = d_sum + d(i)
   END DO; d_sum = d_sum/Xall
   p(:) = d(:) - d_sum
   !
   CALL BRENT_SOLVER(Nall, fft, psi, Qlocal, grid, c_pme, XJZ, asym_mat, &
        & sym_mat, q_i, qmx, p, alpha)
   q_i(:) = q_i(:) + alpha*p(:)
   !   
   Error= MAXVAL(DABS(p(:)))*DABS(alpha)
   !WRITE(*,100) Error, k
   k = k + 1
END DO
WRITE(*,100) Error, k
100 FORMAT("# CG error =", ES9.2, " at ", I3, "th loop")
!
Eself = 0.0D0
q_sum = 0.0D0
Esum = 0.0D0
DO i=1, Nall
   q(i)%q = q_i(i)
   Esum = Esum + XJZ(i,1)*q_i(i) + 0.5D0*XJZ(i,2)*q_i(i)**2
   dq = q_i(i) - qmx(i,1)
   IF (dq < 0.0D0) Esum = Esum + 2.D0*w*dq**2
   dq = q_i(i) - qmx(i,2)
   IF (dq > 0.0D0) Esum = Esum + 2.D0*w*dq**2
   Eself = Eself - q(i)%q**2
   q_sum = q_sum + q_i(i)
END DO
!WRITE(*,150) q_sum
!150 FORMAT("# Charge sum of CTIP loop = ", ES11.4)
sys%Eself = Eself*EPS*a/SQRTPI
!
! Now add constant component of CTIP into Eself
sys%Eself = sys%Eself + Esum
!
RETURN
END SUBROUTINE CTIP_SOLVER
!
! Brent solver
!##############################################################################
SUBROUTINE BRENT_SOLVER(Nall, fft, psi, Qlocal, grid, c_pme, XJZ, asym_mat, &
     & sym_mat, q_i, qmx, p, alpha)
USE DATAFMT, ONLY: GI, DP
IMPLICIT NONE
INTEGER(GI):: Nall, fft(3)
REAL   (DP):: psi(fft(1),fft(2),fft(3)), Qlocal(Nall,4,4,4), c_pme
INTEGER(GI):: grid(Nall,4,4,4,3)
REAL   (DP):: XJZ(Nall,3), asym_mat(Nall,Nall), sym_mat(Nall,Nall)
REAL   (DP):: q_i(Nall), qmx(Nall,2), p(Nall), alpha
!
INTEGER(GI):: i, n
REAL   (DP):: q_l(Nall), g(Nall), phi(Nall)
REAL   (DP):: fa, fb, fc, x, y, z, a, b, c, d, e, xm, s
INTEGER(GI), PARAMETER:: Nmax = 100
REAL   (DP), PARAMETER:: E_crit = 1.D-8, E_end = 1.D-11
!
a = -1.0D-0
q_l(:) = q_i(:) + a*p(:)
CALL PHI_MESH(Nall, fft, psi, Qlocal, q_l, grid, phi); phi(:) = phi(:)*c_pme
CALL dUdq(Nall, XJZ, asym_mat, sym_mat, q_l, qmx, phi, g)
fa = 0.0D0; DO i=1, Nall; fa = fa - g(i)*p(i); END DO
b =  1.0D-0
q_l(:) = q_i(:) + b*p(:)
CALL PHI_MESH(Nall, fft, psi, Qlocal, q_l, grid, phi); phi(:) = phi(:)*c_pme
CALL dUdq(Nall, XJZ, asym_mat, sym_mat, q_l, qmx, phi, g)
fb = 0.0D0; DO i=1, Nall; fb = fb - g(i)*p(i); END DO
c = b
fc = fb
!
DO n=1, Nmax
   IF ( (fb > 0.D0 .AND. fc > 0.D0) .OR. (fb < 0.0D0 .AND. fc < 0.D0)) THEN
      c = a
      fc = fa
      d = b - a
      e = d
   END IF
   IF (DABS(fc) < DABS(fb)) THEN
      a = b
      b = c
      c = a
      fa = fb
      fb = fc
      fc = fa
   END IF
   xm = 0.5D0*(c - b)
   IF (DABS(xm) < E_crit .OR. DABS(fb) < E_end) THEN
      alpha = b
      !WRITE(*,100,ADVANCE = "NO") n;100 FORMAT("# Loop for Brent solver = ",I3)
      RETURN
   END IF
   IF (DABS(e) >= E_crit .AND. DABS(fa) > DABS(fb)) THEN
      s = fb/fa
      IF (a == C) THEN
         x = 2.D0*xm*s
         y = 1.D0 - s
      ELSE
         y = fa/fc
         z = fb/fc
         x = s*(2.*xm*y*(y-z) - (b-a)*(z-1.0D0))
         y = (y - 1.D0)*(z - 1.D0)*(s - 1.D0)
      END IF
      IF (x > 0.D0) y = -y
      x = DABS(x)
      IF (2.D0*x < MIN(3.*xm*y - DABS(E_crit*y), DABS(e*y))) THEN
         e = d
         d = x /y
      ELSE
         d = xm
         e = d
      END IF
   ELSE
      d = xm
      e = d
   END IF
   a = b
   fa = fb
   IF (DABS(d) > E_crit) THEN
      b = b + d
   ELSE
      b = b + SIGN(E_crit, xm)
   END IF
   q_l(:) = q_i(:) + b*p(:)
   CALL PHI_MESH(Nall, fft, psi, Qlocal, q_l, grid, phi); phi(:) = phi(:)*c_pme
   CALL dUdq(Nall, XJZ, asym_mat, sym_mat, q_l, qmx, phi, g)
   fb = 0.0D0; DO i=1, Nall; fb = fb - g(i)*p(i); END DO
END DO
WRITE(*,*) "# **** Brent solver limit breached *** "
STOP
END SUBROUTINE BRENT_SOLVER
!
!##############################################################################
SUBROUTINE dUdq(Nall, XJZ, asym_mat, sym_mat, q_i, qmx, phi, g)
USE DATAFMT, ONLY:GI, DP, w
IMPLICIT NONE
INTEGER(GI):: Nall
REAL   (DP):: XJZ(Nall,3), asym_mat(Nall,Nall), sym_mat(Nall,Nall)
REAL   (DP):: q_i(Nall), qmx(Nall,2), phi(Nall), g(Nall)
!
INTEGER(GI):: i, j
REAL   (DP):: dqmin(Nall), dqmax(Nall)
!
dqmin(:) = q_i(:) - qmx(:,1); dqmax(:) = q_i(:) - qmx(:,2)
DO i=1, Nall
   IF (dqmin(i) >= 0.0D0) dqmin(i) = 0.0D0
   IF (dqmax(i) <= 0.0D0) dqmax(i) = 0.0D0
END DO
g(:) = -(XJZ(:,1) + XJZ(:,2)*q_i(:) + phi(:) + 2.D0*w*( dqmin(:) + dqmax(:)) )
! XJZ(i, 2) = J_i  - 2a/\sqrt{\pi}
DO i=1, Nall-1
   DO j=i+1, Nall
      g(i) = g(i) - XJZ(j,3)*asym_mat(i,j) - q_i(j)*sym_mat(i,j)
      g(j) = g(j) - XJZ(i,3)*asym_mat(j,i) - q_i(i)*sym_mat(j,i)
   END DO
END DO
!
RETURN
! Returns g(:) = -dU/dq_i(:)
END SUBROUTINE dUdq
!
!##############################################################################
SUBROUTINE PHI_MESH(Nall, fft, psi, Qlocal, q_i, grid, phi)
USE DATAFMT
IMPLICIT NONE
!
INTEGER(GI):: Nall, fft(3)
REAL   (DP):: psi(fft(1),fft(2),fft(3)), Qlocal(Nall,4,4,4), q_i(Nall)
REAL   (DP):: phi(Nall)
INTEGER(GI):: grid(Nall,4,4,4,3)
!
REAL   (DP):: E_local, z(3), xfft(3)
INTEGER(GI):: i, j, k, l, ix(3)
COMPLEX(DP):: Qbas(fft(1), fft(2), fft(3)), Qinv(fft(1), fft(2), fft(3)), &
            & Qfin(fft(1), fft(2), fft(3))
!
!
xfft(:) = DBLE(fft(:))
Qbas(:,:,:) = CMPLX(0.0D0, 0.0D0)
DO i=1,Nall
   DO j=1,4
      DO k=1,4
         DO l=1,4
            ix(:) = grid(i,j,k,l,:)  
            Qbas(ix(1),ix(2),ix(3)) = Qbas(ix(1),ix(2),ix(3)) + &
                 & q_i(i)*CMPLX(Qlocal(i,j,k,l),0.0D0)
         END DO
      END DO
   END DO
END DO
!
Qinv(:,:,:) = Qbas(:,:,:)
CALL ZFFT3D(Qinv, fft(1), fft(2), fft(3), 0)
CALL ZFFT3D(Qinv, fft(1), fft(2), fft(3),+1)
DO i=1, fft(1)
   z(1) = (4.D0 + 2.D0*DCOS(2.D0*PI*DBLE(i-1)/xfft(1)))**2 
   DO j=1, fft(2)
      z(2) = (4.D0 + 2.D0*DCOS(2.D0*PI*DBLE(j-1)/xfft(2)))**2 
      DO k=1, fft(3)
         z(3) = (4.D0 + 2.D0*DCOS(2.D0*PI*DBLE(k-1)/xfft(3)))**2 
         Qinv(i,j,k) = Qinv(i,j,k)*PSI(i,j,k)*36.D0**3/z(1)/z(2)/z(3)
      END DO
   END DO
END DO
Qfin(:,:,:) = Qinv(:,:,:)
CALL ZFFT3D(Qfin, fft(1), fft(2), fft(3), 0)
CALL ZFFT3D(Qfin, fft(1), fft(2), fft(3),-1)
DO i=1, Nall
   E_local = 0.0D0
   DO j=1,4
      DO k=1,4
         DO l=1,4
            ix(:) = grid(i,j,k,l,:)
            E_local = E_local + &
                 & Qlocal(i,j,k,l)*DBLE(Qfin(ix(1),ix(2),ix(3)))
         END DO
      END DO
   END DO
   phi(i) = E_local
END DO
!
RETURN
END SUBROUTINE PHI_MESH
