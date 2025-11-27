c=======================================================
c=======================================================
      subroutine wfat_get_ascf0(beta, gamma, n, m, zi, erg, 
     &     gn, ef, fn)
c     DESCRIPTION:
c     gn :
c       $\exp(-\kappa\mu_z) g_\nu$ where $g_\nu$ is structure factor.
      implicit none
      double precision beta       ! in degrees
      double precision gamma      ! in degrees
      integer n
      integer m
      double precision zi
      double precision erg
      double complex gn
      double precision ef
      double complex fn

      double precision kp, bt0, pi_const, fifc
      double complex imnum, a
      parameter ( imnum = cmplx(0.0d0, 1.0d0, kind=kind(1.0d0)) )
      parameter (pi_const = 3.14159265359d0)


      kp = sqrt(2 * abs(erg))
      bt0 = zi - kp*(n + (abs(m)+1)/2.0d0)      
      a = imnum*pi_const/4 + imnum*pi_const*bt0/kp
      call wfat_field_factor(ef, erg, zi, n, m, fifc)

      fn = sqrt(fifc) * exp(a) * gn
      
      end
c=======================================================


c=======================================================
c=======================================================
      subroutine wfat_get_ascf1(beta, gamma, n, m, zi, erg, chmmt, pb,
     &     on, gn, jn, kn, ef, fn)
c     DESCRIPTION:
c     gn :
c       $\exp(-\kappa\mu_z) g_\nu$ where $g_\nu$ is structure factor.
c     jn :
c       $\exp(-\kappa\mu_z) j_\nu$ where $j_\nu$ is one of the two WFAT(1) structure factors.
c     kn :
c       $\exp(-\kappa\mu_z) k_\nu$ where $k_\nu$ is one of the two WFAT(1) structure factors.

      implicit none
      double precision beta       ! in degrees
      double precision gamma      ! in degrees
      integer n
      integer m
      double precision zi
      double precision erg
      double precision chmmt(3)
      double precision pb(3,3)
      double precision on(2)
      double complex gn
      double complex jn
      double complex kn
      double precision ef
      double complex fn

      double precision protmat(3,3), mu(3), pbl(3,3), kp, gmm, bt0
      double precision bt1_1, bt1_2, bt1, aa, bb_, on_, fifc, pi_const
      double complex a, an_, imnum
      parameter ( imnum = cmplx(0.0d0, 1.0d0, kind=kind(1.0d0)) )
      parameter (pi_const = 3.14159265359d0)


      call wfat_passrot(beta, gamma, protmat)
      mu = matmul(transpose(protmat), -chmmt)
      call wfat_represent_sym(3, 3, protmat, pb, pbl)   ! pbl = polarizability tensor in lab frame.

      kp = sqrt(2 * abs(erg))
      gmm = (1.0d0 - m**2) / 4.0d0
      bt0 = zi - kp*(n + (abs(m)+1)/2.0d0)
      bt1_1 = -1.0d0/4/kp**2 * (6*n*(n+abs(m)+1) + m**2 + 3*abs(m) + 2)
      bt1_2 = -1.0d0/2/kp * (2*n+1+abs(m))
      bt1 = bt1_1 + mu(3)*bt1_2
      aa = -2*bt1/kp - (gmm-2*mu(3)*bt0)/kp**3 - 3*bt0**2/kp**5
      bb_ = -kp*pbl(3,3) - mu(3)**2/kp + mu(3)/kp**2 +
     &      4*mu(3)*bt0/kp**3 - (10+18*gmm+3*gmm**2)/24/kp**3 -
     &      (9-6*gmm)*bt0/4/kp**4 - (49+2*gmm)*bt0**2/8/kp**5 +
     &      3*bt0**3/2/kp**6 - bt0**4/8/kp**7
      on_ = on(1) + mu(3)*on(2)
      a = imnum*pi_const/4 + imnum*pi_const*bt0/kp
      an_ = jn + kn
      call wfat_field_factor(ef, erg, zi, n, m, fifc)
            
      fn = sqrt(fifc) * exp(a) * (gn*(1 + aa/2*ef*log(ef/4/kp**2)
     &     + ef/2*(bb_+2*on_-imnum*pi_const*aa)) + 2*an_)

      end
c=======================================================
      
