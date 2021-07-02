! Constitutive law: Non-orthogonal cohesion-friction combined hardening plastic model
! Reference: Zhou X, Lu DC, Su CC, Wang GS, Du XL.  An open source unconstrained 
! A cohesion-friction combined hardening plastic model of concrete based on the non-orthogonal flow rule
! Author: Zhou xin (zhouxin615@126.com)
      
      module global 
       implicit none 
       double precision EK,miu,fc,epd10,cmax,Ac,gampc,phimax,Aphi,
     1 gamphi,u,a1,b1,Ftol,tol,fn(6,1),delta(6,1),Ivol(6,6),
     1 Isym(6,6),PI,Dew(6,6),xc,xphi,xx,Hxc,Hxphi,xgam,qmax
       common EK,miu,fc,epd10,cmax,Ac,gampc,phimax,Aphi,gamphi,
     1 u,a1,b1,Ftol,tol,fn,delta,Ivol,Isym,PI,Dew,xc,xphi,xx,Hxc,Hxphi
      end module

   
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
      use global
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME

      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)
      
! This code is only available for 3D stress state, so NTENS=6
      
      integer NTENS,NSTATV,NPROPS,NDI,NSHR,NOEL,NPT,
     &        LAYER,KSPT,JSTEP,KINC
      double precision SSE,SPD,SCD,RPL,DRPLDE,DTIME,TEMP, 
     &                 DTEMP, PNEWDT,CELENT
      integer K1,K2,Maxit,RN,inittension
      double precision BK,GK,epd,epdw,epd0,c,phi,Rmc,Rmct,pt,qt,thetat,
     1 p,q,theta,pw,qw,thetaw,f,ft,alpha,J2,J2t,J3,J3t,Reps
      double precision II(6,6),IP(6,6),delta1(6,1),sd(6,1),sdt(6,1),
     1 DDp(7,6),Dep(6,6),s(6,1),st(6,1),sw(6,1),sdw(6,1),
     2 x(7,1),dx(7,1),dst(7,1),AIJ(6,6),DDpN(6,7,6),dstra(6,1)

!   The Material parameters of model --------------------------------------------------------------
!   PROPS(1)  - EK      ! Young's modulus
!   PROPS(2)  - miu     ! Poisson's ratio
!   PROPS(3)  - fc      ! Concrete uniaxial compressive strength
!   PROPS(4)  - epd10   ! Equivalent plastic shear strain at the peak of concrete uniaxial compression
!   PROPS(5)  - cmax    ! maximum cohesion in the hardening function
!   PROPS(6)  - Ac      ! shape parameter of hardening curve of  cohesion!   PROPS(7)  - gampc   ! normalized plastic internal variables at the peak of hardening curve of cohesion
!   PROPS(8)  - phimax  ! maximum internal friction angle in the hardening function
!   PROPS(9)  - Aphi    ! shape parameter of hardening curve of internal friction angle!   PROPS(10) - gamphi  ! normalized plastic internal variables at the peak of hardening curve of internal friction angle
!   PROPS(11) - u       ! fractional order

      EK = PROPS(1)
      miu = PROPS(2)
      fc = PROPS(3)
      epd10 = PROPS(4)
      cmax = PROPS(5)
      Ac = PROPS(6)
      gampc = PROPS(7)
      phimax = PROPS(8)
      Aphi = PROPS(9)
      gamphi = PROPS(10)
      u = PROPS(11)

      Ftol = 1.D-6;    ! error tolerance for the NMTR method
      tol = 1.D-14     ! tolerance that is a threshold to determine when to use the L'hospital rule to avoid a zero denominator
      Maxit = 10
      RN = 5
      PI = 3.14159265358979d0
      a1 = 0.001d0
      b1 = 3.03d0
      xx = statev(4)
      xgam = statev(4)
      qmax = statev(3)
      xc = 0.0d0
      xphi = 0.0d0
      Hxphi = 0.0d0
      Hxc = 0.0d0
