c     The codes in this file were taken from density functional repository
c     ftp://ftp.dl.ac.uk/qcg/dft_library/index.html
c     $Revision: 1.10 $, $Date: 2006/02/01 13:48:18
c     The routines are listed below, the rks counterpart are also available
c     uks_x_lda
c     uks_x_b88
c     uks_x_pbe
c     uks_x_pw91
c     uks_c_vwn5
c     uks_c_p86
c     uks_c_lyp
c     uks_c_pw91
c     uks_c_pbe
c     uks_xc_b97
c     uks_xc_hcth407


      real*8 function piecewise(o,a,b)
c
c     Fix for a Maple problem. The Maple Fortran generator is not
c     able to transform its piecewise command to Fortran. Therefore
c     this function was written to take care of this.
c
      implicit none
      logical o
      real*8 a, b
      if (o) then
         piecewise = a
      else
         piecewise = b
      endif
      return
      end

      
      subroutine uks_x_lda
     & (ideriv,npt,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,
     &  zk,vrhoa,vrhob,vsigmaaa,vsigmabb,vsigmaab,
     &  v2rhoa2,v2rhob2,v2rhoab,
     &  v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb,
     &  v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,
     &  v2sigmaaa2,v2sigmaaaab,v2sigmaaabb,
     &  v2sigmaab2,v2sigmaabbb,v2sigmabb2)
c
c     P.A.M. Dirac
c     Proceedings of the Cambridge Philosophical Society, 26 (1930) 376
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit real*8 (a-h,o-z)
      integer ideriv,npt
      real*8 rhoa1(npt),rhob1(npt)
      real*8 sigmaaa1(npt),sigmabb1(npt),sigmaab1(npt)
      real*8 zk(npt),vrhoa(npt),vrhob(npt)
      real*8 vsigmaaa(npt),vsigmabb(npt),vsigmaab(npt)
      real*8 v2rhoa2(npt),v2rhob2(npt),v2rhoab(npt)
      real*8 v2rhoasigmaaa(npt),v2rhoasigmaab(npt)
      real*8 v2rhoasigmabb(npt),v2rhobsigmabb(npt)
      real*8 v2rhobsigmaab(npt),v2rhobsigmaaa(npt)
      real*8 v2sigmaaa2(npt),v2sigmaaaab(npt),v2sigmaaabb(npt)
      real*8 v2sigmaab2(npt),v2sigmaabbb(npt),v2sigmabb2(npt)
      parameter(tol=1.0d-20)
      
      if (ideriv.eq.0) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      t1 = rhob**(1.D0/3.D0)
      zk(i) = -0.9305257363491D0*t1*rhob
      elseif(rhob.lt.tol) then
      rho = rhoa
      t1 = rhoa**(1.D0/3.D0)
      zk(i) = -0.9305257363491D0*t1*rhoa
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      t1 = rhoa**(1.D0/3.D0)
      t4 = rhob**(1.D0/3.D0)
      zk(i) = -0.9305257363491D0*t1*rhoa-0.9305257363491D0*t4*rhob
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      t1 = rhob**(1.D0/3.D0)
      zk(i) = -0.9305257363491D0*t1*rhob
      vrhoa(i) = 0.D0
      vrhob(i) = -0.12407009817988D1*t1
      elseif(rhob.lt.tol) then
      rho = rhoa
      t1 = rhoa**(1.D0/3.D0)
      zk(i) = -0.9305257363491D0*t1*rhoa
      vrhoa(i) = -0.12407009817988D1*t1
      vrhob(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      t1 = rhoa**(1.D0/3.D0)
      t4 = rhob**(1.D0/3.D0)
      zk(i) = -0.9305257363491D0*t1*rhoa-0.9305257363491D0*t4*rhob
      vrhoa(i) = -0.12407009817988D1*t1
      vrhob(i) = -0.12407009817988D1*t4
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vrhob(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      t1 = rhob**(1.D0/3.D0)
      zk(i) = -0.9305257363491D0*t1*rhob
      vrhoa(i) = 0.D0
      vrhob(i) = -0.12407009817988D1*t1
      v2rhoa2(i) = 0.D0
      v2rhoab(i) = 0.D0
      t5 = t1**2
      v2rhob2(i) = -0.4135669939329333D0/t5
      elseif(rhob.lt.tol) then
      rho = rhoa
      t1 = rhoa**(1.D0/3.D0)
      zk(i) = -0.9305257363491D0*t1*rhoa
      vrhoa(i) = -0.12407009817988D1*t1
      vrhob(i) = 0.D0
      t5 = t1**2
      v2rhoa2(i) = -0.4135669939329333D0/t5
      v2rhob2(i) = 0.D0
      v2rhoab(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      t1 = rhoa**(1.D0/3.D0)
      t4 = rhob**(1.D0/3.D0)
      zk(i) = -0.9305257363491D0*t1*rhoa-0.9305257363491D0*t4*rhob
      vrhoa(i) = -0.12407009817988D1*t1
      vrhob(i) = -0.12407009817988D1*t4
      t9 = t1**2
      v2rhoa2(i) = -0.4135669939329333D0/t9
      t12 = t4**2
      v2rhob2(i) = -0.4135669939329333D0/t12
      v2rhoab(i) = 0.D0
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vrhob(i) = 0.0d0
      v2rhoa2(i) = 0.0d0
      v2rhob2(i) = 0.0d0
      v2rhoab(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end
      
      
      subroutine rks_x_lda
     & (ideriv,npt,rhoa1,sigmaaa1,
     &  zk,vrhoa,vsigmaaa,
     &  v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
c
c     P.A.M. Dirac
c     Proceedings of the Cambridge Philosophical Society, 26 (1930) 376
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit real*8 (a-h,o-z)
      integer ideriv,npt
      real*8 rhoa1(npt)
      real*8 sigmaaa1(npt)
      real*8 zk(npt),vrhoa(npt),vsigmaaa(npt)
      real*8 v2rhoa2(npt),v2rhoasigmaaa(npt),v2sigmaaa2(npt)
      parameter(tol=1.0d-20)
      
      if(ideriv.eq.0) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      t1 = rho**(1.D0/3.D0)
      zk(i) = -0.7385587663820224D0*t1*rho
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      t1 = rho**(1.D0/3.D0)
      zk(i) = -0.7385587663820224D0*t1*rho
      vrhoa(i) = -0.9847450218426965D0*t1
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      t1 = rho**(1.D0/3.D0)
      zk(i) = -0.7385587663820224D0*t1*rho
      vrhoa(i) = -0.9847450218426965D0*t1
      t5 = t1**2
      v2rhoa2(i) = -0.6564966812284644D0/t5
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      v2rhoa2(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end


      subroutine uks_x_b88
     & (ideriv,npt,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,
     &  zk,vrhoa,vrhob,vsigmaaa,vsigmabb,vsigmaab,
     &  v2rhoa2,v2rhob2,v2rhoab,
     &  v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb,
     &  v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,
     &  v2sigmaaa2,v2sigmaaaab,v2sigmaaabb,
     &  v2sigmaab2,v2sigmaabbb,v2sigmabb2)
c
c     A.D. Becke
c     Density-functional exchange-energy approximation with correct 
c     asymptotic behaviour
c     Phys. Rev. A38 (1988) 3098-3100
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit real*8 (a-h,o-z)
      integer ideriv,npt
      real*8 rhoa1(npt),rhob1(npt)
      real*8 sigmaaa1(npt),sigmabb1(npt),sigmaab1(npt)
      real*8 zk(npt),vrhoa(npt),vrhob(npt)
      real*8 vsigmaaa(npt),vsigmabb(npt),vsigmaab(npt)
      real*8 v2rhoa2(npt),v2rhob2(npt),v2rhoab(npt)
      real*8 v2rhoasigmaaa(npt),v2rhoasigmaab(npt)
      real*8 v2rhoasigmabb(npt),v2rhobsigmabb(npt)
      real*8 v2rhobsigmaab(npt),v2rhobsigmaaa(npt)
      real*8 v2sigmaaa2(npt),v2sigmaaaab(npt),v2sigmaaabb(npt)
      real*8 v2sigmaab2(npt),v2sigmaabbb(npt),v2sigmabb2(npt)
      parameter(tol=1.0d-20)
      
      if (ideriv.eq.0) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmabb
      t2 = rhob**(1.D0/3.D0)
      t3 = t2*rhob
      t5 = 1/t3
      t7 = dsqrt(sigmabb)
      t8 = t7*t5
      t9 = dlog(t8+dsqrt(1+t8**2))
      zk(i) = -0.9305257363491D0*t3-0.42D-2*t5*sigmabb/(1.D0+0.252D-1
     #*t8*t9)
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = rhoa**(1.D0/3.D0)
      t3 = t2*rhoa
      t5 = 1/t3
      t7 = dsqrt(sigmaaa)
      t8 = t7*t5
      t9 = dlog(t8+dsqrt(1+t8**2))
      zk(i) = -0.9305257363491D0*t3-0.42D-2*t5*sigmaaa/(1.D0+0.252D-1
     #*t8*t9)
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rhoa**(1.D0/3.D0)
      t5 = t4*rhoa
      t7 = 1/t5
      t9 = dsqrt(sigmaaa)
      t10 = t9*t7
      t11 = dlog(t10+dsqrt(1+t10**2))
      t18 = rhob**(1.D0/3.D0)
      t19 = t18*rhob
      t21 = 1/t19
      t23 = dsqrt(sigmabb)
      t24 = t23*t21
      t25 = dlog(t24+dsqrt(1+t24**2))
      zk(i) = -0.9305257363491D0*t5-0.42D-2*t7*sigmaaa/(1.D0+0.252D-1
     #*t10*t11)-0.9305257363491D0*t19-0.42D-2*t21*sigmabb/(1.D0
     #+0.252D-1*t24*t25)
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmabb
      t2 = rhob**(1.D0/3.D0)
      t3 = t2*rhob
      t5 = 1/t3
      t6 = t5*sigmabb
      t7 = dsqrt(sigmabb)
      t8 = t7*t5
      t9 = dlog(t8+dsqrt(1+t8**2))
      t12 = 1.D0+0.252D-1*t8*t9
      t13 = 1/t12
      zk(i) = -0.9305257363491D0*t3-0.42D-2*t6*t13
      vrhoa(i) = 0.D0
      t17 = rhob**2
      t19 = 1/t2/t17
      t23 = t12**2
      t24 = 1/t23
      t29 = t2**2
      t34 = 1/t29/t17
      t37 = dsqrt(1.D0+sigmabb*t34)
      t38 = 1/t37
      vrhob(i) = -0.12407009817988D1*t2+0.56D-2*t19*sigmabb*t13
     #+0.42D-2*t6*t24*(-0.336D-1*t7*t19*t9-0.336D-1*sigmabb/t29/t17
     #/rhob*t38)
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      vsigmabb(i) = -0.42D-2*t5*t13+0.42D-2*t6*t24*(0.126D-1/t7*t5*t9
     #+0.126D-1*t34*t38)
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = rhoa**(1.D0/3.D0)
      t3 = t2*rhoa
      t5 = 1/t3
      t6 = t5*sigmaaa
      t7 = dsqrt(sigmaaa)
      t8 = t7*t5
      t9 = dlog(t8+dsqrt(1+t8**2))
      t12 = 1.D0+0.252D-1*t8*t9
      t13 = 1/t12
      zk(i) = -0.9305257363491D0*t3-0.42D-2*t6*t13
      t17 = rhoa**2
      t19 = 1/t2/t17
      t23 = t12**2
      t24 = 1/t23
      t29 = t2**2
      t34 = 1/t29/t17
      t37 = dsqrt(1.D0+sigmaaa*t34)
      t38 = 1/t37
      vrhoa(i) = -0.12407009817988D1*t2+0.56D-2*t19*sigmaaa*t13
     #+0.42D-2*t6*t24*(-0.336D-1*t7*t19*t9-0.336D-1*sigmaaa/t29/t17
     #/rhoa*t38)
      vrhob(i) = 0.D0
      vsigmaaa(i) = -0.42D-2*t5*t13+0.42D-2*t6*t24*(0.126D-1/t7*t5*t9
     #+0.126D-1*t34*t38)
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rhoa**(1.D0/3.D0)
      t5 = t4*rhoa
      t7 = 1/t5
      t8 = t7*sigmaaa
      t9 = dsqrt(sigmaaa)
      t10 = t9*t7
      t11 = dlog(t10+dsqrt(1+t10**2))
      t14 = 1.D0+0.252D-1*t10*t11
      t15 = 1/t14
      t18 = rhob**(1.D0/3.D0)
      t19 = t18*rhob
      t21 = 1/t19
      t22 = t21*sigmabb
      t23 = dsqrt(sigmabb)
      t24 = t23*t21
      t25 = dlog(t24+dsqrt(1+t24**2))
      t28 = 1.D0+0.252D-1*t24*t25
      t29 = 1/t28
      zk(i) = -0.9305257363491D0*t5-0.42D-2*t8*t15-0.9305257363491D0
     #*t19-0.42D-2*t22*t29
      t33 = rhoa**2
      t35 = 1/t4/t33
      t39 = t14**2
      t40 = 1/t39
      t45 = t4**2
      t50 = 1/t45/t33
      t53 = dsqrt(1.D0+sigmaaa*t50)
      t54 = 1/t53
      vrhoa(i) = -0.12407009817988D1*t4+0.56D-2*t35*sigmaaa*t15
     #+0.42D-2*t8*t40*(-0.336D-1*t9*t35*t11-0.336D-1*sigmaaa/t45/t33
     #/rhoa*t54)
      t62 = rhob**2
      t64 = 1/t18/t62
      t68 = t28**2
      t69 = 1/t68
      t74 = t18**2
      t79 = 1/t74/t62
      t82 = dsqrt(1.D0+sigmabb*t79)
      t83 = 1/t82
      vrhob(i) = -0.12407009817988D1*t18+0.56D-2*t64*sigmabb*t29
     #+0.42D-2*t22*t69*(-0.336D-1*t23*t64*t25-0.336D-1*sigmabb/t74/t62
     #/rhob*t83)
      vsigmaaa(i) = -0.42D-2*t7*t15+0.42D-2*t8*t40*(0.126D-1/t9*t7*t11
     #+0.126D-1*t50*t54)
      vsigmaab(i) = 0.D0
      vsigmabb(i) = -0.42D-2*t21*t29+0.42D-2*t22*t69*(0.126D-1/t23*t21
     #*t25+0.126D-1*t79*t83)
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vrhob(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      vsigmaab(i) = 0.0d0
      vsigmabb(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmabb
      t2 = rhob**(1.D0/3.D0)
      t3 = t2*rhob
      t5 = 1/t3
      t6 = t5*sigmabb
      t7 = dsqrt(sigmabb)
      t8 = t7*t5
      t9 = dlog(t8+dsqrt(1+t8**2))
      t12 = 1.D0+0.252D-1*t8*t9
      t13 = 1/t12
      zk(i) = -0.9305257363491D0*t3-0.42D-2*t6*t13
      vrhoa(i) = 0.D0
      t17 = rhob**2
      t19 = 1/t2/t17
      t20 = t19*sigmabb
      t23 = t12**2
      t24 = 1/t23
      t28 = t17*rhob
      t29 = t2**2
      t34 = 1/t29/t17
      t36 = 1.D0+sigmabb*t34
      t37 = dsqrt(t36)
      t38 = 1/t37
      t41 = -0.336D-1*t7*t19*t9-0.336D-1*sigmabb/t29/t28*t38
      t42 = t24*t41
      vrhob(i) = -0.12407009817988D1*t2+0.56D-2*t20*t13+0.42D-2*t6*t42
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      t53 = 0.126D-1/t7*t5*t9+0.126D-1*t34*t38
      vsigmabb(i) = -0.42D-2*t5*t13+0.42D-2*t6*t24*t53
      v2rhoa2(i) = 0.D0
      v2rhoab(i) = 0.D0
      t60 = 1/t2/t28
      t67 = 1/t23/t12
      t68 = t41**2
      t75 = t17**2
      t81 = sigmabb**2
      t87 = 1/t37/t36
      v2rhob2(i) = -0.4135669939329333D0/t29-0.1306666666666667D-1*t60
     #*sigmabb*t13-0.112D-1*t20*t42-0.84D-2*t6*t67*t68+0.42D-2*t6*t24*
     #(0.784D-1*t7*t60*t9+0.168D0*sigmabb/t29/t75*t38-0.448D-1*t81/t2
     #/t75/t28*t87)
      v2sigmaaa2(i) = 0.D0
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      t97 = t53**2
      v2sigmabb2(i) = 0.84D-2*t5*t24*t53-0.84D-2*t6*t67*t97+0.42D-2*t6
     #*t24*(-0.63D-2/t7/sigmabb*t5*t9+0.63D-2/sigmabb*t34*t38-0.63D-2
     #/t2/t75/rhob*t87)
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = rhoa**(1.D0/3.D0)
      t3 = t2*rhoa
      t5 = 1/t3
      t6 = t5*sigmaaa
      t7 = dsqrt(sigmaaa)
      t8 = t7*t5
      t9 = dlog(t8+dsqrt(1+t8**2))
      t12 = 1.D0+0.252D-1*t8*t9
      t13 = 1/t12
      zk(i) = -0.9305257363491D0*t3-0.42D-2*t6*t13
      t17 = rhoa**2
      t19 = 1/t2/t17
      t20 = t19*sigmaaa
      t23 = t12**2
      t24 = 1/t23
      t28 = t17*rhoa
      t29 = t2**2
      t34 = 1/t29/t17
      t36 = 1.D0+sigmaaa*t34
      t37 = dsqrt(t36)
      t38 = 1/t37
      t41 = -0.336D-1*t7*t19*t9-0.336D-1*sigmaaa/t29/t28*t38
      t42 = t24*t41
      vrhoa(i) = -0.12407009817988D1*t2+0.56D-2*t20*t13+0.42D-2*t6*t42
      vrhob(i) = 0.D0
      t53 = 0.126D-1/t7*t5*t9+0.126D-1*t34*t38
      vsigmaaa(i) = -0.42D-2*t5*t13+0.42D-2*t6*t24*t53
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      t60 = 1/t2/t28
      t67 = 1/t23/t12
      t68 = t41**2
      t75 = t17**2
      t81 = sigmaaa**2
      t87 = 1/t37/t36
      v2rhoa2(i) = -0.4135669939329333D0/t29-0.1306666666666667D-1*t60
     #*sigmaaa*t13-0.112D-1*t20*t42-0.84D-2*t6*t67*t68+0.42D-2*t6*t24*
     #(0.784D-1*t7*t60*t9+0.168D0*sigmaaa/t29/t75*t38-0.448D-1*t81/t2
     #/t75/t28*t87)
      v2rhob2(i) = 0.D0
      v2rhoab(i) = 0.D0
      t97 = t53**2
      v2sigmaaa2(i) = 0.84D-2*t5*t24*t53-0.84D-2*t6*t67*t97+0.42D-2*t6
     #*t24*(-0.63D-2/t7/sigmaaa*t5*t9+0.63D-2/sigmaaa*t34*t38-0.63D-2
     #/t2/t75/rhoa*t87)
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      v2sigmabb2(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rhoa**(1.D0/3.D0)
      t5 = t4*rhoa
      t7 = 1/t5
      t8 = t7*sigmaaa
      t9 = dsqrt(sigmaaa)
      t10 = t9*t7
      t11 = dlog(t10+dsqrt(1+t10**2))
      t14 = 1.D0+0.252D-1*t10*t11
      t15 = 1/t14
      t18 = rhob**(1.D0/3.D0)
      t19 = t18*rhob
      t21 = 1/t19
      t22 = t21*sigmabb
      t23 = dsqrt(sigmabb)
      t24 = t23*t21
      t25 = dlog(t24+dsqrt(1+t24**2))
      t28 = 1.D0+0.252D-1*t24*t25
      t29 = 1/t28
      zk(i) = -0.9305257363491D0*t5-0.42D-2*t8*t15-0.9305257363491D0
     #*t19-0.42D-2*t22*t29
      t33 = rhoa**2
      t35 = 1/t4/t33
      t36 = t35*sigmaaa
      t39 = t14**2
      t40 = 1/t39
      t44 = t33*rhoa
      t45 = t4**2
      t47 = 1/t45/t44
      t50 = 1/t45/t33
      t52 = 1.D0+sigmaaa*t50
      t53 = dsqrt(t52)
      t54 = 1/t53
      t57 = -0.336D-1*t9*t35*t11-0.336D-1*sigmaaa*t47*t54
      t58 = t40*t57
      vrhoa(i) = -0.12407009817988D1*t4+0.56D-2*t36*t15+0.42D-2*t8*t58
      t62 = rhob**2
      t64 = 1/t18/t62
      t65 = t64*sigmabb
      t68 = t28**2
      t69 = 1/t68
      t73 = t62*rhob
      t74 = t18**2
      t76 = 1/t74/t73
      t79 = 1/t74/t62
      t81 = 1.D0+sigmabb*t79
      t82 = dsqrt(t81)
      t83 = 1/t82
      t86 = -0.336D-1*t23*t64*t25-0.336D-1*sigmabb*t76*t83
      t87 = t69*t86
      vrhob(i) = -0.12407009817988D1*t18+0.56D-2*t65*t29+0.42D-2*t22*t87
      t92 = 1/t9
      t98 = 0.126D-1*t92*t7*t11+0.126D-1*t50*t54
      t99 = t40*t98
      vsigmaaa(i) = -0.42D-2*t7*t15+0.42D-2*t8*t99
      vsigmaab(i) = 0.D0
      t104 = 1/t23
      t110 = 0.126D-1*t104*t21*t25+0.126D-1*t79*t83
      t111 = t69*t110
      vsigmabb(i) = -0.42D-2*t21*t29+0.42D-2*t22*t111
      t117 = 1/t4/t44
      t124 = 1/t39/t14
      t125 = t57**2
      t132 = t33**2
      t138 = sigmaaa**2
      t144 = 1/t53/t52
      v2rhoa2(i) = -0.4135669939329333D0/t45-0.1306666666666667D-1
     #*t117*sigmaaa*t15-0.112D-1*t36*t58-0.84D-2*t8*t124*t125+0.42D-2
     #*t8*t40*(0.784D-1*t9*t117*t11+0.168D0*sigmaaa/t45/t132*t54
     #-0.448D-1*t138/t4/t132/t44*t144)
      t154 = 1/t18/t73
      t161 = 1/t68/t28
      t162 = t86**2
      t169 = t62**2
      t175 = sigmabb**2
      t181 = 1/t82/t81
      v2rhob2(i) = -0.4135669939329333D0/t74-0.1306666666666667D-1
     #*t154*sigmabb*t29-0.112D-1*t65*t87-0.84D-2*t22*t161*t162+0.42D-2
     #*t22*t69*(0.784D-1*t23*t154*t25+0.168D0*sigmabb/t74/t169*t83
     #-0.448D-1*t175/t18/t169/t73*t181)
      v2rhoab(i) = 0.D0
      t192 = t7*t40
      v2rhoasigmaaa(i) = 0.56D-2*t35*t15-0.56D-2*t36*t99+0.42D-2*t192
     #*t57-0.84D-2*t8*t124*t57*t98+0.42D-2*t8*t40*(-0.168D-1*t92*t35
     #*t11-0.504D-1*t47*t54+0.168D-1*sigmaaa/t4/t132/t33*t144)
      v2rhoasigmaab(i) = 0.D0
      v2rhoasigmabb(i) = 0.D0
      v2rhobsigmaaa(i) = 0.D0
      v2rhobsigmaab(i) = 0.D0
      t218 = t21*t69
      v2rhobsigmabb(i) = 0.56D-2*t64*t29-0.56D-2*t65*t111+0.42D-2*t218
     #*t86-0.84D-2*t22*t161*t86*t110+0.42D-2*t22*t69*(-0.168D-1*t104
     #*t64*t25-0.504D-1*t76*t83+0.168D-1*sigmabb/t18/t169/t62*t181)
      t242 = t98**2
      v2sigmaaa2(i) = 0.84D-2*t192*t98-0.84D-2*t8*t124*t242+0.42D-2*t8
     #*t40*(-0.63D-2/t9/sigmaaa*t7*t11+0.63D-2/sigmaaa*t50*t54-0.63D-2
     #/t4/t132/rhoa*t144)
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      t266 = t110**2
      v2sigmabb2(i) = 0.84D-2*t218*t110-0.84D-2*t22*t161*t266+0.42D-2
     #*t22*t69*(-0.63D-2/t23/sigmabb*t21*t25+0.63D-2/sigmabb*t79*t83
     #-0.63D-2/t18/t169/rhob*t181)
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vrhob(i) = 0.0d0
      v2rhoa2(i) = 0.0d0
      v2rhob2(i) = 0.0d0
      v2rhoab(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      vsigmaab(i) = 0.0d0
      vsigmabb(i) = 0.0d0
      v2rhoasigmaaa(i) = 0.0d0
      v2rhoasigmaab(i) = 0.0d0
      v2rhoasigmabb(i) = 0.0d0
      v2rhobsigmaaa(i) = 0.0d0
      v2rhobsigmaab(i) = 0.0d0
      v2rhobsigmabb(i) = 0.0d0
      v2sigmaaa2(i) = 0.0d0
      v2sigmaab2(i) = 0.0d0
      v2sigmabb2(i) = 0.0d0
      v2sigmaaaab(i) = 0.0d0
      v2sigmaaabb(i) = 0.0d0
      v2sigmaabbb(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end
      
      
      subroutine rks_x_b88
     & (ideriv,npt,rhoa1,sigmaaa1,
     &  zk,vrhoa,vsigmaaa,
     &  v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
c
c     A.D. Becke
c     Density-functional exchange-energy approximation with correct 
c     asymptotic behaviour
c     Phys. Rev. A38 (1988) 3098-3100
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit real*8 (a-h,o-z)
      integer ideriv,npt
      real*8 rhoa1(npt)
      real*8 sigmaaa1(npt)
      real*8 zk(npt),vrhoa(npt),vsigmaaa(npt)
      real*8 v2rhoa2(npt),v2rhoasigmaaa(npt),v2sigmaaa2(npt)
      parameter(tol=1.0d-20)
      
      if(ideriv.eq.0) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = rho**(1.D0/3.D0)
      t3 = t2*rho
      t5 = 1/t3
      t7 = dsqrt(sigma)
      t8 = t7*t5
      t10 = dlog(0.1259921049894873D1*t8+dsqrt(1+0.1587401051968199D1
     #*t8**2))
      zk(i) = -0.7385587663820224D0*t3-0.5291668409558467D-2*t5*sigma/
     #(1.D0+0.317500104573508D-1*t8*t10)
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = rho**(1.D0/3.D0)
      t3 = t2*rho
      t5 = 1/t3
      t6 = t5*sigma
      t7 = dsqrt(sigma)
      t8 = t7*t5
      t10 = dlog(0.1259921049894873D1*t8+dsqrt(1+0.1587401051968199D1
     #*t8**2))
      t13 = 1.D0+0.317500104573508D-1*t8*t10
      t14 = 1/t13
      zk(i) = -0.7385587663820224D0*t3-0.5291668409558467D-2*t6*t14
      t18 = rho**2
      t20 = 1/t2/t18
      t24 = t13**2
      t25 = 1/t24
      t30 = t2**2
      t35 = 1/t30/t18
      t39 = dsqrt(1.D0+0.1587401051968199D1*sigma*t35)
      t40 = 1/t39
      vrhoa(i) = -0.9847450218426965D0*t2+0.705555787941129D-2*t20
     #*sigma*t14+0.2645834204779234D-2*t6*t25*(-0.8466669455293548D-1
     #*t7*t20*t10-0.106673350692263D0*sigma/t30/t18/rho*t40)
      vsigmaaa(i) = -0.2116667363823387D-1*t5*t14
     #+0.5291668409558467D-2*t6*t25*(0.6350002091470161D-1/t7*t5*t10
     #+0.8000501301919725D-1*t35*t40)
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = rho**(1.D0/3.D0)
      t3 = t2*rho
      t5 = 1/t3
      t6 = t5*sigma
      t7 = dsqrt(sigma)
      t8 = t7*t5
      t10 = dlog(0.1259921049894873D1*t8+dsqrt(1+0.1587401051968199D1
     #*t8**2))
      t13 = 1.D0+0.317500104573508D-1*t8*t10
      t14 = 1/t13
      zk(i) = -0.7385587663820224D0*t3-0.5291668409558467D-2*t6*t14
      t18 = rho**2
      t20 = 1/t2/t18
      t21 = t20*sigma
      t24 = t13**2
      t25 = 1/t24
      t29 = t18*rho
      t30 = t2**2
      t32 = 1/t30/t29
      t35 = 1/t30/t18
      t38 = 1.D0+0.1587401051968199D1*sigma*t35
      t39 = dsqrt(t38)
      t40 = 1/t39
      t43 = -0.8466669455293548D-1*t7*t20*t10-0.106673350692263D0
     #*sigma*t32*t40
      t44 = t25*t43
      vrhoa(i) = -0.9847450218426965D0*t2+0.705555787941129D-2*t21*t14
     #+0.2645834204779234D-2*t6*t44
      t49 = 1/t7
      t55 = 0.6350002091470161D-1*t49*t5*t10+0.8000501301919725D-1*t35
     #*t40
      t56 = t25*t55
      vsigmaaa(i) = -0.2116667363823387D-1*t5*t14
     #+0.5291668409558467D-2*t6*t56
      t62 = 1/t2/t29
      t69 = 1/t24/t13
      t70 = t43**2
      t77 = t18**2
      t83 = sigma**2
      t89 = 1/t39/t38
      v2rhoa2(i) = -0.6564966812284644D0/t30-0.3292593677058602D-1*t62
     #*sigma*t14-0.1411111575882258D-1*t21*t44-0.5291668409558467D-2
     #*t6*t69*t70+0.2645834204779234D-2*t6*t25*(0.3951112412470322D0
     #*t7*t62*t10+0.106673350692263D1*sigma/t30/t77*t40
     #-0.4515557042823225D0*t83/t2/t77/t29*t89)
      t100 = t5*t25
      v2rhoasigmaaa(i) = 0.2822223151764516D-1*t20*t14
     #-0.705555787941129D-2*t21*t56+0.1058333681911693D-1*t100*t43
     #-0.5291668409558467D-2*t6*t69*t43*t55+0.2645834204779234D-2*t6
     #*t25*(-0.169333389105871D0*t49*t20*t10-0.640040104153578D0*t32
     #*t40+0.3386667782117419D0*sigma/t2/t77/t18*t89)
      t124 = t55**2
      v2sigmaaa2(i) = 0.4233334727646774D-1*t100*t55
     #-0.1058333681911693D-1*t6*t69*t124+0.5291668409558467D-2*t6*t25*
     #(-0.1270000418294032D0/t7/sigma*t5*t10+0.1600100260383945D0
     #/sigma*t35*t40-0.2540000836588064D0/t2/t77/rho*t89)
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      v2rhoa2(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      v2rhoasigmaaa(i) = 0.0d0
      v2sigmaaa2(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end


      subroutine uks_x_pbe
     & (ideriv,npt,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,
     &  zk,vrhoa,vrhob,vsigmaaa,vsigmabb,vsigmaab,
     &  v2rhoa2,v2rhob2,v2rhoab,
     &  v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb,
     &  v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,
     &  v2sigmaaa2,v2sigmaaaab,v2sigmaaabb,
     &  v2sigmaab2,v2sigmaabbb,v2sigmabb2)
c
c     J.P. Perdew, K. Burke, and M. Ernzerhof
c     Generalized gradient approximation made simple
c     Phys. Rev. Lett. 77 (1996) 3865-3868
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit real*8 (a-h,o-z)
      integer ideriv,npt
      real*8 rhoa1(npt),rhob1(npt)
      real*8 sigmaaa1(npt),sigmabb1(npt),sigmaab1(npt)
      real*8 zk(npt),vrhoa(npt),vrhob(npt)
      real*8 vsigmaaa(npt),vsigmabb(npt),vsigmaab(npt)
      real*8 v2rhoa2(npt),v2rhob2(npt),v2rhoab(npt)
      real*8 v2rhoasigmaaa(npt),v2rhoasigmaab(npt)
      real*8 v2rhoasigmabb(npt),v2rhobsigmabb(npt)
      real*8 v2rhobsigmaab(npt),v2rhobsigmaaa(npt)
      real*8 v2sigmaaa2(npt),v2sigmaaaab(npt),v2sigmaaabb(npt)
      real*8 v2sigmaab2(npt),v2sigmaabbb(npt),v2sigmabb2(npt)
      parameter(tol=1.0d-20)
      
      if (ideriv.eq.0) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmabb
      t2 = rhob**(1.D0/3.D0)
      t4 = rhob**2
      t5 = t2**2
      zk(i) = -0.9305257363491D0*t2*rhob*(0.1804D1-0.804D0/(1.D0
     #+0.449276922095889D-2*sigmabb/t5/t4))
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = rhoa**(1.D0/3.D0)
      t4 = rhoa**2
      t5 = t2**2
      zk(i) = -0.9305257363491D0*t2*rhoa*(0.1804D1-0.804D0/(1.D0
     #+0.449276922095889D-2*sigmaaa/t5/t4))
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rhoa**(1.D0/3.D0)
      t6 = rhoa**2
      t7 = t4**2
      t18 = rhob**(1.D0/3.D0)
      t20 = rhob**2
      t21 = t18**2
      zk(i) = -0.9305257363491D0*t4*rhoa*(0.1804D1-0.804D0/(1.D0
     #+0.449276922095889D-2*sigmaaa/t7/t6))-0.9305257363491D0*t18*rhob
     #*(0.1804D1-0.804D0/(1.D0+0.449276922095889D-2*sigmabb/t21/t20))
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmabb
      t2 = rhob**(1.D0/3.D0)
      t3 = t2*rhob
      t4 = rhob**2
      t5 = t2**2
      t10 = 1.D0+0.449276922095889D-2*sigmabb/t5/t4
      t13 = 0.1804D1-0.804D0/t10
      zk(i) = -0.9305257363491D0*t3*t13
      vrhoa(i) = 0.D0
      t20 = t10**2
      t21 = 1/t20
      vrhob(i) = -0.12407009817988D1*t2*t13+0.8963286558970112D-2/t2
     #/t4*t21*sigmabb
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      vsigmabb(i) = -0.3361232459613792D-2/t3*t21
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = rhoa**(1.D0/3.D0)
      t3 = t2*rhoa
      t4 = rhoa**2
      t5 = t2**2
      t10 = 1.D0+0.449276922095889D-2*sigmaaa/t5/t4
      t13 = 0.1804D1-0.804D0/t10
      zk(i) = -0.9305257363491D0*t3*t13
      t20 = t10**2
      t21 = 1/t20
      vrhoa(i) = -0.12407009817988D1*t2*t13+0.8963286558970112D-2/t2
     #/t4*t21*sigmaaa
      vrhob(i) = 0.D0
      vsigmaaa(i) = -0.3361232459613792D-2/t3*t21
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rhoa**(1.D0/3.D0)
      t5 = t4*rhoa
      t6 = rhoa**2
      t7 = t4**2
      t12 = 1.D0+0.449276922095889D-2*sigmaaa/t7/t6
      t15 = 0.1804D1-0.804D0/t12
      t18 = rhob**(1.D0/3.D0)
      t19 = t18*rhob
      t20 = rhob**2
      t21 = t18**2
      t26 = 1.D0+0.449276922095889D-2*sigmabb/t21/t20
      t29 = 0.1804D1-0.804D0/t26
      zk(i) = -0.9305257363491D0*t5*t15-0.9305257363491D0*t19*t29
      t36 = t12**2
      t37 = 1/t36
      vrhoa(i) = -0.12407009817988D1*t4*t15+0.8963286558970112D-2/t4
     #/t6*t37*sigmaaa
      t45 = t26**2
      t46 = 1/t45
      vrhob(i) = -0.12407009817988D1*t18*t29+0.8963286558970112D-2/t18
     #/t20*t46*sigmabb
      vsigmaaa(i) = -0.3361232459613792D-2/t5*t37
      vsigmaab(i) = 0.D0
      vsigmabb(i) = -0.3361232459613792D-2/t19*t46
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vrhob(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      vsigmaab(i) = 0.0d0
      vsigmabb(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmabb
      t2 = rhob**(1.D0/3.D0)
      t3 = t2*rhob
      t4 = rhob**2
      t5 = t2**2
      t10 = 1.D0+0.449276922095889D-2*sigmabb/t5/t4
      t13 = 0.1804D1-0.804D0/t10
      zk(i) = -0.9305257363491D0*t3*t13
      vrhoa(i) = 0.D0
      t20 = t10**2
      t21 = 1/t20
      vrhob(i) = -0.12407009817988D1*t2*t13+0.8963286558970112D-2/t2
     #/t4*t21*sigmabb
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      vsigmabb(i) = -0.3361232459613792D-2/t3*t21
      v2rhoa2(i) = 0.D0
      v2rhoab(i) = 0.D0
      t37 = t4**2
      t41 = 1/t20/t10
      t43 = sigmabb**2
      v2rhob2(i) = -0.4135669939329333D0/t5*t13-0.8963286558970112D-2
     #/t2/t4/rhob*t21*sigmabb+0.2147732158441357D-3/t37/t4*t41*t43
      v2sigmaaa2(i) = 0.D0
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      v2sigmabb2(i) = 0.3020248347808158D-4/t37*t41
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = rhoa**(1.D0/3.D0)
      t3 = t2*rhoa
      t4 = rhoa**2
      t5 = t2**2
      t10 = 1.D0+0.449276922095889D-2*sigmaaa/t5/t4
      t13 = 0.1804D1-0.804D0/t10
      zk(i) = -0.9305257363491D0*t3*t13
      t20 = t10**2
      t21 = 1/t20
      vrhoa(i) = -0.12407009817988D1*t2*t13+0.8963286558970112D-2/t2
     #/t4*t21*sigmaaa
      vrhob(i) = 0.D0
      vsigmaaa(i) = -0.3361232459613792D-2/t3*t21
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      t37 = t4**2
      t41 = 1/t20/t10
      t43 = sigmaaa**2
      v2rhoa2(i) = -0.4135669939329333D0/t5*t13-0.8963286558970112D-2
     #/t2/t4/rhoa*t21*sigmaaa+0.2147732158441357D-3/t37/t4*t41*t43
      v2rhob2(i) = 0.D0
      v2rhoab(i) = 0.D0
      v2sigmaaa2(i) = 0.3020248347808158D-4/t37*t41
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      v2sigmabb2(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rhoa**(1.D0/3.D0)
      t5 = t4*rhoa
      t6 = rhoa**2
      t7 = t4**2
      t12 = 1.D0+0.449276922095889D-2*sigmaaa/t7/t6
      t15 = 0.1804D1-0.804D0/t12
      t18 = rhob**(1.D0/3.D0)
      t19 = t18*rhob
      t20 = rhob**2
      t21 = t18**2
      t26 = 1.D0+0.449276922095889D-2*sigmabb/t21/t20
      t29 = 0.1804D1-0.804D0/t26
      zk(i) = -0.9305257363491D0*t5*t15-0.9305257363491D0*t19*t29
      t36 = t12**2
      t37 = 1/t36
      t38 = 1/t4/t6*t37
      vrhoa(i) = -0.12407009817988D1*t4*t15+0.8963286558970112D-2*t38
     #*sigmaaa
      t45 = t26**2
      t46 = 1/t45
      t47 = 1/t18/t20*t46
      vrhob(i) = -0.12407009817988D1*t18*t29+0.8963286558970112D-2*t47
     #*sigmabb
      vsigmaaa(i) = -0.3361232459613792D-2/t5*t37
      vsigmaab(i) = 0.D0
      vsigmabb(i) = -0.3361232459613792D-2/t19*t46
      t65 = t6**2
      t69 = 1/t36/t12
      t71 = sigmaaa**2
      v2rhoa2(i) = -0.4135669939329333D0/t7*t15-0.8963286558970112D-2
     #/t4/t6/rhoa*t37*sigmaaa+0.2147732158441357D-3/t65/t6*t69*t71
      t83 = t20**2
      t87 = 1/t45/t26
      t89 = sigmabb**2
      v2rhob2(i) = -0.4135669939329333D0/t21*t29-0.8963286558970112D-2
     #/t18/t20/rhob*t46*sigmabb+0.2147732158441357D-3/t83/t20*t87*t89
      v2rhoab(i) = 0.D0
      v2rhoasigmaaa(i) = 0.4481643279485056D-2*t38
     #-0.8053995594155087D-4/t65/rhoa*t69*sigmaaa
      v2rhoasigmaab(i) = 0.D0
      v2rhoasigmabb(i) = 0.D0
      v2rhobsigmaaa(i) = 0.D0
      v2rhobsigmaab(i) = 0.D0
      v2rhobsigmabb(i) = 0.4481643279485056D-2*t47
     #-0.8053995594155087D-4/t83/rhob*t87*sigmabb
      v2sigmaaa2(i) = 0.3020248347808158D-4/t65*t69
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      v2sigmabb2(i) = 0.3020248347808158D-4/t83*t87
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vrhob(i) = 0.0d0
      v2rhoa2(i) = 0.0d0
      v2rhob2(i) = 0.0d0
      v2rhoab(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      vsigmaab(i) = 0.0d0
      vsigmabb(i) = 0.0d0
      v2rhoasigmaaa(i) = 0.0d0
      v2rhoasigmaab(i) = 0.0d0
      v2rhoasigmabb(i) = 0.0d0
      v2rhobsigmaaa(i) = 0.0d0
      v2rhobsigmaab(i) = 0.0d0
      v2rhobsigmabb(i) = 0.0d0
      v2sigmaaa2(i) = 0.0d0
      v2sigmaab2(i) = 0.0d0
      v2sigmabb2(i) = 0.0d0
      v2sigmaaaab(i) = 0.0d0
      v2sigmaaabb(i) = 0.0d0
      v2sigmaabbb(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end
      
      
      subroutine rks_x_pbe
     & (ideriv,npt,rhoa1,sigmaaa1,
     &  zk,vrhoa,vsigmaaa,
     &  v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
c
c     J.P. Perdew, K. Burke, and M. Ernzerhof
c     Generalized gradient approximation made simple
c     Phys. Rev. Lett. 77 (1996) 3865-3868
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit real*8 (a-h,o-z)
      integer ideriv,npt
      real*8 rhoa1(npt)
      real*8 sigmaaa1(npt)
      real*8 zk(npt),vrhoa(npt),vsigmaaa(npt)
      real*8 v2rhoa2(npt),v2rhoasigmaaa(npt),v2sigmaaa2(npt)
      parameter(tol=1.0d-20)
      
      if(ideriv.eq.0) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = rho**(1.D0/3.D0)
      t4 = rho**2
      t5 = t2**2
      zk(i) = -0.7385587663820224D0*t2*rho*(0.1804D1-0.804D0/(1.D0
     #+0.7131826587600489D-2*sigma/t5/t4))
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = rho**(1.D0/3.D0)
      t3 = t2*rho
      t4 = rho**2
      t5 = t2**2
      t10 = 1.D0+0.7131826587600489D-2*sigma/t5/t4
      t13 = 0.1804D1-0.804D0/t10
      zk(i) = -0.7385587663820224D0*t3*t13
      t20 = t10**2
      t21 = 1/t20
      vrhoa(i) = -0.9847450218426965D0*t2*t13+0.1129303341188623D-1/t2
     #/t4*t21*sigma
      vsigmaaa(i) = -0.1693955011782934D-1/t3*t21
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = rho**(1.D0/3.D0)
      t3 = t2*rho
      t4 = rho**2
      t5 = t2**2
      t10 = 1.D0+0.7131826587600489D-2*sigma/t5/t4
      t13 = 0.1804D1-0.804D0/t10
      zk(i) = -0.7385587663820224D0*t3*t13
      t20 = t10**2
      t21 = 1/t20
      t22 = 1/t2/t4*t21
      vrhoa(i) = -0.9847450218426965D0*t2*t13+0.1129303341188623D-1
     #*t22*sigma
      vsigmaaa(i) = -0.1693955011782934D-1/t3*t21
      t37 = t4**2
      t41 = 1/t20/t10
      t43 = sigma**2
      v2rhoa2(i) = -0.6564966812284644D0/t5*t13-0.2258606682377246D-1
     #/t2/t4/rho*t21*sigma+0.8590928633765426D-3/t37/t4*t41*t43
      v2rhoasigmaaa(i) = 0.2258606682377246D-1*t22
     #-0.644319647532407D-3/t37/rho*t41*sigma
      v2sigmaaa2(i) = 0.9664794712986105D-3/t37*t41
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      v2rhoa2(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      v2rhoasigmaaa(i) = 0.0d0
      v2sigmaaa2(i) = 0.0d0
      endif ! rho
      enddo
      
      
      endif ! ideriv
      return
      end
            subroutine uks_x_pw91
     & (ideriv,npt,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,
     &  zk,vrhoa,vrhob,vsigmaaa,vsigmabb,vsigmaab,
     &  v2rhoa2,v2rhob2,v2rhoab,
     &  v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb,
     &  v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,
     &  v2sigmaaa2,v2sigmaaaab,v2sigmaaabb,
     &  v2sigmaab2,v2sigmaabbb,v2sigmabb2)
c
c     J.P. Perdew, J.A. Chevary, S.H. Vosko, K.A. Jackson,
c     M.R. Pederson, D.J. Singh, C. Fiolhais
c     Atoms, molecules, solids and surfaces:
c     Applications of the generalized gradient approximation
c     for exchange and correlation
c     Phys. Rev. B 46 (1992) 6671--6687
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit real*8 (a-h,o-z)
      integer ideriv,npt
      real*8 rhoa1(npt),rhob1(npt)
      real*8 sigmaaa1(npt),sigmabb1(npt),sigmaab1(npt)
      real*8 zk(npt),vrhoa(npt),vrhob(npt)
      real*8 vsigmaaa(npt),vsigmabb(npt),vsigmaab(npt)
      real*8 v2rhoa2(npt),v2rhob2(npt),v2rhoab(npt)
      real*8 v2rhoasigmaaa(npt),v2rhoasigmaab(npt)
      real*8 v2rhoasigmabb(npt),v2rhobsigmabb(npt)
      real*8 v2rhobsigmaab(npt),v2rhobsigmaaa(npt)
      real*8 v2sigmaaa2(npt),v2sigmaaaab(npt),v2sigmaaabb(npt)
      real*8 v2sigmaab2(npt),v2sigmaabbb(npt),v2sigmabb2(npt)
      parameter(tol=1.0d-20)
      
      if (ideriv.eq.0) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmabb
      t2 = rhob**(1.D0/3.D0)
      t3 = t2*rhob
      t4 = dsqrt(sigmabb)
      t6 = t4/t3
      t8 = dlog(0.1000005877780776D1*t6+dsqrt(1+0.10000117555961D1*t6*
     #*2))
      t10 = 0.2520026100493014D-1*t6*t8
      t11 = rhob**2
      t12 = t2**2
      t14 = 1/t12/t11
      t17 = dexp(-0.1645530784602056D1*sigmabb*t14)
      t25 = sigmabb**2
      t26 = t11**2
      zk(i) = -0.9305257363491D0*t3*(1.D0+t10+0.1645530784602056D-1*
     #(0.2743D0-0.1508D0*t17)*sigmabb*t14)/(1.D0+t10
     #+0.1083108625229223D-5*t25/t2/t26/rhob)
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = rhoa**(1.D0/3.D0)
      t3 = t2*rhoa
      t4 = dsqrt(sigmaaa)
      t6 = t4/t3
      t8 = dlog(0.1000005877780776D1*t6+dsqrt(1+0.10000117555961D1*t6*
     #*2))
      t10 = 0.2520026100493014D-1*t6*t8
      t11 = rhoa**2
      t12 = t2**2
      t14 = 1/t12/t11
      t17 = dexp(-0.1645530784602056D1*sigmaaa*t14)
      t25 = sigmaaa**2
      t26 = t11**2
      zk(i) = -0.9305257363491D0*t3*(1.D0+t10+0.1645530784602056D-1*
     #(0.2743D0-0.1508D0*t17)*sigmaaa*t14)/(1.D0+t10
     #+0.1083108625229223D-5*t25/t2/t26/rhoa)
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rhoa**(1.D0/3.D0)
      t5 = t4*rhoa
      t6 = dsqrt(sigmaaa)
      t8 = t6/t5
      t10 = dlog(0.1000005877780776D1*t8+dsqrt(1+0.10000117555961D1*t8
     #**2))
      t12 = 0.2520026100493014D-1*t8*t10
      t13 = rhoa**2
      t14 = t4**2
      t16 = 1/t14/t13
      t19 = dexp(-0.1645530784602056D1*sigmaaa*t16)
      t27 = sigmaaa**2
      t28 = t13**2
      t38 = rhob**(1.D0/3.D0)
      t39 = t38*rhob
      t40 = dsqrt(sigmabb)
      t42 = t40/t39
      t44 = dlog(0.1000005877780776D1*t42+dsqrt(1+0.10000117555961D1
     #*t42**2))
      t46 = 0.2520026100493014D-1*t42*t44
      t47 = rhob**2
      t48 = t38**2
      t50 = 1/t48/t47
      t53 = dexp(-0.1645530784602056D1*sigmabb*t50)
      t61 = sigmabb**2
      t62 = t47**2
      zk(i) = -0.9305257363491D0*t5*(1.D0+t12+0.1645530784602056D-1*
     #(0.2743D0-0.1508D0*t19)*sigmaaa*t16)/(1.D0+t12
     #+0.1083108625229223D-5*t27/t4/t28/rhoa)-0.9305257363491D0*t39*
     #(1.D0+t46+0.1645530784602056D-1*(0.2743D0-0.1508D0*t53)*sigmabb
     #*t50)/(1.D0+t46+0.1083108625229223D-5*t61/t38/t62/rhob)
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmabb
      t2 = rhob**(1.D0/3.D0)
      t3 = t2*rhob
      t4 = dsqrt(sigmabb)
      t5 = 1/t3
      t6 = t4*t5
      t8 = dlog(0.1000005877780776D1*t6+dsqrt(1+0.10000117555961D1*t6*
     #*2))
      t10 = 0.2520026100493014D-1*t6*t8
      t11 = rhob**2
      t12 = t2**2
      t14 = 1/t12/t11
      t15 = sigmabb*t14
      t17 = dexp(-0.1645530784602056D1*t15)
      t19 = 0.2743D0-0.1508D0*t17
      t20 = t19*sigmabb
      t23 = 1.D0+t10+0.1645530784602056D-1*t20*t14
      t24 = t3*t23
      t25 = sigmabb**2
      t26 = t11**2
      t29 = 1/t2/t26/rhob
      t32 = 1.D0+t10+0.1083108625229223D-5*t25*t29
      t33 = 1/t32
      zk(i) = -0.9305257363491D0*t24*t33
      vrhoa(i) = 0.D0
      t43 = 0.3360034800657352D-1*t4/t2/t11*t8
      t46 = 1/t12/t11/rhob
      t50 = dsqrt(1.D0+0.10000117555961D1*t15)
      t51 = 1/t50
      t53 = 0.3360054550205309D-1*sigmabb*t46*t51
      t57 = t25/t2/t26/t11
      t66 = t32**2
      t67 = 1/t66
      vrhob(i) = -0.12407009817988D1*t2*t23*t33-0.9305257363491D0*t3*(
     #-t43-t53-0.1088885204563779D-1*t57*t17-0.4388082092272149D-1*t20
     #*t46)*t33+0.9305257363491D0*t24*t67*(-t43-t53
     #-0.5776579334555855D-5*t57)
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      t76 = 0.1260013050246507D-1/t4*t5*t8
      t78 = 0.1260020456326991D-1*t14*t51
      vsigmabb(i) = -0.9305257363491D0*t3*(t76+t78
     #+0.408331951711417D-2*t29*t17*sigmabb+0.1645530784602056D-1*t19
     #*t14)*t33+0.9305257363491D0*t24*t67*(t76+t78
     #+0.2166217250458446D-5*sigmabb*t29)
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = rhoa**(1.D0/3.D0)
      t3 = t2*rhoa
      t4 = dsqrt(sigmaaa)
      t5 = 1/t3
      t6 = t4*t5
      t8 = dlog(0.1000005877780776D1*t6+dsqrt(1+0.10000117555961D1*t6*
     #*2))
      t10 = 0.2520026100493014D-1*t6*t8
      t11 = rhoa**2
      t12 = t2**2
      t14 = 1/t12/t11
      t15 = sigmaaa*t14
      t17 = dexp(-0.1645530784602056D1*t15)
      t19 = 0.2743D0-0.1508D0*t17
      t20 = t19*sigmaaa
      t23 = 1.D0+t10+0.1645530784602056D-1*t20*t14
      t24 = t3*t23
      t25 = sigmaaa**2
      t26 = t11**2
      t29 = 1/t2/t26/rhoa
      t32 = 1.D0+t10+0.1083108625229223D-5*t25*t29
      t33 = 1/t32
      zk(i) = -0.9305257363491D0*t24*t33
      t43 = 0.3360034800657352D-1*t4/t2/t11*t8
      t46 = 1/t12/t11/rhoa
      t50 = dsqrt(1.D0+0.10000117555961D1*t15)
      t51 = 1/t50
      t53 = 0.3360054550205309D-1*sigmaaa*t46*t51
      t57 = t25/t2/t26/t11
      t66 = t32**2
      t67 = 1/t66
      vrhoa(i) = -0.12407009817988D1*t2*t23*t33-0.9305257363491D0*t3*(
     #-t43-t53-0.1088885204563779D-1*t57*t17-0.4388082092272149D-1*t20
     #*t46)*t33+0.9305257363491D0*t24*t67*(-t43-t53
     #-0.5776579334555855D-5*t57)
      vrhob(i) = 0.D0
      t76 = 0.1260013050246507D-1/t4*t5*t8
      t78 = 0.1260020456326991D-1*t14*t51
      vsigmaaa(i) = -0.9305257363491D0*t3*(t76+t78
     #+0.408331951711417D-2*t29*t17*sigmaaa+0.1645530784602056D-1*t19
     #*t14)*t33+0.9305257363491D0*t24*t67*(t76+t78
     #+0.2166217250458446D-5*sigmaaa*t29)
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rhoa**(1.D0/3.D0)
      t5 = t4*rhoa
      t6 = dsqrt(sigmaaa)
      t7 = 1/t5
      t8 = t6*t7
      t10 = dlog(0.1000005877780776D1*t8+dsqrt(1+0.10000117555961D1*t8
     #**2))
      t12 = 0.2520026100493014D-1*t8*t10
      t13 = rhoa**2
      t14 = t4**2
      t16 = 1/t14/t13
      t17 = sigmaaa*t16
      t19 = dexp(-0.1645530784602056D1*t17)
      t21 = 0.2743D0-0.1508D0*t19
      t22 = t21*sigmaaa
      t25 = 1.D0+t12+0.1645530784602056D-1*t22*t16
      t26 = t5*t25
      t27 = sigmaaa**2
      t28 = t13**2
      t31 = 1/t4/t28/rhoa
      t34 = 1.D0+t12+0.1083108625229223D-5*t27*t31
      t35 = 1/t34
      t38 = rhob**(1.D0/3.D0)
      t39 = t38*rhob
      t40 = dsqrt(sigmabb)
      t41 = 1/t39
      t42 = t40*t41
      t44 = dlog(0.1000005877780776D1*t42+dsqrt(1+0.10000117555961D1
     #*t42**2))
      t46 = 0.2520026100493014D-1*t42*t44
      t47 = rhob**2
      t48 = t38**2
      t50 = 1/t48/t47
      t51 = sigmabb*t50
      t53 = dexp(-0.1645530784602056D1*t51)
      t55 = 0.2743D0-0.1508D0*t53
      t56 = t55*sigmabb
      t59 = 1.D0+t46+0.1645530784602056D-1*t56*t50
      t60 = t39*t59
      t61 = sigmabb**2
      t62 = t47**2
      t65 = 1/t38/t62/rhob
      t68 = 1.D0+t46+0.1083108625229223D-5*t61*t65
      t69 = 1/t68
      zk(i) = -0.9305257363491D0*t26*t35-0.9305257363491D0*t60*t69
      t79 = 0.3360034800657352D-1*t6/t4/t13*t10
      t82 = 1/t14/t13/rhoa
      t86 = dsqrt(1.D0+0.10000117555961D1*t17)
      t87 = 1/t86
      t89 = 0.3360054550205309D-1*sigmaaa*t82*t87
      t93 = t27/t4/t28/t13
      t102 = t34**2
      t103 = 1/t102
      vrhoa(i) = -0.12407009817988D1*t4*t25*t35-0.9305257363491D0*t5*(
     #-t79-t89-0.1088885204563779D-1*t93*t19-0.4388082092272149D-1*t22
     #*t82)*t35+0.9305257363491D0*t26*t103*(-t79-t89
     #-0.5776579334555855D-5*t93)
      t116 = 0.3360034800657352D-1*t40/t38/t47*t44
      t119 = 1/t48/t47/rhob
      t123 = dsqrt(1.D0+0.10000117555961D1*t51)
      t124 = 1/t123
      t126 = 0.3360054550205309D-1*sigmabb*t119*t124
      t130 = t61/t38/t62/t47
      t139 = t68**2
      t140 = 1/t139
      vrhob(i) = -0.12407009817988D1*t38*t59*t69-0.9305257363491D0*t39
     #*(-t116-t126-0.1088885204563779D-1*t130*t53
     #-0.4388082092272149D-1*t56*t119)*t69+0.9305257363491D0*t60*t140*
     #(-t116-t126-0.5776579334555855D-5*t130)
      t149 = 0.1260013050246507D-1/t6*t7*t10
      t151 = 0.1260020456326991D-1*t16*t87
      vsigmaaa(i) = -0.9305257363491D0*t5*(t149+t151
     #+0.408331951711417D-2*t31*t19*sigmaaa+0.1645530784602056D-1*t21
     #*t16)*t35+0.9305257363491D0*t26*t103*(t149+t151
     #+0.2166217250458446D-5*sigmaaa*t31)
      vsigmaab(i) = 0.D0
      t170 = 0.1260013050246507D-1/t40*t41*t44
      t172 = 0.1260020456326991D-1*t50*t124
      vsigmabb(i) = -0.9305257363491D0*t39*(t170+t172
     #+0.408331951711417D-2*t65*t53*sigmabb+0.1645530784602056D-1*t55
     #*t50)*t69+0.9305257363491D0*t60*t140*(t170+t172
     #+0.2166217250458446D-5*sigmabb*t65)
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vrhob(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      vsigmaab(i) = 0.0d0
      vsigmabb(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmabb
      t2 = rhob**(1.D0/3.D0)
      t3 = t2*rhob
      t4 = dsqrt(sigmabb)
      t5 = 1/t3
      t6 = t4*t5
      t8 = dlog(0.1000005877780776D1*t6+dsqrt(1+0.10000117555961D1*t6*
     #*2))
      t10 = 0.2520026100493014D-1*t6*t8
      t11 = rhob**2
      t12 = t2**2
      t14 = 1/t12/t11
      t15 = sigmabb*t14
      t17 = dexp(-0.1645530784602056D1*t15)
      t19 = 0.2743D0-0.1508D0*t17
      t20 = t19*sigmabb
      t23 = 1.D0+t10+0.1645530784602056D-1*t20*t14
      t24 = t3*t23
      t25 = sigmabb**2
      t26 = t11**2
      t29 = 1/t2/t26/rhob
      t32 = 1.D0+t10+0.1083108625229223D-5*t25*t29
      t33 = 1/t32
      zk(i) = -0.9305257363491D0*t24*t33
      vrhoa(i) = 0.D0
      t36 = t2*t23
      t43 = 0.3360034800657352D-1*t4/t2/t11*t8
      t44 = t11*rhob
      t46 = 1/t12/t44
      t49 = 1.D0+0.10000117555961D1*t15
      t50 = dsqrt(t49)
      t51 = 1/t50
      t53 = 0.3360054550205309D-1*sigmabb*t46*t51
      t57 = t25/t2/t26/t11
      t62 = -t43-t53-0.1088885204563779D-1*t57*t17
     #-0.4388082092272149D-1*t20*t46
      t63 = t3*t62
      t66 = t32**2
      t67 = 1/t66
      t69 = -t43-t53-0.5776579334555855D-5*t57
      t70 = t67*t69
      vrhob(i) = -0.12407009817988D1*t36*t33-0.9305257363491D0*t63*t33
     #+0.9305257363491D0*t24*t70
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      t76 = 0.1260013050246507D-1/t4*t5*t8
      t78 = 0.1260020456326991D-1*t14*t51
      t79 = t29*t17
      t85 = t3*(t76+t78+0.408331951711417D-2*t79*sigmabb
     #+0.1645530784602056D-1*t19*t14)
      t90 = t76+t78+0.2166217250458446D-5*sigmabb*t29
      t91 = t67*t90
      vsigmabb(i) = -0.9305257363491D0*t85*t33+0.9305257363491D0*t24*t91
      v2rhoa2(i) = 0.D0
      v2rhoab(i) = 0.D0
      t107 = 0.7840081201533821D-1*t4/t2/t44*t8
      t109 = 1/t12/t26
      t112 = 0.1680027275102654D0*sigmabb*t109*t51
      t116 = t25/t2/t26/t44
      t118 = 1/t50/t49
      t120 = 0.4480125399532632D-1*t116*t118
      t124 = t26**2
      t139 = 1/t66/t32
      t140 = t69**2
      v2rhob2(i) = -0.4135669939329333D0/t12*t23*t33
     #-0.24814019635976D1*t2*t62*t33+0.24814019635976D1*t36*t70
     #-0.9305257363491D0*t3*(t107+t112-t120+0.9799966841074009D-1*t116
     #*t17-0.4778117666686413D-1*t25*sigmabb/t124/t11*t17
     #+0.1608963433833121D0*t20*t109)*t33+0.18610514726982D1*t63*t70
     #-0.18610514726982D1*t24*t139*t140+0.9305257363491D0*t24*t67*
     #(t107+t112-t120+0.3658500245218708D-4*t116)
      v2sigmaaa2(i) = 0.D0
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      t153 = 0.6300065251232535D-2/t4/sigmabb*t5*t8
      t157 = 0.6300102281634954D-2/sigmabb*t14*t51
      t159 = 0.6300176343092764D-2*t29*t118
      t171 = t90**2
      v2sigmabb2(i) = -0.9305257363491D0*t3*(-t153+t157-t159
     #-0.6719227968777768D-2/t124*t17*sigmabb+0.8166639034228341D-2
     #*t79)*t33+0.18610514726982D1*t85*t91-0.18610514726982D1*t24*t139
     #*t171+0.9305257363491D0*t24*t67*(-t153+t157-t159
     #+0.2166217250458446D-5*t29)
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = rhoa**(1.D0/3.D0)
      t3 = t2*rhoa
      t4 = dsqrt(sigmaaa)
      t5 = 1/t3
      t6 = t4*t5
      t8 = dlog(0.1000005877780776D1*t6+dsqrt(1+0.10000117555961D1*t6*
     #*2))
      t10 = 0.2520026100493014D-1*t6*t8
      t11 = rhoa**2
      t12 = t2**2
      t14 = 1/t12/t11
      t15 = sigmaaa*t14
      t17 = dexp(-0.1645530784602056D1*t15)
      t19 = 0.2743D0-0.1508D0*t17
      t20 = t19*sigmaaa
      t23 = 1.D0+t10+0.1645530784602056D-1*t20*t14
      t24 = t3*t23
      t25 = sigmaaa**2
      t26 = t11**2
      t29 = 1/t2/t26/rhoa
      t32 = 1.D0+t10+0.1083108625229223D-5*t25*t29
      t33 = 1/t32
      zk(i) = -0.9305257363491D0*t24*t33
      t36 = t2*t23
      t43 = 0.3360034800657352D-1*t4/t2/t11*t8
      t44 = t11*rhoa
      t46 = 1/t12/t44
      t49 = 1.D0+0.10000117555961D1*t15
      t50 = dsqrt(t49)
      t51 = 1/t50
      t53 = 0.3360054550205309D-1*sigmaaa*t46*t51
      t57 = t25/t2/t26/t11
      t62 = -t43-t53-0.1088885204563779D-1*t57*t17
     #-0.4388082092272149D-1*t20*t46
      t63 = t3*t62
      t66 = t32**2
      t67 = 1/t66
      t69 = -t43-t53-0.5776579334555855D-5*t57
      t70 = t67*t69
      vrhoa(i) = -0.12407009817988D1*t36*t33-0.9305257363491D0*t63*t33
     #+0.9305257363491D0*t24*t70
      vrhob(i) = 0.D0
      t76 = 0.1260013050246507D-1/t4*t5*t8
      t78 = 0.1260020456326991D-1*t14*t51
      t79 = t29*t17
      t85 = t3*(t76+t78+0.408331951711417D-2*t79*sigmaaa
     #+0.1645530784602056D-1*t19*t14)
      t90 = t76+t78+0.2166217250458446D-5*sigmaaa*t29
      t91 = t67*t90
      vsigmaaa(i) = -0.9305257363491D0*t85*t33+0.9305257363491D0*t24*t91
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      t107 = 0.7840081201533821D-1*t4/t2/t44*t8
      t109 = 1/t12/t26
      t112 = 0.1680027275102654D0*sigmaaa*t109*t51
      t116 = t25/t2/t26/t44
      t118 = 1/t50/t49
      t120 = 0.4480125399532632D-1*t116*t118
      t124 = t26**2
      t139 = 1/t66/t32
      t140 = t69**2
      v2rhoa2(i) = -0.4135669939329333D0/t12*t23*t33
     #-0.24814019635976D1*t2*t62*t33+0.24814019635976D1*t36*t70
     #-0.9305257363491D0*t3*(t107+t112-t120+0.9799966841074009D-1*t116
     #*t17-0.4778117666686413D-1*t25*sigmaaa/t124/t11*t17
     #+0.1608963433833121D0*t20*t109)*t33+0.18610514726982D1*t63*t70
     #-0.18610514726982D1*t24*t139*t140+0.9305257363491D0*t24*t67*
     #(t107+t112-t120+0.3658500245218708D-4*t116)
      v2rhob2(i) = 0.D0
      v2rhoab(i) = 0.D0
      t153 = 0.6300065251232535D-2/t4/sigmaaa*t5*t8
      t157 = 0.6300102281634954D-2/sigmaaa*t14*t51
      t159 = 0.6300176343092764D-2*t29*t118
      t171 = t90**2
      v2sigmaaa2(i) = -0.9305257363491D0*t3*(-t153+t157-t159
     #-0.6719227968777768D-2/t124*t17*sigmaaa+0.8166639034228341D-2
     #*t79)*t33+0.18610514726982D1*t85*t91-0.18610514726982D1*t24*t139
     #*t171+0.9305257363491D0*t24*t67*(-t153+t157-t159
     #+0.2166217250458446D-5*t29)
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      v2sigmabb2(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rhoa**(1.D0/3.D0)
      t5 = t4*rhoa
      t6 = dsqrt(sigmaaa)
      t7 = 1/t5
      t8 = t6*t7
      t10 = dlog(0.1000005877780776D1*t8+dsqrt(1+0.10000117555961D1*t8
     #**2))
      t12 = 0.2520026100493014D-1*t8*t10
      t13 = rhoa**2
      t14 = t4**2
      t16 = 1/t14/t13
      t17 = sigmaaa*t16
      t19 = dexp(-0.1645530784602056D1*t17)
      t21 = 0.2743D0-0.1508D0*t19
      t22 = t21*sigmaaa
      t25 = 1.D0+t12+0.1645530784602056D-1*t22*t16
      t26 = t5*t25
      t27 = sigmaaa**2
      t28 = t13**2
      t31 = 1/t4/t28/rhoa
      t34 = 1.D0+t12+0.1083108625229223D-5*t27*t31
      t35 = 1/t34
      t38 = rhob**(1.D0/3.D0)
      t39 = t38*rhob
      t40 = dsqrt(sigmabb)
      t41 = 1/t39
      t42 = t40*t41
      t44 = dlog(0.1000005877780776D1*t42+dsqrt(1+0.10000117555961D1
     #*t42**2))
      t46 = 0.2520026100493014D-1*t42*t44
      t47 = rhob**2
      t48 = t38**2
      t50 = 1/t48/t47
      t51 = sigmabb*t50
      t53 = dexp(-0.1645530784602056D1*t51)
      t55 = 0.2743D0-0.1508D0*t53
      t56 = t55*sigmabb
      t59 = 1.D0+t46+0.1645530784602056D-1*t56*t50
      t60 = t39*t59
      t61 = sigmabb**2
      t62 = t47**2
      t65 = 1/t38/t62/rhob
      t68 = 1.D0+t46+0.1083108625229223D-5*t61*t65
      t69 = 1/t68
      zk(i) = -0.9305257363491D0*t26*t35-0.9305257363491D0*t60*t69
      t72 = t4*t25
      t76 = 1/t4/t13
      t79 = 0.3360034800657352D-1*t6*t76*t10
      t80 = t13*rhoa
      t82 = 1/t14/t80
      t85 = 1.D0+0.10000117555961D1*t17
      t86 = dsqrt(t85)
      t87 = 1/t86
      t89 = 0.3360054550205309D-1*sigmaaa*t82*t87
      t92 = 1/t4/t28/t13
      t93 = t27*t92
      t98 = -t79-t89-0.1088885204563779D-1*t93*t19
     #-0.4388082092272149D-1*t22*t82
      t99 = t5*t98
      t102 = t34**2
      t103 = 1/t102
      t105 = -t79-t89-0.5776579334555855D-5*t93
      t106 = t103*t105
      vrhoa(i) = -0.12407009817988D1*t72*t35-0.9305257363491D0*t99*t35
     #+0.9305257363491D0*t26*t106
      t109 = t38*t59
      t113 = 1/t38/t47
      t116 = 0.3360034800657352D-1*t40*t113*t44
      t117 = t47*rhob
      t119 = 1/t48/t117
      t122 = 1.D0+0.10000117555961D1*t51
      t123 = dsqrt(t122)
      t124 = 1/t123
      t126 = 0.3360054550205309D-1*sigmabb*t119*t124
      t129 = 1/t38/t62/t47
      t130 = t61*t129
      t135 = -t116-t126-0.1088885204563779D-1*t130*t53
     #-0.4388082092272149D-1*t56*t119
      t136 = t39*t135
      t139 = t68**2
      t140 = 1/t139
      t142 = -t116-t126-0.5776579334555855D-5*t130
      t143 = t140*t142
      vrhob(i) = -0.12407009817988D1*t109*t69-0.9305257363491D0*t136
     #*t69+0.9305257363491D0*t60*t143
      t146 = 1/t6
      t149 = 0.1260013050246507D-1*t146*t7*t10
      t151 = 0.1260020456326991D-1*t16*t87
      t152 = t31*t19
      t157 = t149+t151+0.408331951711417D-2*t152*sigmaaa
     #+0.1645530784602056D-1*t21*t16
      t158 = t5*t157
      t163 = t149+t151+0.2166217250458446D-5*sigmaaa*t31
      t164 = t103*t163
      vsigmaaa(i) = -0.9305257363491D0*t158*t35+0.9305257363491D0*t26
     #*t164
      vsigmaab(i) = 0.D0
      t167 = 1/t40
      t170 = 0.1260013050246507D-1*t167*t41*t44
      t172 = 0.1260020456326991D-1*t50*t124
      t173 = t65*t53
      t178 = t170+t172+0.408331951711417D-2*t173*sigmabb
     #+0.1645530784602056D-1*t55*t50
      t179 = t39*t178
      t184 = t170+t172+0.2166217250458446D-5*sigmabb*t65
      t185 = t140*t184
      vsigmabb(i) = -0.9305257363491D0*t179*t69+0.9305257363491D0*t60
     #*t185
      t201 = 0.7840081201533821D-1*t6/t4/t80*t10
      t203 = 1/t14/t28
      t206 = 0.1680027275102654D0*sigmaaa*t203*t87
      t210 = t27/t4/t28/t80
      t212 = 1/t86/t85
      t214 = 0.4480125399532632D-1*t210*t212
      t218 = t28**2
      t233 = 1/t102/t34
      t234 = t105**2
      v2rhoa2(i) = -0.4135669939329333D0/t14*t25*t35
     #-0.24814019635976D1*t4*t98*t35+0.24814019635976D1*t72*t106
     #-0.9305257363491D0*t5*(t201+t206-t214+0.9799966841074009D-1*t210
     #*t19-0.4778117666686413D-1*t27*sigmaaa/t218/t13*t19
     #+0.1608963433833121D0*t22*t203)*t35+0.18610514726982D1*t99*t106
     #-0.18610514726982D1*t26*t233*t234+0.9305257363491D0*t26*t103*
     #(t201+t206-t214+0.3658500245218708D-4*t210)
      t256 = 0.7840081201533821D-1*t40/t38/t117*t44
      t258 = 1/t48/t62
      t261 = 0.1680027275102654D0*sigmabb*t258*t124
      t265 = t61/t38/t62/t117
      t267 = 1/t123/t122
      t269 = 0.4480125399532632D-1*t265*t267
      t273 = t62**2
      t288 = 1/t139/t68
      t289 = t142**2
      v2rhob2(i) = -0.4135669939329333D0/t48*t59*t69
     #-0.24814019635976D1*t38*t135*t69+0.24814019635976D1*t109*t143
     #-0.9305257363491D0*t39*(t256+t261-t269+0.9799966841074009D-1
     #*t265*t53-0.4778117666686413D-1*t61*sigmabb/t273/t47*t53
     #+0.1608963433833121D0*t56*t258)*t69+0.18610514726982D1*t136*t143
     #-0.18610514726982D1*t60*t288*t289+0.9305257363491D0*t60*t140*
     #(t256+t261-t269+0.3658500245218708D-4*t265)
      v2rhoab(i) = 0.D0
      t305 = 0.1680017400328676D-1*t146*t76*t10
      t307 = 0.5040081825307963D-1*t82*t87
      t308 = sigmaaa*t92
      t310 = 0.1680047024824737D-1*t308*t212
      v2rhoasigmaaa(i) = -0.12407009817988D1*t4*t157*t35
     #+0.12407009817988D1*t72*t164-0.9305257363491D0*t5*(-t305-t307
     #+t310-0.3266655613691336D-1*t92*t19*sigmaaa
     #+0.1791794125007405D-1*t27/t218/rhoa*t19-0.4388082092272149D-1
     #*t21*t82)*t35+0.9305257363491D0*t99*t164+0.9305257363491D0*t158
     #*t106-0.18610514726982D1*t26*t233*t105*t163+0.9305257363491D0
     #*t26*t103*(-t305-t307+t310-0.1155315866911171D-4*t308)
      v2rhoasigmaab(i) = 0.D0
      v2rhoasigmabb(i) = 0.D0
      v2rhobsigmaaa(i) = 0.D0
      v2rhobsigmaab(i) = 0.D0
      t345 = 0.1680017400328676D-1*t167*t113*t44
      t347 = 0.5040081825307963D-1*t119*t124
      t348 = sigmabb*t129
      t350 = 0.1680047024824737D-1*t348*t267
      v2rhobsigmabb(i) = -0.12407009817988D1*t38*t178*t69
     #+0.12407009817988D1*t109*t185-0.9305257363491D0*t39*(-t345-t347
     #+t350-0.3266655613691336D-1*t129*t53*sigmabb
     #+0.1791794125007405D-1*t61/t273/rhob*t53-0.4388082092272149D-1
     #*t55*t119)*t69+0.9305257363491D0*t136*t185+0.9305257363491D0
     #*t179*t143-0.18610514726982D1*t60*t288*t142*t184
     #+0.9305257363491D0*t60*t140*(-t345-t347+t350
     #-0.1155315866911171D-4*t348)
      t382 = 0.6300065251232535D-2/t6/sigmaaa*t7*t10
      t386 = 0.6300102281634954D-2/sigmaaa*t16*t87
      t388 = 0.6300176343092764D-2*t31*t212
      t400 = t163**2
      v2sigmaaa2(i) = -0.9305257363491D0*t5*(-t382+t386-t388
     #-0.6719227968777768D-2/t218*t19*sigmaaa+0.8166639034228341D-2
     #*t152)*t35+0.18610514726982D1*t158*t164-0.18610514726982D1*t26
     #*t233*t400+0.9305257363491D0*t26*t103*(-t382+t386-t388
     #+0.2166217250458446D-5*t31)
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      t413 = 0.6300065251232535D-2/t40/sigmabb*t41*t44
      t417 = 0.6300102281634954D-2/sigmabb*t50*t124
      t419 = 0.6300176343092764D-2*t65*t267
      t431 = t184**2
      v2sigmabb2(i) = -0.9305257363491D0*t39*(-t413+t417-t419
     #-0.6719227968777768D-2/t273*t53*sigmabb+0.8166639034228341D-2
     #*t173)*t69+0.18610514726982D1*t179*t185-0.18610514726982D1*t60
     #*t288*t431+0.9305257363491D0*t60*t140*(-t413+t417-t419
     #+0.2166217250458446D-5*t65)
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vrhob(i) = 0.0d0
      v2rhoa2(i) = 0.0d0
      v2rhob2(i) = 0.0d0
      v2rhoab(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      vsigmaab(i) = 0.0d0
      vsigmabb(i) = 0.0d0
      v2rhoasigmaaa(i) = 0.0d0
      v2rhoasigmaab(i) = 0.0d0
      v2rhoasigmabb(i) = 0.0d0
      v2rhobsigmaaa(i) = 0.0d0
      v2rhobsigmaab(i) = 0.0d0
      v2rhobsigmabb(i) = 0.0d0
      v2sigmaaa2(i) = 0.0d0
      v2sigmaab2(i) = 0.0d0
      v2sigmabb2(i) = 0.0d0
      v2sigmaaaab(i) = 0.0d0
      v2sigmaaabb(i) = 0.0d0
      v2sigmaabbb(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end
      
      
      subroutine rks_x_pw91
     & (ideriv,npt,rhoa1,sigmaaa1,
     &  zk,vrhoa,vsigmaaa,
     &  v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
c
c     J.P. Perdew, J.A. Chevary, S.H. Vosko, K.A. Jackson,
c     M.R. Pederson, D.J. Singh, C. Fiolhais
c     Atoms, molecules, solids and surfaces:
c     Applications of the generalized gradient approximation
c     for exchange and correlation
c     Phys. Rev. B 46 (1992) 6671--6687
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit real*8 (a-h,o-z)
      integer ideriv,npt
      real*8 rhoa1(npt)
      real*8 sigmaaa1(npt)
      real*8 zk(npt),vrhoa(npt),vsigmaaa(npt)
      real*8 v2rhoa2(npt),v2rhoasigmaaa(npt),v2sigmaaa2(npt)
      parameter(tol=1.0d-20)
      
      if(ideriv.eq.0) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = rho**(1.D0/3.D0)
      t3 = t2*rho
      t4 = dsqrt(sigma)
      t6 = t4/t3
      t8 = dlog(0.1259928455434599D1*t6+dsqrt(1+0.1587419712813814D1
     #*t6**2))
      t10 = 0.3175033930295641D-1*t6*t8
      t11 = rho**2
      t12 = t2**2
      t14 = 1/t12/t11
      t17 = dexp(-0.261211729852336D1*sigma*t14)
      t25 = sigma**2
      t26 = t11**2
      zk(i) = -0.7385587663820224D0*t3*(1.D0+t10+0.261211729852336D-1*
     #(0.2743D0-0.1508D0*t17)*sigma*t14)/(1.D0+t10
     #+0.272926271249799D-5*t25/t2/t26/rho)
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = rho**(1.D0/3.D0)
      t3 = t2*rho
      t4 = dsqrt(sigma)
      t5 = 1/t3
      t6 = t4*t5
      t8 = dlog(0.1259928455434599D1*t6+dsqrt(1+0.1587419712813814D1
     #*t6**2))
      t10 = 0.3175033930295641D-1*t6*t8
      t11 = rho**2
      t12 = t2**2
      t14 = 1/t12/t11
      t15 = sigma*t14
      t17 = dexp(-0.261211729852336D1*t15)
      t19 = 0.2743D0-0.1508D0*t17
      t20 = t19*sigma
      t23 = 1.D0+t10+0.261211729852336D-1*t20*t14
      t24 = t3*t23
      t25 = sigma**2
      t26 = t11**2
      t29 = 1/t2/t26/rho
      t32 = 1.D0+t10+0.272926271249799D-5*t25*t29
      t33 = 1/t32
      zk(i) = -0.7385587663820224D0*t24*t33
      t43 = 0.8466757147455043D-1*t4/t2/t11*t8
      t46 = 1/t12/t11/rho
      t50 = dsqrt(1.D0+0.1587419712813815D1*t15)
      t51 = 1/t50
      t53 = 0.1066750825533289D0*sigma*t46*t51
      t57 = t25/t2/t26/t11
      t66 = t32**2
      t67 = 1/t66
      vrhoa(i) = -0.9847450218426965D0*t2*t23*t33-0.3692793831910112D0
     #*t3*(-t43-t53-0.5487637560595959D-1*t57*t17-0.1393129225879125D0
     #*t20*t46)*t33+0.3692793831910112D0*t24*t67*(-t43-t53
     #-0.2911213559997856D-4*t57)
      t76 = 0.6350067860591282D-1/t4*t5*t8
      t78 = 0.8000631191499664D-1*t14*t51
      vsigmaaa(i) = -0.7385587663820224D0*t3*(t76+t78
     #+0.411572817044697D-1*t29*t17*sigma+0.1044846919409344D0*t19*t14
     #)*t33+0.7385587663820224D0*t24*t67*(t76+t78
     #+0.2183410169998392D-4*sigma*t29)
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = rho**(1.D0/3.D0)
      t3 = t2*rho
      t4 = dsqrt(sigma)
      t5 = 1/t3
      t6 = t4*t5
      t8 = dlog(0.1259928455434599D1*t6+dsqrt(1+0.1587419712813814D1
     #*t6**2))
      t10 = 0.3175033930295641D-1*t6*t8
      t11 = rho**2
      t12 = t2**2
      t14 = 1/t12/t11
      t15 = sigma*t14
      t17 = dexp(-0.261211729852336D1*t15)
      t19 = 0.2743D0-0.1508D0*t17
      t20 = t19*sigma
      t23 = 1.D0+t10+0.261211729852336D-1*t20*t14
      t24 = t3*t23
      t25 = sigma**2
      t26 = t11**2
      t29 = 1/t2/t26/rho
      t32 = 1.D0+t10+0.272926271249799D-5*t25*t29
      t33 = 1/t32
      zk(i) = -0.7385587663820224D0*t24*t33
      t36 = t2*t23
      t40 = 1/t2/t11
      t43 = 0.8466757147455043D-1*t4*t40*t8
      t44 = t11*rho
      t46 = 1/t12/t44
      t49 = 1.D0+0.1587419712813815D1*t15
      t50 = dsqrt(t49)
      t51 = 1/t50
      t53 = 0.1066750825533289D0*sigma*t46*t51
      t56 = 1/t2/t26/t11
      t57 = t25*t56
      t62 = -t43-t53-0.5487637560595959D-1*t57*t17
     #-0.1393129225879125D0*t20*t46
      t63 = t3*t62
      t66 = t32**2
      t67 = 1/t66
      t69 = -t43-t53-0.2911213559997856D-4*t57
      t70 = t67*t69
      vrhoa(i) = -0.9847450218426965D0*t36*t33-0.3692793831910112D0
     #*t63*t33+0.3692793831910112D0*t24*t70
      t73 = 1/t4
      t76 = 0.6350067860591282D-1*t73*t5*t8
      t78 = 0.8000631191499664D-1*t14*t51
      t79 = t29*t17
      t84 = t76+t78+0.411572817044697D-1*t79*sigma
     #+0.1044846919409344D0*t19*t14
      t85 = t3*t84
      t90 = t76+t78+0.2183410169998392D-4*sigma*t29
      t91 = t67*t90
      vsigmaaa(i) = -0.7385587663820224D0*t85*t33+0.7385587663820224D0
     #*t24*t91
      t107 = 0.395115333547902D0*t4/t2/t44*t8
      t109 = 1/t12/t26
      t112 = 0.1066750825533289D1*sigma*t109*t51
      t116 = t25/t2/t26/t44
      t118 = 1/t50/t49
      t120 = 0.4515683437631874D0*t116*t118
      t124 = t26**2
      t139 = 1/t66/t32
      t140 = t69**2
      v2rhoa2(i) = -0.6564966812284644D0/t12*t23*t33
     #-0.1969490043685393D1*t2*t62*t33+0.1969490043685393D1*t36*t70
     #-0.3692793831910112D0*t3*(t107+t112-t120+0.9877747609072727D0
     #*t116*t17-0.764498826669826D0*t25*sigma/t124/t11*t17
     #+0.1021628098978025D1*t20*t109)*t33+0.7385587663820224D0*t63*t70
     #-0.7385587663820224D0*t24*t139*t140+0.3692793831910112D0*t24*t67
     #*(t107+t112-t120+0.3687537175997285D-3*t116)
      t156 = 0.1693351429491009D0*t73*t40*t8
      t158 = 0.6400504953199731D0*t46*t51
      t159 = sigma*t56
      t161 = 0.3386762578223905D0*t159*t118
      v2rhoasigmaaa(i) = -0.9847450218426965D0*t2*t84*t33
     #+0.9847450218426965D0*t36*t91-0.3692793831910112D0*t3*(-t156
     #-t158+t161-0.6585165072715151D0*t56*t17*sigma
     #+0.5733741200023695D0*t25/t124/rho*t17-0.5572516903516501D0*t19
     #*t46)*t33+0.3692793831910112D0*t63*t91+0.3692793831910112D0*t85
     #*t70-0.7385587663820224D0*t24*t139*t69*t90+0.3692793831910112D0
     #*t24*t67*(-t156-t158+t161-0.2328970847998285D-3*t159)
      t193 = 0.1270013572118256D0/t4/sigma*t5*t8
      t197 = 0.1600126238299933D0/sigma*t14*t51
      t199 = 0.2540071933667929D0*t29*t118
      t211 = t90**2
      v2sigmaaa2(i) = -0.7385587663820224D0*t3*(-t193+t197-t199
     #-0.4300305900017772D0/t124*t17*sigma+0.3292582536357576D0*t79)
     #*t33+0.1477117532764045D1*t85*t91-0.1477117532764045D1*t24*t139
     #*t211+0.7385587663820224D0*t24*t67*(-t193+t197-t199
     #+0.8733640679993569D-4*t29)
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      v2rhoa2(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      v2rhoasigmaaa(i) = 0.0d0
      v2sigmaaa2(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end

      
      subroutine uks_c_vwn5
     & (ideriv,npt,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,
     &  zk,vrhoa,vrhob,vsigmaaa,vsigmabb,vsigmaab,
     &  v2rhoa2,v2rhob2,v2rhoab,
     &  v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb,
     &  v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,
     &  v2sigmaaa2,v2sigmaaaab,v2sigmaaabb,
     &  v2sigmaab2,v2sigmaabbb,v2sigmabb2)
c
c     S.H. Vosko, L. Wilk, and M. Nusair
c     Accurate spin-dependent electron liquid correlation energies for
c     local spin density calculations: a critical analysis
c     Can. J. Phys. 58 (1980) 1200-1211
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit real*8 (a-h,o-z)
      integer ideriv,npt
      real*8 rhoa1(npt),rhob1(npt)
      real*8 sigmaaa1(npt),sigmabb1(npt),sigmaab1(npt)
      real*8 zk(npt),vrhoa(npt),vrhob(npt)
      real*8 vsigmaaa(npt),vsigmabb(npt),vsigmaab(npt)
      real*8 v2rhoa2(npt),v2rhob2(npt),v2rhoab(npt)
      real*8 v2rhoasigmaaa(npt),v2rhoasigmaab(npt)
      real*8 v2rhoasigmabb(npt),v2rhobsigmabb(npt)
      real*8 v2rhobsigmaab(npt),v2rhobsigmaaa(npt)
      real*8 v2sigmaaa2(npt),v2sigmaaaab(npt),v2sigmaaabb(npt)
      real*8 v2sigmaab2(npt),v2sigmaabbb(npt),v2sigmabb2(npt)
      parameter(tol=1.0d-20)
      
      if (ideriv.eq.0) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      t1 = 1/rhob
      t2 = t1**(1.D0/6.D0)
      t3 = 0.7876233178997433D0*t2
      t5 = (t3+0.10498D0)**2
      t6 = t1**(1.D0/3.D0)
      t7 = 0.6203504908994D0*t6
      t12 = dlog(t5/(t7+0.2935818660072219D1*t2+0.129352D2))
      t16 = 1/(t7+0.5560951426165705D1*t2+0.180578D2)
      t19 = dlog(0.6203504908994D0*t6*t16)
      t25 = datan(0.4730926909560113D1/(0.1575246635799487D1*t2
     #+0.706042D1))
      t28 = (t3+0.325D0)**2
      t30 = dlog(t28*t16)
      zk(i) = rhob*(0.1D-23*t12+0.1554535D-1*t19+0.5249139316978094D-1
     #*t25+0.2247867095542612D-2*t30)
      elseif(rhob.lt.tol) then
      rho = rhoa
      t1 = 1/rhoa
      t2 = t1**(1.D0/6.D0)
      t3 = 0.7876233178997433D0*t2
      t5 = (t3+0.10498D0)**2
      t6 = t1**(1.D0/3.D0)
      t7 = 0.6203504908994D0*t6
      t12 = dlog(t5/(t7+0.2935818660072219D1*t2+0.129352D2))
      t16 = 1/(t7+0.5560951426165705D1*t2+0.180578D2)
      t19 = dlog(0.6203504908994D0*t6*t16)
      t25 = datan(0.4730926909560113D1/(0.1575246635799487D1*t2
     #+0.706042D1))
      t28 = (t3+0.325D0)**2
      t30 = dlog(t28*t16)
      zk(i) = rhoa*(0.1D-23*t12+0.1554535D-1*t19+0.5249139316978094D-1
     #*t25+0.2247867095542612D-2*t30)
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      t1 = 1/rho
      t2 = t1**(1.D0/3.D0)
      t3 = 0.6203504908994D0*t2
      t4 = t1**(1.D0/6.D0)
      t7 = 1/(t3+0.2935818660072219D1*t4+0.129352D2)
      t10 = dlog(0.6203504908994D0*t2*t7)
      t12 = 0.1575246635799487D1*t4
      t16 = datan(0.615199081975908D1/(t12+0.372744D1))
      t18 = 0.7876233178997433D0*t4
      t20 = (t18+0.10498D0)**2
      t22 = dlog(t20*t7)
      t26 = 1/(t3+0.8908571061768626D0*t4+0.130045D2)
      t29 = dlog(0.6203504908994D0*t2*t26)
      t33 = datan(0.7123108917818118D1/(t12+0.113107D1))
      t36 = (t18+0.47584D-2)**2
      t38 = dlog(t36*t26)
      t40 = t29+0.3177080047439415D0*t33+0.4140337942820627D-3*t38
      t42 = rhoa-1.D0*rhob
      t43 = t42*t1
      t44 = 1.D0+t43
      t45 = t44**(1.D0/3.D0)
      t48 = 1.D0-1.D0*t43
      t49 = t48**(1.D0/3.D0)
      t55 = 1/(t3+0.5560951426165705D1*t4+0.180578D2)
      t58 = dlog(0.6203504908994D0*t2*t55)
      t63 = datan(0.4730926909560113D1/(t12+0.706042D1))
      t66 = (t18+0.325D0)**2
      t68 = dlog(t66*t55)
      t78 = t42**2
      t79 = t78**2
      t81 = rho**2
      t82 = t81**2
      zk(i) = rho*(0.310907D-1*t10+0.3878329487811301D-1*t16
     #+0.9690227711544374D-3*t22-0.1899772193293833D-1*t40*(t45*t44
     #+t49*t48-2.D0)*(1.D0+(-0.101257459063883D3*(0.1554535D-1*t58
     #+0.5249139316978094D-1*t63+0.2247867095542612D-2*t68-0.310907D-1
     #*t10-0.3878329487811301D-1*t16-0.9690227711544374D-3*t22)/t40
     #-1.D0)*t79/t82))
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      t1 = 1/rhob
      t2 = t1**(1.D0/6.D0)
      t3 = 0.7876233178997433D0*t2
      t4 = t3+0.10498D0
      t5 = t4**2
      t6 = t1**(1.D0/3.D0)
      t7 = 0.6203504908994D0*t6
      t9 = t7+0.2935818660072219D1*t2+0.129352D2
      t10 = 1/t9
      t12 = dlog(t5*t10)
      t13 = 0.1D-23*t12
      t15 = t7+0.5560951426165705D1*t2+0.180578D2
      t16 = 1/t15
      t19 = dlog(0.6203504908994D0*t6*t16)
      t20 = 0.1554535D-1*t19
      t21 = 0.1575246635799487D1*t2
      t22 = t21+0.706042D1
      t25 = datan(0.4730926909560113D1/t22)
      t26 = 0.5249139316978094D-1*t25
      t27 = t3+0.325D0
      t28 = t27**2
      t30 = dlog(t28*t16)
      t31 = 0.2247867095542612D-2*t30
      zk(i) = rhob*(t13+t20+t26+t31)
      vrhoa(i) = 0.D0
      t33 = t6**2
      t34 = 1/t33
      t36 = rhob**2
      t37 = 1/t36
      t40 = t9**2
      t41 = 1/t40
      t44 = 0.2067834969664667D0*t34*t37
      t45 = t2**2
      t46 = t45**2
      t48 = 1/t46/t2
      t49 = t48*t37
      t51 = -t44-0.4893031100120365D0*t49
      t55 = 1/t6
      t60 = (t21+0.372744D1)**2
      t61 = 1/t60
      t83 = t15**2
      t84 = 1/t83
      t87 = -t44-0.9268252376942842D0*t49
      t94 = t22**2
      t95 = 1/t94
      vrhob(i) = t13+t20+t26+t31+rhob*(0.2D-21*(-0.2067834969664667D0
     #*t34*t10*t37-0.6203504908994D0*t6*t41*t51)*t55*t9-0.1D-21*t61
     #*t48*t37/(1.D0+0.378469910464D2*t61)+0.1D-23*(
     #-0.2625411059665811D0*t4*t10*t49-1.D0*t5*t41*t51)/t5*t9
     #+0.2505897912236993D-1*(-0.2067834969664667D0*t34*t16*t37
     #-0.6203504908994D0*t6*t84*t87)*t55*t15+0.6519760588836059D-1*t95
     #*t48*t37/(1.D0+0.223816694236D2*t95)+0.2247867095542612D-2*(
     #-0.2625411059665811D0*t27*t16*t49-1.D0*t28*t84*t87)/t28*t15)
      elseif(rhob.lt.tol) then
      rho = rhoa
      t1 = 1/rhoa
      t2 = t1**(1.D0/6.D0)
      t3 = 0.7876233178997433D0*t2
      t4 = t3+0.10498D0
      t5 = t4**2
      t6 = t1**(1.D0/3.D0)
      t7 = 0.6203504908994D0*t6
      t9 = t7+0.2935818660072219D1*t2+0.129352D2
      t10 = 1/t9
      t12 = dlog(t5*t10)
      t13 = 0.1D-23*t12
      t15 = t7+0.5560951426165705D1*t2+0.180578D2
      t16 = 1/t15
      t19 = dlog(0.6203504908994D0*t6*t16)
      t20 = 0.1554535D-1*t19
      t21 = 0.1575246635799487D1*t2
      t22 = t21+0.706042D1
      t25 = datan(0.4730926909560113D1/t22)
      t26 = 0.5249139316978094D-1*t25
      t27 = t3+0.325D0
      t28 = t27**2
      t30 = dlog(t28*t16)
      t31 = 0.2247867095542612D-2*t30
      zk(i) = rhoa*(t13+t20+t26+t31)
      t33 = t6**2
      t34 = 1/t33
      t36 = rhoa**2
      t37 = 1/t36
      t40 = t9**2
      t41 = 1/t40
      t44 = 0.2067834969664667D0*t34*t37
      t45 = t2**2
      t46 = t45**2
      t48 = 1/t46/t2
      t49 = t48*t37
      t51 = -t44-0.4893031100120365D0*t49
      t55 = 1/t6
      t60 = (t21+0.372744D1)**2
      t61 = 1/t60
      t83 = t15**2
      t84 = 1/t83
      t87 = -t44-0.9268252376942842D0*t49
      t94 = t22**2
      t95 = 1/t94
      vrhoa(i) = t13+t20+t26+t31+rhoa*(0.2D-21*(-0.2067834969664667D0
     #*t34*t10*t37-0.6203504908994D0*t6*t41*t51)*t55*t9-0.1D-21*t61
     #*t48*t37/(1.D0+0.378469910464D2*t61)+0.1D-23*(
     #-0.2625411059665811D0*t4*t10*t49-1.D0*t5*t41*t51)/t5*t9
     #+0.2505897912236993D-1*(-0.2067834969664667D0*t34*t16*t37
     #-0.6203504908994D0*t6*t84*t87)*t55*t15+0.6519760588836059D-1*t95
     #*t48*t37/(1.D0+0.223816694236D2*t95)+0.2247867095542612D-2*(
     #-0.2625411059665811D0*t27*t16*t49-1.D0*t28*t84*t87)/t28*t15)
      vrhob(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      t1 = 1/rho
      t2 = t1**(1.D0/3.D0)
      t3 = 0.6203504908994D0*t2
      t4 = t1**(1.D0/6.D0)
      t6 = t3+0.2935818660072219D1*t4+0.129352D2
      t7 = 1/t6
      t10 = dlog(0.6203504908994D0*t2*t7)
      t11 = 0.310907D-1*t10
      t12 = 0.1575246635799487D1*t4
      t13 = t12+0.372744D1
      t16 = datan(0.615199081975908D1/t13)
      t17 = 0.3878329487811301D-1*t16
      t18 = 0.7876233178997433D0*t4
      t19 = t18+0.10498D0
      t20 = t19**2
      t22 = dlog(t20*t7)
      t23 = 0.9690227711544374D-3*t22
      t25 = t3+0.8908571061768626D0*t4+0.130045D2
      t26 = 1/t25
      t29 = dlog(0.6203504908994D0*t2*t26)
      t30 = t12+0.113107D1
      t33 = datan(0.7123108917818118D1/t30)
      t35 = t18+0.47584D-2
      t36 = t35**2
      t38 = dlog(t36*t26)
      t40 = t29+0.3177080047439415D0*t33+0.4140337942820627D-3*t38
      t42 = rhoa-1.D0*rhob
      t43 = t42*t1
      t44 = 1.D0+t43
      t45 = t44**(1.D0/3.D0)
      t48 = 1.D0-1.D0*t43
      t49 = t48**(1.D0/3.D0)
      t51 = t45*t44+t49*t48-2.D0
      t52 = t40*t51
      t54 = t3+0.5560951426165705D1*t4+0.180578D2
      t55 = 1/t54
      t58 = dlog(0.6203504908994D0*t2*t55)
      t60 = t12+0.706042D1
      t63 = datan(0.4730926909560113D1/t60)
      t65 = t18+0.325D0
      t66 = t65**2
      t68 = dlog(t66*t55)
      t73 = 0.1554535D-1*t58+0.5249139316978094D-1*t63
     #+0.2247867095542612D-2*t68-0.310907D-1*t10-0.3878329487811301D-1
     #*t16-0.9690227711544374D-3*t22
      t74 = 1/t40
      t77 = -0.101257459063883D3*t73*t74-1.D0
      t78 = t42**2
      t79 = t78**2
      t80 = t77*t79
      t81 = rho**2
      t82 = t81**2
      t83 = 1/t82
      t85 = 1.D0+t80*t83
      t87 = 0.1899772193293833D-1*t52*t85
      zk(i) = rho*(t11+t17+t23-t87)
      t93 = 0.1333333333333333D1*t45*t1-0.1333333333333333D1*t49*t1
      t100 = t52*t77*t78*t42*t83
      t104 = t2**2
      t105 = 1/t104
      t107 = 1/t81
      t110 = t6**2
      t111 = 1/t110
      t114 = 0.2067834969664667D0*t105*t107
      t115 = t4**2
      t116 = t115**2
      t118 = 1/t116/t4
      t119 = t118*t107
      t121 = -t114-0.4893031100120365D0*t119
      t125 = 1/t2
      t127 = (-0.2067834969664667D0*t105*t7*t107-0.6203504908994D0*t2
     #*t111*t121)*t125*t6
      t129 = t13**2
      t130 = 1/t129
      t136 = t130*t118*t107/(1.D0+0.378469910464D2*t130)
      t147 = (-0.2625411059665811D0*t19*t7*t119-1.D0*t20*t111*t121)
     #/t20*t6
      t152 = t25**2
      t153 = 1/t152
      t156 = -t114-0.1484761843628104D0*t119
      t163 = t30**2
      t164 = 1/t163
      t183 = 0.161199195401647D1*(-0.2067834969664667D0*t105*t26*t107
     #-0.6203504908994D0*t2*t153*t156)*t125*t25+0.5941485651138659D0
     #*t164*t118*t107/(1.D0+0.507386806551D2*t164)
     #+0.4140337942820627D-3*(-0.2625411059665811D0*t35*t26*t119-1.D0
     #*t36*t153*t156)/t36*t25
      t200 = t54**2
      t201 = 1/t200
      t204 = -t114-0.9268252376942842D0*t119
      t211 = t60**2
      t212 = 1/t211
      t237 = t40**2
      t253 = rho*(0.5011795824473985D-1*t127+0.626408570946439D-1*t136
     #+0.9690227711544374D-3*t147-0.1899772193293833D-1*t183*t51*t85
     #-0.1899772193293833D-1*t40*(-0.1333333333333333D1*t45*t42*t107
     #+0.1333333333333333D1*t49*t42*t107)*t85-0.1899772193293833D-1
     #*t52*((-0.101257459063883D3*(0.2505897912236993D-1*(
     #-0.2067834969664667D0*t105*t55*t107-0.6203504908994D0*t2*t201
     #*t204)*t125*t54+0.6519760588836059D-1*t212*t118*t107/(1.D0
     #+0.223816694236D2*t212)+0.2247867095542612D-2*(
     #-0.2625411059665811D0*t65*t55*t119-1.D0*t66*t201*t204)/t66*t54
     #-0.5011795824473985D-1*t127-0.626408570946439D-1*t136
     #-0.9690227711544374D-3*t147)*t74+0.101257459063883D3*t73/t237
     #*t183)*t79*t83-4.D0*t80/t82/rho))
      vrhoa(i) = rho*(-0.1899772193293833D-1*t40*t93*t85
     #-0.7599088773175333D-1*t100)+t11+t17+t23-t87+t253
      vrhob(i) = rho*(0.1899772193293833D-1*t40*t93*t85
     #+0.7599088773175333D-1*t100)+t11+t17+t23-t87+t253
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vrhob(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      t1 = 1/rhob
      t2 = t1**(1.D0/6.D0)
      t3 = 0.7876233178997433D0*t2
      t4 = t3+0.10498D0
      t5 = t4**2
      t6 = t1**(1.D0/3.D0)
      t7 = 0.6203504908994D0*t6
      t9 = t7+0.2935818660072219D1*t2+0.129352D2
      t10 = 1/t9
      t12 = dlog(t5*t10)
      t13 = 0.1D-23*t12
      t15 = t7+0.5560951426165705D1*t2+0.180578D2
      t16 = 1/t15
      t19 = dlog(0.6203504908994D0*t6*t16)
      t20 = 0.1554535D-1*t19
      t21 = 0.1575246635799487D1*t2
      t22 = t21+0.706042D1
      t25 = datan(0.4730926909560113D1/t22)
      t26 = 0.5249139316978094D-1*t25
      t27 = t3+0.325D0
      t28 = t27**2
      t30 = dlog(t28*t16)
      t31 = 0.2247867095542612D-2*t30
      zk(i) = rhob*(t13+t20+t26+t31)
      vrhoa(i) = 0.D0
      t33 = t6**2
      t34 = 1/t33
      t35 = t34*t10
      t36 = rhob**2
      t37 = 1/t36
      t40 = t9**2
      t41 = 1/t40
      t42 = t6*t41
      t44 = 0.2067834969664667D0*t34*t37
      t45 = t2**2
      t46 = t45**2
      t47 = t46*t2
      t48 = 1/t47
      t49 = t48*t37
      t51 = -t44-0.4893031100120365D0*t49
      t55 = 1/t6
      t56 = (-0.2067834969664667D0*t35*t37-0.6203504908994D0*t42*t51)
     #*t55
      t60 = (t21+0.372744D1)**2
      t61 = 1/t60
      t65 = 1/(1.D0+0.378469910464D2*t61)
      t69 = t4*t10
      t72 = t5*t41
      t76 = 1/t5
      t77 = (-0.2625411059665811D0*t69*t49-1.D0*t72*t51)*t76
      t80 = t34*t16
      t83 = t15**2
      t84 = 1/t83
      t85 = t6*t84
      t87 = -t44-0.9268252376942842D0*t49
      t90 = -0.2067834969664667D0*t80*t37-0.6203504908994D0*t85*t87
      t91 = t90*t55
      t92 = t91*t15
      t94 = t22**2
      t95 = 1/t94
      t96 = t95*t48
      t98 = 1.D0+0.223816694236D2*t95
      t99 = 1/t98
      t101 = t96*t37*t99
      t103 = t27*t16
      t106 = t28*t84
      t109 = -0.2625411059665811D0*t103*t49-1.D0*t106*t87
      t110 = 1/t28
      t111 = t109*t110
      t112 = t111*t15
      vrhob(i) = t13+t20+t26+t31+rhob*(0.2D-21*t56*t9-0.1D-21*t61*t48
     #*t37*t65+0.1D-23*t77*t9+0.2505897912236993D-1*t92
     #+0.6519760588836059D-1*t101+0.2247867095542612D-2*t112)
      v2rhoa2(i) = 0.D0
      v2rhoab(i) = 0.D0
      t120 = 1/t47/t1
      t122 = t36**2
      t123 = 1/t122
      t128 = 1/t33/t1
      t130 = t128*t10*t123
      t137 = 1/t36/rhob
      t141 = 1/t40/t9
      t143 = t51**2
      t147 = 0.1378556646443111D0*t128*t123
      t149 = 0.4135669939329333D0*t34*t137
      t150 = t120*t123
      t152 = t48*t137
      t154 = -t147+t149-0.4077525916766971D0*t150+0.978606220024073D0
     #*t152
      t192 = t123*t99
      t201 = t94**2
      t205 = t98**2
      t211 = t128*t16*t123
      t220 = 1/t83/t15
      t222 = t87**2
      t227 = -t147+t149-0.7723543647452368D0*t150+0.1853650475388568D1
     #*t152
      s1 = 0.1D-21*t61*t120*t123*t65+0.2D-21*(-0.1378556646443111D0
     #*t130+0.4135669939329333D0*t34*t41*t37*t51+0.4135669939329333D0
     #*t35*t137+0.12407009817988D1*t6*t141*t143-0.6203504908994D0*t42
     #*t154)*t55*t9+0.2D-21*t56*t51+0.1D-23*(0.3446391616107778D-1
     #*t130+0.5250822119331622D0*t4*t41*t49*t51-0.2187842549721509D0
     #*t69*t150+0.5250822119331622D0*t69*t152+2.D0*t5*t141*t143-1.D0
     #*t72*t154)*t76*t9+0.1D-23*t77*t51+0.8352993040789975D-2*t90/t6
     #/t1*t15*t37+0.3423410311260694D-1/t94/t22*t128*t192
      t263 = s1+0.5433133824030049D-1*t95*t120*t192
     #-0.1303952117767212D0*t96*t137*t99-0.7662163788798043D0/t201/t22
     #*t128*t123/t205+0.2505897912236993D-1*(-0.1378556646443111D0
     #*t211+0.4135669939329333D0*t34*t84*t37*t87+0.4135669939329333D0
     #*t80*t137+0.12407009817988D1*t6*t220*t222-0.6203504908994D0*t85
     #*t227)*t55*t15+0.2505897912236993D-1*t91*t87
     #+0.2247867095542612D-2*(0.3446391616107778D-1*t211
     #+0.5250822119331622D0*t27*t84*t49*t87-0.2187842549721509D0*t103
     #*t150+0.5250822119331622D0*t103*t152+2.D0*t28*t220*t222-1.D0
     #*t106*t227)*t110*t15+0.5901575133296437D-3*t109/t28/t27*t15*t48
     #*t37+0.2247867095542612D-2*t111*t87
      v2rhob2(i) = 0.5011795824473985D-1*t92+0.1303952117767212D0*t101
     #+0.4495734191085224D-2*t112+rhob*t263
      elseif(rhob.lt.tol) then
      rho = rhoa
      t1 = 1/rhoa
      t2 = t1**(1.D0/6.D0)
      t3 = 0.7876233178997433D0*t2
      t4 = t3+0.10498D0
      t5 = t4**2
      t6 = t1**(1.D0/3.D0)
      t7 = 0.6203504908994D0*t6
      t9 = t7+0.2935818660072219D1*t2+0.129352D2
      t10 = 1/t9
      t12 = dlog(t5*t10)
      t13 = 0.1D-23*t12
      t15 = t7+0.5560951426165705D1*t2+0.180578D2
      t16 = 1/t15
      t19 = dlog(0.6203504908994D0*t6*t16)
      t20 = 0.1554535D-1*t19
      t21 = 0.1575246635799487D1*t2
      t22 = t21+0.706042D1
      t25 = datan(0.4730926909560113D1/t22)
      t26 = 0.5249139316978094D-1*t25
      t27 = t3+0.325D0
      t28 = t27**2
      t30 = dlog(t28*t16)
      t31 = 0.2247867095542612D-2*t30
      zk(i) = rhoa*(t13+t20+t26+t31)
      t33 = t6**2
      t34 = 1/t33
      t35 = t34*t10
      t36 = rhoa**2
      t37 = 1/t36
      t40 = t9**2
      t41 = 1/t40
      t42 = t6*t41
      t44 = 0.2067834969664667D0*t34*t37
      t45 = t2**2
      t46 = t45**2
      t47 = t46*t2
      t48 = 1/t47
      t49 = t48*t37
      t51 = -t44-0.4893031100120365D0*t49
      t55 = 1/t6
      t56 = (-0.2067834969664667D0*t35*t37-0.6203504908994D0*t42*t51)
     #*t55
      t60 = (t21+0.372744D1)**2
      t61 = 1/t60
      t65 = 1/(1.D0+0.378469910464D2*t61)
      t69 = t4*t10
      t72 = t5*t41
      t76 = 1/t5
      t77 = (-0.2625411059665811D0*t69*t49-1.D0*t72*t51)*t76
      t80 = t34*t16
      t83 = t15**2
      t84 = 1/t83
      t85 = t6*t84
      t87 = -t44-0.9268252376942842D0*t49
      t90 = -0.2067834969664667D0*t80*t37-0.6203504908994D0*t85*t87
      t91 = t90*t55
      t92 = t91*t15
      t94 = t22**2
      t95 = 1/t94
      t96 = t95*t48
      t98 = 1.D0+0.223816694236D2*t95
      t99 = 1/t98
      t101 = t96*t37*t99
      t103 = t27*t16
      t106 = t28*t84
      t109 = -0.2625411059665811D0*t103*t49-1.D0*t106*t87
      t110 = 1/t28
      t111 = t109*t110
      t112 = t111*t15
      vrhoa(i) = t13+t20+t26+t31+rhoa*(0.2D-21*t56*t9-0.1D-21*t61*t48
     #*t37*t65+0.1D-23*t77*t9+0.2505897912236993D-1*t92
     #+0.6519760588836059D-1*t101+0.2247867095542612D-2*t112)
      vrhob(i) = 0.D0
      t120 = 1/t47/t1
      t122 = t36**2
      t123 = 1/t122
      t128 = 1/t33/t1
      t130 = t128*t10*t123
      t137 = 1/t36/rhoa
      t141 = 1/t40/t9
      t143 = t51**2
      t147 = 0.1378556646443111D0*t128*t123
      t149 = 0.4135669939329333D0*t34*t137
      t150 = t120*t123
      t152 = t48*t137
      t154 = -t147+t149-0.4077525916766971D0*t150+0.978606220024073D0
     #*t152
      t192 = t123*t99
      t201 = t94**2
      t205 = t98**2
      t211 = t128*t16*t123
      t220 = 1/t83/t15
      t222 = t87**2
      t227 = -t147+t149-0.7723543647452368D0*t150+0.1853650475388568D1
     #*t152
      s1 = 0.1D-21*t61*t120*t123*t65+0.2D-21*(-0.1378556646443111D0
     #*t130+0.4135669939329333D0*t34*t41*t37*t51+0.4135669939329333D0
     #*t35*t137+0.12407009817988D1*t6*t141*t143-0.6203504908994D0*t42
     #*t154)*t55*t9+0.2D-21*t56*t51+0.1D-23*(0.3446391616107778D-1
     #*t130+0.5250822119331622D0*t4*t41*t49*t51-0.2187842549721509D0
     #*t69*t150+0.5250822119331622D0*t69*t152+2.D0*t5*t141*t143-1.D0
     #*t72*t154)*t76*t9+0.1D-23*t77*t51+0.8352993040789975D-2*t90/t6
     #/t1*t15*t37+0.3423410311260694D-1/t94/t22*t128*t192
      t263 = s1+0.5433133824030049D-1*t95*t120*t192
     #-0.1303952117767212D0*t96*t137*t99-0.7662163788798043D0/t201/t22
     #*t128*t123/t205+0.2505897912236993D-1*(-0.1378556646443111D0
     #*t211+0.4135669939329333D0*t34*t84*t37*t87+0.4135669939329333D0
     #*t80*t137+0.12407009817988D1*t6*t220*t222-0.6203504908994D0*t85
     #*t227)*t55*t15+0.2505897912236993D-1*t91*t87
     #+0.2247867095542612D-2*(0.3446391616107778D-1*t211
     #+0.5250822119331622D0*t27*t84*t49*t87-0.2187842549721509D0*t103
     #*t150+0.5250822119331622D0*t103*t152+2.D0*t28*t220*t222-1.D0
     #*t106*t227)*t110*t15+0.5901575133296437D-3*t109/t28/t27*t15*t48
     #*t37+0.2247867095542612D-2*t111*t87
      v2rhoa2(i) = 0.5011795824473985D-1*t92+0.1303952117767212D0*t101
     #+0.4495734191085224D-2*t112+rhoa*t263
      v2rhob2(i) = 0.D0
      v2rhoab(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      t1 = 1/rho
      t2 = t1**(1.D0/3.D0)
      t3 = 0.6203504908994D0*t2
      t4 = t1**(1.D0/6.D0)
      t6 = t3+0.2935818660072219D1*t4+0.129352D2
      t7 = 1/t6
      t10 = dlog(0.6203504908994D0*t2*t7)
      t11 = 0.310907D-1*t10
      t12 = 0.1575246635799487D1*t4
      t13 = t12+0.372744D1
      t16 = datan(0.615199081975908D1/t13)
      t17 = 0.3878329487811301D-1*t16
      t18 = 0.7876233178997433D0*t4
      t19 = t18+0.10498D0
      t20 = t19**2
      t22 = dlog(t20*t7)
      t23 = 0.9690227711544374D-3*t22
      t25 = t3+0.8908571061768626D0*t4+0.130045D2
      t26 = 1/t25
      t29 = dlog(0.6203504908994D0*t2*t26)
      t30 = t12+0.113107D1
      t33 = datan(0.7123108917818118D1/t30)
      t35 = t18+0.47584D-2
      t36 = t35**2
      t38 = dlog(t36*t26)
      t40 = t29+0.3177080047439415D0*t33+0.4140337942820627D-3*t38
      t42 = rhoa-1.D0*rhob
      t43 = t42*t1
      t44 = 1.D0+t43
      t45 = t44**(1.D0/3.D0)
      t48 = 1.D0-1.D0*t43
      t49 = t48**(1.D0/3.D0)
      t51 = t45*t44+t49*t48-2.D0
      t52 = t40*t51
      t54 = t3+0.5560951426165705D1*t4+0.180578D2
      t55 = 1/t54
      t58 = dlog(0.6203504908994D0*t2*t55)
      t60 = t12+0.706042D1
      t63 = datan(0.4730926909560113D1/t60)
      t65 = t18+0.325D0
      t66 = t65**2
      t68 = dlog(t66*t55)
      t73 = 0.1554535D-1*t58+0.5249139316978094D-1*t63
     #+0.2247867095542612D-2*t68-0.310907D-1*t10-0.3878329487811301D-1
     #*t16-0.9690227711544374D-3*t22
      t74 = 1/t40
      t77 = -0.101257459063883D3*t73*t74-1.D0
      t78 = t42**2
      t79 = t78**2
      t80 = t77*t79
      t81 = rho**2
      t82 = t81**2
      t83 = 1/t82
      t85 = 1.D0+t80*t83
      t87 = 0.1899772193293833D-1*t52*t85
      zk(i) = rho*(t11+t17+t23-t87)
      t93 = 0.1333333333333333D1*t45*t1-0.1333333333333333D1*t49*t1
      t94 = t40*t93
      t95 = t94*t85
      t96 = 0.1899772193293833D-1*t95
      t97 = t78*t42
      t98 = t77*t97
      t99 = t98*t83
      t100 = t52*t99
      t104 = t2**2
      t105 = 1/t104
      t106 = t105*t7
      t107 = 1/t81
      t110 = t6**2
      t111 = 1/t110
      t112 = t2*t111
      t114 = 0.2067834969664667D0*t105*t107
      t115 = t4**2
      t116 = t115**2
      t117 = t116*t4
      t118 = 1/t117
      t119 = t118*t107
      t121 = -t114-0.4893031100120365D0*t119
      t124 = -0.2067834969664667D0*t106*t107-0.6203504908994D0*t112*t121
      t125 = 1/t2
      t126 = t124*t125
      t127 = t126*t6
      t129 = t13**2
      t130 = 1/t129
      t131 = t130*t118
      t133 = 1.D0+0.378469910464D2*t130
      t134 = 1/t133
      t136 = t131*t107*t134
      t138 = t19*t7
      t141 = t20*t111
      t144 = -0.2625411059665811D0*t138*t119-1.D0*t141*t121
      t145 = 1/t20
      t146 = t144*t145
      t147 = t146*t6
      t149 = t105*t26
      t152 = t25**2
      t153 = 1/t152
      t154 = t2*t153
      t156 = -t114-0.1484761843628104D0*t119
      t159 = -0.2067834969664667D0*t149*t107-0.6203504908994D0*t154*t156
      t160 = t159*t125
      t163 = t30**2
      t164 = 1/t163
      t165 = t164*t118
      t167 = 1.D0+0.507386806551D2*t164
      t168 = 1/t167
      t172 = t35*t26
      t175 = t36*t153
      t178 = -0.2625411059665811D0*t172*t119-1.D0*t175*t156
      t179 = 1/t36
      t180 = t178*t179
      t183 = 0.161199195401647D1*t160*t25+0.5941485651138659D0*t165
     #*t107*t168+0.4140337942820627D-3*t180*t25
      t184 = t183*t51
      t185 = t184*t85
      t187 = t45*t42
      t190 = t49*t42
      t193 = -0.1333333333333333D1*t187*t107+0.1333333333333333D1*t190
     #*t107
      t194 = t40*t193
      t195 = t194*t85
      t197 = t105*t55
      t200 = t54**2
      t201 = 1/t200
      t202 = t2*t201
      t204 = -t114-0.9268252376942842D0*t119
      t207 = -0.2067834969664667D0*t197*t107-0.6203504908994D0*t202*t204
      t208 = t207*t125
      t211 = t60**2
      t212 = 1/t211
      t213 = t212*t118
      t215 = 1.D0+0.223816694236D2*t212
      t216 = 1/t215
      t220 = t65*t55
      t223 = t66*t201
      t226 = -0.2625411059665811D0*t220*t119-1.D0*t223*t204
      t227 = 1/t66
      t228 = t226*t227
      t234 = 0.2505897912236993D-1*t208*t54+0.6519760588836059D-1*t213
     #*t107*t216+0.2247867095542612D-2*t228*t54-0.5011795824473985D-1
     #*t127-0.626408570946439D-1*t136-0.9690227711544374D-3*t147
      t237 = t40**2
      t238 = 1/t237
      t239 = t73*t238
      t242 = -0.101257459063883D3*t234*t74+0.101257459063883D3*t239*t183
      t243 = t242*t79
      t246 = 1/t82/rho
      t249 = t243*t83-4.D0*t80*t246
      t250 = t52*t249
      t253 = rho*(0.5011795824473985D-1*t127+0.626408570946439D-1*t136
     #+0.9690227711544374D-3*t147-0.1899772193293833D-1*t185
     #-0.1899772193293833D-1*t195-0.1899772193293833D-1*t250)
      vrhoa(i) = rho*(-t96-0.7599088773175333D-1*t100)+t11+t17+t23-t87
     #+t253
      t254 = -t93
      t255 = t40*t254
      t256 = t255*t85
      t257 = 0.1899772193293833D-1*t256
      vrhob(i) = rho*(-t257+0.7599088773175333D-1*t100)+t11+t17+t23
     #-t87+t253
      t261 = t45**2
      t262 = 1/t261
      t265 = t49**2
      t266 = 1/t265
      t269 = 0.4444444444444444D0*t262*t107+0.4444444444444444D0*t266
     #*t107
      t272 = 0.1899772193293833D-1*t40*t269*t85
      t273 = t94*t99
      t277 = t52*t77*t78*t83
      t278 = 0.22797266319526D0*t277
      t285 = 0.1899772193293833D-1*t183*t93*t85
      t286 = t184*t99
      t287 = 0.7599088773175333D-1*t286
      t290 = 1/t81/rho
      t300 = -0.4444444444444444D0*t262*t42*t290-0.1333333333333333D1
     #*t45*t107-0.4444444444444444D0*t266*t42*t290
     #+0.1333333333333333D1*t49*t107
      t303 = 0.1899772193293833D-1*t40*t300*t85
      t304 = t194*t99
      t305 = 0.7599088773175333D-1*t304
      t307 = 0.1899772193293833D-1*t94*t249
      t309 = t242*t97*t83
      t311 = t98*t246
      t313 = 4.D0*t309-16.D0*t311
      t317 = rho*(-t285-t287-t303-t305-t307-0.1899772193293833D-1*t52
     #*t313)
      t318 = t52*t309
      t320 = t52*t311
      t323 = rho*(-t285-t303-t307-t287-t305-0.7599088773175333D-1*t318
     #+0.3039635509270133D0*t320)
      t324 = 0.1002359164894797D0*t127
      t325 = 0.1252817141892878D0*t136
      t326 = 0.1938045542308875D-2*t147
      t327 = 0.3799544386587666D-1*t185
      t328 = 0.3799544386587666D-1*t195
      t329 = 0.3799544386587666D-1*t250
      t331 = 1/t104/t1
      t333 = t331*t7*t83
      t342 = 1/t110/t6
      t344 = t121**2
      t348 = 0.1378556646443111D0*t331*t83
      t350 = 0.4135669939329333D0*t105*t290
      t352 = 1/t117/t1
      t353 = t352*t83
      t355 = t118*t290
      t357 = -t348+t350-0.4077525916766971D0*t353+0.978606220024073D0
     #*t355
      t362 = (-0.1378556646443111D0*t333+0.4135669939329333D0*t105
     #*t111*t107*t121+0.4135669939329333D0*t106*t290
     #+0.12407009817988D1*t2*t342*t344-0.6203504908994D0*t112*t357)
     #*t125*t6
      t365 = 1/t2/t1
      t368 = t124*t365*t6*t107
      t370 = t126*t121
      t375 = t83*t134
      t376 = 1/t129/t13*t331*t375
      t379 = t130*t352*t375
      t382 = t131*t290*t134
      t384 = t129**2
      t388 = t133**2
      t391 = 1/t384/t13*t331*t83/t388
      t409 = (0.3446391616107778D-1*t333+0.5250822119331622D0*t19*t111
     #*t119*t121-0.2187842549721509D0*t138*t353+0.5250822119331622D0
     #*t138*t355+2.D0*t20*t342*t344-1.D0*t141*t357)*t145*t6
      t416 = t144/t20/t19*t6*t118*t107
      t418 = t146*t121
      t421 = t331*t26*t83
      t430 = 1/t152/t25
      t432 = t156**2
      t437 = -t348+t350-0.1237301536356754D0*t353+0.2969523687256209D0
     #*t355
      t453 = t83*t168
      t462 = t163**2
      t466 = t167**2
      t498 = 0.161199195401647D1*(-0.1378556646443111D0*t421
     #+0.4135669939329333D0*t105*t153*t107*t156+0.4135669939329333D0
     #*t149*t290+0.12407009817988D1*t2*t430*t432-0.6203504908994D0
     #*t154*t437)*t125*t25+0.5373306513388232D0*t159*t365*t25*t107
     #+0.161199195401647D1*t160*t156+0.3119768427869032D0/t163/t30
     #*t331*t453+0.495123804261555D0*t164*t352*t453
     #-0.1188297130227732D1*t165*t290*t168-0.1582929339795102D2/t462
     #/t30*t331*t83/t466+0.4140337942820627D-3*(0.3446391616107778D-1
     #*t421+0.5250822119331622D0*t35*t153*t119*t156
     #-0.2187842549721509D0*t172*t353+0.5250822119331622D0*t172*t355
     #+2.D0*t36*t430*t432-1.D0*t175*t437)*t179*t25
     #+0.1087008902583527D-3*t178/t36/t35*t25*t118*t107
     #+0.4140337942820627D-3*t180*t156
      t524 = t331*t55*t83
      t533 = 1/t200/t54
      t535 = t204**2
      t540 = -t348+t350-0.7723543647452368D0*t353+0.1853650475388568D1
     #*t355
      t556 = t83*t216
      t565 = t211**2
      t569 = t215**2
      s1 = 0.2505897912236993D-1*(-0.1378556646443111D0*t524
     #+0.4135669939329333D0*t105*t201*t107*t204+0.4135669939329333D0
     #*t197*t290+0.12407009817988D1*t2*t533*t535-0.6203504908994D0
     #*t202*t540)*t125*t54+0.8352993040789975D-2*t207*t365*t54*t107
     #+0.2505897912236993D-1*t208*t204+0.3423410311260694D-1/t211/t60
     #*t331*t556+0.5433133824030049D-1*t212*t352*t556
     #-0.1303952117767212D0*t213*t290*t216-0.7662163788798043D0/t565
     #/t60*t331*t83/t569+0.2247867095542612D-2*(0.3446391616107778D-1
     #*t524+0.5250822119331622D0*t65*t201*t119*t204
     #-0.2187842549721509D0*t220*t353+0.5250822119331622D0*t220*t355
     #+2.D0*t66*t533*t535-1.D0*t223*t540)*t227*t54
     #+0.5901575133296437D-3*t226/t66/t65*t54*t118*t107
     #+0.2247867095542612D-2*t228*t204
      t611 = s1-0.5011795824473985D-1*t362-0.1670598608157995D-1*t368
     #-0.5011795824473985D-1*t370-0.3289159980064473D-1*t376
     #-0.5220071424553658D-1*t379+0.1252817141892878D0*t382
     #+0.1244848083156773D1*t391-0.9690227711544374D-3*t409
     #-0.2544083100456872D-3*t416-0.9690227711544374D-3*t418
      t620 = t183**2
      t637 = 0.5011795824473985D-1*t362+0.1670598608157995D-1*t368
     #+0.5011795824473985D-1*t370+0.3289159980064473D-1*t376
     #+0.5220071424553658D-1*t379-0.1252817141892878D0*t382
     #-0.1244848083156773D1*t391+0.9690227711544374D-3*t409
     #+0.2544083100456872D-3*t416+0.9690227711544374D-3*t418
     #-0.1899772193293833D-1*t498*t51*t85-0.3799544386587666D-1*t183
     #*t193*t85-0.3799544386587666D-1*t184*t249-0.1899772193293833D-1
     #*t40*(0.4444444444444444D0*t262*t78*t83+0.2666666666666667D1
     #*t187*t290+0.4444444444444444D0*t266*t78*t83
     #-0.2666666666666667D1*t190*t290)*t85-0.3799544386587666D-1*t194
     #*t249-0.1899772193293833D-1*t52*((-0.101257459063883D3*t611*t74
     #+0.2025149181277661D3*t234*t238*t183-0.2025149181277661D3*t73
     #/t237/t40*t620+0.101257459063883D3*t239*t498)*t79*t83-8.D0*t243
     #*t246+20.D0*t80/t82/t81)
      t638 = rho*t637
      v2rhoa2(i) = rho*(-t272-0.1519817754635067D0*t273-t278)
     #-0.3799544386587666D-1*t95-0.1519817754635067D0*t100+t317+t323
     #+t324+t325+t326-t327-t328-t329+t638
      t639 = t255*t99
      t647 = 0.1899772193293833D-1*t183*t254*t85
      t648 = 0.7599088773175333D-1*t286
      t652 = -0.1899772193293833D-1*t40*t300*t85
      t653 = 0.7599088773175333D-1*t304
      t655 = 0.1899772193293833D-1*t255*t249
      t660 = rho*(-t647+t648-t652+t653-t655+0.1899772193293833D-1*t52
     #*t313)
      t664 = rho*(-t647-t652-t655+t648+t653+0.7599088773175333D-1*t318
     #-0.3039635509270133D0*t320)
      v2rhob2(i) = rho*(-t272+0.1519817754635067D0*t639-t278)
     #-0.3799544386587666D-1*t256+0.1519817754635067D0*t100+t660+t664
     #+t324+t325+t326-t327-t328-t329+t638
      v2rhoab(i) = rho*(0.1899772193293833D-1*t40*t269*t85
     #+0.7599088773175333D-1*t273-0.7599088773175333D-1*t639
     #+0.22797266319526D0*t277)-t257+0.5D0*t660-t96+0.5D0*t323+t324
     #+t325+t326-t327-t328-t329+t638+0.5D0*t317+0.5D0*t664
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vrhob(i) = 0.0d0
      v2rhoa2(i) = 0.0d0
      v2rhob2(i) = 0.0d0
      v2rhoab(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end
      
      
      subroutine rks_c_vwn5
     & (ideriv,npt,rhoa1,sigmaaa1,
     &  zk,vrhoa,vsigmaaa,
     &  v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
c
c     S.H. Vosko, L. Wilk, and M. Nusair
c     Accurate spin-dependent electron liquid correlation energies for
c     local spin density calculations: a critical analysis
c     Can. J. Phys. 58 (1980) 1200-1211
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit real*8 (a-h,o-z)
      integer ideriv,npt
      real*8 rhoa1(npt)
      real*8 sigmaaa1(npt)
      real*8 zk(npt),vrhoa(npt),vsigmaaa(npt)
      real*8 v2rhoa2(npt),v2rhoasigmaaa(npt),v2sigmaaa2(npt)
      parameter(tol=1.0d-20)
      
      if(ideriv.eq.0) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      t1 = 1/rho
      t2 = t1**(1.D0/3.D0)
      t4 = t1**(1.D0/6.D0)
      t7 = 1/(0.6203504908994D0*t2+0.2935818660072219D1*t4+0.129352D2)
      t10 = dlog(0.6203504908994D0*t2*t7)
      t16 = datan(0.615199081975908D1/(0.1575246635799487D1*t4
     #+0.372744D1))
      t20 = (0.7876233178997433D0*t4+0.10498D0)**2
      t22 = dlog(t20*t7)
      zk(i) = rho*(0.310907D-1*t10+0.3878329487811301D-1*t16
     #+0.9690227711544374D-3*t22)
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      t1 = 1/rho
      t2 = t1**(1.D0/3.D0)
      t4 = t1**(1.D0/6.D0)
      t6 = 0.6203504908994D0*t2+0.2935818660072219D1*t4+0.129352D2
      t7 = 1/t6
      t10 = dlog(0.6203504908994D0*t2*t7)
      t11 = 0.310907D-1*t10
      t13 = 0.1575246635799487D1*t4+0.372744D1
      t16 = datan(0.615199081975908D1/t13)
      t17 = 0.3878329487811301D-1*t16
      t19 = 0.7876233178997433D0*t4+0.10498D0
      t20 = t19**2
      t22 = dlog(t20*t7)
      t23 = 0.9690227711544374D-3*t22
      zk(i) = rho*(t11+t17+t23)
      t25 = t2**2
      t26 = 1/t25
      t28 = rho**2
      t29 = 1/t28
      t32 = t6**2
      t33 = 1/t32
      t37 = t4**2
      t38 = t37**2
      t40 = 1/t38/t4
      t41 = t40*t29
      t43 = -0.2067834969664667D0*t26*t29-0.4893031100120365D0*t41
      t51 = t13**2
      t52 = 1/t51
      vrhoa(i) = t11+t17+t23+rho*(0.5011795824473985D-1*(
     #-0.2067834969664667D0*t26*t7*t29-0.6203504908994D0*t2*t33*t43)
     #/t2*t6+0.626408570946439D-1*t52*t40*t29/(1.D0+0.378469910464D2
     #*t52)+0.9690227711544374D-3*(-0.2625411059665811D0*t19*t7*t41
     #-1.D0*t20*t33*t43)/t20*t6)
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      t1 = 1/rho
      t2 = t1**(1.D0/3.D0)
      t4 = t1**(1.D0/6.D0)
      t6 = 0.6203504908994D0*t2+0.2935818660072219D1*t4+0.129352D2
      t7 = 1/t6
      t10 = dlog(0.6203504908994D0*t2*t7)
      t11 = 0.310907D-1*t10
      t13 = 0.1575246635799487D1*t4+0.372744D1
      t16 = datan(0.615199081975908D1/t13)
      t17 = 0.3878329487811301D-1*t16
      t19 = 0.7876233178997433D0*t4+0.10498D0
      t20 = t19**2
      t22 = dlog(t20*t7)
      t23 = 0.9690227711544374D-3*t22
      zk(i) = rho*(t11+t17+t23)
      t25 = t2**2
      t26 = 1/t25
      t27 = t26*t7
      t28 = rho**2
      t29 = 1/t28
      t32 = t6**2
      t33 = 1/t32
      t34 = t2*t33
      t37 = t4**2
      t38 = t37**2
      t39 = t38*t4
      t40 = 1/t39
      t41 = t40*t29
      t43 = -0.2067834969664667D0*t26*t29-0.4893031100120365D0*t41
      t46 = -0.2067834969664667D0*t27*t29-0.6203504908994D0*t34*t43
      t47 = 1/t2
      t48 = t46*t47
      t49 = t48*t6
      t51 = t13**2
      t52 = 1/t51
      t53 = t52*t40
      t55 = 1.D0+0.378469910464D2*t52
      t56 = 1/t55
      t58 = t53*t29*t56
      t60 = t19*t7
      t63 = t20*t33
      t66 = -0.2625411059665811D0*t60*t41-1.D0*t63*t43
      t67 = 1/t20
      t68 = t66*t67
      t69 = t68*t6
      vrhoa(i) = t11+t17+t23+rho*(0.5011795824473985D-1*t49
     #+0.626408570946439D-1*t58+0.9690227711544374D-3*t69)
      t77 = 1/t25/t1
      t79 = t28**2
      t80 = 1/t79
      t81 = t77*t7*t80
      t88 = 1/t28/rho
      t92 = 1/t32/t6
      t94 = t43**2
      t102 = 1/t39/t1
      t103 = t102*t80
      t105 = t40*t88
      t107 = -0.1378556646443111D0*t77*t80+0.4135669939329333D0*t26
     #*t88-0.4077525916766971D0*t103+0.978606220024073D0*t105
      t125 = t80*t56
      t134 = t51**2
      t138 = t55**2
      s1 = 0.2004718329789594D0*t49+0.2505634283785756D0*t58
      v2rhoa2(i) = s1+0.387609108461775D-2*t69+2.D0*rho*
     #(0.5011795824473985D-1*(-0.1378556646443111D0*t81
     #+0.4135669939329333D0*t26*t33*t29*t43+0.4135669939329333D0*t27
     #*t88+0.12407009817988D1*t2*t92*t94-0.6203504908994D0*t34*t107)
     #*t47*t6+0.1670598608157995D-1*t46/t2/t1*t6*t29
     #+0.5011795824473985D-1*t48*t43+0.3289159980064473D-1/t51/t13*t77
     #*t125+0.5220071424553658D-1*t52*t102*t125-0.1252817141892878D0
     #*t53*t88*t56-0.1244848083156773D1/t134/t13*t77*t80/t138
     #+0.9690227711544374D-3*(0.3446391616107778D-1*t81
     #+0.5250822119331622D0*t19*t33*t41*t43-0.2187842549721509D0*t60
     #*t103+0.5250822119331622D0*t60*t105+2.D0*t20*t92*t94-1.D0*t63
     #*t107)*t67*t6+0.2544083100456872D-3*t66/t20/t19*t6*t40*t29
     #+0.9690227711544374D-3*t68*t43)
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      v2rhoa2(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end


      subroutine uks_c_p86
     & (ideriv,npt,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,
     &  zk,vrhoa,vrhob,vsigmaaa,vsigmabb,vsigmaab,
     &  v2rhoa2,v2rhob2,v2rhoab,
     &  v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb,
     &  v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,
     &  v2sigmaaa2,v2sigmaaaab,v2sigmaaabb,
     &  v2sigmaab2,v2sigmaabbb,v2sigmabb2)
c
c     J.P. Perdew
c     Density-functional approximation for the correlation energy of
c     the inhomogeneous electron gas
c     Phys. Rev. B33 (1986) 8822-8824
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit real*8 (a-h,o-z)
      integer ideriv,npt
      real*8 rhoa1(npt),rhob1(npt)
      real*8 sigmaaa1(npt),sigmabb1(npt),sigmaab1(npt)
      real*8 zk(npt),vrhoa(npt),vrhob(npt)
      real*8 vsigmaaa(npt),vsigmabb(npt),vsigmaab(npt)
      real*8 v2rhoa2(npt),v2rhob2(npt),v2rhoab(npt)
      real*8 v2rhoasigmaaa(npt),v2rhoasigmaab(npt)
      real*8 v2rhoasigmabb(npt),v2rhobsigmabb(npt)
      real*8 v2rhobsigmaab(npt),v2rhobsigmaaa(npt)
      real*8 v2sigmaaa2(npt),v2sigmaaaab(npt),v2sigmaaabb(npt)
      real*8 v2sigmaab2(npt),v2sigmaabbb(npt),v2sigmabb2(npt)
      parameter(tol=1.0d-20)
      logical t7
      
      if (ideriv.eq.0) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmabb
      t2 = 1/rhob
      t3 = t2**(1.D0/3.D0)
      t4 = 0.6203504908994D0*t3
      t6 = t2**(1.D0/6.D0)
      t12 = dlog(t4)
      t18 = piecewise(1.D0 .le. t4,-0.843D-1/(1.D0
     #+0.1101176160755631D1*t6+0.1619735131738333D0*t3),0.1555D-1*t12
     #-0.269D-1+0.43424534362958D-3*t3*t12-0.297768235631712D-2*t3)
      t20 = dsqrt(sigmabb)
      t22 = t3**2
      t31 = 0.1667D-2+(0.2568D-2+0.1443307452126544D-1*t3
     #+0.2843543831490386D-5*t22)/(1.D0+0.5411317332115466D1*t3
     #+0.1816419932959077D0*t22+0.1763993811759022D-1*t2)
      t34 = rhob**(1.D0/6.D0)
      t39 = dexp(-0.81290825D-3*t20/t31/t34/rhob)
      t41 = rhob**(1.D0/3.D0)
      zk(i) = rhob*t18+0.7937005259840997D0*t39*t31*sigmabb/t41/rhob
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = 1/rhoa
      t3 = t2**(1.D0/3.D0)
      t4 = 0.6203504908994D0*t3
      t6 = t2**(1.D0/6.D0)
      t12 = dlog(t4)
      t18 = piecewise(1.D0 .le. t4,-0.843D-1/(1.D0
     #+0.1101176160755631D1*t6+0.1619735131738333D0*t3),0.1555D-1*t12
     #-0.269D-1+0.43424534362958D-3*t3*t12-0.297768235631712D-2*t3)
      t20 = dsqrt(sigmaaa)
      t22 = t3**2
      t31 = 0.1667D-2+(0.2568D-2+0.1443307452126544D-1*t3
     #+0.2843543831490386D-5*t22)/(1.D0+0.5411317332115466D1*t3
     #+0.1816419932959077D0*t22+0.1763993811759022D-1*t2)
      t34 = rhoa**(1.D0/6.D0)
      t39 = dexp(-0.81290825D-3*t20/t31/t34/rhoa)
      t41 = rhoa**(1.D0/3.D0)
      zk(i) = rhoa*t18+0.7937005259840997D0*t39*t31*sigmaaa/t41/rhoa
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = 1/rho
      t5 = t4**(1.D0/3.D0)
      t6 = 0.6203504908994D0*t5
      t8 = t4**(1.D0/6.D0)
      t13 = 0.1423D0/(1.D0+0.8292885914166397D0*t8+0.20682485366586D0
     #*t5)
      t22 = (rhoa-1.D0*rhob)*t4
      t23 = 1.D0+t22
      t24 = t23**(1.D0/3.D0)
      t27 = 1.D0-1.D0*t22
      t28 = t27**(1.D0/3.D0)
      t30 = t24*t23+t28*t27-2.D0
      t34 = dlog(t6)
      t36 = t5*t34
      t46 = piecewise(1.D0 .le. t6,-t13+0.1923661050931536D1*(
     #-0.843D-1/(1.D0+0.1101176160755631D1*t8+0.1619735131738333D0*t5)
     #+t13)*t30,0.311D-1*t34-0.48D-1+0.12407009817988D-2*t36
     #-0.719606569443304D-2*t5+0.1923661050931536D1*(-0.1555D-1*t34
     #+0.211D-1-0.80645563816922D-3*t36+0.421838333811592D-2*t5)*t30)
      t48 = dsqrt(sigma)
      t50 = t5**2
      t59 = 0.1667D-2+(0.2568D-2+0.1443307452126544D-1*t5
     #+0.2843543831490386D-5*t50)/(1.D0+0.5411317332115466D1*t5
     #+0.1816419932959077D0*t50+0.1763993811759022D-1*t4)
      t62 = rho**(1.D0/6.D0)
      t67 = dexp(-0.81290825D-3*t48/t59/t62/rho)
      t69 = rho**(1.D0/3.D0)
      t74 = 0.5D0+0.5D0*t22
      t75 = t74**(1.D0/3.D0)
      t76 = t75**2
      t79 = 0.5D0-0.5D0*t22
      t80 = t79**(1.D0/3.D0)
      t81 = t80**2
      t84 = dsqrt(t76*t74+t81*t79)
      zk(i) = rho*t46+0.7937005259840997D0*t67*t59*sigma/t69/rho/t84
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmabb
      t2 = 1/rhob
      t3 = t2**(1.D0/3.D0)
      t4 = 0.6203504908994D0*t3
      t7 = 1.D0 .le. t4
      t6 = t2**(1.D0/6.D0)
      t9 = 1.D0+0.1101176160755631D1*t6+0.1619735131738333D0*t3
      t12 = dlog(t4)
      t18 = piecewise(t7,-0.843D-1/t9,0.1555D-1*t12-0.269D-1
     #+0.43424534362958D-3*t3*t12-0.297768235631712D-2*t3)
      t20 = dsqrt(sigmabb)
      t22 = t3**2
      t24 = 0.2568D-2+0.1443307452126544D-1*t3+0.2843543831490386D-5*t22
      t28 = 1.D0+0.5411317332115466D1*t3+0.1816419932959077D0*t22
     #+0.1763993811759022D-1*t2
      t29 = 1/t28
      t31 = 0.1667D-2+t24*t29
      t33 = t20/t31
      t34 = rhob**(1.D0/6.D0)
      t36 = 1/t34/rhob
      t39 = dexp(-0.81290825D-3*t33*t36)
      t40 = t39*t31
      t41 = rhob**(1.D0/3.D0)
      t43 = 1/t41/rhob
      t44 = sigmabb*t43
      zk(i) = rhob*t18+0.7937005259840997D0*t40*t44
      vrhoa(i) = 0.D0
      t47 = t9**2
      t49 = t6**2
      t50 = t49**2
      t53 = rhob**2
      t54 = 1/t53
      t57 = 1/t22
      t58 = t57*t54
      t71 = piecewise(t7,0.843D-1/t47*(-0.1835293601259385D0/t50/t6
     #*t54-0.5399117105794445D-1*t58),-0.5183333333333333D-2*t2
     #-0.1447484478765267D-3*t57*t12*t54-0.1447484478765267D-3*t3*t2
     #+0.99256078543904D-3*t58)
      t73 = t31**2
      t78 = 1/t3*t54
      t82 = t28**2
      t91 = (-0.4811024840421814D-2*t58-0.1895695887660258D-5*t78)*t29
     #-1.D0*t24/t82*(-0.1803772444038489D1*t58-0.1210946621972718D0
     #*t78-0.1763993811759022D-1*t54)
      vrhob(i) = t18+rhob*t71+0.7937005259840997D0*(0.81290825D-3*t20
     #/t73*t36*t91+0.9483929583333333D-3*t33/t34/t53)*t39*t31*sigmabb
     #*t43+0.7937005259840997D0*t39*t91*t44-0.10582673679788D1*t40
     #*sigmabb/t41/t53
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      t113 = dsqrt(rhob)
      vsigmabb(i) = -0.322602852800907D-3*t20/t113/t53*t39
     #+0.7937005259840997D0*t40*t43
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = 1/rhoa
      t3 = t2**(1.D0/3.D0)
      t4 = 0.6203504908994D0*t3
      t7 = 1.D0 .le. t4
      t6 = t2**(1.D0/6.D0)
      t9 = 1.D0+0.1101176160755631D1*t6+0.1619735131738333D0*t3
      t12 = dlog(t4)
      t18 = piecewise(t7,-0.843D-1/t9,0.1555D-1*t12-0.269D-1
     #+0.43424534362958D-3*t3*t12-0.297768235631712D-2*t3)
      t20 = dsqrt(sigmaaa)
      t22 = t3**2
      t24 = 0.2568D-2+0.1443307452126544D-1*t3+0.2843543831490386D-5*t22
      t28 = 1.D0+0.5411317332115466D1*t3+0.1816419932959077D0*t22
     #+0.1763993811759022D-1*t2
      t29 = 1/t28
      t31 = 0.1667D-2+t24*t29
      t33 = t20/t31
      t34 = rhoa**(1.D0/6.D0)
      t36 = 1/t34/rhoa
      t39 = dexp(-0.81290825D-3*t33*t36)
      t40 = t39*t31
      t41 = rhoa**(1.D0/3.D0)
      t43 = 1/t41/rhoa
      t44 = sigmaaa*t43
      zk(i) = rhoa*t18+0.7937005259840997D0*t40*t44
      t47 = t9**2
      t49 = t6**2
      t50 = t49**2
      t53 = rhoa**2
      t54 = 1/t53
      t57 = 1/t22
      t58 = t57*t54
      t71 = piecewise(t7,0.843D-1/t47*(-0.1835293601259385D0/t50/t6
     #*t54-0.5399117105794445D-1*t58),-0.5183333333333333D-2*t2
     #-0.1447484478765267D-3*t57*t12*t54-0.1447484478765267D-3*t3*t2
     #+0.99256078543904D-3*t58)
      t73 = t31**2
      t78 = 1/t3*t54
      t82 = t28**2
      t91 = (-0.4811024840421814D-2*t58-0.1895695887660258D-5*t78)*t29
     #-1.D0*t24/t82*(-0.1803772444038489D1*t58-0.1210946621972718D0
     #*t78-0.1763993811759022D-1*t54)
      vrhoa(i) = t18+rhoa*t71+0.7937005259840997D0*(0.81290825D-3*t20
     #/t73*t36*t91+0.9483929583333333D-3*t33/t34/t53)*t39*t31*sigmaaa
     #*t43+0.7937005259840997D0*t39*t91*t44-0.10582673679788D1*t40
     #*sigmaaa/t41/t53
      vrhob(i) = 0.D0
      t113 = dsqrt(rhoa)
      vsigmaaa(i) = -0.322602852800907D-3*t20/t113/t53*t39
     #+0.7937005259840997D0*t40*t43
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = 1/rho
      t5 = t4**(1.D0/3.D0)
      t6 = 0.6203504908994D0*t5
      t7 = 1.D0 .le. t6
      t8 = t4**(1.D0/6.D0)
      t11 = 1.D0+0.8292885914166397D0*t8+0.20682485366586D0*t5
      t13 = 0.1423D0/t11
      t16 = 1.D0+0.1101176160755631D1*t8+0.1619735131738333D0*t5
      t19 = -0.843D-1/t16+t13
      t21 = rhoa-1.D0*rhob
      t22 = t21*t4
      t23 = 1.D0+t22
      t24 = t23**(1.D0/3.D0)
      t27 = 1.D0-1.D0*t22
      t28 = t27**(1.D0/3.D0)
      t30 = t24*t23+t28*t27-2.D0
      t34 = dlog(t6)
      t36 = t5*t34
      t42 = -0.1555D-1*t34+0.211D-1-0.80645563816922D-3*t36
     #+0.421838333811592D-2*t5
      t46 = piecewise(t7,-t13+0.1923661050931536D1*t19*t30,0.311D-1
     #*t34-0.48D-1+0.12407009817988D-2*t36-0.719606569443304D-2*t5
     #+0.1923661050931536D1*t42*t30)
      t48 = dsqrt(sigma)
      t50 = t5**2
      t52 = 0.2568D-2+0.1443307452126544D-1*t5+0.2843543831490386D-5*t50
      t56 = 1.D0+0.5411317332115466D1*t5+0.1816419932959077D0*t50
     #+0.1763993811759022D-1*t4
      t57 = 1/t56
      t59 = 0.1667D-2+t52*t57
      t61 = t48/t59
      t62 = rho**(1.D0/6.D0)
      t64 = 1/t62/rho
      t67 = dexp(-0.81290825D-3*t61*t64)
      t68 = t67*t59
      t69 = rho**(1.D0/3.D0)
      t71 = 1/t69/rho
      t74 = 0.5D0+0.5D0*t22
      t75 = t74**(1.D0/3.D0)
      t76 = t75**2
      t79 = 0.5D0-0.5D0*t22
      t80 = t79**(1.D0/3.D0)
      t81 = t80**2
      t83 = t76*t74+t81*t79
      t84 = dsqrt(t83)
      t85 = 1/t84
      t86 = sigma*t71*t85
      zk(i) = rho*t46+0.7937005259840997D0*t68*t86
      t93 = 0.1333333333333333D1*t24*t4-0.1333333333333333D1*t28*t4
      t98 = piecewise(t7,0.1923661050931536D1*t19
     #*t93,0.1923661050931536D1*t42*t93)
      t100 = t68*sigma
      t103 = t71/t84/t83
      t108 = 0.8333333333333333D0*t76*t4-0.8333333333333333D0*t81*t4
      t112 = t11**2
      t114 = t8**2
      t115 = t114**2
      t118 = rho**2
      t119 = 1/t118
      t120 = 1/t115/t8*t119
      t122 = 1/t50
      t123 = t122*t119
      t127 = 0.1423D0/t112*(-0.1382147652361066D0*t120
     #-0.6894161788861999D-1*t123)
      t128 = t16**2
      t144 = -0.1333333333333333D1*t24*t21*t119+0.1333333333333333D1
     #*t28*t21*t119
      t150 = t122*t34*t119
      t152 = t5*t4
      t165 = piecewise(t7,t127+0.1923661050931536D1*(0.843D-1/t128*(
     #-0.1835293601259385D0*t120-0.5399117105794445D-1*t123)-t127)*t30
     #+0.1923661050931536D1*t19*t144,-0.1036666666666667D-1*t4
     #-0.4135669939329333D-3*t150-0.4135669939329333D-3*t152
     #+0.2398688564811013D-2*t123+0.1923661050931536D1*
     #(0.5183333333333333D-2*t4+0.2688185460564067D-3*t150
     #+0.2688185460564067D-3*t152-0.1406127779371973D-2*t123)*t30
     #+0.1923661050931536D1*t42*t144)
      t166 = rho*t165
      t167 = t59**2
      t172 = 1/t5*t119
      t176 = t56**2
      t185 = (-0.4811024840421814D-2*t123-0.1895695887660258D-5*t172)
     #*t57-1.D0*t52/t176*(-0.1803772444038489D1*t123
     #-0.1210946621972718D0*t172-0.1763993811759022D-1*t119)
      t197 = 0.7937005259840997D0*(0.81290825D-3*t48/t167*t64*t185
     #+0.9483929583333333D-3*t61/t62/t118)*t67*t59*t86
      t200 = 0.7937005259840997D0*t67*t185*t86
      t206 = 0.10582673679788D1*t68*sigma/t69/t118*t85
      t216 = 0.3968502629920499D0*t100*t103*(-0.8333333333333333D0*t76
     #*t21*t119+0.8333333333333333D0*t81*t21*t119)
      vrhoa(i) = rho*t98-0.3968502629920499D0*t100*t103*t108+t46+t166
     #+t197+t200-t206-t216
      t217 = -t93
      t222 = piecewise(t7,0.1923661050931536D1*t19
     #*t217,0.1923661050931536D1*t42*t217)
      vrhob(i) = rho*t222+0.3968502629920499D0*t100*t103*t108+t46+t166
     #+t197+t200-t206-t216
      t228 = dsqrt(rho)
      t233 = t48/t228/t118*t67*t85
      t236 = t68*t71*t85
      vsigmaaa(i) = -0.322602852800907D-3*t233+0.7937005259840997D0*t236
      vsigmaab(i) = -0.645205705601814D-3*t233+0.1587401051968199D1*t236
      vsigmabb(i) = vsigmaaa(i)
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vrhob(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      vsigmaab(i) = 0.0d0
      vsigmabb(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmabb
      t2 = 1/rhob
      t3 = t2**(1.D0/3.D0)
      t4 = 0.6203504908994D0*t3
      t7 = 1.D0 .le. t4
      t6 = t2**(1.D0/6.D0)
      t9 = 1.D0+0.1101176160755631D1*t6+0.1619735131738333D0*t3
      t12 = dlog(t4)
      t18 = piecewise(t7,-0.843D-1/t9,0.1555D-1*t12-0.269D-1
     #+0.43424534362958D-3*t3*t12-0.297768235631712D-2*t3)
      t20 = dsqrt(sigmabb)
      t22 = t3**2
      t24 = 0.2568D-2+0.1443307452126544D-1*t3+0.2843543831490386D-5*t22
      t28 = 1.D0+0.5411317332115466D1*t3+0.1816419932959077D0*t22
     #+0.1763993811759022D-1*t2
      t29 = 1/t28
      t31 = 0.1667D-2+t24*t29
      t32 = 1/t31
      t33 = t20*t32
      t34 = rhob**(1.D0/6.D0)
      t36 = 1/t34/rhob
      t39 = dexp(-0.81290825D-3*t33*t36)
      t40 = t39*t31
      t41 = rhob**(1.D0/3.D0)
      t43 = 1/t41/rhob
      t44 = sigmabb*t43
      zk(i) = rhob*t18+0.7937005259840997D0*t40*t44
      vrhoa(i) = 0.D0
      t47 = t9**2
      t48 = 1/t47
      t49 = t6**2
      t50 = t49**2
      t51 = t50*t6
      t52 = 1/t51
      t53 = rhob**2
      t54 = 1/t53
      t57 = 1/t22
      t58 = t57*t54
      t60 = -0.1835293601259385D0*t52*t54-0.5399117105794445D-1*t58
      t64 = t57*t12
      t67 = t3*t2
      t71 = piecewise(t7,0.843D-1*t48*t60,-0.5183333333333333D-2*t2
     #-0.1447484478765267D-3*t64*t54-0.1447484478765267D-3*t67
     #+0.99256078543904D-3*t58)
      t73 = t31**2
      t75 = t20/t73
      t77 = 1/t3
      t78 = t77*t54
      t80 = -0.4811024840421814D-2*t58-0.1895695887660258D-5*t78
      t82 = t28**2
      t83 = 1/t82
      t84 = t24*t83
      t88 = -0.1803772444038489D1*t58-0.1210946621972718D0*t78
     #-0.1763993811759022D-1*t54
      t91 = t80*t29-1.D0*t84*t88
      t96 = 1/t34/t53
      t99 = 0.81290825D-3*t75*t36*t91+0.9483929583333333D-3*t33*t96
      t100 = t99*t39
      t101 = t31*sigmabb
      t102 = t101*t43
      t105 = t39*t91
      t109 = 1/t41/t53
      t110 = sigmabb*t109
      vrhob(i) = t18+rhob*t71+0.7937005259840997D0*t100*t102
     #+0.7937005259840997D0*t105*t44-0.10582673679788D1*t40*t110
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      t113 = dsqrt(rhob)
      t115 = 1/t113/t53
      vsigmabb(i) = -0.322602852800907D-3*t20*t115*t39
     #+0.7937005259840997D0*t40*t43
      v2rhoa2(i) = 0.D0
      v2rhoab(i) = 0.D0
      t124 = t60**2
      t129 = t53**2
      t130 = 1/t129
      t133 = t53*rhob
      t134 = 1/t133
      t138 = 1/t22/t2
      t139 = t138*t130
      t141 = t57*t134
      t158 = piecewise(t7,-0.1686D0/t47/t9*t124+0.843D-1*t48*(
     #-0.1529411334382821D0/t51/t2*t130+0.367058720251877D0*t52*t134
     #-0.3599411403862963D-1*t139+0.1079823421158889D0*t141
     #),0.5183333333333333D-2*t54-0.9649896525101778D-4*t138*t12*t130
     #-0.1888622605627062D-2*t141+0.2894968957530533D-3*t64*t134
     #+0.1447484478765267D-3*t3*t54+0.6617071902926934D-3*t139)
      t163 = t91**2
      t173 = 1/t67*t130
      t175 = t77*t134
      t185 = t88**2
      t196 = (-0.3207349893614542D-2*t139+0.9622049680843627D-2*t141
     #-0.6318986292200858D-6*t173+0.3791391775320515D-5*t175)*t29-2.D0
     #*t80*t83*t88+2.D0*t24/t82/t28*t185-1.D0*t84*(
     #-0.1202514962692326D1*t139+0.3607544888076978D1*t141
     #-0.4036488739909061D-1*t173+0.2421893243945437D0*t175
     #+0.3527987623518044D-1*t134)
      t208 = t99**2
      v2rhob2(i) = 2.D0*t71+rhob*t158+0.7937005259840997D0*(
     #-0.16258165D-2*t20/t73/t31*t36*t163-0.1896785916666667D-2*t75
     #*t96*t91+0.81290825D-3*t75*t36*t196-0.2054851409722222D-2*t33
     #/t34/t133)*t39*t102+0.7937005259840997D0*t208*t39*t102
     #+0.1587401051968199D1*t100*t91*sigmabb*t43-0.2116534735957599D1
     #*t100*t101*t109+0.7937005259840997D0*t39*t196*t44
     #-0.2116534735957599D1*t105*t110+0.2469290525283866D1*t40*sigmabb
     #/t41/t133
      v2sigmaaa2(i) = 0.D0
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      t233 = t41**2
      v2sigmabb2(i) = -0.4839042792013605D-3/t20*t115*t39
     #+0.1311232602576965D-6/t233/t133*t32*t39
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = 1/rhoa
      t3 = t2**(1.D0/3.D0)
      t4 = 0.6203504908994D0*t3
      t7 = 1.D0 .le. t4
      t6 = t2**(1.D0/6.D0)
      t9 = 1.D0+0.1101176160755631D1*t6+0.1619735131738333D0*t3
      t12 = dlog(t4)
      t18 = piecewise(t7,-0.843D-1/t9,0.1555D-1*t12-0.269D-1
     #+0.43424534362958D-3*t3*t12-0.297768235631712D-2*t3)
      t20 = dsqrt(sigmaaa)
      t22 = t3**2
      t24 = 0.2568D-2+0.1443307452126544D-1*t3+0.2843543831490386D-5*t22
      t28 = 1.D0+0.5411317332115466D1*t3+0.1816419932959077D0*t22
     #+0.1763993811759022D-1*t2
      t29 = 1/t28
      t31 = 0.1667D-2+t24*t29
      t32 = 1/t31
      t33 = t20*t32
      t34 = rhoa**(1.D0/6.D0)
      t36 = 1/t34/rhoa
      t39 = dexp(-0.81290825D-3*t33*t36)
      t40 = t39*t31
      t41 = rhoa**(1.D0/3.D0)
      t43 = 1/t41/rhoa
      t44 = sigmaaa*t43
      zk(i) = rhoa*t18+0.7937005259840997D0*t40*t44
      t47 = t9**2
      t48 = 1/t47
      t49 = t6**2
      t50 = t49**2
      t51 = t50*t6
      t52 = 1/t51
      t53 = rhoa**2
      t54 = 1/t53
      t57 = 1/t22
      t58 = t57*t54
      t60 = -0.1835293601259385D0*t52*t54-0.5399117105794445D-1*t58
      t64 = t57*t12
      t67 = t3*t2
      t71 = piecewise(t7,0.843D-1*t48*t60,-0.5183333333333333D-2*t2
     #-0.1447484478765267D-3*t64*t54-0.1447484478765267D-3*t67
     #+0.99256078543904D-3*t58)
      t73 = t31**2
      t75 = t20/t73
      t77 = 1/t3
      t78 = t77*t54
      t80 = -0.4811024840421814D-2*t58-0.1895695887660258D-5*t78
      t82 = t28**2
      t83 = 1/t82
      t84 = t24*t83
      t88 = -0.1803772444038489D1*t58-0.1210946621972718D0*t78
     #-0.1763993811759022D-1*t54
      t91 = t80*t29-1.D0*t84*t88
      t96 = 1/t34/t53
      t99 = 0.81290825D-3*t75*t36*t91+0.9483929583333333D-3*t33*t96
      t100 = t99*t39
      t101 = t31*sigmaaa
      t102 = t101*t43
      t105 = t39*t91
      t109 = 1/t41/t53
      t110 = sigmaaa*t109
      vrhoa(i) = t18+rhoa*t71+0.7937005259840997D0*t100*t102
     #+0.7937005259840997D0*t105*t44-0.10582673679788D1*t40*t110
      vrhob(i) = 0.D0
      t113 = dsqrt(rhoa)
      t115 = 1/t113/t53
      vsigmaaa(i) = -0.322602852800907D-3*t20*t115*t39
     #+0.7937005259840997D0*t40*t43
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      t124 = t60**2
      t129 = t53**2
      t130 = 1/t129
      t133 = t53*rhoa
      t134 = 1/t133
      t138 = 1/t22/t2
      t139 = t138*t130
      t141 = t57*t134
      t158 = piecewise(t7,-0.1686D0/t47/t9*t124+0.843D-1*t48*(
     #-0.1529411334382821D0/t51/t2*t130+0.367058720251877D0*t52*t134
     #-0.3599411403862963D-1*t139+0.1079823421158889D0*t141
     #),0.5183333333333333D-2*t54-0.9649896525101778D-4*t138*t12*t130
     #-0.1888622605627062D-2*t141+0.2894968957530533D-3*t64*t134
     #+0.1447484478765267D-3*t3*t54+0.6617071902926934D-3*t139)
      t163 = t91**2
      t173 = 1/t67*t130
      t175 = t77*t134
      t185 = t88**2
      t196 = (-0.3207349893614542D-2*t139+0.9622049680843627D-2*t141
     #-0.6318986292200858D-6*t173+0.3791391775320515D-5*t175)*t29-2.D0
     #*t80*t83*t88+2.D0*t24/t82/t28*t185-1.D0*t84*(
     #-0.1202514962692326D1*t139+0.3607544888076978D1*t141
     #-0.4036488739909061D-1*t173+0.2421893243945437D0*t175
     #+0.3527987623518044D-1*t134)
      t208 = t99**2
      v2rhoa2(i) = 2.D0*t71+rhoa*t158+0.7937005259840997D0*(
     #-0.16258165D-2*t20/t73/t31*t36*t163-0.1896785916666667D-2*t75
     #*t96*t91+0.81290825D-3*t75*t36*t196-0.2054851409722222D-2*t33
     #/t34/t133)*t39*t102+0.7937005259840997D0*t208*t39*t102
     #+0.1587401051968199D1*t100*t91*sigmaaa*t43-0.2116534735957599D1
     #*t100*t101*t109+0.7937005259840997D0*t39*t196*t44
     #-0.2116534735957599D1*t105*t110+0.2469290525283866D1*t40*sigmaaa
     #/t41/t133
      v2rhob2(i) = 0.D0
      v2rhoab(i) = 0.D0
      t233 = t41**2
      v2sigmaaa2(i) = -0.4839042792013605D-3/t20*t115*t39
     #+0.1311232602576965D-6/t233/t133*t32*t39
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      v2sigmabb2(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = 1/rho
      t5 = t4**(1.D0/3.D0)
      t6 = 0.6203504908994D0*t5
      t7 = 1.D0 .le. t6
      t8 = t4**(1.D0/6.D0)
      t11 = 1.D0+0.8292885914166397D0*t8+0.20682485366586D0*t5
      t13 = 0.1423D0/t11
      t16 = 1.D0+0.1101176160755631D1*t8+0.1619735131738333D0*t5
      t19 = -0.843D-1/t16+t13
      t21 = rhoa-1.D0*rhob
      t22 = t21*t4
      t23 = 1.D0+t22
      t24 = t23**(1.D0/3.D0)
      t27 = 1.D0-1.D0*t22
      t28 = t27**(1.D0/3.D0)
      t30 = t24*t23+t28*t27-2.D0
      t34 = dlog(t6)
      t36 = t5*t34
      t42 = -0.1555D-1*t34+0.211D-1-0.80645563816922D-3*t36
     #+0.421838333811592D-2*t5
      t46 = piecewise(t7,-t13+0.1923661050931536D1*t19*t30,0.311D-1
     #*t34-0.48D-1+0.12407009817988D-2*t36-0.719606569443304D-2*t5
     #+0.1923661050931536D1*t42*t30)
      t48 = dsqrt(sigma)
      t50 = t5**2
      t52 = 0.2568D-2+0.1443307452126544D-1*t5+0.2843543831490386D-5*t50
      t56 = 1.D0+0.5411317332115466D1*t5+0.1816419932959077D0*t50
     #+0.1763993811759022D-1*t4
      t57 = 1/t56
      t59 = 0.1667D-2+t52*t57
      t60 = 1/t59
      t61 = t48*t60
      t62 = rho**(1.D0/6.D0)
      t64 = 1/t62/rho
      t67 = dexp(-0.81290825D-3*t61*t64)
      t68 = t67*t59
      t69 = rho**(1.D0/3.D0)
      t71 = 1/t69/rho
      t72 = sigma*t71
      t74 = 0.5D0+0.5D0*t22
      t75 = t74**(1.D0/3.D0)
      t76 = t75**2
      t79 = 0.5D0-0.5D0*t22
      t80 = t79**(1.D0/3.D0)
      t81 = t80**2
      t83 = t76*t74+t81*t79
      t84 = dsqrt(t83)
      t85 = 1/t84
      t86 = t72*t85
      zk(i) = rho*t46+0.7937005259840997D0*t68*t86
      t93 = 0.1333333333333333D1*t24*t4-0.1333333333333333D1*t28*t4
      t98 = piecewise(t7,0.1923661050931536D1*t19
     #*t93,0.1923661050931536D1*t42*t93)
      t100 = t68*sigma
      t102 = 1/t84/t83
      t103 = t71*t102
      t108 = 0.8333333333333333D0*t76*t4-0.8333333333333333D0*t81*t4
      t109 = t103*t108
      t112 = t11**2
      t113 = 1/t112
      t114 = t8**2
      t115 = t114**2
      t116 = t115*t8
      t117 = 1/t116
      t118 = rho**2
      t119 = 1/t118
      t120 = t117*t119
      t122 = 1/t50
      t123 = t122*t119
      t125 = -0.1382147652361066D0*t120-0.6894161788861999D-1*t123
      t127 = 0.1423D0*t113*t125
      t128 = t16**2
      t129 = 1/t128
      t132 = -0.1835293601259385D0*t120-0.5399117105794445D-1*t123
      t135 = 0.843D-1*t129*t132-t127
      t138 = t24*t21
      t141 = t28*t21
      t144 = -0.1333333333333333D1*t138*t119+0.1333333333333333D1*t141
     #*t119
      t149 = t122*t34
      t150 = t149*t119
      t152 = t5*t4
      t159 = 0.5183333333333333D-2*t4+0.2688185460564067D-3*t150
     #+0.2688185460564067D-3*t152-0.1406127779371973D-2*t123
      t165 = piecewise(t7,t127+0.1923661050931536D1*t135*t30
     #+0.1923661050931536D1*t19*t144,-0.1036666666666667D-1*t4
     #-0.4135669939329333D-3*t150-0.4135669939329333D-3*t152
     #+0.2398688564811013D-2*t123+0.1923661050931536D1*t159*t30
     #+0.1923661050931536D1*t42*t144)
      t166 = rho*t165
      t167 = t59**2
      t168 = 1/t167
      t169 = t48*t168
      t171 = 1/t5
      t172 = t171*t119
      t174 = -0.4811024840421814D-2*t123-0.1895695887660258D-5*t172
      t176 = t56**2
      t177 = 1/t176
      t178 = t52*t177
      t182 = -0.1803772444038489D1*t123-0.1210946621972718D0*t172
     #-0.1763993811759022D-1*t119
      t185 = t174*t57-1.D0*t178*t182
      t186 = t64*t185
      t190 = 1/t62/t118
      t193 = 0.81290825D-3*t169*t186+0.9483929583333333D-3*t61*t190
      t194 = t193*t67
      t195 = t194*t59
      t197 = 0.7937005259840997D0*t195*t86
      t198 = t67*t185
      t200 = 0.7937005259840997D0*t198*t86
      t202 = 1/t69/t118
      t204 = sigma*t202*t85
      t206 = 0.10582673679788D1*t68*t204
      t207 = t76*t21
      t210 = t81*t21
      t213 = -0.8333333333333333D0*t207*t119+0.8333333333333333D0*t210
     #*t119
      t214 = t103*t213
      t216 = 0.3968502629920499D0*t100*t214
      vrhoa(i) = rho*t98-0.3968502629920499D0*t100*t109+t46+t166+t197
     #+t200-t206-t216
      t217 = -t93
      t222 = piecewise(t7,0.1923661050931536D1*t19
     #*t217,0.1923661050931536D1*t42*t217)
      t224 = -t108
      t225 = t103*t224
      vrhob(i) = rho*t222-0.3968502629920499D0*t100*t225+t46+t166+t197
     #+t200-t206-t216
      t228 = dsqrt(rho)
      t230 = 1/t228/t118
      t231 = t48*t230
      t232 = t67*t85
      t233 = t231*t232
      t235 = t71*t85
      t236 = t68*t235
      vsigmaaa(i) = -0.322602852800907D-3*t233+0.7937005259840997D0*t236
      vsigmaab(i) = -0.645205705601814D-3*t233+0.1587401051968199D1*t236
      vsigmabb(i) = vsigmaaa(i)
      t240 = 2.D0*t165
      t244 = t185**2
      t252 = 1/t50/t4
      t253 = t118**2
      t254 = 1/t253
      t255 = t252*t254
      t257 = t118*rho
      t258 = 1/t257
      t259 = t122*t258
      t262 = 1/t152*t254
      t264 = t171*t258
      t274 = t182**2
      t285 = (-0.3207349893614542D-2*t255+0.9622049680843627D-2*t259
     #-0.6318986292200858D-6*t262+0.3791391775320515D-5*t264)*t57-2.D0
     #*t174*t177*t182+2.D0*t52/t176/t56*t274-1.D0*t178*(
     #-0.1202514962692326D1*t255+0.3607544888076978D1*t259
     #-0.4036488739909061D-1*t262+0.2421893243945437D0*t264
     #+0.3527987623518044D-1*t258)
      t297 = 0.7937005259840997D0*(-0.16258165D-2*t48/t167/t59*t64
     #*t244-0.1896785916666667D-2*t169*t190*t185+0.81290825D-3*t169
     #*t64*t285-0.2054851409722222D-2*t61/t62/t257)*t67*t59*t86
      t298 = t193**2
      t302 = 0.7937005259840997D0*t298*t67*t59*t86
      t306 = 0.7937005259840997D0*t195*t72*t102*t213
      t308 = 0.2116534735957599D1*t195*t204
      t311 = 0.1587401051968199D1*t194*t185*t86
      t314 = 0.7937005259840997D0*t67*t285*t86
      t315 = t198*sigma
      t317 = 0.7937005259840997D0*t315*t214
      t319 = 0.2116534735957599D1*t198*t204
      t320 = t24**2
      t321 = 1/t320
      t324 = t28**2
      t325 = 1/t324
      t328 = 0.4444444444444444D0*t321*t119+0.4444444444444444D0*t325
     #*t119
      t333 = piecewise(t7,0.1923661050931536D1*t19
     #*t328,0.1923661050931536D1*t42*t328)
      t334 = rho*t333
      t335 = t83**2
      t338 = t71/t84/t335
      t339 = t108**2
      t343 = 1/t75
      t346 = 1/t80
      t349 = 0.2777777777777778D0*t343*t119+0.2777777777777778D0*t346
     #*t119
      t352 = 0.3968502629920499D0*t100*t103*t349
      t353 = t240+t297+t302-t306-t308+t311+t314-t317-t319+t334
     #+0.5952753944880748D0*t100*t338*t339-t352
      t355 = t202*t102
      t358 = 0.10582673679788D1*t100*t355*t213
      t364 = 0.2469290525283866D1*t68*sigma/t69/t257*t85
      t375 = -0.2777777777777778D0*t343*t21*t258-0.8333333333333333D0
     #*t76*t119-0.2777777777777778D0*t346*t21*t258
     #+0.8333333333333333D0*t81*t119
      t377 = t100*t103*t375
      t379 = t315*t109
      t383 = t195*t72*t102*t108
      t397 = -0.4444444444444444D0*t321*t21*t258-0.1333333333333333D1
     #*t24*t119-0.4444444444444444D0*t325*t21*t258
     #+0.1333333333333333D1*t28*t119
      t406 = piecewise(t7,0.1923661050931536D1*t135*t93
     #+0.1923661050931536D1*t19*t397,0.1923661050931536D1*t159*t93
     #+0.1923661050931536D1*t42*t397)
      t407 = rho*t406
      t411 = t100*t338*t213*t108
      t414 = t100*t355*t108
      t416 = t213**2
      t419 = 0.5952753944880748D0*t100*t338*t416
      t420 = t21**2
      t434 = 0.3968502629920499D0*t100*t103*(0.2777777777777778D0*t343
     #*t420*t254+0.1666666666666667D1*t207*t258+0.2777777777777778D0
     #*t346*t420*t254-0.1666666666666667D1*t210*t258)
      t437 = t125**2
      t439 = 0.2846D0/t112/t11*t437
      t442 = 1/t116/t4*t254
      t444 = t117*t258
      t450 = 0.1423D0*t113*(-0.1151789710300888D0*t442
     #+0.2764295304722132D0*t444-0.4596107859241333D-1*t255
     #+0.13788323577724D0*t259)
      t453 = t132**2
      t478 = 0.4444444444444444D0*t321*t420*t254+0.2666666666666667D1
     #*t138*t258+0.4444444444444444D0*t325*t420*t254
     #-0.2666666666666667D1*t141*t258
      t484 = t252*t34*t254
      t487 = t149*t258
      t489 = t5*t119
      t506 = piecewise(t7,-t439+t450+0.1923661050931536D1*(-0.1686D0
     #/t128/t16*t453+0.843D-1*t129*(-0.1529411334382821D0*t442
     #+0.367058720251877D0*t444-0.3599411403862963D-1*t255
     #+0.1079823421158889D0*t259)+t439-t450)*t30+0.3847322101863073D1
     #*t135*t144+0.1923661050931536D1*t19*t478,0.1036666666666667D-1
     #*t119-0.2757113292886222D-3*t484-0.4521665800333405D-2*t259
     #+0.8271339878658667D-3*t487+0.4135669939329333D-3*t489
     #+0.1599125709874009D-2*t255+0.1923661050931536D1*(
     #-0.5183333333333333D-2*t119+0.1792123640376044D-3*t484
     #+0.2633043194706342D-2*t259-0.5376370921128133D-3*t487
     #-0.2688185460564067D-3*t489-0.9374185195813156D-3*t255)*t30
     #+0.3847322101863073D1*t159*t144+0.1923661050931536D1*t42*t478)
      t507 = rho*t506
      t508 = 2.D0*t98+t358+t364-0.7937005259840997D0*t377
     #-0.7937005259840997D0*t379-0.7937005259840997D0*t383+2.D0*t407
     #+0.119055078897615D1*t411+0.10582673679788D1*t414+t419-t434+t507
      v2rhoa2(i) = t353+t508
      t510 = t224**2
      t516 = -t397
      t525 = piecewise(t7,0.1923661050931536D1*t135*t217
     #+0.1923661050931536D1*t19*t516,0.1923661050931536D1*t159*t217
     #+0.1923661050931536D1*t42*t516)
      t526 = rho*t525
      t530 = -t100*t103*t375
      t532 = t315*t225
      t534 = t240+t297+2.D0*t222+0.5952753944880748D0*t100*t338*t510
     #+t302-t306-t308+t311+t314+2.D0*t526-0.7937005259840997D0*t530
     #-0.7937005259840997D0*t532
      t536 = t100*t355*t224
      t540 = t100*t338*t213*t224
      t544 = t195*t72*t102*t224
      t546 = 0.10582673679788D1*t536+0.119055078897615D1*t540
     #-0.7937005259840997D0*t544-t317-t319+t334-t352+t358+t364+t419
     #-t434+t507
      v2rhob2(i) = t534+t546
      t560 = -t328
      t565 = piecewise(t7,0.1923661050931536D1*t19
     #*t560,0.1923661050931536D1*t42*t560)
      t567 = t526-0.3968502629920499D0*t530-0.3968502629920499D0*t532
     #+0.5291336839893998D0*t536+0.5952753944880748D0*t540
     #-0.3968502629920499D0*t544+0.3968502629920499D0*t100*t103*t349
     #+0.5952753944880748D0*t100*t338*t108*t224+t98+rho*t565+t240+t222
     #+t297+t302-t306
      t573 = -t308+t311+t314-t319-t317+t358+t364-0.3968502629920499D0
     #*t377-0.3968502629920499D0*t379-0.3968502629920499D0*t383+t407
     #+0.5952753944880748D0*t411+t507+0.5291336839893998D0*t414-t434
     #+t419
      v2rhoab(i) = t567+t573
      t574 = t67*t102
      t576 = t231*t574*t108
      t578 = t68*t109
      t580 = 1/t48
      t590 = (0.406454125D-3*t580*t168*t186+0.4741964791666667D-3*t580
     #*t60*t190)*t67*t59*t86
      t591 = 0.7937005259840997D0*t590
      t595 = t193*t48*t230*t67*t85
      t596 = 0.322602852800907D-3*t595
      t599 = t194*t59*t71*t85
      t600 = 0.7937005259840997D0*t599
      t603 = t61*t230*t198*t85
      t604 = 0.322602852800907D-3*t603
      t605 = t198*t235
      t606 = 0.7937005259840997D0*t605
      t610 = t48/t228/t257*t232
      t611 = 0.430137137067876D-3*t610
      t613 = t68*t202*t85
      t614 = 0.10582673679788D1*t613
      t616 = t231*t574*t213
      t617 = 0.1613014264004535D-3*t616
      t618 = t68*t214
      t619 = 0.3968502629920499D0*t618
      v2rhoasigmaaa(i) = 0.1613014264004535D-3*t576
     #-0.3968502629920499D0*t578+t591-t596+t600-t604+t606+t611-t614
     #+t617-t619
      t622 = 0.1587401051968199D1*t590
      t623 = 0.645205705601814D-3*t595
      t624 = 0.1587401051968199D1*t599
      t625 = 0.645205705601814D-3*t603
      t626 = 0.1587401051968199D1*t605
      t627 = 0.8602742741357521D-3*t610
      t628 = 0.2116534735957599D1*t613
      t629 = 0.322602852800907D-3*t616
      t630 = 0.7937005259840997D0*t618
      v2rhoasigmaab(i) = 0.322602852800907D-3*t576
     #-0.7937005259840997D0*t578+t622-t623+t624-t625+t626+t627-t628
     #+t629-t630
      v2rhoasigmabb(i) = v2rhoasigmaaa(i)
      t632 = t231*t574*t224
      t634 = t68*t225
      v2rhobsigmaaa(i) = 0.1613014264004535D-3*t632
     #-0.3968502629920499D0*t634+t591-t596+t600-t604+t606+t611-t614
     #+t617-t619
      v2rhobsigmaab(i) = 0.322602852800907D-3*t632
     #-0.7937005259840997D0*t634+t622-t623+t624-t625+t626+t627-t628
     #+t629-t630
      v2rhobsigmabb(i) = v2rhobsigmaaa(i)
      t639 = t580*t230*t232
      t641 = t69**2
      t645 = 1/t641/t257*t60*t232
      v2sigmaaa2(i) = -0.4839042792013605D-3*t639
     #+0.1311232602576965D-6*t645
      v2sigmaaaab(i) = -0.9678085584027211D-3*t639
     #+0.2622465205153929D-6*t645
      v2sigmaaabb(i) = v2sigmaaa2(i)
      v2sigmaab2(i) = -0.1935617116805442D-2*t639
     #+0.5244930410307859D-6*t645
      v2sigmaabbb(i) = v2sigmaaaab(i)
      v2sigmabb2(i) = v2sigmaaabb(i)
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vrhob(i) = 0.0d0
      v2rhoa2(i) = 0.0d0
      v2rhob2(i) = 0.0d0
      v2rhoab(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      vsigmaab(i) = 0.0d0
      vsigmabb(i) = 0.0d0
      v2rhoasigmaaa(i) = 0.0d0
      v2rhoasigmaab(i) = 0.0d0
      v2rhoasigmabb(i) = 0.0d0
      v2rhobsigmaaa(i) = 0.0d0
      v2rhobsigmaab(i) = 0.0d0
      v2rhobsigmabb(i) = 0.0d0
      v2sigmaaa2(i) = 0.0d0
      v2sigmaab2(i) = 0.0d0
      v2sigmabb2(i) = 0.0d0
      v2sigmaaaab(i) = 0.0d0
      v2sigmaaabb(i) = 0.0d0
      v2sigmaabbb(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end
      
      
      subroutine rks_c_p86
     & (ideriv,npt,rhoa1,sigmaaa1,
     &  zk,vrhoa,vsigmaaa,
     &  v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
c
c     J.P. Perdew
c     Density-functional approximation for the correlation energy of
c     the inhomogeneous electron gas
c     Phys. Rev. B33 (1986) 8822-8824
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit real*8 (a-h,o-z)
      integer ideriv,npt
      real*8 rhoa1(npt)
      real*8 sigmaaa1(npt)
      real*8 zk(npt),vrhoa(npt),vsigmaaa(npt)
      real*8 v2rhoa2(npt),v2rhoasigmaaa(npt),v2sigmaaa2(npt)
      parameter(tol=1.0d-20)
      logical t5
      
      if(ideriv.eq.0) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = 1/rho
      t3 = t2**(1.D0/3.D0)
      t4 = 0.6203504908994D0*t3
      t6 = t2**(1.D0/6.D0)
      t12 = dlog(t4)
      t18 = piecewise(1.D0 .le. t4,-0.1423D0/(1.D0
     #+0.8292885914166397D0*t6+0.20682485366586D0*t3),0.311D-1*t12
     #-0.48D-1+0.12407009817988D-2*t3*t12-0.719606569443304D-2*t3)
      t20 = dsqrt(sigma)
      t22 = t3**2
      t31 = 0.1667D-2+(0.2568D-2+0.1443307452126544D-1*t3
     #+0.2843543831490386D-5*t22)/(1.D0+0.5411317332115466D1*t3
     #+0.1816419932959077D0*t22+0.1763993811759022D-1*t2)
      t34 = rho**(1.D0/6.D0)
      t39 = dexp(-0.81290825D-3*t20/t31/t34/rho)
      t41 = rho**(1.D0/3.D0)
      zk(i) = rho*t18+0.1D1*t39*t31*sigma/t41/rho
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = 1/rho
      t3 = t2**(1.D0/3.D0)
      t4 = 0.6203504908994D0*t3
      t5 = 1.D0 .le. t4
      t6 = t2**(1.D0/6.D0)
      t9 = 1.D0+0.8292885914166397D0*t6+0.20682485366586D0*t3
      t12 = dlog(t4)
      t18 = piecewise(t5,-0.1423D0/t9,0.311D-1*t12-0.48D-1
     #+0.12407009817988D-2*t3*t12-0.719606569443304D-2*t3)
      t20 = dsqrt(sigma)
      t22 = t3**2
      t24 = 0.2568D-2+0.1443307452126544D-1*t3+0.2843543831490386D-5*t22
      t28 = 1.D0+0.5411317332115466D1*t3+0.1816419932959077D0*t22
     #+0.1763993811759022D-1*t2
      t29 = 1/t28
      t31 = 0.1667D-2+t24*t29
      t33 = t20/t31
      t34 = rho**(1.D0/6.D0)
      t36 = 1/t34/rho
      t39 = dexp(-0.81290825D-3*t33*t36)
      t40 = t39*t31
      t41 = rho**(1.D0/3.D0)
      t43 = 1/t41/rho
      t44 = sigma*t43
      zk(i) = rho*t18+0.1D1*t40*t44
      t47 = piecewise(t5,0.D0,0.D0)
      t49 = t9**2
      t51 = t6**2
      t52 = t51**2
      t55 = rho**2
      t56 = 1/t55
      t59 = 1/t22
      t60 = t59*t56
      t73 = piecewise(t5,0.1423D0/t49*(-0.1382147652361066D0/t52/t6
     #*t56-0.6894161788861999D-1*t60),-0.1036666666666667D-1*t2
     #-0.4135669939329333D-3*t59*t12*t56-0.4135669939329333D-3*t3*t2
     #+0.2398688564811013D-2*t60)
      t75 = t31**2
      t80 = 1/t3*t56
      t84 = t28**2
      t93 = (-0.4811024840421814D-2*t60-0.1895695887660258D-5*t80)*t29
     #-1.D0*t24/t84*(-0.1803772444038489D1*t60-0.1210946621972718D0
     #*t80-0.1763993811759022D-1*t56)
      vrhoa(i) = rho*t47+t18+rho*t73+0.1D1*(0.81290825D-3*t20/t75*t36
     #*t93+0.9483929583333333D-3*t33/t34/t55)*t39*t31*sigma*t43+0.1D1
     #*t39*t93*t44-0.1333333333333333D1*t40*sigma/t41/t55
      t115 = dsqrt(rho)
      vsigmaaa(i) = -0.16258165D-2*t20/t115/t55*t39+0.4D1*t40*t43
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = 1/rho
      t3 = t2**(1.D0/3.D0)
      t4 = 0.6203504908994D0*t3
      t5 = 1.D0 .le. t4
      t6 = t2**(1.D0/6.D0)
      t9 = 1.D0+0.8292885914166397D0*t6+0.20682485366586D0*t3
      t11 = 0.1423D0/t9
      t12 = dlog(t4)
      t14 = t3*t12
      t18 = piecewise(t5,-t11,0.311D-1*t12-0.48D-1+0.12407009817988D-2
     #*t14-0.719606569443304D-2*t3)
      t20 = dsqrt(sigma)
      t22 = t3**2
      t24 = 0.2568D-2+0.1443307452126544D-1*t3+0.2843543831490386D-5*t22
      t28 = 1.D0+0.5411317332115466D1*t3+0.1816419932959077D0*t22
     #+0.1763993811759022D-1*t2
      t29 = 1/t28
      t31 = 0.1667D-2+t24*t29
      t32 = 1/t31
      t33 = t20*t32
      t34 = rho**(1.D0/6.D0)
      t36 = 1/t34/rho
      t39 = dexp(-0.81290825D-3*t33*t36)
      t40 = t39*t31
      t41 = rho**(1.D0/3.D0)
      t43 = 1/t41/rho
      t44 = sigma*t43
      zk(i) = rho*t18+0.1D1*t40*t44
      t47 = piecewise(t5,0.D0,0.D0)
      t48 = rho*t47
      t49 = t9**2
      t50 = 1/t49
      t51 = t6**2
      t52 = t51**2
      t53 = t52*t6
      t54 = 1/t53
      t55 = rho**2
      t56 = 1/t55
      t59 = 1/t22
      t60 = t59*t56
      t62 = -0.1382147652361066D0*t54*t56-0.6894161788861999D-1*t60
      t66 = t59*t12
      t69 = t3*t2
      t73 = piecewise(t5,0.1423D0*t50*t62,-0.1036666666666667D-1*t2
     #-0.4135669939329333D-3*t66*t56-0.4135669939329333D-3*t69
     #+0.2398688564811013D-2*t60)
      t75 = t31**2
      t76 = 1/t75
      t77 = t20*t76
      t79 = 1/t3
      t80 = t79*t56
      t82 = -0.4811024840421814D-2*t60-0.1895695887660258D-5*t80
      t84 = t28**2
      t85 = 1/t84
      t86 = t24*t85
      t90 = -0.1803772444038489D1*t60-0.1210946621972718D0*t80
     #-0.1763993811759022D-1*t56
      t93 = t82*t29-1.D0*t86*t90
      t94 = t36*t93
      t98 = 1/t34/t55
      t101 = 0.81290825D-3*t77*t94+0.9483929583333333D-3*t33*t98
      t102 = t101*t39
      t103 = t31*sigma
      t104 = t103*t43
      t107 = t39*t93
      t111 = 1/t41/t55
      t112 = sigma*t111
      vrhoa(i) = t48+t18+rho*t73+0.1D1*t102*t104+0.1D1*t107*t44
     #-0.1333333333333333D1*t40*t112
      t115 = dsqrt(rho)
      t117 = 1/t115/t55
      vsigmaaa(i) = -0.16258165D-2*t20*t117*t39+0.4D1*t40*t43
      t129 = (-0.843D-1/(1.D0+0.1101176160755631D1*t6
     #+0.1619735131738333D0*t3)+t11)*t56
      t135 = (-0.1555D-1*t12+0.211D-1-0.80645563816922D-3*t14
     #+0.421838333811592D-2*t3)*t56
      t137 = piecewise(t5,0.1709920934161366D1
     #*t129,0.1709920934161366D1*t135)
      t139 = t55*rho
      t147 = t62**2
      t152 = t55**2
      t153 = 1/t152
      t156 = 1/t139
      t160 = 1/t22/t2
      t161 = t160*t153
      t163 = t59*t156
      t180 = piecewise(t5,-0.2846D0/t49/t9*t147+0.1423D0*t50*(
     #-0.1151789710300888D0/t53/t2*t153+0.2764295304722132D0*t54*t156
     #-0.4596107859241333D-1*t161+0.13788323577724D0*t163
     #),0.1036666666666667D-1*t56-0.2757113292886222D-3*t160*t12*t153
     #-0.4521665800333405D-2*t163+0.8271339878658667D-3*t66*t156
     #+0.4135669939329333D-3*t3*t56+0.1599125709874009D-2*t161)
      t185 = piecewise(t5,-0.1709920934161366D1*t129,
     #-0.1709920934161366D1*t135)
      t193 = t93**2
      t203 = 1/t69*t153
      t205 = t79*t156
      t215 = t90**2
      t226 = (-0.3207349893614542D-2*t161+0.9622049680843627D-2*t163
     #-0.6318986292200858D-6*t203+0.3791391775320515D-5*t205)*t29-2.D0
     #*t82*t85*t90+2.D0*t24/t84/t28*t215-1.D0*t86*(
     #-0.1202514962692326D1*t161+0.3607544888076978D1*t163
     #-0.4036488739909061D-1*t203+0.2421893243945437D0*t205
     #+0.3527987623518044D-1*t156)
      t247 = t101**2
      v2rhoa2(i) = rho*t137+0.6222222222222222D1*t40*sigma/t41/t139
     #+2.D0*rho*t180+rho*t185-0.5333333333333333D1*t102*t103*t111
     #+0.2D1*(-0.16258165D-2*t20/t75/t31*t36*t193
     #-0.1896785916666667D-2*t77*t98*t93+0.81290825D-3*t77*t36*t226
     #-0.2054851409722222D-2*t33/t34/t139)*t39*t104+0.4D1*t102*t93
     #*sigma*t43+0.2D1*t39*t226*t44-0.5333333333333333D1*t107*t112
     #+0.2D1*t247*t39*t104+4.D0*t48+4.D0*t47+4.D0*t73
      t254 = 1/t20
      t266 = t117*t39
      v2rhoasigmaaa(i) = 0.4D1*(0.406454125D-3*t254*t76*t94
     #+0.4741964791666667D-3*t254*t32*t98)*t39*t104-0.16258165D-2*t101
     #*t20*t266+0.4D1*t102*t31*t43-0.16258165D-2*t33*t266*t93+0.4D1
     #*t107*t43+0.2167755333333333D-2*t20/t115/t139*t39
     #-0.5333333333333333D1*t40*t111
      t287 = t41**2
      v2sigmaaa2(i) = -0.9754899D-2*t254*t117*t39+0.264327929167225D-5
     #/t287/t139*t32*t39
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      v2rhoa2(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      v2rhoasigmaaa(i) = 0.0d0
      v2sigmaaa2(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end


      subroutine uks_c_lyp
     & (ideriv,npt,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,
     &  zk,vrhoa,vrhob,vsigmaaa,vsigmabb,vsigmaab,
     &  v2rhoa2,v2rhob2,v2rhoab,
     &  v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb,
     &  v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,
     &  v2sigmaaa2,v2sigmaaaab,v2sigmaaabb,
     &  v2sigmaab2,v2sigmaabbb,v2sigmabb2)
c
c     C. Lee, W. Yang, and R.G. Parr
c     Development of the Colle-Salvetti correlation-energy formula into
c     a functional of the electron density
c     Phys. Rev. B37 (1988) 785-789
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit real*8 (a-h,o-z)
      integer ideriv,npt
      real*8 rhoa1(npt),rhob1(npt)
      real*8 sigmaaa1(npt),sigmabb1(npt),sigmaab1(npt)
      real*8 zk(npt),vrhoa(npt),vrhob(npt)
      real*8 vsigmaaa(npt),vsigmabb(npt),vsigmaab(npt)
      real*8 v2rhoa2(npt),v2rhob2(npt),v2rhoab(npt)
      real*8 v2rhoasigmaaa(npt),v2rhoasigmaab(npt)
      real*8 v2rhoasigmabb(npt),v2rhobsigmabb(npt)
      real*8 v2rhobsigmaab(npt),v2rhobsigmaaa(npt)
      real*8 v2sigmaaa2(npt),v2sigmaaaab(npt),v2sigmaaabb(npt)
      real*8 v2sigmaab2(npt),v2sigmaabbb(npt),v2sigmabb2(npt)
      parameter(tol=1.0d-20)
      
      if (ideriv.eq.0) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmabb
      zk(i) = 0.D0
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      zk(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rho**(1.D0/3.D0)
      t5 = 1/t4
      t8 = 1/(1.D0+0.349D0*t5)
      t10 = 1/rho
      t11 = rhob*t10
      t14 = 0.2533D0*t5
      t15 = dexp(-t14)
      t17 = rho**2
      t19 = t4**2
      t23 = rhoa**2
      t24 = rhoa**(1.D0/3.D0)
      t25 = t24**2
      t28 = rhob**2
      t29 = rhob**(1.D0/3.D0)
      t30 = t29**2
      t34 = t5*t8
      t56 = 0.6666666666666667D0*t17
      zk(i) = -0.19672D0*t8*rhoa*t11-0.649176D-2*t15*t8/t19/t17/rho*
     #(rhoa*rhob*(0.3646239897876478D2*t25*t23+0.3646239897876478D2
     #*t30*t28+(0.2611111111111111D1-0.9850555555555556D-1*t5
     #-0.1357222222222222D0*t34)*sigma-1.D0*(0.25D1
     #-0.1407222222222222D-1*t5-0.1938888888888889D-1*t34)*(sigmaaa
     #+sigmabb)-0.1111111111111111D0*(t14+0.349D0*t34-11.D0)*(rhoa*t10
     #*sigmaaa+t11*sigmabb))-0.6666666666666667D0*t17*sigma+(t56-1.D0
     #*t23)*sigmabb+(t56-1.D0*t28)*sigmaaa)
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmabb
      zk(i) = 0.D0
      vrhoa(i) = 0.D0
      vrhob(i) = 0.D0
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      zk(i) = 0.D0
      vrhoa(i) = 0.D0
      vrhob(i) = 0.D0
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rho**(1.D0/3.D0)
      t5 = 1/t4
      t7 = 1.D0+0.349D0*t5
      t8 = 1/t7
      t9 = t8*rhoa
      t10 = 1/rho
      t11 = rhob*t10
      t14 = 0.2533D0*t5
      t15 = dexp(-t14)
      t16 = t15*t8
      t17 = rho**2
      t19 = t4**2
      t21 = 1/t19/t17/rho
      t22 = rhoa*rhob
      t23 = rhoa**2
      t24 = rhoa**(1.D0/3.D0)
      t25 = t24**2
      t28 = rhob**2
      t29 = rhob**(1.D0/3.D0)
      t30 = t29**2
      t34 = t5*t8
      t36 = 0.2611111111111111D1-0.9850555555555556D-1*t5
     #-0.1357222222222222D0*t34
      t41 = sigmaaa+sigmabb
      t45 = t14+0.349D0*t34-11.D0
      t49 = rhoa*t10*sigmaaa+t11*sigmabb
      t52 = 0.3646239897876478D2*t25*t23+0.3646239897876478D2*t30*t28
     #+t36*sigma-1.D0*(0.25D1-0.1407222222222222D-1*t5
     #-0.1938888888888889D-1*t34)*t41-0.1111111111111111D0*t45*t49
      t56 = 0.6666666666666667D0*t17
      t57 = 1.D0*t23
      t60 = 1.D0*t28
      t63 = t22*t52-0.6666666666666667D0*t17*sigma+(t56-t57)*sigmabb+
     #(t56-t60)*sigmaaa
      zk(i) = -0.19672D0*t9*t11-0.649176D-2*t16*t21*t63
      t73 = t45*t10
      t84 = t7**2
      t85 = 1/t84
      t91 = 0.2288509333333333D-1*t85*rhoa*rhob/t4/t17
      t92 = 1/t17
      t93 = rhob*t92
      t95 = 0.19672D0*t9*t93
      t96 = t17**2
      t98 = 1/t96/rho
      t102 = 0.548120936D-3*t98*t15*t8*t63
      t106 = 0.75520808D-3*t15*t85*t98*t63
      t111 = 0.2380312D-1*t16/t19/t96*t63
      t113 = 1/t4/rho
      t115 = t113*t8
      t119 = 1/t19/rho*t85
      t154 = 0.649176D-2*t16*t21*(t22*((0.3283518518518519D-1*t113
     #+0.4524074074074074D-1*t115-0.1578901851851852D-1*t119)*sigma
     #-1.D0*(0.4690740740740741D-2*t113+0.6462962962962963D-2*t115
     #-0.2255574074074074D-2*t119)*t41-0.1111111111111111D0*(
     #-0.8443333333333333D-1*t113-0.1163333333333333D0*t115
     #+0.4060033333333333D-1*t119)*t49-0.1111111111111111D0*t45*(-1.D0
     #*rhoa*t92*sigmaaa-1.D0*t93*sigmabb))-0.1333333333333333D1*rho
     #*sigma+0.1333333333333333D1*rho*sigmabb+0.1333333333333333D1*rho
     #*sigmaaa)
      vrhoa(i) = -0.19672D0*t8*rhob*t10-0.649176D-2*t16*t21*(rhob*t52
     #+t22*(0.9723306394337274D2*t25*rhoa-0.1111111111111111D0*t73
     #*sigmaaa)-2.D0*rhoa*sigmabb)-t91+t95-t102-t106+t111-t154
      vrhob(i) = -0.19672D0*t9*t10-0.649176D-2*t16*t21*(rhoa*t52+t22*
     #(0.9723306394337274D2*t30*rhob-0.1111111111111111D0*t73*sigmabb)
     #-2.D0*rhob*sigmaaa)-t91+t95-t102-t106+t111-t154
      t170 = 0.1407222222222222D-1*t5
      t171 = 0.1938888888888889D-1*t34
      t185 = t16*t21*(t22*t36-0.6666666666666667D0*t17)
      t186 = 0.649176D-2*t185
      vsigmaaa(i) = -0.649176D-2*t16*t21*(t22*(-0.25D1+t170+t171
     #-0.1111111111111111D0*t45*rhoa*t10)+t56-t60)-t186
      vsigmaab(i) = -0.1298352D-1*t185
      vsigmabb(i) = -0.649176D-2*t16*t21*(t22*(-0.25D1+t170+t171
     #-0.1111111111111111D0*t45*rhob*t10)+t56-t57)-t186
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vrhob(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      vsigmaab(i) = 0.0d0
      vsigmabb(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmabb
      zk(i) = 0.D0
      vrhoa(i) = 0.D0
      vrhob(i) = 0.D0
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      v2rhoa2(i) = 0.D0
      v2rhoab(i) = 0.D0
      v2rhob2(i) = 0.D0
      v2sigmaaa2(i) = 0.D0
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      v2sigmabb2(i) = 0.D0
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      zk(i) = 0.D0
      vrhoa(i) = 0.D0
      vrhob(i) = 0.D0
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      v2rhoa2(i) = 0.D0
      v2rhob2(i) = 0.D0
      v2rhoab(i) = 0.D0
      v2sigmaaa2(i) = 0.D0
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      v2sigmabb2(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rho**(1.D0/3.D0)
      t5 = 1/t4
      t7 = 1.D0+0.349D0*t5
      t8 = 1/t7
      t9 = t8*rhoa
      t10 = 1/rho
      t11 = rhob*t10
      t14 = 0.2533D0*t5
      t15 = dexp(-t14)
      t16 = t15*t8
      t17 = rho**2
      t18 = t17*rho
      t19 = t4**2
      t21 = 1/t19/t18
      t22 = rhoa*rhob
      t23 = rhoa**2
      t24 = rhoa**(1.D0/3.D0)
      t25 = t24**2
      t27 = 0.3646239897876478D2*t25*t23
      t28 = rhob**2
      t29 = rhob**(1.D0/3.D0)
      t30 = t29**2
      t32 = 0.3646239897876478D2*t30*t28
      t34 = t5*t8
      t36 = 0.2611111111111111D1-0.9850555555555556D-1*t5
     #-0.1357222222222222D0*t34
      t37 = t36*sigma
      t41 = sigmaaa+sigmabb
      t43 = 1.D0*(0.25D1-0.1407222222222222D-1*t5
     #-0.1938888888888889D-1*t34)*t41
      t45 = t14+0.349D0*t34-11.D0
      t49 = rhoa*t10*sigmaaa+t11*sigmabb
      t51 = 0.1111111111111111D0*t45*t49
      t52 = t27+t32+t37-t43-t51
      t56 = 0.6666666666666667D0*t17
      t57 = 1.D0*t23
      t60 = 1.D0*t28
      t63 = t22*t52-0.6666666666666667D0*t17*sigma+(t56-t57)*sigmabb+
     #(t56-t60)*sigmaaa
      zk(i) = -0.19672D0*t9*t11-0.649176D-2*t16*t21*t63
      t67 = t8*rhob
      t71 = t25*rhoa
      t73 = t45*t10
      t76 = 0.9723306394337274D2*t71-0.1111111111111111D0*t73*sigmaaa
      t80 = rhob*t52+t22*t76-2.D0*rhoa*sigmabb
      t84 = t7**2
      t85 = 1/t84
      t86 = t85*rhoa
      t88 = 1/t4/t17
      t91 = 0.2288509333333333D-1*t86*rhob*t88
      t92 = 1/t17
      t93 = rhob*t92
      t95 = 0.19672D0*t9*t93
      t96 = t17**2
      t97 = t96*rho
      t98 = 1/t97
      t99 = t98*t15
      t100 = t8*t63
      t102 = 0.548120936D-3*t99*t100
      t103 = t15*t85
      t106 = 0.75520808D-3*t103*t98*t63
      t108 = 1/t19/t96
      t111 = 0.2380312D-1*t16*t108*t63
      t113 = 1/t4/rho
      t115 = t113*t8
      t119 = 1/t19/rho*t85
      t121 = 0.3283518518518519D-1*t113+0.4524074074074074D-1*t115
     #-0.1578901851851852D-1*t119
      t132 = -0.8443333333333333D-1*t113-0.1163333333333333D0*t115
     #+0.4060033333333333D-1*t119
      t140 = -1.D0*rhoa*t92*sigmaaa-1.D0*t93*sigmabb
      t143 = t121*sigma-1.D0*(0.4690740740740741D-2*t113
     #+0.6462962962962963D-2*t115-0.2255574074074074D-2*t119)*t41
     #-0.1111111111111111D0*t132*t49-0.1111111111111111D0*t45*t140
      t151 = t22*t143-0.1333333333333333D1*rho*sigma
     #+0.1333333333333333D1*rho*sigmabb+0.1333333333333333D1*rho*sigmaaa
      t154 = 0.649176D-2*t16*t21*t151
      vrhoa(i) = -0.19672D0*t67*t10-0.649176D-2*t16*t21*t80-t91+t95
     #-t102-t106+t111-t154
      t158 = t30*rhob
      t162 = 0.9723306394337274D2*t158-0.1111111111111111D0*t73*sigmabb
      t166 = rhoa*t52+t22*t162-2.D0*rhob*sigmaaa
      vrhob(i) = -0.19672D0*t9*t10-0.649176D-2*t16*t21*t166-t91+t95
     #-t102-t106+t111-t154
      t170 = 0.1407222222222222D-1*t5
      t171 = 0.1938888888888889D-1*t34
      t172 = t45*rhoa
      t175 = -0.25D1+t170+t171-0.1111111111111111D0*t172*t10
      t177 = t22*t175+t56-t60
      t183 = t22*t36-0.6666666666666667D0*t17
      t185 = t16*t21*t183
      t186 = 0.649176D-2*t185
      vsigmaaa(i) = -0.649176D-2*t16*t21*t177-t186
      vsigmaab(i) = -0.1298352D-1*t185
      t188 = t45*rhob
      t191 = -0.25D1+t170+t171-0.1111111111111111D0*t188*t10
      t193 = t22*t191+t56-t57
      vsigmabb(i) = -0.649176D-2*t16*t21*t193-t186
      t198 = t85*rhob*t88
      t200 = t67*t92
      t203 = t99*t8*t80
      t206 = t103*t98*t80
      t209 = t16*t108*t80
      t212 = t132*t10
      t215 = t45*t92
      t222 = t16*t21*(rhob*t143+t22*(-0.1111111111111111D0*t212
     #*sigmaaa+0.1111111111111111D0*t215*sigmaaa))
      t237 = 0.1110812266666667D0*t16/t19/t97*t63
      t239 = t88*t8
      t243 = 1/t19/t17*t85
      t245 = 1/t18
      t247 = 1/t84/t7
      t248 = t245*t247
      t271 = rhob*t245
      t285 = 0.649176D-2*t16*t21*(t22*((-0.4378024691358025D-1*t88
     #-0.6032098765432099D-1*t239+0.3157803703703704D-1*t243
     #-0.3673578308641975D-2*t248)*sigma-1.D0*(-0.6254320987654321D-2
     #*t88-0.8617283950617284D-2*t239+0.4511148148148148D-2*t243
     #-0.5247969012345679D-3*t248)*t41-0.1111111111111111D0*
     #(0.1125777777777778D0*t88+0.1551111111111111D0*t239
     #-0.8120066666666667D-1*t243+0.9446344222222222D-2*t248)*t49
     #-0.2222222222222222D0*t132*t140-0.1111111111111111D0*t45*(2.D0
     #*rhoa*t245*sigmaaa+2.D0*t271*sigmabb))-0.1333333333333333D1
     #*sigma+0.1333333333333333D1*sigmabb+0.1333333333333333D1*sigmaaa)
      t286 = t96*t17
      t287 = 1/t286
      t290 = 0.6545136693333333D-2*t103*t287*t63
      t293 = 0.4750381445333333D-2*t287*t15*t100
      t298 = 0.7628364444444444D-1*t86*rhob/t4/t18
      t300 = rhob*t21
      t302 = 0.5324598382222222D-2*t247*rhoa*t300
      t305 = 0.1096241872D-2*t99*t8*t151
      t307 = 1/t4/t286
      t308 = t307*t15
      t310 = 0.4627967769626667D-4*t308*t100
      t313 = 0.1275294711093333D-3*t308*t85*t63
      t316 = 0.151041616D-2*t103*t98*t151
      t320 = 0.1757117466133333D-3*t15*t247*t307*t63
      t323 = 0.4760624D-1*t16*t108*t151
      t325 = 0.39344D0*t9*t271
      v2rhoa2(i) = -0.4577018666666667D-1*t198+0.39344D0*t200
     #-0.1096241872D-2*t203-0.151041616D-2*t206+0.4760624D-1*t209
     #-0.1298352D-1*t222-0.649176D-2*t16*t21*(2.D0*rhob*t76
     #+0.1620551065722879D3*t71*rhob-2.D0*sigmabb)-t237-t285+t290+t293
     #+t298-t302-t305-t310-t313-t316-t320+t323-t325
      t335 = t86*t88
      t337 = t9*t92
      t340 = t99*t8*t166
      t343 = t103*t98*t166
      t346 = t16*t108*t166
      t357 = t16*t21*(rhoa*t143+t22*(-0.1111111111111111D0*t212
     #*sigmabb+0.1111111111111111D0*t215*sigmabb))
      v2rhob2(i) = -t237-t285+t290+t293+t298-t302-t305-t310-t313-t316
     #-t320+t323-0.649176D-2*t16*t21*(2.D0*rhoa*t162
     #+0.1620551065722879D3*rhoa*t158-2.D0*sigmaaa)
     #-0.4577018666666667D-1*t335+0.39344D0*t337-0.1096241872D-2*t340
     #-0.151041616D-2*t343+0.4760624D-1*t346-0.1298352D-1*t357-t325
      t377 = 0.19672D0*t337-0.2288509333333333D-1*t335-0.548120936D-3
     #*t340-0.75520808D-3*t343+0.2380312D-1*t346-0.19672D0*t8*t10
     #-0.649176D-2*t16*t21*(t27+t32+t37-t43-t51+rhob*t162+rhoa*t76)
     #-0.649176D-2*t357-0.548120936D-3*t203-0.75520808D-3*t206
     #+0.2380312D-1*t209-0.649176D-2*t222-t237
      t380 = -t285+t290+t293+t298-t302-t305-t310-t313-t316-t320+t323
     #-t325-0.2288509333333333D-1*t198+0.19672D0*t200
      v2rhoab(i) = t377+t380
      t383 = 0.1111111111111111D0*t22*t73
      t390 = 0.548120936D-3*t99*t8*t177
      t393 = 0.75520808D-3*t103*t98*t177
      t396 = 0.2380312D-1*t16*t108*t177
      t397 = 0.4690740740740741D-2*t113
      t398 = 0.6462962962962963D-2*t115
      t399 = 0.2255574074074074D-2*t119
      t407 = 0.1333333333333333D1*rho
      t411 = 0.649176D-2*t16*t21*(t22*(-t397-t398+t399
     #-0.1111111111111111D0*t132*rhoa*t10+0.1111111111111111D0*t172
     #*t92)+t407)
      t413 = t16*t300*t36
      t414 = 0.649176D-2*t413
      t416 = t99*t8*t183
      t417 = 0.548120936D-3*t416
      t419 = t103*t98*t183
      t420 = 0.75520808D-3*t419
      t422 = t16*t108*t183
      t423 = 0.2380312D-1*t422
      t428 = t16*t21*(t22*t121-0.1333333333333333D1*rho)
      t429 = 0.649176D-2*t428
      v2rhoasigmaaa(i) = -0.649176D-2*t16*t21*(rhob*t175-t383)-t390
     #-t393+t396-t411-t414-t417-t420+t423-t429
      t431 = 0.1096241872D-2*t416
      t432 = 0.151041616D-2*t419
      t433 = 0.4760624D-1*t422
      t434 = 0.1298352D-1*t428
      v2rhoasigmaab(i) = -0.1298352D-1*t413-t431-t432+t433-t434
      t443 = 0.548120936D-3*t99*t8*t193
      t446 = 0.75520808D-3*t103*t98*t193
      t449 = 0.2380312D-1*t16*t108*t193
      t460 = 0.649176D-2*t16*t21*(t22*(-t397-t398+t399
     #-0.1111111111111111D0*t132*rhob*t10+0.1111111111111111D0*t188
     #*t92)+t407)
      v2rhoasigmabb(i) = -0.649176D-2*t16*t21*(rhob*t191-2.D0*rhoa)
     #-t443-t446+t449-t460-t414-t417-t420+t423-t429
      t469 = t16*t21*rhoa*t36
      t470 = 0.649176D-2*t469
      v2rhobsigmaaa(i) = -0.649176D-2*t16*t21*(rhoa*t175-2.D0*rhob)
     #-t390-t393+t396-t411-t470-t417-t420+t423-t429
      v2rhobsigmaab(i) = -0.1298352D-1*t469-t431-t432+t433-t434
      v2rhobsigmabb(i) = -0.649176D-2*t16*t21*(rhoa*t191-t383)-t443
     #-t446+t449-t460-t470-t417-t420+t423-t429
      v2sigmaaa2(i) = 0.D0
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      v2sigmabb2(i) = 0.D0
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vrhob(i) = 0.0d0
      v2rhoa2(i) = 0.0d0
      v2rhob2(i) = 0.0d0
      v2rhoab(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      vsigmaab(i) = 0.0d0
      vsigmabb(i) = 0.0d0
      v2rhoasigmaaa(i) = 0.0d0
      v2rhoasigmaab(i) = 0.0d0
      v2rhoasigmabb(i) = 0.0d0
      v2rhobsigmaaa(i) = 0.0d0
      v2rhobsigmaab(i) = 0.0d0
      v2rhobsigmabb(i) = 0.0d0
      v2sigmaaa2(i) = 0.0d0
      v2sigmaab2(i) = 0.0d0
      v2sigmabb2(i) = 0.0d0
      v2sigmaaaab(i) = 0.0d0
      v2sigmaaabb(i) = 0.0d0
      v2sigmaabbb(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end
      
      
      subroutine rks_c_lyp
     & (ideriv,npt,rhoa1,sigmaaa1,
     &  zk,vrhoa,vsigmaaa,
     &  v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
c
c     C. Lee, W. Yang, and R.G. Parr
c     Development of the Colle-Salvetti correlation-energy formula into
c     a functional of the electron density
c     Phys. Rev. B37 (1988) 785-789
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit real*8 (a-h,o-z)
      integer ideriv,npt
      real*8 rhoa1(npt)
      real*8 sigmaaa1(npt)
      real*8 zk(npt),vrhoa(npt),vsigmaaa(npt)
      real*8 v2rhoa2(npt),v2rhoasigmaaa(npt),v2sigmaaa2(npt)
      parameter(tol=1.0d-20)
      
      if(ideriv.eq.0) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = rho**(1.D0/3.D0)
      t3 = 1/t2
      t6 = 1/(1.D0+0.349D0*t3)
      t9 = 0.2533D0*t3
      t10 = dexp(-t9)
      t12 = rho**2
      t14 = t2**2
      t20 = t3*t6
      zk(i) = -0.4918D-1*t6*rho-0.649176D-2*t10*t6/t14/t12/rho*(0.25D0
     #*t12*(0.1148493600075277D2*t14*t12+(0.2611111111111111D1
     #-0.9850555555555556D-1*t3-0.1357222222222222D0*t20)*sigma-0.5D0*
     #(0.25D1-0.1407222222222222D-1*t3-0.1938888888888889D-1*t20)
     #*sigma-0.2777777777777778D-1*(t9+0.349D0*t20-11.D0)*sigma)
     #-0.4583333333333333D0*t12*sigma)
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = rho**(1.D0/3.D0)
      t3 = 1/t2
      t5 = 1.D0+0.349D0*t3
      t6 = 1/t5
      t9 = 0.2533D0*t3
      t10 = dexp(-t9)
      t11 = t10*t6
      t12 = rho**2
      t14 = t2**2
      t16 = 1/t14/t12/rho
      t20 = t3*t6
      t22 = 0.2611111111111111D1-0.9850555555555556D-1*t3
     #-0.1357222222222222D0*t20
      t30 = t9+0.349D0*t20-11.D0
      t33 = 0.1148493600075277D2*t14*t12+t22*sigma-0.5D0*(0.25D1
     #-0.1407222222222222D-1*t3-0.1938888888888889D-1*t20)*sigma
     #-0.2777777777777778D-1*t30*sigma
      t38 = 0.25D0*t12*t33-0.4583333333333333D0*t12*sigma
      zk(i) = -0.4918D-1*t6*rho-0.649176D-2*t11*t16*t38
      t45 = t14*rho
      t49 = t30/rho*sigma
      t54 = rho*sigma
      t60 = t5**2
      t61 = 1/t60
      t64 = t12**2
      t66 = 1/t64/rho
      t81 = 1/t2/rho
      t83 = t81*t6
      t85 = 1/t45
      t86 = t85*t61
      vrhoa(i) = -0.4918D-1*t6-0.649176D-2*t11*t16*(0.5D0*rho*t33
     #+0.25D0*t12*(0.3062649600200738D2*t45-0.2777777777777778D-1*t49)
     #-0.25D0*t54)-0.5721273333333333D-2*t61*t3-0.548120936D-3*t66*t10
     #*t6*t38-0.75520808D-3*t10*t61*t66*t38+0.2380312D-1*t11/t14/t64
     #*t38-0.649176D-2*t11*t16*(0.25D0*t12*((0.3283518518518519D-1*t81
     #+0.4524074074074074D-1*t83-0.1578901851851852D-1*t86)*sigma
     #-0.5D0*(0.4690740740740741D-2*t81+0.6462962962962963D-2*t83
     #-0.2255574074074074D-2*t86)*sigma-0.2777777777777778D-1*(
     #-0.8443333333333333D-1*t81-0.1163333333333333D0*t83
     #+0.4060033333333333D-1*t86)*sigma+0.2777777777777778D-1*t49)
     #-0.6666666666666667D0*t54)
      vsigmaaa(i) = 0.7213066666666667D-3*t11*t85-0.2596704D-1*t11*t16
     #*(0.25D0*t12*t22-0.6666666666666667D0*t12)
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = rho**(1.D0/3.D0)
      t3 = 1/t2
      t5 = 1.D0+0.349D0*t3
      t6 = 1/t5
      t9 = 0.2533D0*t3
      t10 = dexp(-t9)
      t11 = t10*t6
      t12 = rho**2
      t13 = t12*rho
      t14 = t2**2
      t16 = 1/t14/t13
      t17 = t14*t12
      t18 = 0.1148493600075277D2*t17
      t20 = t3*t6
      t22 = 0.2611111111111111D1-0.9850555555555556D-1*t3
     #-0.1357222222222222D0*t20
      t23 = t22*sigma
      t28 = 0.5D0*(0.25D1-0.1407222222222222D-1*t3
     #-0.1938888888888889D-1*t20)*sigma
      t30 = t9+0.349D0*t20-11.D0
      t32 = 0.2777777777777778D-1*t30*sigma
      t33 = t18+t23-t28-t32
      t38 = 0.25D0*t12*t33-0.4583333333333333D0*t12*sigma
      zk(i) = -0.4918D-1*t6*rho-0.649176D-2*t11*t16*t38
      t45 = t14*rho
      t47 = 1/rho
      t48 = t30*t47
      t49 = t48*sigma
      t51 = 0.3062649600200738D2*t45-0.2777777777777778D-1*t49
      t54 = rho*sigma
      t56 = 0.5D0*rho*t33+0.25D0*t12*t51-0.25D0*t54
      t60 = t5**2
      t61 = 1/t60
      t64 = t12**2
      t65 = t64*rho
      t66 = 1/t65
      t67 = t66*t10
      t68 = t6*t38
      t71 = t10*t61
      t76 = 1/t14/t64
      t81 = 1/t2/rho
      t83 = t81*t6
      t85 = 1/t45
      t86 = t85*t61
      t88 = 0.3283518518518519D-1*t81+0.4524074074074074D-1*t83
     #-0.1578901851851852D-1*t86
      t99 = -0.8443333333333333D-1*t81-0.1163333333333333D0*t83
     #+0.4060033333333333D-1*t86
      t103 = t88*sigma-0.5D0*(0.4690740740740741D-2*t81
     #+0.6462962962962963D-2*t83-0.2255574074074074D-2*t86)*sigma
     #-0.2777777777777778D-1*t99*sigma+0.2777777777777778D-1*t49
      t107 = 0.25D0*t12*t103-0.6666666666666667D0*t54
      vrhoa(i) = -0.4918D-1*t6-0.649176D-2*t11*t16*t56
     #-0.5721273333333333D-2*t61*t3-0.548120936D-3*t67*t68
     #-0.75520808D-3*t71*t66*t38+0.2380312D-1*t11*t76*t38-0.649176D-2
     #*t11*t16*t107
      t116 = 0.25D0*t12*t22-0.6666666666666667D0*t12
      vsigmaaa(i) = 0.7213066666666667D-3*t11*t85-0.2596704D-1*t11*t16
     #*t116
      t121 = 1/t60/t5
      t126 = t64*t12
      t127 = 1/t126
      t137 = t99*t47*sigma
      t141 = t30/t12*sigma
      t150 = rho*t51
      t165 = 1/t2/t126
      t166 = t165*t10
      t188 = 1/t2/t12
      t190 = t188*t6
      t192 = 1/t17
      t193 = t192*t61
      t195 = 1/t13
      t196 = t195*t121
      s1 = -0.2662299191111111D-2*t121*t85-0.7628364444444444D-2*t61
     #*t81+0.9500762890666667D-2*t127*t10*t68+0.9521248D-1*t11*t76*t56
     #-0.2596704D-1*t11*t16*(0.5D0*rho*t103+0.25D0*t12*(
     #-0.2777777777777778D-1*t137+0.2777777777777778D-1*t141))
     #-0.649176D-2*t11*t16*(t18+t23-t28-t32+t150)+0.9521248D-1*t11*t76
     #*t107-0.2192483744D-2*t67*t6*t56-0.302083232D-2*t71*t66*t56
      s2 = s1-0.2550589422186667D-3*t166*t61*t38-0.302083232D-2*t71
     #*t66*t107-0.3514234932266667D-3*t10*t121*t165*t38
     #-0.2192483744D-2*t67*t6*t107
      v2rhoa2(i) = s2-0.9255935539253333D-4*t166*t68
     #-0.2221624533333333D0*t11/t14/t65*t38-0.1298352D-1*t11*t16*
     #(0.25D0*t12*((-0.4378024691358025D-1*t188-0.6032098765432099D-1
     #*t190+0.3157803703703704D-1*t193-0.3673578308641975D-2*t196)
     #*sigma-0.5D0*(-0.6254320987654321D-2*t188-0.8617283950617284D-2
     #*t190+0.4511148148148148D-2*t193-0.5247969012345679D-3*t196)
     #*sigma-0.2777777777777778D-1*(0.1125777777777778D0*t188
     #+0.1551111111111111D0*t190-0.8120066666666667D-1*t193
     #+0.9446344222222222D-2*t196)*sigma+0.5555555555555556D-1*t137
     #-0.5555555555555556D-1*t141)-0.6666666666666667D0*sigma)
     #-0.649176D-2*t11*t16*(0.1D1*t150+0.2552208000167282D2*t17-0.5D0
     #*sigma)+0.1309027338666667D-1*t71*t127*t38
      v2rhoasigmaaa(i) = -0.649176D-2*t11*t16*(-0.9444444444444444D0
     #*rho-0.2777777777777778D-1*rho*t30)+0.6090232622222222D-4*t195
     #*t10*t6+0.8391200888888889D-4*t71*t195+0.9978075555555556D-2*t11
     #*t192-0.1298352D-1*t11*t16*(0.25D0*t12*(-0.3D-22*t83+0.1D-22*t86
     #+0.5555555555555556D-1*t48)+0.1333333333333333D1*rho)
     #-0.1298352D-1*t11*t192*t22-0.2192483744D-2*t67*t6*t116
     #-0.302083232D-2*t71*t66*t116+0.9521248D-1*t11*t76*t116
     #-0.2596704D-1*t11*t16*(0.25D0*t12*t88-0.1333333333333333D1*rho)
      v2sigmaaa2(i) = 0.D0
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      v2rhoa2(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      v2rhoasigmaaa(i) = 0.0d0
      v2sigmaaa2(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end
      

      subroutine uks_c_pw91
     & (ideriv,npt,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,
     &  zk,vrhoa,vrhob,vsigmaaa,vsigmabb,vsigmaab,
     &  v2rhoa2,v2rhob2,v2rhoab,
     &  v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb,
     &  v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,
     &  v2sigmaaa2,v2sigmaaaab,v2sigmaaabb,
     &  v2sigmaab2,v2sigmaabbb,v2sigmabb2)
c
c     J.P. Perdew, J.A. Chevary, S.H. Vosko, K.A. Jackson, 
c     M.R. Pederson, D.J. Singh, C. Fiolhais
c     Atoms, molecules, solids and surfaces: 
c     Applications of the generalized gradient approximation 
c     for exchange and correlation
c     Phys. Rev. B 46 (1992) 6671--6687
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit real*8 (a-h,o-z)
      integer ideriv,npt
      real*8 rhoa1(npt),rhob1(npt)
      real*8 sigmaaa1(npt),sigmabb1(npt),sigmaab1(npt)
      real*8 zk(npt),vrhoa(npt),vrhob(npt)
      real*8 vsigmaaa(npt),vsigmabb(npt),vsigmaab(npt)
      real*8 v2rhoa2(npt),v2rhob2(npt),v2rhoab(npt)
      real*8 v2rhoasigmaaa(npt),v2rhoasigmaab(npt)
      real*8 v2rhoasigmabb(npt),v2rhobsigmabb(npt)
      real*8 v2rhobsigmaab(npt),v2rhobsigmaaa(npt)
      real*8 v2sigmaaa2(npt),v2sigmaaaab(npt),v2sigmaaabb(npt)
      real*8 v2sigmaab2(npt),v2sigmaabbb(npt),v2sigmabb2(npt)
      parameter(tol=1.0d-20)
      
      if (ideriv.eq.0) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmabb
      t2 = 1/rhob
      t3 = t2**(1.D0/3.D0)
      t6 = t2**(1.D0/6.D0)
      t9 = dsqrt(t2)
      t11 = t3**2
      t17 = dlog(1.D0+0.321646831778707D2/(0.1112037486309468D2*t6
     #+0.3844746237447211D1*t3+0.1644733775567609D1*t9
     #+0.2405871291288192D0*t11))
      t18 = (1.D0+0.1274696188700087D0*t3)*t17
      t20 = rhob**2
      t21 = rhob**(1.D0/3.D0)
      t23 = 1/t21/t20
      t27 = dexp(0.2513789108010112D1*t18)
      t28 = t27-1.D0
      t29 = 1/t28
      t30 = sigmabb**2
      t32 = t20**2
      t33 = t21**2
      t35 = 1/t33/t32
      t42 = t28**2
      t52 = dlog(1.D0+0.2697586091519874D1*(0.1007494971260294D0
     #*sigmabb*t23+0.2738174287780083D-1*t29*t30*t35)/(1.D0
     #+0.2717804421747983D0*t29*sigmabb*t23+0.7386460874872889D-1/t42
     #*t30*t35))
      t70 = dexp(-0.1645530784602056D1/t33/t20*sigmabb)
      zk(i) = rhob*(-0.3109D-1*t18+0.1236778371778789D-1*t52
     #+0.7937005259840997D0*(0.1D-2*(0.2568D1+0.1443307452126544D2*t3
     #+0.2843543831490386D-2*t11)/(1.D0+0.5411317332115466D1*t3
     #+0.1816419932959077D0*t11+0.1763993811759022D-1*t2)
     #-0.1853571428571429D-2)*sigmabb*t23*t70)
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = 1/rhoa
      t3 = t2**(1.D0/3.D0)
      t6 = t2**(1.D0/6.D0)
      t9 = dsqrt(t2)
      t11 = t3**2
      t17 = dlog(1.D0+0.321646831778707D2/(0.1112037486309468D2*t6
     #+0.3844746237447211D1*t3+0.1644733775567609D1*t9
     #+0.2405871291288192D0*t11))
      t18 = (1.D0+0.1274696188700087D0*t3)*t17
      t20 = rhoa**2
      t21 = rhoa**(1.D0/3.D0)
      t23 = 1/t21/t20
      t27 = dexp(0.2513789108010112D1*t18)
      t28 = t27-1.D0
      t29 = 1/t28
      t30 = sigmaaa**2
      t32 = t20**2
      t33 = t21**2
      t35 = 1/t33/t32
      t42 = t28**2
      t52 = dlog(1.D0+0.2697586091519874D1*(0.1007494971260294D0
     #*sigmaaa*t23+0.2738174287780083D-1*t29*t30*t35)/(1.D0
     #+0.2717804421747983D0*t29*sigmaaa*t23+0.7386460874872889D-1/t42
     #*t30*t35))
      t70 = dexp(-0.1645530784602056D1/t33/t20*sigmaaa)
      zk(i) = rhoa*(-0.3109D-1*t18+0.1236778371778789D-1*t52
     #+0.7937005259840997D0*(0.1D-2*(0.2568D1+0.1443307452126544D2*t3
     #+0.2843543831490386D-2*t11)/(1.D0+0.5411317332115466D1*t3
     #+0.1816419932959077D0*t11+0.1763993811759022D-1*t2)
     #-0.1853571428571429D-2)*sigmaaa*t23*t70)
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = 1/rho
      t5 = t4**(1.D0/3.D0)
      t8 = t4**(1.D0/6.D0)
      t11 = dsqrt(t4)
      t13 = t5**2
      t19 = dlog(1.D0+0.160818243221511D2/(0.598255043577108D1*t8
     #+0.2225569421150687D1*t5+0.8004286349993634D0*t11
     #+0.1897004325747559D0*t13))
      t21 = 0.62182D-1*(1.D0+0.1325688999052018D0*t5)*t19
      t32 = dlog(1.D0+0.2960857464321668D2/(0.8157414703487641D1*t8
     #+0.2247591863577616D1*t5+0.4300972471276643D0*t11
     #+0.1911512595127338D0*t13))
      t35 = rhoa-1.D0*rhob
      t36 = t35*t4
      t37 = 1.D0+t36
      t38 = t37**(1.D0/3.D0)
      t41 = 1.D0-1.D0*t36
      t42 = t41**(1.D0/3.D0)
      t44 = t38*t37+t42*t41-2.D0
      t45 = t35**2
      t46 = t45**2
      t47 = rho**2
      t48 = t47**2
      t50 = t46/t48
      t55 = 0.3799574853701528D-1*(1.D0+0.6901399211255825D-1*t5)*t32
     #*t44*(1.D0-1.D0*t50)
      t66 = dlog(1.D0+0.321646831778707D2/(0.1112037486309468D2*t8
     #+0.3844746237447211D1*t5+0.1644733775567609D1*t11
     #+0.2405871291288192D0*t13))
      t72 = 0.1923661050931536D1*(-0.3109D-1*(1.D0
     #+0.1274696188700087D0*t5)*t66+t21)*t44*t50
      t73 = t38**2
      t75 = t42**2
      t77 = 0.5D0*t73+0.5D0*t75
      t78 = t77**2
      t79 = t78*t77
      t80 = 1/t78
      t82 = rho**(1.D0/3.D0)
      t84 = 1/t82/t47
      t91 = dexp(-0.4042761511756372D2*(-t21+t55+t72)/t79)
      t92 = t91-1.D0
      t93 = 1/t92
      t94 = sigma**2
      t96 = t78**2
      t98 = t82**2
      t101 = 1/t96/t98/t48
      t109 = t92**2
      t119 = dlog(1.D0+0.2697586091519874D1*(0.634682060977037D-1
     #*sigma*t80*t84+0.1086645186223595D-1*t93*t94*t101)/(1.D0
     #+0.1712109500228824D0*t93*sigma*t80*t84+0.2931318940773793D-1
     #/t109*t94*t101))
      t140 = dexp(-0.261211729852336D1*t78/t98/t47*sigma)
      zk(i) = rho*(-t21+t55+t72+0.2473556743557577D-1*t79*t119+0.1D1*
     #(0.1D-2*(0.2568D1+0.1443307452126544D2*t5+0.2843543831490386D-2
     #*t13)/(1.D0+0.5411317332115466D1*t5+0.1816419932959077D0*t13
     #+0.1763993811759022D-1*t4)-0.1853571428571429D-2)*t77*sigma*t84
     #*t140)
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmabb
      t2 = 1/rhob
      t3 = t2**(1.D0/3.D0)
      t5 = 1.D0+0.1274696188700087D0*t3
      t6 = t2**(1.D0/6.D0)
      t9 = dsqrt(t2)
      t11 = t3**2
      t13 = 0.1112037486309468D2*t6+0.3844746237447211D1*t3
     #+0.1644733775567609D1*t9+0.2405871291288192D0*t11
      t16 = 1.D0+0.321646831778707D2/t13
      t17 = dlog(t16)
      t18 = t5*t17
      t19 = 0.3109D-1*t18
      t20 = rhob**2
      t21 = rhob**(1.D0/3.D0)
      t23 = 1/t21/t20
      t27 = dexp(0.2513789108010112D1*t18)
      t28 = t27-1.D0
      t29 = 1/t28
      t30 = sigmabb**2
      t31 = t29*t30
      t32 = t20**2
      t33 = t21**2
      t35 = 1/t33/t32
      t38 = 0.1007494971260294D0*sigmabb*t23+0.2738174287780083D-1*t31
     #*t35
      t39 = t29*sigmabb
      t42 = t28**2
      t43 = 1/t42
      t44 = t43*t30
      t47 = 1.D0+0.2717804421747983D0*t39*t23+0.7386460874872889D-1
     #*t44*t35
      t48 = 1/t47
      t51 = 1.D0+0.2697586091519874D1*t38*t48
      t52 = dlog(t51)
      t53 = 0.1236778371778789D-1*t52
      t56 = 0.2568D1+0.1443307452126544D2*t3+0.2843543831490386D-2*t11
      t60 = 1.D0+0.5411317332115466D1*t3+0.1816419932959077D0*t11
     #+0.1763993811759022D-1*t2
      t61 = 1/t60
      t64 = 0.1D-2*t56*t61-0.1853571428571429D-2
      t65 = t64*sigmabb
      t70 = dexp(-0.1645530784602056D1/t33/t20*sigmabb)
      t71 = t23*t70
      t73 = 0.7937005259840997D0*t65*t71
      zk(i) = rhob*(-t19+t53+t73)
      vrhoa(i) = 0.D0
      t76 = 1/t20
      t77 = 1/t11*t76
      t79 = 0.1321010150222857D-2*t77*t17
      t80 = t13**2
      t83 = t6**2
      t84 = t83**2
      t94 = 1/t3*t76
      t100 = 0.1D1*t5/t80*(-0.1853395810515781D1/t84/t6*t76
     #-0.128158207914907D1*t77-0.8223668877838045D0/t9*t76
     #-0.1603914194192128D0*t94)/t16
      t103 = 1/t21/t20/rhob
      t106 = t79+t100
      t108 = t35*t106*t27
      t111 = t32*rhob
      t113 = 1/t33/t111
      t119 = t47**2
      t121 = t38/t119
      t122 = t43*sigmabb
      t140 = 1/t51
      t148 = t60**2
      vrhob(i) = -t19+t53+t73+rhob*(t79+t100+0.1236778371778789D-1*
     #(0.2697586091519874D1*(-0.2350821599607352D0*sigmabb*t103
     #+0.2213957124623647D1*t44*t108-0.1277814667630705D0*t31*t113)
     #*t48-0.2697586091519874D1*t121*(0.2197487022544806D2*t122*t23
     #*t106*t27-0.6341543650745294D0*t39*t103+0.1194467989321217D2/t42
     #/t28*t30*t108-0.3447015074940681D0*t44*t113))*t140
     #+0.7937005259840997D0*(0.1D-2*(-0.4811024840421814D1*t77
     #-0.1895695887660258D-2*t94)*t61-0.1D-2*t56/t148*(
     #-0.1803772444038489D1*t77-0.1210946621972718D0*t94
     #-0.1763993811759022D-1*t76))*sigmabb*t71-0.1851967893962899D1
     #*t65*t103*t70+0.3482823064697813D1*t64*t30/t32/t20*t70)
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      vsigmabb(i) = rhob*(0.1236778371778789D-1*(0.2697586091519874D1*
     #(0.1007494971260294D0*t23+0.5476348575560166D-1*t39*t35)*t48
     #-0.2697586091519874D1*t121*(0.2717804421747983D0*t29*t23
     #+0.1477292174974578D0*t122*t35))*t140+0.7937005259840997D0*t64
     #*t23*t70-0.130605864926168D1*t65/t111*t70)
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = 1/rhoa
      t3 = t2**(1.D0/3.D0)
      t5 = 1.D0+0.1274696188700087D0*t3
      t6 = t2**(1.D0/6.D0)
      t9 = dsqrt(t2)
      t11 = t3**2
      t13 = 0.1112037486309468D2*t6+0.3844746237447211D1*t3
     #+0.1644733775567609D1*t9+0.2405871291288192D0*t11
      t16 = 1.D0+0.321646831778707D2/t13
      t17 = dlog(t16)
      t18 = t5*t17
      t19 = 0.3109D-1*t18
      t20 = rhoa**2
      t21 = rhoa**(1.D0/3.D0)
      t23 = 1/t21/t20
      t27 = dexp(0.2513789108010112D1*t18)
      t28 = t27-1.D0
      t29 = 1/t28
      t30 = sigmaaa**2
      t31 = t29*t30
      t32 = t20**2
      t33 = t21**2
      t35 = 1/t33/t32
      t38 = 0.1007494971260294D0*sigmaaa*t23+0.2738174287780083D-1*t31
     #*t35
      t39 = t29*sigmaaa
      t42 = t28**2
      t43 = 1/t42
      t44 = t43*t30
      t47 = 1.D0+0.2717804421747983D0*t39*t23+0.7386460874872889D-1
     #*t44*t35
      t48 = 1/t47
      t51 = 1.D0+0.2697586091519874D1*t38*t48
      t52 = dlog(t51)
      t53 = 0.1236778371778789D-1*t52
      t56 = 0.2568D1+0.1443307452126544D2*t3+0.2843543831490386D-2*t11
      t60 = 1.D0+0.5411317332115466D1*t3+0.1816419932959077D0*t11
     #+0.1763993811759022D-1*t2
      t61 = 1/t60
      t64 = 0.1D-2*t56*t61-0.1853571428571429D-2
      t65 = t64*sigmaaa
      t70 = dexp(-0.1645530784602056D1/t33/t20*sigmaaa)
      t71 = t23*t70
      t73 = 0.7937005259840997D0*t65*t71
      zk(i) = rhoa*(-t19+t53+t73)
      t76 = 1/t20
      t77 = 1/t11*t76
      t79 = 0.1321010150222857D-2*t77*t17
      t80 = t13**2
      t83 = t6**2
      t84 = t83**2
      t94 = 1/t3*t76
      t100 = 0.1D1*t5/t80*(-0.1853395810515781D1/t84/t6*t76
     #-0.128158207914907D1*t77-0.8223668877838045D0/t9*t76
     #-0.1603914194192128D0*t94)/t16
      t103 = 1/t21/t20/rhoa
      t106 = t79+t100
      t108 = t35*t106*t27
      t111 = t32*rhoa
      t113 = 1/t33/t111
      t119 = t47**2
      t121 = t38/t119
      t122 = t43*sigmaaa
      t140 = 1/t51
      t148 = t60**2
      vrhoa(i) = -t19+t53+t73+rhoa*(t79+t100+0.1236778371778789D-1*
     #(0.2697586091519874D1*(-0.2350821599607352D0*sigmaaa*t103
     #+0.2213957124623647D1*t44*t108-0.1277814667630705D0*t31*t113)
     #*t48-0.2697586091519874D1*t121*(0.2197487022544806D2*t122*t23
     #*t106*t27-0.6341543650745294D0*t39*t103+0.1194467989321217D2/t42
     #/t28*t30*t108-0.3447015074940681D0*t44*t113))*t140
     #+0.7937005259840997D0*(0.1D-2*(-0.4811024840421814D1*t77
     #-0.1895695887660258D-2*t94)*t61-0.1D-2*t56/t148*(
     #-0.1803772444038489D1*t77-0.1210946621972718D0*t94
     #-0.1763993811759022D-1*t76))*sigmaaa*t71-0.1851967893962899D1
     #*t65*t103*t70+0.3482823064697813D1*t64*t30/t32/t20*t70)
      vrhob(i) = 0.D0
      vsigmaaa(i) = rhoa*(0.1236778371778789D-1*(0.2697586091519874D1*
     #(0.1007494971260294D0*t23+0.5476348575560166D-1*t39*t35)*t48
     #-0.2697586091519874D1*t121*(0.2717804421747983D0*t29*t23
     #+0.1477292174974578D0*t122*t35))*t140+0.7937005259840997D0*t64
     #*t23*t70-0.130605864926168D1*t65/t111*t70)
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = 1/rho
      t5 = t4**(1.D0/3.D0)
      t7 = 1.D0+0.1325688999052018D0*t5
      t8 = t4**(1.D0/6.D0)
      t11 = dsqrt(t4)
      t13 = t5**2
      t15 = 0.598255043577108D1*t8+0.2225569421150687D1*t5
     #+0.8004286349993634D0*t11+0.1897004325747559D0*t13
      t18 = 1.D0+0.160818243221511D2/t15
      t19 = dlog(t18)
      t21 = 0.62182D-1*t7*t19
      t23 = 1.D0+0.6901399211255825D-1*t5
      t28 = 0.8157414703487641D1*t8+0.2247591863577616D1*t5
     #+0.4300972471276643D0*t11+0.1911512595127338D0*t13
      t31 = 1.D0+0.2960857464321668D2/t28
      t32 = dlog(t31)
      t33 = t23*t32
      t35 = rhoa-1.D0*rhob
      t36 = t35*t4
      t37 = 1.D0+t36
      t38 = t37**(1.D0/3.D0)
      t41 = 1.D0-1.D0*t36
      t42 = t41**(1.D0/3.D0)
      t44 = t38*t37+t42*t41-2.D0
      t45 = t35**2
      t46 = t45**2
      t47 = rho**2
      t48 = t47**2
      t49 = 1/t48
      t50 = t46*t49
      t52 = 1.D0-1.D0*t50
      t55 = 0.3799574853701528D-1*t33*t44*t52
      t57 = 1.D0+0.1274696188700087D0*t5
      t62 = 0.1112037486309468D2*t8+0.3844746237447211D1*t5
     #+0.1644733775567609D1*t11+0.2405871291288192D0*t13
      t65 = 1.D0+0.321646831778707D2/t62
      t66 = dlog(t65)
      t69 = -0.3109D-1*t57*t66+t21
      t70 = t69*t44
      t72 = 0.1923661050931536D1*t70*t50
      t73 = t38**2
      t75 = t42**2
      t77 = 0.5D0*t73+0.5D0*t75
      t78 = t77**2
      t79 = t78*t77
      t80 = 1/t78
      t81 = sigma*t80
      t82 = rho**(1.D0/3.D0)
      t84 = 1/t82/t47
      t87 = -t21+t55+t72
      t88 = 1/t79
      t91 = dexp(-0.4042761511756372D2*t87*t88)
      t92 = t91-1.D0
      t93 = 1/t92
      t94 = sigma**2
      t95 = t93*t94
      t96 = t78**2
      t97 = 1/t96
      t98 = t82**2
      t100 = 1/t98/t48
      t101 = t97*t100
      t104 = 0.634682060977037D-1*t81*t84+0.1086645186223595D-1*t95*t101
      t105 = t93*sigma
      t106 = t80*t84
      t109 = t92**2
      t110 = 1/t109
      t111 = t110*t94
      t114 = 1.D0+0.1712109500228824D0*t105*t106+0.2931318940773793D-1
     #*t111*t101
      t115 = 1/t114
      t118 = 1.D0+0.2697586091519874D1*t104*t115
      t119 = dlog(t118)
      t121 = 0.2473556743557577D-1*t79*t119
      t124 = 0.2568D1+0.1443307452126544D2*t5+0.2843543831490386D-2*t13
      t128 = 1.D0+0.5411317332115466D1*t5+0.1816419932959077D0*t13
     #+0.1763993811759022D-1*t4
      t129 = 1/t128
      t132 = 0.1D-2*t124*t129-0.1853571428571429D-2
      t133 = t132*t77
      t136 = 1/t98/t47
      t140 = dexp(-0.261211729852336D1*t78*t136*sigma)
      t141 = sigma*t84*t140
      t143 = 0.1D1*t133*t141
      zk(i) = rho*(-t21+t55+t72+t121+t143)
      t149 = 0.1333333333333333D1*t38*t4-0.1333333333333333D1*t42*t4
      t152 = 0.3799574853701528D-1*t33*t149*t52
      t153 = t45*t35
      t156 = t33*t44*t153*t49
      t157 = 0.1519829941480611D0*t156
      t160 = 0.1923661050931536D1*t69*t149*t50
      t162 = t70*t153*t49
      t163 = 0.7694644203726145D1*t162
      t164 = t78*t119
      t165 = 1/t38
      t168 = 1/t42
      t171 = 0.3333333333333333D0*t165*t4-0.3333333333333333D0*t168*t4
      t174 = sigma*t88
      t178 = t111*t97
      t182 = t87*t97
      t185 = -0.4042761511756372D2*(t152-t157+t160+t163)*t88
     #+0.1212828453526912D3*t182*t171
      t187 = t100*t185*t91
      t192 = 1/t96/t77*t100
      t193 = t192*t171
      t199 = t114**2
      t201 = t104/t199
      t202 = t110*sigma
      t203 = t202*t80
      t208 = t88*t84
      t215 = 1/t109/t92*t94*t97
      t225 = 1/t118
      t232 = t132*t78*t94
      t233 = t48*rho
      t234 = 1/t233
      t242 = 1/t47
      t243 = 1/t13*t242
      t244 = t243*t19
      t245 = 0.2747799777968419D-2*t244
      t246 = t15**2
      t249 = t8**2
      t250 = t249**2
      t253 = 1/t250/t8*t242
      t257 = 1/t11*t242
      t260 = 1/t5*t242
      t265 = t7/t246*(-0.99709173929518D0*t253-0.7418564737168958D0
     #*t243-0.4002143174996817D0*t257-0.1264669550498372D0*t260)/t18
      t266 = 0.1D1*t265
      t270 = 0.8740794299481065D-3*t243*t32*t44*t52
      t271 = t28**2
      t284 = 0.1124999956683108D1*t23/t271*(-0.135956911724794D1*t253
     #-0.7491972878592054D0*t243-0.2150486235638321D0*t257
     #-0.1274341730084892D0*t260)/t31*t44*t52
      t291 = -0.1333333333333333D1*t38*t35*t242+0.1333333333333333D1
     #*t42*t35*t242
      t294 = 0.3799574853701528D-1*t33*t291*t52
      t298 = 0.1519829941480611D0*t33*t44*t46*t234
      t301 = t62**2
      t318 = 0.1923661050931536D1*(0.1321010150222857D-2*t243*t66
     #+0.1D1*t57/t301*(-0.1853395810515781D1*t253-0.128158207914907D1
     #*t243-0.8223668877838045D0*t257-0.1603914194192128D0*t260)/t65
     #-0.2747799777968419D-2*t244-0.1D1*t265)*t44*t50
      t321 = 0.1923661050931536D1*t69*t291*t50
      t324 = 0.7694644203726145D1*t70*t46*t234
      t331 = -0.3333333333333333D0*t165*t35*t242+0.3333333333333333D0
     #*t168*t35*t242
      t337 = t47*rho
      t339 = 1/t82/t337
      t347 = -0.4042761511756372D2*(t245+t266-t270-t284+t294+t298+t318
     #+t321-t324)*t88+0.1212828453526912D3*t182*t331
      t349 = t100*t347*t91
      t352 = t192*t331
      t357 = t97/t98/t233
      t391 = t128**2
      s1 = t245+t266-t270-t284+t294+t298+t318
      s2 = s1+t321-t324+0.7420670230672731D-1*t164*t331
      s3 = s2+0.2473556743557577D-1*t79*(0.2697586091519874D1*(
     #-0.1269364121954074D0*t174*t84*t331-0.148092480894642D0*t81*t339
     #-0.1086645186223595D-1*t178*t349-0.4346580744894379D-1*t95*t352
     #-0.5071010869043442D-1*t95*t357)*t115-0.2697586091519874D1*t201*
     #(-0.1712109500228824D0*t203*t84*t347*t91-0.3424219000457648D0
     #*t105*t208*t331-0.3994922167200589D0*t105*t80*t339
     #-0.5862637881547585D-1*t215*t349-0.1172527576309517D0*t111*t352
     #-0.136794883902777D0*t111*t357))*t225+0.1D1*(0.1D-2*(
     #-0.4811024840421814D1*t243-0.1895695887660258D-2*t260)*t129
     #-0.1D-2*t124/t391*(-0.1803772444038489D1*t243
     #-0.1210946621972718D0*t260-0.1763993811759022D-1*t242))*t77*t141
      t426 = s3+0.1D1*t132*t331*t141-0.2333333333333333D1*t133*sigma
     #*t339*t140+0.1D1*t133*sigma*t84*(-0.522423459704672D1*t77*t136
     #*sigma*t331+0.6965646129395627D1*t78/t98/t337*sigma)*t140
      t427 = rho*t426
      vrhoa(i) = rho*(t152-t157+t160+t163+0.7420670230672731D-1*t164
     #*t171+0.2473556743557577D-1*t79*(0.2697586091519874D1*(
     #-0.1269364121954074D0*t174*t84*t171-0.1086645186223595D-1*t178
     #*t187-0.4346580744894379D-1*t95*t193)*t115-0.2697586091519874D1
     #*t201*(-0.1712109500228824D0*t203*t84*t185*t91
     #-0.3424219000457648D0*t105*t208*t171-0.5862637881547585D-1*t215
     #*t187-0.1172527576309517D0*t111*t193))*t225+0.1D1*t132*t171*t141
     #-0.522423459704672D1*t232*t234*t171*t140)-t21+t55+t72+t121+t143
     #+t427
      t428 = -t149
      t431 = 0.3799574853701528D-1*t33*t428*t52
      t432 = 0.1519829941480611D0*t156
      t435 = 0.1923661050931536D1*t69*t428*t50
      t436 = 0.7694644203726145D1*t162
      t437 = -t171
      t448 = -0.4042761511756372D2*(t431+t432+t435-t436)*t88
     #+0.1212828453526912D3*t182*t437
      t450 = t100*t448*t91
      t453 = t192*t437
      vrhob(i) = rho*(t431+t432+t435-t436+0.7420670230672731D-1*t164
     #*t437+0.2473556743557577D-1*t79*(0.2697586091519874D1*(
     #-0.1269364121954074D0*t174*t84*t437-0.1086645186223595D-1*t178
     #*t450-0.4346580744894379D-1*t95*t453)*t115-0.2697586091519874D1
     #*t201*(-0.1712109500228824D0*t203*t84*t448*t91
     #-0.3424219000457648D0*t105*t208*t437-0.5862637881547585D-1*t215
     #*t450-0.1172527576309517D0*t111*t453))*t225+0.1D1*t132*t437*t141
     #-0.522423459704672D1*t232*t234*t437*t140)-t21+t55+t72+t121+t143
     #+t427
      vsigmaaa(i) = rho*(0.2473556743557577D-1*t79*
     #(0.2697586091519874D1*(0.634682060977037D-1*t106
     #+0.2173290372447189D-1*t105*t101)*t115-0.2697586091519874D1*t201
     #*(0.1712109500228824D0*t93*t80*t84+0.5862637881547585D-1*t202
     #*t101))*t225+0.1D1*t133*t84*t140-0.261211729852336D1*t132*t79
     #*sigma*t234*t140)
      vsigmaab(i) = 2.D0*vsigmaaa(i)
      vsigmabb(i) = vsigmaaa(i)
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vrhob(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      vsigmaab(i) = 0.0d0
      vsigmabb(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmabb
      t2 = 1/rhob
      t3 = t2**(1.D0/3.D0)
      t5 = 1.D0+0.1274696188700087D0*t3
      t6 = t2**(1.D0/6.D0)
      t9 = dsqrt(t2)
      t11 = t3**2
      t13 = 0.1112037486309468D2*t6+0.3844746237447211D1*t3
     #+0.1644733775567609D1*t9+0.2405871291288192D0*t11
      t16 = 1.D0+0.321646831778707D2/t13
      t17 = dlog(t16)
      t18 = t5*t17
      t19 = 0.3109D-1*t18
      t20 = rhob**2
      t21 = rhob**(1.D0/3.D0)
      t23 = 1/t21/t20
      t27 = dexp(0.2513789108010112D1*t18)
      t28 = t27-1.D0
      t29 = 1/t28
      t30 = sigmabb**2
      t31 = t29*t30
      t32 = t20**2
      t33 = t21**2
      t35 = 1/t33/t32
      t38 = 0.1007494971260294D0*sigmabb*t23+0.2738174287780083D-1*t31
     #*t35
      t39 = t29*sigmabb
      t42 = t28**2
      t43 = 1/t42
      t44 = t43*t30
      t47 = 1.D0+0.2717804421747983D0*t39*t23+0.7386460874872889D-1
     #*t44*t35
      t48 = 1/t47
      t51 = 1.D0+0.2697586091519874D1*t38*t48
      t52 = dlog(t51)
      t53 = 0.1236778371778789D-1*t52
      t56 = 0.2568D1+0.1443307452126544D2*t3+0.2843543831490386D-2*t11
      t60 = 1.D0+0.5411317332115466D1*t3+0.1816419932959077D0*t11
     #+0.1763993811759022D-1*t2
      t61 = 1/t60
      t64 = 0.1D-2*t56*t61-0.1853571428571429D-2
      t65 = t64*sigmabb
      t70 = dexp(-0.1645530784602056D1/t33/t20*sigmabb)
      t71 = t23*t70
      t73 = 0.7937005259840997D0*t65*t71
      zk(i) = rhob*(-t19+t53+t73)
      vrhoa(i) = 0.D0
      t75 = 1/t11
      t76 = 1/t20
      t77 = t75*t76
      t78 = t77*t17
      t79 = 0.1321010150222857D-2*t78
      t80 = t13**2
      t81 = 1/t80
      t82 = t5*t81
      t83 = t6**2
      t84 = t83**2
      t85 = t84*t6
      t86 = 1/t85
      t90 = 1/t9
      t93 = 1/t3
      t94 = t93*t76
      t96 = -0.1853395810515781D1*t86*t76-0.128158207914907D1*t77
     #-0.8223668877838045D0*t90*t76-0.1603914194192128D0*t94
      t97 = 1/t16
      t99 = t82*t96*t97
      t100 = 0.1D1*t99
      t101 = t20*rhob
      t103 = 1/t21/t101
      t106 = t79+t100
      t108 = t35*t106*t27
      t111 = t32*rhob
      t113 = 1/t33/t111
      t116 = -0.2350821599607352D0*sigmabb*t103+0.2213957124623647D1
     #*t44*t108-0.1277814667630705D0*t31*t113
      t119 = t47**2
      t120 = 1/t119
      t121 = t38*t120
      t122 = t43*sigmabb
      t130 = 1/t42/t28
      t131 = t130*t30
      t136 = 0.2197487022544806D2*t122*t23*t106*t27
     #-0.6341543650745294D0*t39*t103+0.1194467989321217D2*t131*t108
     #-0.3447015074940681D0*t44*t113
      t139 = 0.2697586091519874D1*t116*t48-0.2697586091519874D1*t121
     #*t136
      t140 = 1/t51
      t141 = t139*t140
      t145 = -0.4811024840421814D1*t77-0.1895695887660258D-2*t94
      t148 = t60**2
      t149 = 1/t148
      t150 = t56*t149
      t154 = -0.1803772444038489D1*t77-0.1210946621972718D0*t94
     #-0.1763993811759022D-1*t76
      t157 = 0.1D-2*t145*t61-0.1D-2*t150*t154
      t158 = t157*sigmabb
      t159 = t158*t71
      t161 = t103*t70
      t162 = t65*t161
      t164 = t64*t30
      t165 = t32*t20
      t167 = 1/t165*t70
      t168 = t164*t167
      vrhob(i) = -t19+t53+t73+rhob*(t79+t100+0.1236778371778789D-1
     #*t141+0.7937005259840997D0*t159-0.1851967893962899D1*t162
     #+0.3482823064697813D1*t168)
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      t175 = 0.1007494971260294D0*t23+0.5476348575560166D-1*t39*t35
      t182 = 0.2717804421747983D0*t29*t23+0.1477292174974578D0*t122*t35
      t185 = 0.2697586091519874D1*t175*t48-0.2697586091519874D1*t121
     #*t182
      t191 = 1/t111
      vsigmabb(i) = rhob*(0.1236778371778789D-1*t185*t140
     #+0.7937005259840997D0*t64*t23*t70-0.130605864926168D1*t65*t191
     #*t70)
      v2rhoa2(i) = 0.D0
      v2rhoab(i) = 0.D0
      t204 = dlog(1.D0+0.160818243221511D2/(0.598255043577108D1*t6
     #+0.2225569421150687D1*t3+0.8004286349993634D0*t9
     #+0.1897004325747559D0*t11))
      t215 = 1/t32
      t216 = 1/t11/t2*t215
      t218 = 1/t101
      t219 = t75*t218
      t223 = 1/t3/t2*t215
      t225 = t93*t218
      t236 = t154**2
      t253 = t32**2
      t260 = t32*t101
      t266 = 1/t21/t32
      t275 = t139**2
      t276 = t51**2
      t277 = 1/t276
      t282 = t106**2
      t283 = t35*t282
      t285 = dexp(0.5027578216020225D1*t18)
      t286 = t283*t285
      t290 = t113*t106*t27
      t294 = 0.1D-22*t219*t204
      t296 = 0.8806734334819047D-3*t216*t17
      t298 = 0.2642020300445714D-2*t219*t17
      t302 = 0.8497974591333914D-1*t77*t81*t96*t97
      t306 = t96**2
      t309 = 0.2D1*t5/t80/t13*t306*t97
      t329 = 0.1D1*t82*(-0.1544496508763151D1/t85/t2*t215
     #+0.3706791621031562D1*t86*t218-0.854388052766047D0*t216
     #+0.2563164158298141D1*t219-0.4111834438919023D0/t9/t2*t215
     #+0.1644733775567609D1*t90*t218-0.5346380647307093D-1*t223
     #+0.3207828388384256D0*t225)*t97
      t330 = t80**2
      t333 = t16**2
      t337 = 0.321646831778707D2*t5/t330*t306/t333
      t338 = -t294+t296-t298-t302-t309+t329+t337
      t340 = t35*t338*t27
      t343 = t283*t27
      t347 = 1/t33/t165
      t358 = t38/t119/t47
      t359 = t136**2
      t363 = t23*t282
      t380 = t42**2
      s1 = 0.7937005259840997D0*(0.1D-2*(-0.3207349893614542D1*t216
     #+0.9622049680843627D1*t219-0.6318986292200858D-3*t223
     #+0.3791391775320515D-2*t225)*t61-0.2D-2*t145*t149*t154+0.2D-2
     #*t56/t148/t60*t236-0.1D-2*t150*(-0.1202514962692326D1*t216
     #+0.3607544888076978D1*t219-0.4036488739909061D-1*t223
     #+0.2421893243945437D0*t225+0.3527987623518044D-1*t218))*sigmabb
     #*t71+0.1528291352075288D2*t64*t30*sigmabb/t33/t253/rhob*t70
     #-0.2902352553914844D2*t164/t260*t70+0.6173226313209665D1*t65
     #*t266*t70-0.3703935787925799D1*t158*t161+0.6965646129395627D1
     #*t157*t30*t167-0.1236778371778789D-1*t275*t277
      t399 = s1+0.1236778371778789D-1*(0.2697586091519874D1*
     #(0.7836071998691172D0*sigmabb*t266+0.3580200260842915D3*t131
     #*t286-0.2066359982982071D2*t44*t290+0.2213957124623647D1*t44
     #*t340-0.1790100130421458D3*t44*t343+0.7240949783240664D0*t31
     #*t347)*t48-0.5395172183039748D1*t116*t120*t136
     #+0.5395172183039748D1*t358*t359-0.2697586091519874D1*t121*
     #(0.35535663829313D4*t130*sigmabb*t363*t285-0.1025493943854243D3
     #*t122*t103*t106*t27+0.2197487022544806D2*t122*t23*t338*t27
     #-0.177678319146565D4*t122*t363*t27+0.2113847883581765D1*t39*t266
     #+0.2897369528551702D4/t380*t30*t286-0.1114836790033136D3*t131
     #*t290+0.1194467989321217D2*t131*t340-0.9657898428505673D3*t131
     #*t343+0.1953308542466386D1*t44*t347))*t140+t337+t329-t294+t296
     #-t298-t309-t302
      v2rhob2(i) = -0.1D-22*t77*t204+0.2642020300445714D-2*t78+0.2D1
     #*t99+0.2473556743557577D-1*t141+0.1587401051968199D1*t159
     #-0.3703935787925799D1*t162+0.6965646129395627D1*t168+rhob*t399
      v2sigmaaa2(i) = 0.D0
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      t407 = t182**2
      t416 = t185**2
      v2sigmabb2(i) = rhob*(0.1236778371778789D-1*
     #(0.1477292174974578D0*t29*t35*t48-0.5395172183039748D1*t175*t120
     #*t182+0.5395172183039748D1*t358*t407-0.3985122824322565D0*t121
     #*t43*t35)*t140-0.1236778371778789D-1*t416*t277
     #-0.261211729852336D1*t64*t191*t70+0.2149159713855873D1*t65/t33
     #/t260*t70)
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = 1/rhoa
      t3 = t2**(1.D0/3.D0)
      t5 = 1.D0+0.1274696188700087D0*t3
      t6 = t2**(1.D0/6.D0)
      t9 = dsqrt(t2)
      t11 = t3**2
      t13 = 0.1112037486309468D2*t6+0.3844746237447211D1*t3
     #+0.1644733775567609D1*t9+0.2405871291288192D0*t11
      t16 = 1.D0+0.321646831778707D2/t13
      t17 = dlog(t16)
      t18 = t5*t17
      t19 = 0.3109D-1*t18
      t20 = rhoa**2
      t21 = rhoa**(1.D0/3.D0)
      t23 = 1/t21/t20
      t27 = dexp(0.2513789108010112D1*t18)
      t28 = t27-1.D0
      t29 = 1/t28
      t30 = sigmaaa**2
      t31 = t29*t30
      t32 = t20**2
      t33 = t21**2
      t35 = 1/t33/t32
      t38 = 0.1007494971260294D0*sigmaaa*t23+0.2738174287780083D-1*t31
     #*t35
      t39 = t29*sigmaaa
      t42 = t28**2
      t43 = 1/t42
      t44 = t43*t30
      t47 = 1.D0+0.2717804421747983D0*t39*t23+0.7386460874872889D-1
     #*t44*t35
      t48 = 1/t47
      t51 = 1.D0+0.2697586091519874D1*t38*t48
      t52 = dlog(t51)
      t53 = 0.1236778371778789D-1*t52
      t56 = 0.2568D1+0.1443307452126544D2*t3+0.2843543831490386D-2*t11
      t60 = 1.D0+0.5411317332115466D1*t3+0.1816419932959077D0*t11
     #+0.1763993811759022D-1*t2
      t61 = 1/t60
      t64 = 0.1D-2*t56*t61-0.1853571428571429D-2
      t65 = t64*sigmaaa
      t70 = dexp(-0.1645530784602056D1/t33/t20*sigmaaa)
      t71 = t23*t70
      t73 = 0.7937005259840997D0*t65*t71
      zk(i) = rhoa*(-t19+t53+t73)
      t75 = 1/t11
      t76 = 1/t20
      t77 = t75*t76
      t78 = t77*t17
      t79 = 0.1321010150222857D-2*t78
      t80 = t13**2
      t81 = 1/t80
      t82 = t5*t81
      t83 = t6**2
      t84 = t83**2
      t85 = t84*t6
      t86 = 1/t85
      t90 = 1/t9
      t93 = 1/t3
      t94 = t93*t76
      t96 = -0.1853395810515781D1*t86*t76-0.128158207914907D1*t77
     #-0.8223668877838045D0*t90*t76-0.1603914194192128D0*t94
      t97 = 1/t16
      t99 = t82*t96*t97
      t100 = 0.1D1*t99
      t101 = t20*rhoa
      t103 = 1/t21/t101
      t106 = t79+t100
      t108 = t35*t106*t27
      t111 = t32*rhoa
      t113 = 1/t33/t111
      t116 = -0.2350821599607352D0*sigmaaa*t103+0.2213957124623647D1
     #*t44*t108-0.1277814667630705D0*t31*t113
      t119 = t47**2
      t120 = 1/t119
      t121 = t38*t120
      t122 = t43*sigmaaa
      t130 = 1/t42/t28
      t131 = t130*t30
      t136 = 0.2197487022544806D2*t122*t23*t106*t27
     #-0.6341543650745294D0*t39*t103+0.1194467989321217D2*t131*t108
     #-0.3447015074940681D0*t44*t113
      t139 = 0.2697586091519874D1*t116*t48-0.2697586091519874D1*t121
     #*t136
      t140 = 1/t51
      t141 = t139*t140
      t145 = -0.4811024840421814D1*t77-0.1895695887660258D-2*t94
      t148 = t60**2
      t149 = 1/t148
      t150 = t56*t149
      t154 = -0.1803772444038489D1*t77-0.1210946621972718D0*t94
     #-0.1763993811759022D-1*t76
      t157 = 0.1D-2*t145*t61-0.1D-2*t150*t154
      t158 = t157*sigmaaa
      t159 = t158*t71
      t161 = t103*t70
      t162 = t65*t161
      t164 = t64*t30
      t165 = t32*t20
      t167 = 1/t165*t70
      t168 = t164*t167
      vrhoa(i) = -t19+t53+t73+rhoa*(t79+t100+0.1236778371778789D-1
     #*t141+0.7937005259840997D0*t159-0.1851967893962899D1*t162
     #+0.3482823064697813D1*t168)
      vrhob(i) = 0.D0
      t175 = 0.1007494971260294D0*t23+0.5476348575560166D-1*t39*t35
      t182 = 0.2717804421747983D0*t29*t23+0.1477292174974578D0*t122*t35
      t185 = 0.2697586091519874D1*t175*t48-0.2697586091519874D1*t121
     #*t182
      t191 = 1/t111
      vsigmaaa(i) = rhoa*(0.1236778371778789D-1*t185*t140
     #+0.7937005259840997D0*t64*t23*t70-0.130605864926168D1*t65*t191
     #*t70)
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      t204 = dlog(1.D0+0.160818243221511D2/(0.598255043577108D1*t6
     #+0.2225569421150687D1*t3+0.8004286349993634D0*t9
     #+0.1897004325747559D0*t11))
      t215 = 1/t32
      t216 = 1/t11/t2*t215
      t218 = 1/t101
      t219 = t75*t218
      t223 = 1/t3/t2*t215
      t225 = t93*t218
      t236 = t154**2
      t253 = t32**2
      t260 = t32*t101
      t266 = 1/t21/t32
      t275 = t139**2
      t276 = t51**2
      t277 = 1/t276
      t282 = t106**2
      t283 = t35*t282
      t285 = dexp(0.5027578216020225D1*t18)
      t286 = t283*t285
      t290 = t113*t106*t27
      t294 = 0.1D-22*t219*t204
      t296 = 0.8806734334819047D-3*t216*t17
      t298 = 0.2642020300445714D-2*t219*t17
      t302 = 0.8497974591333914D-1*t77*t81*t96*t97
      t306 = t96**2
      t309 = 0.2D1*t5/t80/t13*t306*t97
      t329 = 0.1D1*t82*(-0.1544496508763151D1/t85/t2*t215
     #+0.3706791621031562D1*t86*t218-0.854388052766047D0*t216
     #+0.2563164158298141D1*t219-0.4111834438919023D0/t9/t2*t215
     #+0.1644733775567609D1*t90*t218-0.5346380647307093D-1*t223
     #+0.3207828388384256D0*t225)*t97
      t330 = t80**2
      t333 = t16**2
      t337 = 0.321646831778707D2*t5/t330*t306/t333
      t338 = -t294+t296-t298-t302-t309+t329+t337
      t340 = t35*t338*t27
      t343 = t283*t27
      t347 = 1/t33/t165
      t358 = t38/t119/t47
      t359 = t136**2
      t363 = t23*t282
      t380 = t42**2
      s1 = 0.7937005259840997D0*(0.1D-2*(-0.3207349893614542D1*t216
     #+0.9622049680843627D1*t219-0.6318986292200858D-3*t223
     #+0.3791391775320515D-2*t225)*t61-0.2D-2*t145*t149*t154+0.2D-2
     #*t56/t148/t60*t236-0.1D-2*t150*(-0.1202514962692326D1*t216
     #+0.3607544888076978D1*t219-0.4036488739909061D-1*t223
     #+0.2421893243945437D0*t225+0.3527987623518044D-1*t218))*sigmaaa
     #*t71+0.1528291352075288D2*t64*t30*sigmaaa/t33/t253/rhoa*t70
     #-0.2902352553914844D2*t164/t260*t70+0.6173226313209665D1*t65
     #*t266*t70-0.3703935787925799D1*t158*t161+0.6965646129395627D1
     #*t157*t30*t167-0.1236778371778789D-1*t275*t277
      t399 = s1+0.1236778371778789D-1*(0.2697586091519874D1*
     #(0.7836071998691172D0*sigmaaa*t266+0.3580200260842915D3*t131
     #*t286-0.2066359982982071D2*t44*t290+0.2213957124623647D1*t44
     #*t340-0.1790100130421458D3*t44*t343+0.7240949783240664D0*t31
     #*t347)*t48-0.5395172183039748D1*t116*t120*t136
     #+0.5395172183039748D1*t358*t359-0.2697586091519874D1*t121*
     #(0.35535663829313D4*t130*sigmaaa*t363*t285-0.1025493943854243D3
     #*t122*t103*t106*t27+0.2197487022544806D2*t122*t23*t338*t27
     #-0.177678319146565D4*t122*t363*t27+0.2113847883581765D1*t39*t266
     #+0.2897369528551702D4/t380*t30*t286-0.1114836790033136D3*t131
     #*t290+0.1194467989321217D2*t131*t340-0.9657898428505673D3*t131
     #*t343+0.1953308542466386D1*t44*t347))*t140+t337+t329-t294+t296
     #-t298-t309-t302
      v2rhoa2(i) = -0.1D-22*t77*t204+0.2642020300445714D-2*t78+0.2D1
     #*t99+0.2473556743557577D-1*t141+0.1587401051968199D1*t159
     #-0.3703935787925799D1*t162+0.6965646129395627D1*t168+rhoa*t399
      v2rhob2(i) = 0.D0
      v2rhoab(i) = 0.D0
      t407 = t182**2
      t416 = t185**2
      v2sigmaaa2(i) = rhoa*(0.1236778371778789D-1*
     #(0.1477292174974578D0*t29*t35*t48-0.5395172183039748D1*t175*t120
     #*t182+0.5395172183039748D1*t358*t407-0.3985122824322565D0*t121
     #*t43*t35)*t140-0.1236778371778789D-1*t416*t277
     #-0.261211729852336D1*t64*t191*t70+0.2149159713855873D1*t65/t33
     #/t260*t70)
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      v2sigmabb2(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = 1/rho
      t5 = t4**(1.D0/3.D0)
      t7 = 1.D0+0.1325688999052018D0*t5
      t8 = t4**(1.D0/6.D0)
      t11 = dsqrt(t4)
      t13 = t5**2
      t15 = 0.598255043577108D1*t8+0.2225569421150687D1*t5
     #+0.8004286349993634D0*t11+0.1897004325747559D0*t13
      t18 = 1.D0+0.160818243221511D2/t15
      t19 = dlog(t18)
      t21 = 0.62182D-1*t7*t19
      t23 = 1.D0+0.6901399211255825D-1*t5
      t28 = 0.8157414703487641D1*t8+0.2247591863577616D1*t5
     #+0.4300972471276643D0*t11+0.1911512595127338D0*t13
      t31 = 1.D0+0.2960857464321668D2/t28
      t32 = dlog(t31)
      t33 = t23*t32
      t35 = rhoa-1.D0*rhob
      t36 = t35*t4
      t37 = 1.D0+t36
      t38 = t37**(1.D0/3.D0)
      t39 = t38*t37
      t41 = 1.D0-1.D0*t36
      t42 = t41**(1.D0/3.D0)
      t43 = t42*t41
      t44 = t39+t43-2.D0
      t45 = t35**2
      t46 = t45**2
      t47 = rho**2
      t48 = t47**2
      t49 = 1/t48
      t50 = t46*t49
      t52 = 1.D0-1.D0*t50
      t53 = t44*t52
      t55 = 0.3799574853701528D-1*t33*t53
      t57 = 1.D0+0.1274696188700087D0*t5
      t62 = 0.1112037486309468D2*t8+0.3844746237447211D1*t5
     #+0.1644733775567609D1*t11+0.2405871291288192D0*t13
      t65 = 1.D0+0.321646831778707D2/t62
      t66 = dlog(t65)
      t69 = -0.3109D-1*t57*t66+t21
      t70 = t69*t44
      t72 = 0.1923661050931536D1*t70*t50
      t73 = t38**2
      t75 = t42**2
      t77 = 0.5D0*t73+0.5D0*t75
      t78 = t77**2
      t79 = t78*t77
      t80 = 1/t78
      t81 = sigma*t80
      t82 = rho**(1.D0/3.D0)
      t84 = 1/t82/t47
      t87 = -t21+t55+t72
      t88 = 1/t79
      t89 = t87*t88
      t91 = dexp(-0.4042761511756372D2*t89)
      t92 = t91-1.D0
      t93 = 1/t92
      t94 = sigma**2
      t95 = t93*t94
      t96 = t78**2
      t97 = 1/t96
      t98 = t82**2
      t100 = 1/t98/t48
      t101 = t97*t100
      t104 = 0.634682060977037D-1*t81*t84+0.1086645186223595D-1*t95*t101
      t105 = t93*sigma
      t106 = t80*t84
      t109 = t92**2
      t110 = 1/t109
      t111 = t110*t94
      t114 = 1.D0+0.1712109500228824D0*t105*t106+0.2931318940773793D-1
     #*t111*t101
      t115 = 1/t114
      t118 = 1.D0+0.2697586091519874D1*t104*t115
      t119 = dlog(t118)
      t121 = 0.2473556743557577D-1*t79*t119
      t124 = 0.2568D1+0.1443307452126544D2*t5+0.2843543831490386D-2*t13
      t128 = 1.D0+0.5411317332115466D1*t5+0.1816419932959077D0*t13
     #+0.1763993811759022D-1*t4
      t129 = 1/t128
      t132 = 0.1D-2*t124*t129-0.1853571428571429D-2
      t133 = t132*t77
      t136 = 1/t98/t47
      t140 = dexp(-0.261211729852336D1*t78*t136*sigma)
      t141 = sigma*t84*t140
      t143 = 0.1D1*t133*t141
      zk(i) = rho*(-t21+t55+t72+t121+t143)
      t149 = 0.1333333333333333D1*t38*t4-0.1333333333333333D1*t42*t4
      t151 = t33*t149*t52
      t152 = 0.3799574853701528D-1*t151
      t153 = t45*t35
      t154 = t44*t153
      t156 = t33*t154*t49
      t157 = 0.1519829941480611D0*t156
      t158 = t69*t149
      t159 = t158*t50
      t160 = 0.1923661050931536D1*t159
      t161 = t153*t49
      t162 = t70*t161
      t163 = 0.7694644203726145D1*t162
      t164 = t78*t119
      t165 = 1/t38
      t168 = 1/t42
      t171 = 0.3333333333333333D0*t165*t4-0.3333333333333333D0*t168*t4
      t172 = t164*t171
      t173 = 0.7420670230672731D-1*t172
      t174 = sigma*t88
      t175 = t84*t171
      t178 = t111*t97
      t179 = t152-t157+t160+t163
      t182 = t87*t97
      t185 = -0.4042761511756372D2*t179*t88+0.1212828453526912D3*t182
     #*t171
      t186 = t100*t185
      t187 = t186*t91
      t190 = t96*t77
      t191 = 1/t190
      t192 = t191*t100
      t193 = t192*t171
      t196 = -0.1269364121954074D0*t174*t175-0.1086645186223595D-1
     #*t178*t187-0.4346580744894379D-1*t95*t193
      t199 = t114**2
      t200 = 1/t199
      t201 = t104*t200
      t202 = t110*sigma
      t203 = t202*t80
      t204 = t84*t185
      t205 = t204*t91
      t208 = t88*t84
      t209 = t208*t171
      t213 = 1/t109/t92
      t214 = t213*t94
      t215 = t214*t97
      t220 = -0.1712109500228824D0*t203*t205-0.3424219000457648D0*t105
     #*t209-0.5862637881547585D-1*t215*t187-0.1172527576309517D0*t111
     #*t193
      t223 = 0.2697586091519874D1*t196*t115-0.2697586091519874D1*t201
     #*t220
      t224 = t79*t223
      t225 = 1/t118
      t226 = t224*t225
      t227 = 0.2473556743557577D-1*t226
      t228 = t132*t171
      t229 = t228*t141
      t230 = 0.1D1*t229
      t232 = t132*t78*t94
      t233 = t48*rho
      t234 = 1/t233
      t236 = t234*t171*t140
      t237 = t232*t236
      t238 = 0.522423459704672D1*t237
      t241 = 1/t13
      t242 = 1/t47
      t243 = t241*t242
      t244 = t243*t19
      t245 = 0.2747799777968419D-2*t244
      t246 = t15**2
      t247 = 1/t246
      t248 = t7*t247
      t249 = t8**2
      t250 = t249**2
      t251 = t250*t8
      t252 = 1/t251
      t253 = t252*t242
      t256 = 1/t11
      t257 = t256*t242
      t259 = 1/t5
      t260 = t259*t242
      t262 = -0.99709173929518D0*t253-0.7418564737168958D0*t243
     #-0.4002143174996817D0*t257-0.1264669550498372D0*t260
      t263 = 1/t18
      t265 = t248*t262*t263
      t266 = 0.1D1*t265
      t267 = t32*t44
      t268 = t267*t52
      t269 = t243*t268
      t270 = 0.8740794299481065D-3*t269
      t271 = t28**2
      t272 = 1/t271
      t273 = t23*t272
      t278 = -0.135956911724794D1*t253-0.7491972878592054D0*t243
     #-0.2150486235638321D0*t257-0.1274341730084892D0*t260
      t279 = t273*t278
      t280 = 1/t31
      t281 = t280*t44
      t282 = t281*t52
      t283 = t279*t282
      t284 = 0.1124999956683108D1*t283
      t285 = t38*t35
      t288 = t42*t35
      t291 = -0.1333333333333333D1*t285*t242+0.1333333333333333D1*t288
     #*t242
      t293 = t33*t291*t52
      t294 = 0.3799574853701528D-1*t293
      t295 = t44*t46
      t297 = t33*t295*t234
      t298 = 0.1519829941480611D0*t297
      t301 = t62**2
      t302 = 1/t301
      t303 = t57*t302
      t308 = -0.1853395810515781D1*t253-0.128158207914907D1*t243
     #-0.8223668877838045D0*t257-0.1603914194192128D0*t260
      t309 = 1/t65
      t315 = 0.1321010150222857D-2*t243*t66+0.1D1*t303*t308*t309
     #-0.2747799777968419D-2*t244-0.1D1*t265
      t316 = t315*t44
      t317 = t316*t50
      t318 = 0.1923661050931536D1*t317
      t319 = t69*t291
      t320 = t319*t50
      t321 = 0.1923661050931536D1*t320
      t322 = t46*t234
      t323 = t70*t322
      t324 = 0.7694644203726145D1*t323
      t325 = t165*t35
      t328 = t168*t35
      t331 = -0.3333333333333333D0*t325*t242+0.3333333333333333D0*t328
     #*t242
      t332 = t164*t331
      t334 = t84*t331
      t337 = t47*rho
      t339 = 1/t82/t337
      t342 = t245+t266-t270-t284+t294+t298+t318+t321-t324
      t347 = -0.4042761511756372D2*t342*t88+0.1212828453526912D3*t182
     #*t331
      t348 = t100*t347
      t349 = t348*t91
      t352 = t192*t331
      t356 = 1/t98/t233
      t357 = t97*t356
      t360 = -0.1269364121954074D0*t174*t334-0.148092480894642D0*t81
     #*t339-0.1086645186223595D-1*t178*t349-0.4346580744894379D-1*t95
     #*t352-0.5071010869043442D-1*t95*t357
      t363 = t84*t347
      t364 = t363*t91
      t367 = t208*t331
      t370 = t80*t339
      t379 = -0.1712109500228824D0*t203*t364-0.3424219000457648D0*t105
     #*t367-0.3994922167200589D0*t105*t370-0.5862637881547585D-1*t215
     #*t349-0.1172527576309517D0*t111*t352-0.136794883902777D0*t111*t357
      t382 = 0.2697586091519874D1*t360*t115-0.2697586091519874D1*t201
     #*t379
      t383 = t79*t382
      t384 = t383*t225
      t388 = -0.4811024840421814D1*t243-0.1895695887660258D-2*t260
      t391 = t128**2
      t392 = 1/t391
      t393 = t124*t392
      t397 = -0.1803772444038489D1*t243-0.1210946621972718D0*t260
     #-0.1763993811759022D-1*t242
      t400 = 0.1D-2*t388*t129-0.1D-2*t393*t397
      t401 = t400*t77
      t402 = t401*t141
      t404 = t132*t331
      t405 = t404*t141
      t408 = sigma*t339*t140
      t409 = t133*t408
      t411 = t133*sigma
      t412 = t77*t136
      t413 = sigma*t331
      t417 = 1/t98/t337
      t418 = t78*t417
      t421 = -0.522423459704672D1*t412*t413+0.6965646129395627D1*t418
     #*sigma
      t423 = t84*t421*t140
      t424 = t411*t423
      t426 = t245+t266-t270-t284+t294+t298+t318+t321-t324
     #+0.7420670230672731D-1*t332+0.2473556743557577D-1*t384+0.1D1
     #*t402+0.1D1*t405-0.2333333333333333D1*t409+0.1D1*t424
      t427 = rho*t426
      vrhoa(i) = rho*(t152-t157+t160+t163+t173+t227+t230-t238)-t21+t55
     #+t72+t121+t143+t427
      t428 = -t149
      t430 = t33*t428*t52
      t431 = 0.3799574853701528D-1*t430
      t432 = 0.1519829941480611D0*t156
      t433 = t69*t428
      t434 = t433*t50
      t435 = 0.1923661050931536D1*t434
      t436 = 0.7694644203726145D1*t162
      t437 = -t171
      t438 = t164*t437
      t439 = 0.7420670230672731D-1*t438
      t440 = t84*t437
      t443 = t431+t432+t435-t436
      t448 = -0.4042761511756372D2*t443*t88+0.1212828453526912D3*t182
     #*t437
      t449 = t100*t448
      t450 = t449*t91
      t453 = t192*t437
      t456 = -0.1269364121954074D0*t174*t440-0.1086645186223595D-1
     #*t178*t450-0.4346580744894379D-1*t95*t453
      t459 = t84*t448
      t460 = t459*t91
      t463 = t208*t437
      t470 = -0.1712109500228824D0*t203*t460-0.3424219000457648D0*t105
     #*t463-0.5862637881547585D-1*t215*t450-0.1172527576309517D0*t111
     #*t453
      t473 = 0.2697586091519874D1*t456*t115-0.2697586091519874D1*t201
     #*t470
      t474 = t79*t473
      t475 = t474*t225
      t476 = 0.2473556743557577D-1*t475
      t477 = t132*t437
      t478 = t477*t141
      t479 = 0.1D1*t478
      t481 = t234*t437*t140
      t482 = t232*t481
      t483 = 0.522423459704672D1*t482
      vrhob(i) = rho*(t431+t432+t435-t436+t439+t476+t479-t483)-t21+t55
     #+t72+t121+t143+t427
      t489 = 0.634682060977037D-1*t106+0.2173290372447189D-1*t105*t101
      t492 = t93*t80
      t497 = 0.1712109500228824D0*t492*t84+0.5862637881547585D-1*t202
     #*t101
      t500 = 0.2697586091519874D1*t489*t115-0.2697586091519874D1*t201
     #*t497
      t502 = t79*t500*t225
      t503 = 0.2473556743557577D-1*t502
      t504 = t84*t140
      t505 = t133*t504
      t506 = 0.1D1*t505
      t507 = t132*t79
      t509 = sigma*t234*t140
      t510 = t507*t509
      t511 = 0.261211729852336D1*t510
      vsigmaaa(i) = rho*(t503+t506-t511)
      vsigmaab(i) = 2.D0*vsigmaaa(i)
      vsigmabb(i) = vsigmaaa(i)
      t516 = 0.2D1*t265
      t517 = 0.7599149707403056D-1*t293
      t520 = 1/t73
      t523 = 1/t75
      t526 = 0.4444444444444444D0*t520*t242+0.4444444444444444D0*t523
     #*t242
      t529 = 0.3799574853701528D-1*t33*t526*t52
      t532 = t33*t149*t153*t49
      t533 = 0.3039659882961223D0*t532
      t536 = t33*t44*t45*t49
      t537 = 0.4559489824441834D0*t536
      t540 = 0.1923661050931536D1*t69*t526*t50
      t541 = t158*t161
      t542 = 0.1538928840745229D2*t541
      t544 = t70*t45*t49
      t545 = 0.2308393261117844D2*t544
      t546 = t77*t119
      t547 = t171**2
      t550 = t78*t223
      t551 = t225*t171
      t554 = 1/t39
      t557 = 1/t43
      t560 = -0.1111111111111111D0*t554*t242-0.1111111111111111D0*t557
     #*t242
      t562 = 0.7420670230672731D-1*t164*t560
      t563 = sigma*t97
      t569 = 0.1269364121954074D0*t174*t84*t560
      t570 = t185**2
      t571 = t100*t570
      t573 = dexp(-0.8085523023512745D2*t89)
      t574 = t571*t573
      t577 = t111*t191
      t578 = t91*t171
      t579 = t186*t578
      t585 = t179*t97
      t588 = t87*t191
      t592 = 0.1212828453526912D3*t182*t560
      t593 = -0.4042761511756372D2*(t529-t533-t537+t540+t542+t545)*t88
     #+0.2425656907053823D3*t585*t171-0.4851313814107647D3*t588*t547
     #+t592
      t595 = t100*t593*t91
      t598 = t571*t91
      t602 = 1/t96/t78
      t603 = t602*t100
      t604 = t603*t547
      t607 = t192*t560
      t609 = 0.4346580744894379D-1*t95*t607
      t613 = t196*t200
      t618 = t104/t199/t114
      t619 = t220**2
      t622 = t213*sigma
      t623 = t622*t80
      t624 = t84*t570
      t628 = t202*t88
      t639 = t97*t84
      t645 = 0.3424219000457648D0*t105*t208*t560
      t646 = t109**2
      t649 = 1/t646*t94*t97
      t652 = t214*t191
      t662 = 0.1172527576309517D0*t111*t607
      t663 = 0.3424219000457648D0*t623*t624*t573+0.6848438000915295D0
     #*t628*t204*t578-0.1712109500228824D0*t203*t84*t593*t91
     #-0.1712109500228824D0*t203*t624*t91+0.1027265700137294D1*t105
     #*t639*t547-t645+0.1758791364464276D0*t649*t574
     #+0.4690110305238068D0*t652*t579-0.5862637881547585D-1*t215*t595
     #-0.5862637881547585D-1*t215*t598+0.5862637881547585D0*t111*t604
     #-t662
      t670 = t223**2
      t672 = t118**2
      t673 = 1/t672
      t678 = 0.1D1*t132*t560*t141
      t681 = t234*t77
      t682 = t681*t140
      t688 = 0.522423459704672D1*t232*t234*t560*t140
      t690 = t507*t94*sigma
      t691 = t48*t337
      t693 = 1/t98/t691
      t698 = t529-t533-t537+t540+t542+t545+0.1484134046134546D0*t546
     #*t547+0.1484134046134546D0*t550*t551+t562+0.2473556743557577D-1
     #*t79*(0.2697586091519874D1*(0.3808092365862222D0*t563*t84*t547
     #-t569+0.2173290372447189D-1*t215*t574+0.8693161489788758D-1*t577
     #*t579-0.1086645186223595D-1*t178*t595-0.1086645186223595D-1*t178
     #*t598+0.2173290372447189D0*t95*t604-t609)*t115
     #-0.5395172183039748D1*t613*t220+0.5395172183039748D1*t618*t619
     #-0.2697586091519874D1*t201*t663)*t225-0.2473556743557577D-1*t79
     #*t670*t673+t678-0.1567270379114016D2*t132*t547*t94*t682-t688
     #+0.272926271249799D2*t690*t693*t547*t140
      t701 = 0.5495599555936838D-2*t244
      t702 = 0.1484134046134546D0*t332
      t705 = 0.7599149707403056D-1*t151+0.3847322101863073D1*t159
     #+0.1538928840745229D2*t162+t516+t517+0.1484134046134546D0*t172
     #-0.3039659882961223D0*t156+rho*t698+0.4947113487115154D-1*t226
     #+t701+t702+0.2D1*t229-0.1044846919409344D2*t237
      t706 = 0.3039659882961223D0*t297
      t707 = 0.3847322101863073D1*t317
      t708 = 0.3847322101863073D1*t320
      t709 = 0.1748158859896213D-2*t269
      t710 = 0.2249999913366216D1*t283
      t711 = 0.1538928840745229D2*t323
      t712 = 0.4947113487115154D-1*t384
      t713 = 0.4666666666666667D1*t409
      t714 = 0.2D1*t424
      t715 = 0.2D1*t405
      t716 = 0.2D1*t402
      t717 = t78*t382
      t718 = t225*t331
      t724 = t278**2
      t727 = 0.2249999913366216D1*t23/t271/t28*t724*t282
      t728 = t382**2
      t738 = 1/t82/t48
      t743 = t421**2
      t752 = t331**2
      t756 = t77*t417
      t762 = 1/t337
      t770 = -0.1111111111111111D0*t554*t45*t49+0.6666666666666667D0
     #*t325*t762-0.1111111111111111D0*t557*t45*t49
     #-0.6666666666666667D0*t328*t762
      t782 = t404*sigma
      t792 = 1/t13/t4*t49
      t794 = t241*t762
      t798 = 1/t5/t4*t49
      t800 = t259*t762
      t811 = t397**2
      t841 = t356*t347*t91
      t844 = t91*t331
      t845 = t348*t844
      t851 = t262**2
      t853 = t7/t246/t15*t851*t263
      t854 = 0.2D1*t853
      t857 = 1/t251/t4*t49
      t859 = t252*t762
      t865 = 1/t11/t4*t49
      t867 = t256*t762
      t873 = t248*(-0.8309097827459833D0*t857+0.199418347859036D1*t859
     #-0.4945709824779306D0*t792+0.1483712947433792D1*t794
     #-0.2001071587498409D0*t865+0.8004286349993634D0*t867
     #-0.4215565168327908D-1*t798+0.2529339100996745D0*t800)*t263
      t874 = 0.1D1*t873
      t875 = t246**2
      t878 = t18**2
      t881 = t7/t875*t851/t878
      t882 = 0.160818243221511D2*t881
      t894 = t308**2
      t910 = t301**2
      t913 = t65**2
      t918 = t792*t19
      t920 = t794*t19
      t924 = t243*t247*t262*t263
      t929 = 0.8806734334819047D-3*t792*t66-0.2642020300445714D-2*t794
     #*t66-0.8497974591333914D-1*t243*t302*t308*t309-0.2D1*t57/t301
     #/t62*t894*t309+0.1D1*t303*(-0.1544496508763151D1*t857
     #+0.3706791621031562D1*t859-0.854388052766047D0*t792
     #+0.2563164158298141D1*t794-0.4111834438919023D0*t865
     #+0.1644733775567609D1*t867-0.5346380647307093D-1*t798
     #+0.3207828388384256D0*t800)*t309+0.321646831778707D2*t57/t910
     #*t894/t913-0.1831866518645613D-2*t918+0.5495599555936838D-2*t920
     #+0.8837926660346786D-1*t924+0.2D1*t853-0.1D1*t873
     #-0.160818243221511D2*t881
      t932 = 0.1923661050931536D1*t929*t44*t50
      t943 = 0.4444444444444444D0*t520*t45*t49+0.2666666666666667D1
     #*t285*t762+0.4444444444444444D0*t523*t45*t49
     #-0.2666666666666667D1*t288*t762
      t946 = 0.3799574853701528D-1*t33*t943*t52
      t947 = 0.8837926660346786D-1*t924
      t948 = t271**2
      t952 = t31**2
      t957 = 0.3330964519106732D2*t23/t948*t724/t952*t44*t52
      t969 = 0.1124999956683108D1*t273*(-0.1132974264373283D1*t857
     #+0.271913823449588D1*t859-0.4994648585728036D0*t792
     #+0.1498394575718411D1*t794-0.1075243117819161D0*t865
     #+0.4300972471276643D0*t867-0.4247805766949639D-1*t798
     #+0.2548683460169784D0*t800)*t282
      t970 = 0.1831866518645613D-2*t918
      t972 = 0.1748158859896213D-2*t794*t268
      t976 = 0.2249999913366216D1*t279*t280*t291*t52
      t977 = t727-t854+t874+t882+t932+t946-t947-t957-t969+t970+t972-t976
      t978 = 0.5495599555936838D-2*t920
      t982 = 0.3039659882961223D0*t33*t291*t46*t234
      t985 = 0.8999999653464862D1*t279*t281*t322
      t986 = t48*t47
      t987 = 1/t986
      t990 = 0.3847322101863073D2*t70*t46*t987
      t992 = 0.1538928840745229D2*t319*t322
      t995 = 0.3847322101863073D1*t315*t291*t50
      t997 = 0.1538928840745229D2*t316*t322
      t1000 = 0.1923661050931536D1*t69*t943*t50
      t1003 = 0.7599149707403056D0*t33*t295*t987
      t1007 = 0.1748158859896213D-2*t243*t32*t291*t52
      t1009 = 0.5827196199654043D-3*t792*t268
      t1014 = 0.5176049209143758D-1*t243*t272*t278*t280*t53
      t1019 = 0.6992635439584852D-2*t241/t691*t267*t46
      t1020 = -t978+t982-t985+t990-t992+t995-t997+t1000-t1003-t1007
     #-t1009+t1014-t1019
      t1024 = t342*t97
      t1031 = -0.4042761511756372D2*(t977+t1020)*t88
     #+0.2425656907053823D3*t1024*t331-0.4851313814107647D3*t588*t752
     #+0.1212828453526912D3*t182*t770
      t1033 = t100*t1031*t91
      t1036 = t347**2
      t1037 = t100*t1036
      t1038 = t1037*t573
      t1041 = t1037*t91
      t1046 = t97/t98/t986
      t1049 = t603*t752
      t1052 = t191*t356
      t1053 = t1052*t331
      t1056 = t192*t770
      t1059 = 0.3808092365862222D0*t563*t84*t752+0.5923699235785679D0
     #*t174*t339*t331-0.1269364121954074D0*t174*t84*t770
     #+0.4936416029821399D0*t81*t738+0.1014202173808688D0*t178*t841
     #+0.8693161489788758D-1*t577*t845-0.1086645186223595D-1*t178
     #*t1033+0.2173290372447189D-1*t215*t1038-0.1086645186223595D-1
     #*t178*t1041+0.2873572825791284D0*t95*t1046+0.2173290372447189D0
     #*t95*t1049+0.4056808695234754D0*t95*t1053-0.4346580744894379D-1
     #*t95*t1056
      t1062 = t360*t200
      t1065 = t379**2
      t1085 = t84*t1036
      t1100 = t88*t339
      t1117 = -0.5862637881547585D-1*t215*t1041+0.6848438000915295D0
     #*t628*t363*t844+0.7989844334401178D0*t203*t339*t347*t91
     #+0.4690110305238068D0*t652*t845-0.1712109500228824D0*t203*t84
     #*t1031*t91-0.5862637881547585D-1*t215*t1033+0.3424219000457648D0
     #*t623*t1085*t573+0.1331640722400196D1*t105*t80*t738
     #-0.1712109500228824D0*t203*t1085*t91+0.547179535611108D0*t215
     #*t841+0.1027265700137294D1*t105*t639*t752+0.1597968866880236D1
     #*t105*t1100*t331-0.3424219000457648D0*t105*t208*t770
     #+0.1758791364464276D0*t649*t1038+0.7751710087824029D0*t111*t1046
     #+0.5862637881547585D0*t111*t1049+0.1094359071222216D1*t111*t1053
     #-0.1172527576309517D0*t111*t1056
      s1 = 0.1484134046134546D0*t717*t718+t727-0.2473556743557577D-1
     #*t79*t728*t673-0.4666666666666667D1*t401*t408+0.2D1*t400*t331
     #*t141+0.7777777777777778D1*t133*sigma*t738*t140+0.1D1*t411*t84
     #*t743*t140-0.4666666666666667D1*t411*t339*t421*t140+0.1D1*t411
     #*t84*(-0.522423459704672D1*t752*t136*sigma+0.2786258451758251D2
     #*t756*t413-0.522423459704672D1*t412*sigma*t770
     #-0.2554070247445063D2*t78*t100*sigma)*t140+0.2D1*t782*t423
      t1126 = s1-0.4666666666666667D1*t404*t408+0.2D1*t401*sigma*t423
     #+0.1D1*(0.1D-2*(-0.3207349893614542D1*t792+0.9622049680843627D1
     #*t794-0.6318986292200858D-3*t798+0.3791391775320515D-2*t800)
     #*t129-0.2D-2*t388*t392*t397+0.2D-2*t124/t391/t128*t811-0.1D-2
     #*t393*(-0.1202514962692326D1*t792+0.3607544888076978D1*t794
     #-0.4036488739909061D-1*t798+0.2421893243945437D0*t800
     #+0.3527987623518044D-1*t762))*t77*t141+0.1D1*t132*t770*t141
     #+0.2473556743557577D-1*t79*(0.2697586091519874D1*t1059*t115
     #-0.5395172183039748D1*t1062*t379+0.5395172183039748D1*t618*t1065
     #-0.2697586091519874D1*t201*t1117)*t225+0.1484134046134546D0*t546
     #*t752-t854+t874+t882+t932
      t1130 = -t985+t990-t992+t995-t997+t1000-t1003-t1007-t1009+t1014
     #-t1019
      t1133 = rho*(t1126+t946-t947-t957+0.7420670230672731D-1*t164
     #*t770-t969+t970+t972-t976-t978+t982+t1130)
      t1137 = 0.1519829941480611D0*t33*t149*t46*t234
      t1139 = t33*t154*t234
      t1140 = 0.6079319765922445D0*t1139
      t1151 = -0.4444444444444444D0*t520*t35*t762-0.1333333333333333D1
     #*t38*t242-0.4444444444444444D0*t523*t35*t762
     #+0.1333333333333333D1*t42*t242
      t1154 = 0.3799574853701528D-1*t33*t1151*t52
      t1157 = 0.1923661050931536D1*t315*t149*t50
      t1158 = t316*t161
      t1159 = 0.7694644203726145D1*t1158
      t1162 = 0.1923661050931536D1*t69*t1151*t50
      t1164 = t70*t153*t234
      t1165 = 0.3077857681490458D2*t1164
      t1169 = 0.8740794299481065D-3*t243*t32*t149*t52
      t1172 = t241*t987*t267*t153
      t1173 = 0.3496317719792426D-2*t1172
      t1175 = t279*t281*t161
      t1176 = 0.4499999826732431D1*t1175
      t1180 = 0.1124999956683108D1*t279*t280*t149*t52
      t1183 = t33*t291*t153*t49
      t1184 = 0.1519829941480611D0*t1183
      t1186 = 0.7694644203726145D1*t158*t322
      t1189 = t1137+t1140+t1154+t1157+t1159+t1162-t1165-t1169+t1173
     #+t1176-t1180-t1184-t1186+0.7420670230672731D-1*t717*t551
      t1200 = 0.1111111111111111D0*t554*t35*t762-0.3333333333333333D0
     #*t165*t242+0.1111111111111111D0*t557*t35*t762
     #+0.3333333333333333D0*t168*t242
      t1202 = 0.7420670230672731D-1*t164*t1200
      t1204 = 0.7420670230672731D-1*t550*t718
      t1205 = t319*t161
      t1206 = 0.7694644203726145D1*t1205
      t1207 = t331*t171
      t1209 = 0.1484134046134546D0*t546*t1207
      t1212 = 0.2473556743557577D-1*t383*t673*t223
      t1229 = 0.1D1*t400*t171*t141
      t1231 = 0.2333333333333333D1*t228*t408
      t1234 = 0.1D1*t132*t1200*t141
      t1235 = t404*t94
      t1236 = t171*t140
      t1238 = t1235*t681*t1236
      t1242 = t232*t987*t171*t140
      t1245 = t400*t78*t94
      t1247 = 0.522423459704672D1*t1245*t236
      t1248 = t228*sigma
      t1250 = 0.1D1*t1248*t423
      t1251 = t234*t421
      t1254 = 0.522423459704672D1*t232*t1251*t1236
      t1255 = t334*t171
      t1264 = t573*t185
      t1265 = t348*t1264
      t1268 = t348*t578
      t1271 = -t1169+t1173-t1180+t1176+t1154-t1184+t1137+t1140+t1157
     #+t1159+t1162+t1206-t1186-t1165
      t1282 = -0.4042761511756372D2*t1271*t88+0.1212828453526912D3
     #*t1024*t171+0.1212828453526912D3*t585*t331-0.4851313814107647D3
     #*t588*t1207+0.1212828453526912D3*t182*t1200
      t1284 = t100*t1282*t91
      t1287 = t185*t91
      t1288 = t348*t1287
      t1291 = t100*t331
      t1292 = t1291*t1287
      t1295 = t95*t602
      t1296 = t1291*t171
      t1299 = t192*t1200
      t1303 = t356*t185*t91
      t1306 = t1052*t171
      t1309 = 0.3808092365862222D0*t563*t1255-0.1269364121954074D0
     #*t174*t84*t1200+0.296184961789284D0*t174*t339*t171
     #+0.2173290372447189D-1*t215*t1265+0.4346580744894379D-1*t577
     #*t1268-0.1086645186223595D-1*t178*t1284-0.1086645186223595D-1
     #*t178*t1288+0.4346580744894379D-1*t577*t1292
     #+0.2173290372447189D0*t1295*t1296-0.4346580744894379D-1*t95
     #*t1299+0.5071010869043442D-1*t178*t1303+0.2028404347617377D0*t95
     #*t1306
      t1335 = t105*t97
      t1358 = t111*t602
      t1367 = 0.3424219000457648D0*t623*t363*t1264
     #+0.3424219000457648D0*t628*t363*t578-0.1712109500228824D0*t203
     #*t84*t1282*t91-0.1712109500228824D0*t203*t363*t1287
     #+0.3424219000457648D0*t628*t334*t1287+0.1027265700137294D1*t1335
     #*t1255-0.3424219000457648D0*t105*t208*t1200+0.3994922167200589D0
     #*t203*t339*t185*t91+0.7989844334401178D0*t105*t1100*t171
     #+0.1758791364464276D0*t649*t1265+0.2345055152619034D0*t652*t1268
     #-0.5862637881547585D-1*t215*t1284-0.5862637881547585D-1*t215
     #*t1288+0.2345055152619034D0*t652*t1292+0.5862637881547585D0
     #*t1358*t1296-0.1172527576309517D0*t111*t1299+0.273589767805554D0
     #*t215*t1303+0.547179535611108D0*t111*t1306
      t1373 = 0.2473556743557577D-1*t79*(0.2697586091519874D1*t1309
     #*t115-0.2697586091519874D1*t1062*t220-0.2697586091519874D1*t613
     #*t379+0.5395172183039748D1*t618*t379*t220-0.2697586091519874D1
     #*t201*t1367)*t225
      t1374 = t1202+t1204+t1206+t1209-t1212+0.1D1*t411*t84*(
     #-0.522423459704672D1*t171*t136*t413-0.522423459704672D1*t412
     #*sigma*t1200+0.1393129225879125D2*t756*sigma*t171)*t140+t1229
     #-t1231+t1234-0.522423459704672D1*t1238+0.1218988072644235D2
     #*t1242-t1247+t1250-t1254+t1373
      t1376 = rho*(t1189+t1374)
      t1383 = t1202+t1204+t1206+t1209-t1212+t1229-t1231+t1234
     #-0.1044846919409344D2*t1238+0.261211729852336D2*t1242-t1247
     #+t1250-t1254-0.522423459704672D1*t232*t234*t1200*t140+t1373
      t1385 = rho*(t1189+t1383)
      t1386 = t706+t707+t708-t709-t710-t711+t712-t713+t714+t715+t716
     #+t1133+t1376+t1385
      v2rhoa2(i) = t705+t1386
      t1392 = -0.1538928840745229D2*t162+t516+t517
     #+0.4947113487115154D-1*t475+0.3039659882961223D0*t156+0.2D1*t478
     #-0.1044846919409344D2*t482+t701+t702+t706+t707+t708-t709
      t1394 = t225*t437
      t1396 = 0.7420670230672731D-1*t717*t1394
      t1397 = 0.6079319765922445D0*t1139
      t1398 = 0.7694644203726145D1*t1158
      t1399 = 0.3077857681490458D2*t1164
      t1400 = t673*t473
      t1402 = 0.2473556743557577D-1*t383*t1400
      t1403 = 0.3496317719792426D-2*t1172
      t1404 = 0.4499999826732431D1*t1175
      t1405 = 0.1519829941480611D0*t1183
      t1406 = -t1200
      t1412 = 0.7420670230672731D-1*t164*t1406
      t1415 = t232*t987*t437*t140
      t1418 = 0.2333333333333333D1*t477*t408
      t1419 = t437*t140
      t1420 = t681*t1419
      t1421 = t1235*t1420
      t1423 = t477*sigma
      t1425 = 0.1D1*t1423*t423
      t1426 = t1396-t1397-t1398+t1399-t1402-t1403-t1404+t1405
     #-0.522423459704672D1*t232*t234*t1406*t140+t1412
     #+0.261211729852336D2*t1415-t1418-0.1044846919409344D2*t1421+t1425
      t1430 = 0.7694644203726145D1*t1205
      t1433 = t334*t437
      t1442 = t573*t448
      t1443 = t348*t1442
      t1446 = t91*t437
      t1447 = t348*t1446
      t1453 = 0.8740794299481065D-3*t243*t32*t428*t52
      t1457 = 0.1124999956683108D1*t279*t280*t428*t52
      t1458 = -t1151
      t1461 = 0.3799574853701528D-1*t33*t1458*t52
      t1465 = 0.1519829941480611D0*t33*t428*t46*t234
      t1468 = 0.1923661050931536D1*t315*t428*t50
      t1471 = 0.1923661050931536D1*t69*t1458*t50
      t1473 = 0.7694644203726145D1*t433*t322
      t1474 = -t1453-t1403-t1457-t1404+t1461+t1405+t1465-t1397+t1468
     #-t1398+t1471-t1430-t1473+t1399
      t1479 = t443*t97
      t1482 = t331*t437
      t1487 = -0.4042761511756372D2*t1474*t88+0.1212828453526912D3
     #*t1024*t437+0.1212828453526912D3*t1479*t331-0.4851313814107647D3
     #*t588*t1482+0.1212828453526912D3*t182*t1406
      t1489 = t100*t1487*t91
      t1492 = t448*t91
      t1493 = t348*t1492
      t1496 = t1291*t1492
      t1499 = t1291*t437
      t1502 = t192*t1406
      t1506 = t356*t448*t91
      t1509 = t1052*t437
      t1512 = 0.3808092365862222D0*t563*t1433-0.1269364121954074D0
     #*t174*t84*t1406+0.296184961789284D0*t174*t339*t437
     #+0.2173290372447189D-1*t215*t1443+0.4346580744894379D-1*t577
     #*t1447-0.1086645186223595D-1*t178*t1489-0.1086645186223595D-1
     #*t178*t1493+0.4346580744894379D-1*t577*t1496
     #+0.2173290372447189D0*t1295*t1499-0.4346580744894379D-1*t95
     #*t1502+0.5071010869043442D-1*t178*t1506+0.2028404347617377D0*t95
     #*t1509
      t1517 = t456*t200
      t1569 = 0.3424219000457648D0*t623*t363*t1442
     #+0.3424219000457648D0*t628*t363*t1446-0.1712109500228824D0*t203
     #*t84*t1487*t91-0.1712109500228824D0*t203*t363*t1492
     #+0.3424219000457648D0*t628*t334*t1492+0.1027265700137294D1*t1335
     #*t1433-0.3424219000457648D0*t105*t208*t1406+0.3994922167200589D0
     #*t203*t339*t448*t91+0.7989844334401178D0*t105*t1100*t437
     #+0.1758791364464276D0*t649*t1443+0.2345055152619034D0*t652*t1447
     #-0.5862637881547585D-1*t215*t1489-0.5862637881547585D-1*t215
     #*t1493+0.2345055152619034D0*t652*t1496+0.5862637881547585D0
     #*t1358*t1499-0.1172527576309517D0*t111*t1502+0.273589767805554D0
     #*t215*t1506+0.547179535611108D0*t111*t1509
      t1578 = t78*t473
      t1587 = 0.1D1*t132*t1406*t141-t1430-0.522423459704672D1*t1245
     #*t481+0.2473556743557577D-1*t79*(0.2697586091519874D1*t1512*t115
     #-0.2697586091519874D1*t1062*t470-0.2697586091519874D1*t1517*t379
     #+0.5395172183039748D1*t618*t379*t470-0.2697586091519874D1*t201
     #*t1569)*t225+0.1484134046134546D0*t546*t1482
     #+0.7420670230672731D-1*t1578*t718-0.522423459704672D1*t232*t1251
     #*t1419+t1461+t1468+t1471-t1473+0.1D1*t400*t437*t141-t1453-t1457
     #+t1465
      t1589 = rho*(t1426+t1587)
      t1608 = t1396-t1397-t1398+t1399-t1402-t1403-t1404+t1405+t1412
     #+0.1218988072644235D2*t1415-t1418-0.522423459704672D1*t1421
     #+0.1D1*t411*t84*(-0.522423459704672D1*t437*t136*t413
     #-0.522423459704672D1*t412*sigma*t1406+0.1393129225879125D2*t756
     #*sigma*t437)*t140+t1425
      t1610 = rho*(t1608+t1587)
      t1613 = t33*t428*t153*t49
      t1614 = 0.3039659882961223D0*t1613
      t1615 = t433*t161
      t1616 = 0.1538928840745229D2*t1615
      t1617 = t437**2
      t1625 = t448**2
      t1626 = t100*t1625
      t1627 = t1626*t573
      t1630 = t449*t1446
      t1640 = -0.4042761511756372D2*(t529+t1614-t537+t540-t1616+t545)
     #*t88+0.2425656907053823D3*t1479*t437-0.4851313814107647D3*t588
     #*t1617+t592
      t1642 = t100*t1640*t91
      t1645 = t1626*t91
      t1648 = t603*t1617
      t1656 = t470**2
      t1659 = t84*t1625
      t1686 = 0.3424219000457648D0*t623*t1659*t573
     #+0.6848438000915295D0*t628*t459*t1446-0.1712109500228824D0*t203
     #*t84*t1640*t91-0.1712109500228824D0*t203*t1659*t91
     #+0.1027265700137294D1*t105*t639*t1617-t645+0.1758791364464276D0
     #*t649*t1627+0.4690110305238068D0*t652*t1630
     #-0.5862637881547585D-1*t215*t1642-0.5862637881547585D-1*t215
     #*t1645+0.5862637881547585D0*t111*t1648-t662
      t1693 = t473**2
      t1705 = t529+t1614-t537+t540-t1616+t545+0.1484134046134546D0
     #*t546*t1617+0.1484134046134546D0*t1578*t1394+t562
     #+0.2473556743557577D-1*t79*(0.2697586091519874D1*
     #(0.3808092365862222D0*t563*t84*t1617-t569+0.2173290372447189D-1
     #*t215*t1627+0.8693161489788758D-1*t577*t1630
     #-0.1086645186223595D-1*t178*t1642-0.1086645186223595D-1*t178
     #*t1645+0.2173290372447189D0*t95*t1648-t609)*t115
     #-0.5395172183039748D1*t1517*t470+0.5395172183039748D1*t618*t1656
     #-0.2697586091519874D1*t201*t1686)*t225-0.2473556743557577D-1*t79
     #*t1693*t673+t678-0.1567270379114016D2*t132*t1617*t94*t682-t688
     #+0.272926271249799D2*t690*t693*t1617*t140
      t1707 = -t710-t711+t712+0.1484134046134546D0*t438+t1589-t713
     #+0.7599149707403056D-1*t430+0.3847322101863073D1*t434+t714+t715
     #+t716+t1610+t1133+rho*t1705
      v2rhob2(i) = t1392+t1707
      t1708 = -t526
      t1711 = 0.3799574853701528D-1*t33*t1708*t52
      t1712 = 0.1519829941480611D0*t532
      t1713 = 0.1519829941480611D0*t1613
      t1714 = 0.4559489824441834D0*t536
      t1717 = 0.1923661050931536D1*t69*t1708*t50
      t1718 = 0.7694644203726145D1*t541
      t1719 = 0.7694644203726145D1*t1615
      t1720 = 0.2308393261117844D2*t544
      t1721 = t171*t437
      t1726 = -t560
      t1731 = t175*t437
      t1737 = t186*t1442
      t1740 = t186*t1446
      t1754 = -0.4042761511756372D2*(t1711+t1712-t1713+t1714+t1717
     #-t1718+t1719-t1720)*t88+0.1212828453526912D3*t585*t437
     #+0.1212828453526912D3*t1479*t171-0.4851313814107647D3*t588*t1721
     #+0.1212828453526912D3*t182*t1726
      t1756 = t100*t1754*t91
      t1759 = t186*t1492
      t1762 = t100*t171
      t1763 = t1762*t1492
      t1766 = t1762*t437
      t1769 = t192*t1726
      t1817 = 0.3424219000457648D0*t623*t204*t1442
     #+0.3424219000457648D0*t628*t204*t1446-0.1712109500228824D0*t203
     #*t84*t1754*t91-0.1712109500228824D0*t203*t204*t1492
     #+0.3424219000457648D0*t628*t175*t1492+0.1027265700137294D1*t1335
     #*t1731-0.3424219000457648D0*t105*t208*t1726+0.1758791364464276D0
     #*t649*t1737+0.2345055152619034D0*t652*t1740
     #-0.5862637881547585D-1*t215*t1756-0.5862637881547585D-1*t215
     #*t1759+0.2345055152619034D0*t652*t1763+0.5862637881547585D0
     #*t1358*t1766-0.1172527576309517D0*t111*t1769
      t1836 = t693*t171
      s1 = t1711+t1712-t1713+t1714+t1717-t1718+t1719-t1720
     #+0.1484134046134546D0*t546*t1721
      s2 = s1+0.7420670230672731D-1*t1578*t551+0.7420670230672731D-1
     #*t164*t1726+0.7420670230672731D-1*t550*t1394
     #+0.2473556743557577D-1*t79*(0.2697586091519874D1*
     #(0.3808092365862222D0*t563*t1731-0.1269364121954074D0*t174*t84
     #*t1726+0.2173290372447189D-1*t215*t1737+0.4346580744894379D-1
     #*t577*t1740-0.1086645186223595D-1*t178*t1756
     #-0.1086645186223595D-1*t178*t1759+0.4346580744894379D-1*t577
     #*t1763+0.2173290372447189D0*t1295*t1766-0.4346580744894379D-1
     #*t95*t1769)*t115-0.2697586091519874D1*t613*t470
     #-0.2697586091519874D1*t1517*t220+0.5395172183039748D1*t618*t220
     #*t470-0.2697586091519874D1*t201*t1817)*t225
      t1840 = s2-0.2473556743557577D-1*t224*t1400+0.1D1*t132*t1726
     #*t141-0.1567270379114016D2*t228*t94*t1420-0.522423459704672D1
     #*t232*t234*t1726*t140+0.272926271249799D2*t690*t1836*t1419
      t1842 = t152+t160+t516+t517+t173+rho*t1840+t476+t479-t483+t227
     #+t701+t702+t230-t238+t706+t707
      t1847 = t708-t709-t710-t711+t712+t439+0.5D0*t1589-t713+t431+t435
     #+t714+t715+t716+0.5D0*t1610+t1133+0.5D0*t1376+0.5D0*t1385
      v2rhoab(i) = t1842+t1847
      t1848 = t78*t500
      t1852 = t202*t97
      t1862 = t489*t200
      t1868 = t110*t80
      t1871 = t93*t88
      t1874 = t622*t97
      t1886 = t673*t500
      t1892 = t234*t78*t140
      t1896 = t132*t96*t94
      t1901 = rho*(0.7420670230672731D-1*t1848*t551
     #+0.2473556743557577D-1*t79*(0.2697586091519874D1*(
     #-0.1269364121954074D0*t209-0.2173290372447189D-1*t1852*t187
     #-0.8693161489788758D-1*t105*t193)*t115-0.2697586091519874D1*t613
     #*t497-0.2697586091519874D1*t1862*t220+0.5395172183039748D1*t618
     #*t220*t497-0.2697586091519874D1*t201*(-0.1712109500228824D0
     #*t1868*t205-0.3424219000457648D0*t1871*t175-0.1172527576309517D0
     #*t1874*t187-0.2345055152619034D0*t202*t193))*t225
     #-0.2473556743557577D-1*t224*t1886+0.1D1*t228*t504
     #-0.130605864926168D2*t1248*t1892+0.1364631356248995D2*t1896
     #*t1836*t140)
      s1 = 0.7420670230672731D-1*t1848*t718+0.2473556743557577D-1*t79*
     #(0.2697586091519874D1*(-0.1269364121954074D0*t367
     #-0.148092480894642D0*t370-0.2173290372447189D-1*t1852*t349
     #-0.8693161489788758D-1*t105*t352-0.1014202173808688D0*t105*t357)
     #*t115-0.2697586091519874D1*t1062*t497-0.2697586091519874D1*t1862
     #*t379+0.5395172183039748D1*t618*t379*t497-0.2697586091519874D1
     #*t201*(-0.1712109500228824D0*t1868*t364-0.3424219000457648D0
     #*t1871*t334-0.3994922167200589D0*t492*t339-0.1172527576309517D0
     #*t1874*t349-0.2345055152619034D0*t202*t352-0.273589767805554D0
     #*t202*t357))*t225-0.2473556743557577D-1*t383*t1886+0.1D1*t401
     #*t504-0.261211729852336D1*t400*t79*t509+0.1D1*t404*t504
      t1973 = s1-0.261211729852336D1*t782*t1892-0.2333333333333333D1
     #*t133*t339*t140+0.6094940363221173D1*t507*sigma*t987*t140+0.1D1
     #*t133*t423+0.1D1*t411*t84*(-0.522423459704672D1*t412*t331
     #+0.6965646129395627D1*t418)*t140-0.261211729852336D1*t507*sigma
     #*t1251*t140
      t1974 = rho*t1973
      v2rhoasigmaaa(i) = t1901+t503+t506-t511+t1974
      t1976 = 0.4947113487115154D-1*t502
      t1977 = 0.2D1*t505
      t1978 = 0.522423459704672D1*t510
      t1979 = 2.D0*t1974
      v2rhoasigmaab(i) = 2.D0*t1901+t1976+t1977-t1978+t1979
      v2rhoasigmabb(i) = v2rhoasigmaaa(i)
      t2023 = rho*(0.7420670230672731D-1*t1848*t1394
     #+0.2473556743557577D-1*t79*(0.2697586091519874D1*(
     #-0.1269364121954074D0*t463-0.2173290372447189D-1*t1852*t450
     #-0.8693161489788758D-1*t105*t453)*t115-0.2697586091519874D1
     #*t1517*t497-0.2697586091519874D1*t1862*t470+0.5395172183039748D1
     #*t618*t470*t497-0.2697586091519874D1*t201*(-0.1712109500228824D0
     #*t1868*t460-0.3424219000457648D0*t1871*t440-0.1172527576309517D0
     #*t1874*t450-0.2345055152619034D0*t202*t453))*t225
     #-0.2473556743557577D-1*t474*t1886+0.1D1*t477*t504
     #-0.130605864926168D2*t1423*t1892+0.1364631356248995D2*t1896*t693
     #*t437*t140)
      v2rhobsigmaaa(i) = t2023+t503+t506-t511+t1974
      v2rhobsigmaab(i) = 2.D0*t2023+t1976+t1977-t1978+t1979
      v2rhobsigmabb(i) = v2rhobsigmaaa(i)
      t2031 = t497**2
      t2042 = t500**2
      v2sigmaaa2(i) = rho*(0.2473556743557577D-1*t79*
     #(0.5862637881547585D-1*t93*t97*t100*t115-0.5395172183039748D1
     #*t1862*t497+0.5395172183039748D1*t618*t2031-0.1581497040888031D0
     #*t201*t110*t97*t100)*t225-0.2473556743557577D-1*t79*t2042*t673
     #-0.522423459704672D1*t507*t234*t140+0.6823156781244976D1*t132
     #*t190*sigma*t693*t140)
      v2sigmaaaab(i) = 2.D0*v2sigmaaa2(i)
      v2sigmaaabb(i) = v2sigmaaa2(i)
      v2sigmaab2(i) = 4.D0*v2sigmaaabb(i)
      v2sigmaabbb(i) = v2sigmaaaab(i)
      v2sigmabb2(i) = v2sigmaaabb(i)
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vrhob(i) = 0.0d0
      v2rhoa2(i) = 0.0d0
      v2rhob2(i) = 0.0d0
      v2rhoab(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      vsigmaab(i) = 0.0d0
      vsigmabb(i) = 0.0d0
      v2rhoasigmaaa(i) = 0.0d0
      v2rhoasigmaab(i) = 0.0d0
      v2rhoasigmabb(i) = 0.0d0
      v2rhobsigmaaa(i) = 0.0d0
      v2rhobsigmaab(i) = 0.0d0
      v2rhobsigmabb(i) = 0.0d0
      v2sigmaaa2(i) = 0.0d0
      v2sigmaab2(i) = 0.0d0
      v2sigmabb2(i) = 0.0d0
      v2sigmaaaab(i) = 0.0d0
      v2sigmaaabb(i) = 0.0d0
      v2sigmaabbb(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end
      
      
      subroutine rks_c_pw91
     & (ideriv,npt,rhoa1,sigmaaa1,
     &  zk,vrhoa,vsigmaaa,
     &  v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
c
c     J.P. Perdew, J.A. Chevary, S.H. Vosko, K.A. Jackson, 
c     M.R. Pederson, D.J. Singh, C. Fiolhais
c     Atoms, molecules, solids and surfaces: 
c     Applications of the generalized gradient approximation 
c     for exchange and correlation
c     Phys. Rev. B 46 (1992) 6671--6687
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit real*8 (a-h,o-z)
      integer ideriv,npt
      real*8 rhoa1(npt)
      real*8 sigmaaa1(npt)
      real*8 zk(npt),vrhoa(npt),vsigmaaa(npt)
      real*8 v2rhoa2(npt),v2rhoasigmaaa(npt),v2sigmaaa2(npt)
      parameter(tol=1.0d-20)
      
      if(ideriv.eq.0) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = 1/rho
      t3 = t2**(1.D0/3.D0)
      t6 = t2**(1.D0/6.D0)
      t9 = dsqrt(t2)
      t11 = t3**2
      t17 = dlog(1.D0+0.160818243221511D2/(0.598255043577108D1*t6
     #+0.2225569421150687D1*t3+0.8004286349993634D0*t9
     #+0.1897004325747559D0*t11))
      t18 = (1.D0+0.1325688999052018D0*t3)*t17
      t20 = rho**2
      t21 = rho**(1.D0/3.D0)
      t23 = 1/t21/t20
      t27 = dexp(0.2513869963240348D1*t18)
      t28 = t27-1.D0
      t29 = 1/t28
      t30 = sigma**2
      t32 = t20**2
      t33 = t21**2
      t35 = 1/t33/t32
      t42 = t28**2
      t52 = dlog(1.D0+0.2697586091519874D1*(0.634682060977037D-1*sigma
     #*t23+0.1086645186223595D-1*t29*t30*t35)/(1.D0
     #+0.1712109500228824D0*t29*sigma*t23+0.2931318940773793D-1/t42
     #*t30*t35))
      t70 = dexp(-0.261211729852336D1/t33/t20*sigma)
      zk(i) = rho*(-0.62182D-1*t18+0.2473556743557577D-1*t52+0.1D1*
     #(0.1D-2*(0.2568D1+0.1443307452126544D2*t3+0.2843543831490386D-2
     #*t11)/(1.D0+0.5411317332115466D1*t3+0.1816419932959077D0*t11
     #+0.1763993811759022D-1*t2)-0.1853571428571429D-2)*sigma*t23*t70)
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = 1/rho
      t3 = t2**(1.D0/3.D0)
      t5 = 1.D0+0.1325688999052018D0*t3
      t6 = t2**(1.D0/6.D0)
      t9 = dsqrt(t2)
      t11 = t3**2
      t13 = 0.598255043577108D1*t6+0.2225569421150687D1*t3
     #+0.8004286349993634D0*t9+0.1897004325747559D0*t11
      t16 = 1.D0+0.160818243221511D2/t13
      t17 = dlog(t16)
      t18 = t5*t17
      t19 = 0.62182D-1*t18
      t20 = rho**2
      t21 = rho**(1.D0/3.D0)
      t23 = 1/t21/t20
      t27 = dexp(0.2513869963240348D1*t18)
      t28 = t27-1.D0
      t29 = 1/t28
      t30 = sigma**2
      t31 = t29*t30
      t32 = t20**2
      t33 = t21**2
      t35 = 1/t33/t32
      t38 = 0.634682060977037D-1*sigma*t23+0.1086645186223595D-1*t31*t35
      t39 = t29*sigma
      t42 = t28**2
      t43 = 1/t42
      t44 = t43*t30
      t47 = 1.D0+0.1712109500228824D0*t39*t23+0.2931318940773793D-1
     #*t44*t35
      t48 = 1/t47
      t51 = 1.D0+0.2697586091519874D1*t38*t48
      t52 = dlog(t51)
      t53 = 0.2473556743557577D-1*t52
      t56 = 0.2568D1+0.1443307452126544D2*t3+0.2843543831490386D-2*t11
      t60 = 1.D0+0.5411317332115466D1*t3+0.1816419932959077D0*t11
     #+0.1763993811759022D-1*t2
      t61 = 1/t60
      t64 = 0.1D-2*t56*t61-0.1853571428571429D-2
      t65 = t64*sigma
      t70 = dexp(-0.261211729852336D1/t33/t20*sigma)
      t71 = t23*t70
      t73 = 0.1D1*t65*t71
      zk(i) = rho*(-t19+t53+t73)
      t76 = 1/t20
      t77 = 1/t11*t76
      t78 = t77*t17
      t80 = t13**2
      t83 = t6**2
      t84 = t83**2
      t94 = 1/t3*t76
      t99 = t5/t80*(-0.99709173929518D0/t84/t6*t76
     #-0.7418564737168958D0*t77-0.4002143174996817D0/t9*t76
     #-0.1264669550498372D0*t94)/t16
      t103 = 1/t21/t20/rho
      t108 = -0.1110869918438343D0*t78-0.4042761511756372D2*t99
      t110 = t35*t108*t27
      t113 = t32*rho
      t115 = 1/t33/t113
      t121 = t47**2
      t123 = t38/t121
      t124 = t43*sigma
      t142 = 1/t51
      t150 = t60**2
      vrhoa(i) = -t19+t53+t73+rho*(0.2747799777968419D-2*t78+0.1D1*t99
     #+0.2473556743557577D-1*(0.2697586091519874D1*(
     #-0.148092480894642D0*sigma*t103-0.1086645186223595D-1*t44*t110
     #-0.5071010869043442D-1*t31*t115)*t48-0.2697586091519874D1*t123*(
     #-0.1712109500228824D0*t124*t23*t108*t27-0.3994922167200589D0*t39
     #*t103-0.5862637881547585D-1/t42/t28*t30*t110-0.136794883902777D0
     #*t44*t115))*t142+0.1D1*(0.1D-2*(-0.4811024840421814D1*t77
     #-0.1895695887660258D-2*t94)*t61-0.1D-2*t56/t150*(
     #-0.1803772444038489D1*t77-0.1210946621972718D0*t94
     #-0.1763993811759022D-1*t76))*sigma*t71-0.2333333333333333D1*t65
     #*t103*t70+0.6965646129395627D1*t64*t30/t32/t20*t70)
      vsigmaaa(i) = 4.D0*rho*(0.2473556743557577D-1*
     #(0.2697586091519874D1*(0.634682060977037D-1*t23
     #+0.2173290372447189D-1*t39*t35)*t48-0.2697586091519874D1*t123*
     #(0.1712109500228824D0*t29*t23+0.5862637881547585D-1*t124*t35))
     #*t142+0.1D1*t64*t23*t70-0.261211729852336D1*t65/t113*t70)
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = 1/rho
      t3 = t2**(1.D0/3.D0)
      t5 = 1.D0+0.1325688999052018D0*t3
      t6 = t2**(1.D0/6.D0)
      t9 = dsqrt(t2)
      t11 = t3**2
      t13 = 0.598255043577108D1*t6+0.2225569421150687D1*t3
     #+0.8004286349993634D0*t9+0.1897004325747559D0*t11
      t16 = 1.D0+0.160818243221511D2/t13
      t17 = dlog(t16)
      t18 = t5*t17
      t19 = 0.62182D-1*t18
      t20 = rho**2
      t21 = rho**(1.D0/3.D0)
      t23 = 1/t21/t20
      t27 = dexp(0.2513869963240348D1*t18)
      t28 = t27-1.D0
      t29 = 1/t28
      t30 = sigma**2
      t31 = t29*t30
      t32 = t20**2
      t33 = t21**2
      t35 = 1/t33/t32
      t38 = 0.634682060977037D-1*sigma*t23+0.1086645186223595D-1*t31*t35
      t39 = t29*sigma
      t42 = t28**2
      t43 = 1/t42
      t44 = t43*t30
      t47 = 1.D0+0.1712109500228824D0*t39*t23+0.2931318940773793D-1
     #*t44*t35
      t48 = 1/t47
      t51 = 1.D0+0.2697586091519874D1*t38*t48
      t52 = dlog(t51)
      t53 = 0.2473556743557577D-1*t52
      t56 = 0.2568D1+0.1443307452126544D2*t3+0.2843543831490386D-2*t11
      t60 = 1.D0+0.5411317332115466D1*t3+0.1816419932959077D0*t11
     #+0.1763993811759022D-1*t2
      t61 = 1/t60
      t64 = 0.1D-2*t56*t61-0.1853571428571429D-2
      t65 = t64*sigma
      t70 = dexp(-0.261211729852336D1/t33/t20*sigma)
      t71 = t23*t70
      t73 = 0.1D1*t65*t71
      zk(i) = rho*(-t19+t53+t73)
      t75 = 1/t11
      t76 = 1/t20
      t77 = t75*t76
      t78 = t77*t17
      t80 = t13**2
      t81 = 1/t80
      t82 = t5*t81
      t83 = t6**2
      t84 = t83**2
      t85 = t84*t6
      t86 = 1/t85
      t90 = 1/t9
      t93 = 1/t3
      t94 = t93*t76
      t96 = -0.99709173929518D0*t86*t76-0.7418564737168958D0*t77
     #-0.4002143174996817D0*t90*t76-0.1264669550498372D0*t94
      t97 = 1/t16
      t99 = t82*t96*t97
      t101 = t20*rho
      t103 = 1/t21/t101
      t108 = -0.1110869918438343D0*t78-0.4042761511756372D2*t99
      t110 = t35*t108*t27
      t113 = t32*rho
      t115 = 1/t33/t113
      t118 = -0.148092480894642D0*sigma*t103-0.1086645186223595D-1*t44
     #*t110-0.5071010869043442D-1*t31*t115
      t121 = t47**2
      t122 = 1/t121
      t123 = t38*t122
      t124 = t43*sigma
      t132 = 1/t42/t28
      t133 = t132*t30
      t138 = -0.1712109500228824D0*t124*t23*t108*t27
     #-0.3994922167200589D0*t39*t103-0.5862637881547585D-1*t133*t110
     #-0.136794883902777D0*t44*t115
      t141 = 0.2697586091519874D1*t118*t48-0.2697586091519874D1*t123
     #*t138
      t142 = 1/t51
      t143 = t141*t142
      t147 = -0.4811024840421814D1*t77-0.1895695887660258D-2*t94
      t150 = t60**2
      t151 = 1/t150
      t152 = t56*t151
      t156 = -0.1803772444038489D1*t77-0.1210946621972718D0*t94
     #-0.1763993811759022D-1*t76
      t159 = 0.1D-2*t147*t61-0.1D-2*t152*t156
      t160 = t159*sigma
      t161 = t160*t71
      t163 = t103*t70
      t164 = t65*t163
      t166 = t64*t30
      t167 = t32*t20
      t169 = 1/t167*t70
      t170 = t166*t169
      vrhoa(i) = -t19+t53+t73+rho*(0.2747799777968419D-2*t78+0.1D1*t99
     #+0.2473556743557577D-1*t143+0.1D1*t161-0.2333333333333333D1*t164
     #+0.6965646129395627D1*t170)
      t177 = 0.634682060977037D-1*t23+0.2173290372447189D-1*t39*t35
      t184 = 0.1712109500228824D0*t29*t23+0.5862637881547585D-1*t124*t35
      t187 = 0.2697586091519874D1*t177*t48-0.2697586091519874D1*t123
     #*t184
      t188 = t187*t142
      t191 = t64*t23*t70
      t193 = 1/t113
      t194 = t193*t70
      t195 = t65*t194
      vsigmaaa(i) = 4.D0*rho*(0.2473556743557577D-1*t188+0.1D1*t191
     #-0.261211729852336D1*t195)
      t215 = dlog(1.D0+0.2960857464321668D2/(0.8157414703487641D1*t6
     #+0.2247591863577616D1*t3+0.4300972471276643D0*t9
     #+0.1911512595127338D0*t11))
      t217 = (1.D0+0.6901399211255825D-1*t3)*t215*t76
      t219 = t52*t76
      t222 = 1/t21/t32
      t223 = sigma*t222
      t228 = 0.1365402220407279D1*t217-0.1675913308826898D1*t18*t76
      t230 = t35*t228*t27
      t234 = 1/t33/t167
      t235 = t31*t234
      t244 = t39*t222
      t248 = t44*t234
      t258 = t64*t222*sigma*t70
      t260 = t32*t101
      t263 = t166/t260*t70
      t270 = -t228
      t272 = t35*t270*t27
      t299 = 1/t32
      t302 = 1/t101
      t307 = 1/t11/t2*t299
      t309 = t75*t302
      t319 = 1/t3/t2*t299
      t321 = t93*t302
      t325 = t82*(-0.8309097827459833D0/t85/t2*t299
     #+0.199418347859036D1*t86*t302-0.4945709824779306D0*t307
     #+0.1483712947433792D1*t309-0.2001071587498409D0/t9/t2*t299
     #+0.8004286349993634D0*t90*t302-0.4215565168327908D-1*t319
     #+0.2529339100996745D0*t321)*t97
      t327 = t80**2
      t330 = t96**2
      t331 = t16**2
      t334 = t5/t327*t330/t331
      t336 = t309*t17
      t338 = t307*t17
      t342 = t77*t81*t96*t97
      t348 = t5/t80/t13*t330*t97
      t352 = t32**2
      t374 = t156**2
      t391 = t141**2
      t392 = t51**2
      t393 = 1/t392
      t397 = t115*t108*t27
      t406 = -0.4042761511756372D2*t325-0.6501498040842D3*t334
     #+0.2221739836876686D0*t336-0.740579945625562D-1*t338
     #+0.3572962974617552D1*t342+0.8085523023512745D2*t348
      t408 = t35*t406*t27
      t411 = t108**2
      t412 = t35*t411
      t414 = dexp(0.5027739926480695D1*t18)
      t415 = t412*t414
      t418 = t412*t27
      t426 = t118*t122
      t431 = t38/t121/t47
      t432 = t138**2
      t443 = t132*sigma
      t444 = t23*t411
      t459 = t42**2
      s1 = 0.1D1*t325+0.160818243221511D2*t334-0.5495599555936838D-2
     #*t336+0.1831866518645613D-2*t338-0.8837926660346786D-1*t342
     #-0.2D1*t348+0.4852022599996427D2*t64*t30*sigma/t33/t352/rho*t70
      s2 = s1-0.5804705107829689D2*t263+0.7777777777777778D1*t258
     #+0.1D1*(0.1D-2*(-0.3207349893614542D1*t307+0.9622049680843627D1
     #*t309-0.6318986292200858D-3*t319+0.3791391775320515D-2*t321)*t61
     #-0.2D-2*t147*t151*t156+0.2D-2*t56/t150/t60*t374-0.1D-2*t152*(
     #-0.1202514962692326D1*t307+0.3607544888076978D1*t309
     #-0.4036488739909061D-1*t319+0.2421893243945437D0*t321
     #+0.3527987623518044D-1*t302))*sigma*t71
      t473 = s2-0.4666666666666667D1*t160*t163-0.2473556743557577D-1
     #*t391*t393+0.2473556743557577D-1*(0.2697586091519874D1*
     #(0.1014202173808688D0*t44*t397-0.1086645186223595D-1*t44*t408
     #+0.2173290372447189D-1*t133*t415-0.1086645186223595D-1*t44*t418
     #+0.2873572825791284D0*t235+0.4936416029821399D0*t223)*t48
     #-0.5395172183039748D1*t426*t138+0.5395172183039748D1*t431*t432
     #-0.2697586091519874D1*t123*(-0.1712109500228824D0*t124*t23*t406
     #*t27+0.7989844334401178D0*t124*t103*t108*t27
     #+0.3424219000457648D0*t443*t444*t414-0.1712109500228824D0*t124
     #*t444*t27+0.1331640722400196D1*t244+0.547179535611108D0*t133
     #*t397-0.5862637881547585D-1*t133*t408-0.5862637881547585D-1*t133
     #*t418+0.7751710087824029D0*t248+0.1758791364464276D0/t459*t30
     #*t415))*t142+0.1393129225879125D2*t159*t30*t169
      s1 = 0.2786258451758251D2*t170-0.9333333333333333D1*t164
     #+0.9894226974230308D-1*t143+0.4D1*t161
      s2 = s1+0.1099119911187368D-1*t78+0.4D1*t99
      s3 = s2+rho*(-0.3377399869956914D-1*t217+0.1649037829038385D-1
     #*t219+0.2473556743557577D-1*(0.2697586091519874D1*(
     #-0.2820809159897942D-1*t223-0.1086645186223595D-1*t44*t230
     #-0.9659068321987509D-2*t235)*t48-0.2697586091519874D1*t123*(
     #-0.1712109500228824D0*t124*t23*t228*t27-0.760937555657255D-1
     #*t244-0.5862637881547585D-1*t133*t230-0.2605616836243371D-1*t248
     #))*t142+0.2222222222222222D0*t258-0.1160941021565938D1*t263)
      v2rhoa2(i) = s3+rho*(0.3377399869956914D-1*t217
     #-0.1649037829038385D-1*t219+0.2473556743557577D-1*
     #(0.2697586091519874D1*(0.2820809159897942D-1*t223
     #-0.1086645186223595D-1*t44*t272+0.9659068321987509D-2*t235)*t48
     #-0.2697586091519874D1*t123*(-0.1712109500228824D0*t124*t23*t270
     #*t27+0.760937555657255D-1*t244-0.5862637881547585D-1*t133*t272
     #+0.2605616836243371D-1*t248))*t142-0.2222222222222222D0*t258
     #+0.1160941021565938D1*t263)+2.D0*rho*t473
      t489 = t177*t122
      v2rhoasigmaaa(i) = 0.9894226974230308D-1*t188+0.4D1*t191
     #-0.1044846919409344D2*t195+4.D0*rho*(0.2473556743557577D-1*
     #(0.2697586091519874D1*(-0.148092480894642D0*t103
     #-0.2173290372447189D-1*t124*t110-0.1014202173808688D0*t39*t115)
     #*t48-0.2697586091519874D1*t426*t184-0.2697586091519874D1*t489
     #*t138+0.5395172183039748D1*t431*t138*t184-0.2697586091519874D1
     #*t123*(-0.1712109500228824D0*t43*t23*t108*t27
     #-0.3994922167200589D0*t29*t103-0.1172527576309517D0*t443*t110
     #-0.273589767805554D0*t124*t115))*t142-0.2473556743557577D-1*t141
     #*t393*t187+0.1D1*t159*t23*t70-0.261211729852336D1*t160*t194
     #-0.2333333333333333D1*t64*t103*t70+0.2002623262201243D2*t65*t169
     #-0.181950847499866D2*t166/t33/t352*t70)
      t537 = t184**2
      t546 = t187**2
      v2sigmaaa2(i) = 16.D0*rho*(0.2473556743557577D-1*
     #(0.5862637881547585D-1*t29*t35*t48-0.5395172183039748D1*t489
     #*t184+0.5395172183039748D1*t431*t537-0.1581497040888031D0*t123
     #*t43*t35)*t142-0.2473556743557577D-1*t546*t393
     #-0.522423459704672D1*t64*t193*t70+0.6823156781244976D1*t65/t33
     #/t260*t70)
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      v2rhoa2(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      v2rhoasigmaaa(i) = 0.0d0
      v2sigmaaa2(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end


      subroutine uks_c_pbe
     & (ideriv,npt,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,
     &  zk,vrhoa,vrhob,vsigmaaa,vsigmabb,vsigmaab,
     &  v2rhoa2,v2rhob2,v2rhoab,
     &  v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb,
     &  v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,
     &  v2sigmaaa2,v2sigmaaaab,v2sigmaaabb,
     &  v2sigmaab2,v2sigmaabbb,v2sigmabb2)
c
c     J.P. Perdew, K. Burke, and M. Ernzerhof
c     Generalized gradient approximation made simple
c     Phys. Rev. Lett. 77 (1996) 3865-3868
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit real*8 (a-h,o-z)
      integer ideriv,npt
      real*8 rhoa1(npt),rhob1(npt)
      real*8 sigmaaa1(npt),sigmabb1(npt),sigmaab1(npt)
      real*8 zk(npt),vrhoa(npt),vrhob(npt)
      real*8 vsigmaaa(npt),vsigmabb(npt),vsigmaab(npt)
      real*8 v2rhoa2(npt),v2rhob2(npt),v2rhoab(npt)
      real*8 v2rhoasigmaaa(npt),v2rhoasigmaab(npt)
      real*8 v2rhoasigmabb(npt),v2rhobsigmabb(npt)
      real*8 v2rhobsigmaab(npt),v2rhobsigmaaa(npt)
      real*8 v2sigmaaa2(npt),v2sigmaaaab(npt),v2sigmaaabb(npt)
      real*8 v2sigmaab2(npt),v2sigmaabbb(npt),v2sigmabb2(npt)
      parameter(tol=1.0d-20)
      
      if (ideriv.eq.0) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmabb
      t2 = 1/rhob
      t3 = t2**(1.D0/3.D0)
      t6 = t2**(1.D0/6.D0)
      t9 = dsqrt(t2)
      t11 = t3**2
      t17 = dlog(1.D0+0.3216395899738507D2/(0.1112037486309468D2*t6
     #+0.3844746237447211D1*t3+0.1644733775567609D1*t9
     #+0.2405871291288192D0*t11))
      t18 = (1.D0+0.1274696188700087D0*t3)*t17
      t20 = rhob**2
      t21 = rhob**(1.D0/3.D0)
      t23 = 1/t21/t20
      t26 = dexp(0.2000000587336264D1*t18)
      t27 = t26-1.D0
      t31 = 0.2162211495206379D0/t27*sigmabb*t23
      t33 = t27**2
      t35 = sigmabb**2
      t37 = t20**2
      t38 = t21**2
      t49 = dlog(1.D0+0.2162211495206379D0*sigmabb*t23*(1.D0+t31)/
     #(1.D0+t31+0.4675158550002605D-1/t33*t35/t38/t37))
      zk(i) = rhob*(-0.310907D-1*t18+0.1554534543482745D-1*t49)
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = 1/rhoa
      t3 = t2**(1.D0/3.D0)
      t6 = t2**(1.D0/6.D0)
      t9 = dsqrt(t2)
      t11 = t3**2
      t17 = dlog(1.D0+0.3216395899738507D2/(0.1112037486309468D2*t6
     #+0.3844746237447211D1*t3+0.1644733775567609D1*t9
     #+0.2405871291288192D0*t11))
      t18 = (1.D0+0.1274696188700087D0*t3)*t17
      t20 = rhoa**2
      t21 = rhoa**(1.D0/3.D0)
      t23 = 1/t21/t20
      t26 = dexp(0.2000000587336264D1*t18)
      t27 = t26-1.D0
      t31 = 0.2162211495206379D0/t27*sigmaaa*t23
      t33 = t27**2
      t35 = sigmaaa**2
      t37 = t20**2
      t38 = t21**2
      t49 = dlog(1.D0+0.2162211495206379D0*sigmaaa*t23*(1.D0+t31)/
     #(1.D0+t31+0.4675158550002605D-1/t33*t35/t38/t37))
      zk(i) = rhoa*(-0.310907D-1*t18+0.1554534543482745D-1*t49)
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = 1/rho
      t5 = t4**(1.D0/3.D0)
      t8 = t4**(1.D0/6.D0)
      t11 = dsqrt(t4)
      t13 = t5**2
      t19 = dlog(1.D0+0.1608197949869254D2/(0.598255043577108D1*t8
     #+0.2225569421150687D1*t5+0.8004286349993634D0*t11
     #+0.1897004325747559D0*t13))
      t21 = 0.621814D-1*(1.D0+0.1325688999052018D0*t5)*t19
      t32 = dlog(1.D0+0.2960874997779344D2/(0.8157414703487641D1*t8
     #+0.2247591863577616D1*t5+0.4300972471276643D0*t11
     #+0.1911512595127338D0*t13))
      t35 = rhoa-1.D0*rhob
      t36 = t35*t4
      t37 = 1.D0+t36
      t38 = t37**(1.D0/3.D0)
      t41 = 1.D0-1.D0*t36
      t42 = t41**(1.D0/3.D0)
      t44 = t38*t37+t42*t41-2.D0
      t45 = t35**2
      t46 = t45**2
      t47 = rho**2
      t48 = t47**2
      t50 = t46/t48
      t55 = 0.37995525D-1*(1.D0+0.6901399211255825D-1*t5)*t32*t44*
     #(1.D0-1.D0*t50)
      t66 = dlog(1.D0+0.3216395899738507D2/(0.1112037486309468D2*t8
     #+0.3844746237447211D1*t5+0.1644733775567609D1*t11
     #+0.2405871291288192D0*t13))
      t72 = 0.1923661050931536D1*(-0.310907D-1*(1.D0
     #+0.1274696188700087D0*t5)*t66+t21)*t44*t50
      t73 = t38**2
      t75 = t42**2
      t77 = 0.5D0*t73+0.5D0*t75
      t78 = t77**2
      t79 = t78*t77
      t80 = 1/t78
      t82 = rho**(1.D0/3.D0)
      t84 = 1/t82/t47
      t89 = dexp(-0.3216396844291482D2*(-t21+t55+t72)/t79)
      t90 = t89-1.D0
      t95 = 0.1362107888567592D0/t90*sigma*t80*t84
      t98 = t90**2
      t100 = sigma**2
      t102 = t78**2
      t104 = t82**2
      t116 = dlog(1.D0+0.1362107888567592D0*sigma*t80*t84*(1.D0+t95)/
     #(1.D0+t95+0.1855337900098064D-1/t98*t100/t102/t104/t48))
      zk(i) = rho*(-t21+t55+t72+0.310906908696549D-1*t79*t116)
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmabb
      t2 = 1/rhob
      t3 = t2**(1.D0/3.D0)
      t5 = 1.D0+0.1274696188700087D0*t3
      t6 = t2**(1.D0/6.D0)
      t9 = dsqrt(t2)
      t11 = t3**2
      t13 = 0.1112037486309468D2*t6+0.3844746237447211D1*t3
     #+0.1644733775567609D1*t9+0.2405871291288192D0*t11
      t16 = 1.D0+0.3216395899738507D2/t13
      t17 = dlog(t16)
      t18 = t5*t17
      t19 = 0.310907D-1*t18
      t20 = rhob**2
      t21 = rhob**(1.D0/3.D0)
      t23 = 1/t21/t20
      t24 = sigmabb*t23
      t26 = dexp(0.2000000587336264D1*t18)
      t27 = t26-1.D0
      t28 = 1/t27
      t29 = t28*sigmabb
      t31 = 0.2162211495206379D0*t29*t23
      t32 = 1.D0+t31
      t33 = t27**2
      t34 = 1/t33
      t35 = sigmabb**2
      t36 = t34*t35
      t37 = t20**2
      t38 = t21**2
      t40 = 1/t38/t37
      t41 = t36*t40
      t43 = 1.D0+t31+0.4675158550002605D-1*t41
      t44 = 1/t43
      t45 = t32*t44
      t48 = 1.D0+0.2162211495206379D0*t24*t45
      t49 = dlog(t48)
      t50 = 0.1554534543482745D-1*t49
      zk(i) = rhob*(-t19+t50)
      vrhoa(i) = 0.D0
      t54 = 1/t21/t20/rhob
      t60 = 1/t38/t37/rhob
      t62 = t28*t44
      t65 = t43**2
      t67 = t32/t65
      t78 = 1/t48
      t82 = 1/t20
      t83 = 1/t11*t82
      t84 = t83*t17
      t86 = t13**2
      t89 = t6**2
      t90 = t89**2
      t105 = t5/t86*(-0.1853395810515781D1/t90/t6*t82
     #-0.128158207914907D1*t83-0.8223668877838045D0/t9*t82
     #-0.1603914194192128D0/t3*t82)/t16
      t109 = -0.8497977086918237D-1*t84-0.6432793688582964D2*t105
      t114 = t34*sigmabb
      vrhob(i) = -t19+t50+0.1554534543482745D-1*rhob*(
     #-0.5045160155481551D0*sigmabb*t54*t45-0.1090870328333941D0*t35
     #*t60*t62-0.2162211495206379D0*t24*t67*(-0.5045160155481551D0*t29
     #*t54-0.2181740656667882D0*t36*t60))*t78+rhob*
     #(0.1321039893133927D-2*t84+0.1D1*t105+0.1554534543482745D-1*(
     #-0.4675158550002605D-1*t41*t109*t26*t44-0.2162211495206379D0*t24
     #*t67*(-0.2162211495206379D0*t114*t23*t109*t26
     #-0.935031710000521D-1/t33/t27*t35*t40*t109*t26))*t78)
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.1554534543482745D-1*rhob*(0.2162211495206379D0
     #*t23*t32*t44+0.4675158550002605D-1*sigmabb*t40*t62
     #-0.2162211495206379D0*t24*t67*(0.2162211495206379D0*t28*t23
     #+0.935031710000521D-1*t114*t40))*t78
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = 1/rhoa
      t3 = t2**(1.D0/3.D0)
      t5 = 1.D0+0.1274696188700087D0*t3
      t6 = t2**(1.D0/6.D0)
      t9 = dsqrt(t2)
      t11 = t3**2
      t13 = 0.1112037486309468D2*t6+0.3844746237447211D1*t3
     #+0.1644733775567609D1*t9+0.2405871291288192D0*t11
      t16 = 1.D0+0.3216395899738507D2/t13
      t17 = dlog(t16)
      t18 = t5*t17
      t19 = 0.310907D-1*t18
      t20 = rhoa**2
      t21 = rhoa**(1.D0/3.D0)
      t23 = 1/t21/t20
      t24 = sigmaaa*t23
      t26 = dexp(0.2000000587336264D1*t18)
      t27 = t26-1.D0
      t28 = 1/t27
      t29 = t28*sigmaaa
      t31 = 0.2162211495206379D0*t29*t23
      t32 = 1.D0+t31
      t33 = t27**2
      t34 = 1/t33
      t35 = sigmaaa**2
      t36 = t34*t35
      t37 = t20**2
      t38 = t21**2
      t40 = 1/t38/t37
      t41 = t36*t40
      t43 = 1.D0+t31+0.4675158550002605D-1*t41
      t44 = 1/t43
      t45 = t32*t44
      t48 = 1.D0+0.2162211495206379D0*t24*t45
      t49 = dlog(t48)
      t50 = 0.1554534543482745D-1*t49
      zk(i) = rhoa*(-t19+t50)
      t54 = 1/t21/t20/rhoa
      t60 = 1/t38/t37/rhoa
      t62 = t28*t44
      t65 = t43**2
      t67 = t32/t65
      t78 = 1/t48
      t82 = 1/t20
      t83 = 1/t11*t82
      t84 = t83*t17
      t86 = t13**2
      t89 = t6**2
      t90 = t89**2
      t105 = t5/t86*(-0.1853395810515781D1/t90/t6*t82
     #-0.128158207914907D1*t83-0.8223668877838045D0/t9*t82
     #-0.1603914194192128D0/t3*t82)/t16
      t109 = -0.8497977086918237D-1*t84-0.6432793688582964D2*t105
      t114 = t34*sigmaaa
      vrhoa(i) = -t19+t50+0.1554534543482745D-1*rhoa*(
     #-0.5045160155481551D0*sigmaaa*t54*t45-0.1090870328333941D0*t35
     #*t60*t62-0.2162211495206379D0*t24*t67*(-0.5045160155481551D0*t29
     #*t54-0.2181740656667882D0*t36*t60))*t78+rhoa*
     #(0.1321039893133927D-2*t84+0.1D1*t105+0.1554534543482745D-1*(
     #-0.4675158550002605D-1*t41*t109*t26*t44-0.2162211495206379D0*t24
     #*t67*(-0.2162211495206379D0*t114*t23*t109*t26
     #-0.935031710000521D-1/t33/t27*t35*t40*t109*t26))*t78)
      vrhob(i) = 0.D0
      vsigmaaa(i) = 0.1554534543482745D-1*rhoa*(0.2162211495206379D0
     #*t23*t32*t44+0.4675158550002605D-1*sigmaaa*t40*t62
     #-0.2162211495206379D0*t24*t67*(0.2162211495206379D0*t28*t23
     #+0.935031710000521D-1*t114*t40))*t78
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = 1/rho
      t5 = t4**(1.D0/3.D0)
      t7 = 1.D0+0.1325688999052018D0*t5
      t8 = t4**(1.D0/6.D0)
      t11 = dsqrt(t4)
      t13 = t5**2
      t15 = 0.598255043577108D1*t8+0.2225569421150687D1*t5
     #+0.8004286349993634D0*t11+0.1897004325747559D0*t13
      t18 = 1.D0+0.1608197949869254D2/t15
      t19 = dlog(t18)
      t21 = 0.621814D-1*t7*t19
      t23 = 1.D0+0.6901399211255825D-1*t5
      t28 = 0.8157414703487641D1*t8+0.2247591863577616D1*t5
     #+0.4300972471276643D0*t11+0.1911512595127338D0*t13
      t31 = 1.D0+0.2960874997779344D2/t28
      t32 = dlog(t31)
      t33 = t23*t32
      t35 = rhoa-1.D0*rhob
      t36 = t35*t4
      t37 = 1.D0+t36
      t38 = t37**(1.D0/3.D0)
      t41 = 1.D0-1.D0*t36
      t42 = t41**(1.D0/3.D0)
      t44 = t38*t37+t42*t41-2.D0
      t45 = t35**2
      t46 = t45**2
      t47 = rho**2
      t48 = t47**2
      t49 = 1/t48
      t50 = t46*t49
      t52 = 1.D0-1.D0*t50
      t55 = 0.37995525D-1*t33*t44*t52
      t57 = 1.D0+0.1274696188700087D0*t5
      t62 = 0.1112037486309468D2*t8+0.3844746237447211D1*t5
     #+0.1644733775567609D1*t11+0.2405871291288192D0*t13
      t65 = 1.D0+0.3216395899738507D2/t62
      t66 = dlog(t65)
      t69 = -0.310907D-1*t57*t66+t21
      t70 = t69*t44
      t72 = 0.1923661050931536D1*t70*t50
      t73 = t38**2
      t75 = t42**2
      t77 = 0.5D0*t73+0.5D0*t75
      t78 = t77**2
      t79 = t78*t77
      t80 = 1/t78
      t81 = sigma*t80
      t82 = rho**(1.D0/3.D0)
      t84 = 1/t82/t47
      t85 = -t21+t55+t72
      t86 = 1/t79
      t89 = dexp(-0.3216396844291482D2*t85*t86)
      t90 = t89-1.D0
      t91 = 1/t90
      t92 = t91*sigma
      t93 = t80*t84
      t95 = 0.1362107888567592D0*t92*t93
      t96 = 1.D0+t95
      t98 = t90**2
      t99 = 1/t98
      t100 = sigma**2
      t101 = t99*t100
      t102 = t78**2
      t103 = 1/t102
      t104 = t82**2
      t106 = 1/t104/t48
      t107 = t103*t106
      t110 = 1.D0+t95+0.1855337900098064D-1*t101*t107
      t111 = 1/t110
      t115 = 1.D0+0.1362107888567592D0*t81*t84*t96*t111
      t116 = dlog(t115)
      t118 = 0.310906908696549D-1*t79*t116
      zk(i) = rho*(-t21+t55+t72+t118)
      t124 = 0.1333333333333333D1*t38*t4-0.1333333333333333D1*t42*t4
      t127 = 0.37995525D-1*t33*t124*t52
      t128 = t45*t35
      t131 = t33*t44*t128*t49
      t132 = 0.1519821D0*t131
      t135 = 0.1923661050931536D1*t69*t124*t50
      t137 = t70*t128*t49
      t138 = 0.7694644203726145D1*t137
      t139 = t78*t116
      t140 = 1/t38
      t143 = 1/t42
      t146 = 0.3333333333333333D0*t140*t4-0.3333333333333333D0*t143*t4
      t150 = sigma*t86*t84
      t151 = t96*t111
      t155 = t99*sigma
      t156 = t155*t80
      t160 = t85*t103
      t163 = -0.3216396844291482D2*(t127-t132+t135+t138)*t86
     #+0.9649190532874446D2*t160*t146
      t167 = 0.1362107888567592D0*t156*t84*t163*t89
      t168 = t86*t84
      t171 = 0.2724215777135184D0*t92*t168*t146
      t177 = t81*t84
      t178 = t110**2
      t180 = t96/t178
      t184 = 1/t98/t90*t100*t103
      t191 = 1/t102/t77*t106
      t201 = 1/t115
      t207 = 1/t47
      t208 = 1/t13*t207
      t209 = t208*t19
      t210 = 0.2747773264188438D-2*t209
      t211 = t15**2
      t214 = t8**2
      t215 = t214**2
      t218 = 1/t215/t8*t207
      t222 = 1/t11*t207
      t225 = 1/t5*t207
      t230 = t7/t211*(-0.99709173929518D0*t218-0.7418564737168958D0
     #*t208-0.4002143174996817D0*t222-0.1264669550498372D0*t225)/t18
      t231 = 0.1D1*t230
      t235 = 0.87407428755417D-3*t208*t32*t44*t52
      t236 = t28**2
      t249 = 0.1125D1*t23/t236*(-0.135956911724794D1*t218
     #-0.7491972878592054D0*t208-0.2150486235638321D0*t222
     #-0.1274341730084892D0*t225)/t31*t44*t52
      t256 = -0.1333333333333333D1*t38*t35*t207+0.1333333333333333D1
     #*t42*t35*t207
      t259 = 0.37995525D-1*t33*t256*t52
      t261 = t48*rho
      t262 = 1/t261
      t265 = 0.1519821D0*t33*t44*t46*t262
      t268 = t62**2
      t285 = 0.1923661050931536D1*(0.1321039893133927D-2*t208*t66
     #+0.1D1*t57/t268*(-0.1853395810515781D1*t218-0.128158207914907D1
     #*t208-0.8223668877838045D0*t222-0.1603914194192128D0*t225)/t65
     #-0.2747773264188438D-2*t209-0.1D1*t230)*t44*t50
      t288 = 0.1923661050931536D1*t69*t256*t50
      t291 = 0.7694644203726145D1*t70*t46*t262
      t298 = -0.3333333333333333D0*t140*t35*t207+0.3333333333333333D0
     #*t143*t35*t207
      t306 = 1/t82/t47/rho
      t316 = -0.3216396844291482D2*(t210+t231-t235-t249+t259+t265+t285
     #+t288-t291)*t86+0.9649190532874446D2*t160*t298
      t320 = 0.1362107888567592D0*t156*t84*t316*t89
      t323 = 0.2724215777135184D0*t92*t168*t298
      t326 = 0.3178251739991049D0*t92*t80*t306
      t352 = t210+t231-t235-t249+t259+t265+t285+t288-t291
     #+0.9327207260896469D-1*t139*t298+0.310906908696549D-1*t79*(
     #-0.2724215777135184D0*t150*t151*t298-0.3178251739991049D0*t81
     #*t306*t96*t111+0.1362107888567592D0*t81*t84*(-t320-t323-t326)
     #*t111-0.1362107888567592D0*t177*t180*(-t320-t323-t326
     #-0.3710675800196129D-1*t184*t106*t316*t89-0.7421351600392257D-1
     #*t101*t191*t298-0.8658243533790967D-1*t101*t103/t104/t261))*t201
      t353 = rho*t352
      vrhoa(i) = rho*(t127-t132+t135+t138+0.9327207260896469D-1*t139
     #*t146+0.310906908696549D-1*t79*(-0.2724215777135184D0*t150*t151
     #*t146+0.1362107888567592D0*t81*t84*(-t167-t171)*t111
     #-0.1362107888567592D0*t177*t180*(-t167-t171
     #-0.3710675800196129D-1*t184*t106*t163*t89-0.7421351600392257D-1
     #*t101*t191*t146))*t201)-t21+t55+t72+t118+t353
      t354 = -t124
      t357 = 0.37995525D-1*t33*t354*t52
      t358 = 0.1519821D0*t131
      t361 = 0.1923661050931536D1*t69*t354*t50
      t362 = 0.7694644203726145D1*t137
      t363 = -t146
      t374 = -0.3216396844291482D2*(t357+t358+t361-t362)*t86
     #+0.9649190532874446D2*t160*t363
      t378 = 0.1362107888567592D0*t156*t84*t374*t89
      t381 = 0.2724215777135184D0*t92*t168*t363
      vrhob(i) = rho*(t357+t358+t361-t362+0.9327207260896469D-1*t139
     #*t363+0.310906908696549D-1*t79*(-0.2724215777135184D0*t150*t151
     #*t363+0.1362107888567592D0*t81*t84*(-t378-t381)*t111
     #-0.1362107888567592D0*t177*t180*(-t378-t381
     #-0.3710675800196129D-1*t184*t106*t374*t89-0.7421351600392257D-1
     #*t101*t191*t363))*t201)-t21+t55+t72+t118+t353
      t423 = rho*t79*(0.1362107888567592D0*t93*t151
     #+0.1855337900098064D-1*sigma*t103*t106*t91*t111
     #-0.1362107888567592D0*t177*t180*(0.1362107888567592D0*t91*t80
     #*t84+0.3710675800196129D-1*t155*t107))*t201
      vsigmaaa(i) = 0.310906908696549D-1*t423
      vsigmaab(i) = 0.6218138173930979D-1*t423
      vsigmabb(i) = vsigmaaa(i)
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vrhob(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      vsigmaab(i) = 0.0d0
      vsigmabb(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmabb
      t2 = 1/rhob
      t3 = t2**(1.D0/3.D0)
      t5 = 1.D0+0.1274696188700087D0*t3
      t6 = t2**(1.D0/6.D0)
      t9 = dsqrt(t2)
      t11 = t3**2
      t13 = 0.1112037486309468D2*t6+0.3844746237447211D1*t3
     #+0.1644733775567609D1*t9+0.2405871291288192D0*t11
      t16 = 1.D0+0.3216395899738507D2/t13
      t17 = dlog(t16)
      t18 = t5*t17
      t19 = 0.310907D-1*t18
      t20 = rhob**2
      t21 = rhob**(1.D0/3.D0)
      t23 = 1/t21/t20
      t24 = sigmabb*t23
      t26 = dexp(0.2000000587336264D1*t18)
      t27 = t26-1.D0
      t28 = 1/t27
      t29 = t28*sigmabb
      t31 = 0.2162211495206379D0*t29*t23
      t32 = 1.D0+t31
      t33 = t27**2
      t34 = 1/t33
      t35 = sigmabb**2
      t36 = t34*t35
      t37 = t20**2
      t38 = t21**2
      t40 = 1/t38/t37
      t41 = t36*t40
      t43 = 1.D0+t31+0.4675158550002605D-1*t41
      t44 = 1/t43
      t45 = t32*t44
      t48 = 1.D0+0.2162211495206379D0*t24*t45
      t49 = dlog(t48)
      t50 = 0.1554534543482745D-1*t49
      zk(i) = rhob*(-t19+t50)
      vrhoa(i) = 0.D0
      t52 = t20*rhob
      t54 = 1/t21/t52
      t55 = sigmabb*t54
      t60 = 1/t38/t37/rhob
      t61 = t35*t60
      t62 = t28*t44
      t65 = t43**2
      t66 = 1/t65
      t67 = t32*t66
      t70 = t36*t60
      t72 = -0.5045160155481551D0*t29*t54-0.2181740656667882D0*t70
      t73 = t67*t72
      t76 = -0.5045160155481551D0*t55*t45-0.1090870328333941D0*t61*t62
     #-0.2162211495206379D0*t24*t73
      t77 = rhob*t76
      t78 = 1/t48
      t81 = 1/t11
      t82 = 1/t20
      t83 = t81*t82
      t84 = t83*t17
      t86 = t13**2
      t87 = 1/t86
      t88 = t5*t87
      t89 = t6**2
      t90 = t89**2
      t91 = t90*t6
      t92 = 1/t91
      t96 = 1/t9
      t99 = 1/t3
      t102 = -0.1853395810515781D1*t92*t82-0.128158207914907D1*t83
     #-0.8223668877838045D0*t96*t82-0.1603914194192128D0*t99*t82
      t103 = 1/t16
      t105 = t88*t102*t103
      t109 = -0.8497977086918237D-1*t84-0.6432793688582964D2*t105
      t110 = t109*t26
      t111 = t110*t44
      t114 = t34*sigmabb
      t120 = 1/t33/t27
      t121 = t120*t35
      t126 = -0.2162211495206379D0*t114*t23*t109*t26
     #-0.935031710000521D-1*t121*t40*t109*t26
      t127 = t67*t126
      t130 = -0.4675158550002605D-1*t41*t111-0.2162211495206379D0*t24
     #*t127
      t131 = t130*t78
      vrhob(i) = -t19+t50+0.1554534543482745D-1*t77*t78+rhob*
     #(0.1321039893133927D-2*t84+0.1D1*t105+0.1554534543482745D-1*t131)
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      t135 = t23*t32
      t138 = sigmabb*t40
      t145 = 0.2162211495206379D0*t28*t23+0.935031710000521D-1*t114*t40
      t149 = 0.2162211495206379D0*t135*t44+0.4675158550002605D-1*t138
     #*t62-0.2162211495206379D0*t24*t67*t145
      vsigmabb(i) = 0.1554534543482745D-1*rhob*t149*t78
      v2rhoa2(i) = 0.D0
      v2rhoab(i) = 0.D0
      t155 = 1/t21/t37
      t161 = 1/t38/t37/t20
      t167 = t28*t66
      t172 = 1/t65/t43
      t173 = t32*t172
      t174 = t72**2
      t190 = t76**2
      t192 = t48**2
      t193 = 1/t192
      t199 = t70*t111
      t204 = 0.4675158550002605D-1*t41*t110*t66*t72
      t206 = 0.5045160155481551D0*t55*t127
      t209 = 0.1090870328333941D0*t61*t167*t126
      t214 = 0.4324422990412758D0*t24*t32*t172*t126*t72
      t226 = 0.2162211495206379D0*t24*t67*(0.5045160155481551D0*t114
     #*t54*t109*t26+0.4363481313335765D0*t121*t60*t109*t26)
      t230 = t130*t193
      t244 = 1/t37
      t245 = 1/t11/t2*t244
      t246 = t245*t17
      t248 = 1/t52
      t249 = t81*t248
      t250 = t249*t17
      t254 = t83*t87*t102*t103
      t259 = t102**2
      t261 = t5/t86/t13*t259*t103
      t285 = t88*(-0.1544496508763151D1/t91/t2*t244
     #+0.3706791621031562D1*t92*t248-0.854388052766047D0*t245
     #+0.2563164158298141D1*t249-0.4111834438919023D0/t9/t2*t244
     #+0.1644733775567609D1*t96*t248-0.5346380647307093D-1/t3/t2*t244
     #+0.3207828388384256D0*t99*t248)*t103
      t287 = t86**2
      t290 = t16**2
      t293 = t5/t287*t259/t290
      t297 = t109**2
      t298 = t26**2
      t309 = -0.5665318057945491D-1*t246+0.1699595417383647D0*t250
     #+0.546657173168712D1*t254+0.1286558737716593D3*t261
     #-0.6432793688582964D2*t285-0.2069041124382199D4*t293
      t322 = t126**2
      t327 = t23*t297
      t338 = t33**2
      t341 = t40*t297
      t359 = t130**2
      s1 = 0.310906908696549D-1*t76*t78+0.1554534543482745D-1*rhob*
     #(0.1681720051827184D1*sigmabb*t155*t45+0.8726962626671529D0*t35
     #*t161*t62+0.100903203109631D1*t55*t73+0.2181740656667882D0*t61
     #*t167*t72+0.4324422990412758D0*t24*t173*t174
     #-0.2162211495206379D0*t24*t67*(0.1681720051827184D1*t29*t155
     #+0.1236319705445133D1*t36*t161))*t78-0.1554534543482745D-1*rhob
     #*t190*t193+0.2642079786267853D-2*t84+0.2D1*t105
      s2 = s1+0.310906908696549D-1*t131+rhob*(0.1554534543482745D-1*
     #(0.2181740656667882D0*t199+t204+t206+t209+t214-t226)*t78
     #-0.1554534543482745D-1*t230*t76)
      s3 = s2+0.1554534543482745D-1*rhob*(0.2181740656667882D0*t199
     #+t206+t209+t204+t214-t226)*t78
      v2rhob2(i) = s3-0.1554534543482745D-1*t77*t230+rhob*
     #(0.8806932620892844D-3*t246-0.2642079786267853D-2*t250
     #-0.8497974591333914D-1*t254-0.2D1*t261+0.1D1*t285
     #+0.3216395899738507D2*t293+0.1554534543482745D-1*
     #(0.935031710000521D-1*t35*t40*t120*t297*t298*t44
     #-0.4675158550002605D-1*t41*t309*t26*t44-0.4675158550002605D-1
     #*t41*t297*t26*t44+0.935031710000521D-1*t41*t110*t66*t126
     #+0.4324422990412758D0*t24*t173*t322-0.2162211495206379D0*t24*t67
     #*(0.4324422990412758D0*t120*sigmabb*t327*t298
     #-0.2162211495206379D0*t114*t23*t309*t26-0.2162211495206379D0
     #*t114*t327*t26+0.2805095130001563D0/t338*t35*t341*t298
     #-0.935031710000521D-1*t121*t40*t309*t26-0.935031710000521D-1
     #*t121*t341*t26))*t78-0.1554534543482745D-1*t359*t193)
      v2sigmaaa2(i) = 0.D0
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      t373 = t145**2
      t387 = t149**2
      v2sigmabb2(i) = 0.1554534543482745D-1*rhob*(0.935031710000521D-1
     #*t40*t28*t44-0.4324422990412758D0*t135*t66*t145
     #-0.935031710000521D-1*t138*t167*t145+0.4324422990412758D0*t24
     #*t173*t373-0.2021736311745604D-1*sigmabb/t37/t52*t67*t34)*t78
     #-0.1554534543482745D-1*rhob*t387*t193
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = 1/rhoa
      t3 = t2**(1.D0/3.D0)
      t5 = 1.D0+0.1274696188700087D0*t3
      t6 = t2**(1.D0/6.D0)
      t9 = dsqrt(t2)
      t11 = t3**2
      t13 = 0.1112037486309468D2*t6+0.3844746237447211D1*t3
     #+0.1644733775567609D1*t9+0.2405871291288192D0*t11
      t16 = 1.D0+0.3216395899738507D2/t13
      t17 = dlog(t16)
      t18 = t5*t17
      t19 = 0.310907D-1*t18
      t20 = rhoa**2
      t21 = rhoa**(1.D0/3.D0)
      t23 = 1/t21/t20
      t24 = sigmaaa*t23
      t26 = dexp(0.2000000587336264D1*t18)
      t27 = t26-1.D0
      t28 = 1/t27
      t29 = t28*sigmaaa
      t31 = 0.2162211495206379D0*t29*t23
      t32 = 1.D0+t31
      t33 = t27**2
      t34 = 1/t33
      t35 = sigmaaa**2
      t36 = t34*t35
      t37 = t20**2
      t38 = t21**2
      t40 = 1/t38/t37
      t41 = t36*t40
      t43 = 1.D0+t31+0.4675158550002605D-1*t41
      t44 = 1/t43
      t45 = t32*t44
      t48 = 1.D0+0.2162211495206379D0*t24*t45
      t49 = dlog(t48)
      t50 = 0.1554534543482745D-1*t49
      zk(i) = rhoa*(-t19+t50)
      t52 = t20*rhoa
      t54 = 1/t21/t52
      t55 = sigmaaa*t54
      t60 = 1/t38/t37/rhoa
      t61 = t35*t60
      t62 = t28*t44
      t65 = t43**2
      t66 = 1/t65
      t67 = t32*t66
      t70 = t36*t60
      t72 = -0.5045160155481551D0*t29*t54-0.2181740656667882D0*t70
      t73 = t67*t72
      t76 = -0.5045160155481551D0*t55*t45-0.1090870328333941D0*t61*t62
     #-0.2162211495206379D0*t24*t73
      t77 = rhoa*t76
      t78 = 1/t48
      t81 = 1/t11
      t82 = 1/t20
      t83 = t81*t82
      t84 = t83*t17
      t86 = t13**2
      t87 = 1/t86
      t88 = t5*t87
      t89 = t6**2
      t90 = t89**2
      t91 = t90*t6
      t92 = 1/t91
      t96 = 1/t9
      t99 = 1/t3
      t102 = -0.1853395810515781D1*t92*t82-0.128158207914907D1*t83
     #-0.8223668877838045D0*t96*t82-0.1603914194192128D0*t99*t82
      t103 = 1/t16
      t105 = t88*t102*t103
      t109 = -0.8497977086918237D-1*t84-0.6432793688582964D2*t105
      t110 = t109*t26
      t111 = t110*t44
      t114 = t34*sigmaaa
      t120 = 1/t33/t27
      t121 = t120*t35
      t126 = -0.2162211495206379D0*t114*t23*t109*t26
     #-0.935031710000521D-1*t121*t40*t109*t26
      t127 = t67*t126
      t130 = -0.4675158550002605D-1*t41*t111-0.2162211495206379D0*t24
     #*t127
      t131 = t130*t78
      vrhoa(i) = -t19+t50+0.1554534543482745D-1*t77*t78+rhoa*
     #(0.1321039893133927D-2*t84+0.1D1*t105+0.1554534543482745D-1*t131)
      vrhob(i) = 0.D0
      t135 = t23*t32
      t138 = sigmaaa*t40
      t145 = 0.2162211495206379D0*t28*t23+0.935031710000521D-1*t114*t40
      t149 = 0.2162211495206379D0*t135*t44+0.4675158550002605D-1*t138
     #*t62-0.2162211495206379D0*t24*t67*t145
      vsigmaaa(i) = 0.1554534543482745D-1*rhoa*t149*t78
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      t155 = 1/t21/t37
      t161 = 1/t38/t37/t20
      t167 = t28*t66
      t172 = 1/t65/t43
      t173 = t32*t172
      t174 = t72**2
      t190 = t76**2
      t192 = t48**2
      t193 = 1/t192
      t199 = t70*t111
      t204 = 0.4675158550002605D-1*t41*t110*t66*t72
      t206 = 0.5045160155481551D0*t55*t127
      t209 = 0.1090870328333941D0*t61*t167*t126
      t214 = 0.4324422990412758D0*t24*t32*t172*t126*t72
      t226 = 0.2162211495206379D0*t24*t67*(0.5045160155481551D0*t114
     #*t54*t109*t26+0.4363481313335765D0*t121*t60*t109*t26)
      t230 = t130*t193
      t244 = 1/t37
      t245 = 1/t11/t2*t244
      t246 = t245*t17
      t248 = 1/t52
      t249 = t81*t248
      t250 = t249*t17
      t254 = t83*t87*t102*t103
      t259 = t102**2
      t261 = t5/t86/t13*t259*t103
      t285 = t88*(-0.1544496508763151D1/t91/t2*t244
     #+0.3706791621031562D1*t92*t248-0.854388052766047D0*t245
     #+0.2563164158298141D1*t249-0.4111834438919023D0/t9/t2*t244
     #+0.1644733775567609D1*t96*t248-0.5346380647307093D-1/t3/t2*t244
     #+0.3207828388384256D0*t99*t248)*t103
      t287 = t86**2
      t290 = t16**2
      t293 = t5/t287*t259/t290
      t297 = t109**2
      t298 = t26**2
      t309 = -0.5665318057945491D-1*t246+0.1699595417383647D0*t250
     #+0.546657173168712D1*t254+0.1286558737716593D3*t261
     #-0.6432793688582964D2*t285-0.2069041124382199D4*t293
      t322 = t126**2
      t327 = t23*t297
      t338 = t33**2
      t341 = t40*t297
      t359 = t130**2
      s1 = 0.310906908696549D-1*t76*t78+0.1554534543482745D-1*rhoa*
     #(0.1681720051827184D1*sigmaaa*t155*t45+0.8726962626671529D0*t35
     #*t161*t62+0.100903203109631D1*t55*t73+0.2181740656667882D0*t61
     #*t167*t72+0.4324422990412758D0*t24*t173*t174
     #-0.2162211495206379D0*t24*t67*(0.1681720051827184D1*t29*t155
     #+0.1236319705445133D1*t36*t161))*t78-0.1554534543482745D-1*rhoa
     #*t190*t193+0.2642079786267853D-2*t84+0.2D1*t105
      s2 = s1+0.310906908696549D-1*t131+rhoa*(0.1554534543482745D-1*
     #(0.2181740656667882D0*t199+t204+t206+t209+t214-t226)*t78
     #-0.1554534543482745D-1*t230*t76)
      s3 = s2+0.1554534543482745D-1*rhoa*(0.2181740656667882D0*t199
     #+t206+t209+t204+t214-t226)*t78
      v2rhoa2(i) = s3-0.1554534543482745D-1*t77*t230+rhoa*
     #(0.8806932620892844D-3*t246-0.2642079786267853D-2*t250
     #-0.8497974591333914D-1*t254-0.2D1*t261+0.1D1*t285
     #+0.3216395899738507D2*t293+0.1554534543482745D-1*
     #(0.935031710000521D-1*t35*t40*t120*t297*t298*t44
     #-0.4675158550002605D-1*t41*t309*t26*t44-0.4675158550002605D-1
     #*t41*t297*t26*t44+0.935031710000521D-1*t41*t110*t66*t126
     #+0.4324422990412758D0*t24*t173*t322-0.2162211495206379D0*t24*t67
     #*(0.4324422990412758D0*t120*sigmaaa*t327*t298
     #-0.2162211495206379D0*t114*t23*t309*t26-0.2162211495206379D0
     #*t114*t327*t26+0.2805095130001563D0/t338*t35*t341*t298
     #-0.935031710000521D-1*t121*t40*t309*t26-0.935031710000521D-1
     #*t121*t341*t26))*t78-0.1554534543482745D-1*t359*t193)
      v2rhob2(i) = 0.D0
      v2rhoab(i) = 0.D0
      t373 = t145**2
      t387 = t149**2
      v2sigmaaa2(i) = 0.1554534543482745D-1*rhoa*(0.935031710000521D-1
     #*t40*t28*t44-0.4324422990412758D0*t135*t66*t145
     #-0.935031710000521D-1*t138*t167*t145+0.4324422990412758D0*t24
     #*t173*t373-0.2021736311745604D-1*sigmaaa/t37/t52*t67*t34)*t78
     #-0.1554534543482745D-1*rhoa*t387*t193
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      v2sigmabb2(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = 1/rho
      t5 = t4**(1.D0/3.D0)
      t7 = 1.D0+0.1325688999052018D0*t5
      t8 = t4**(1.D0/6.D0)
      t11 = dsqrt(t4)
      t13 = t5**2
      t15 = 0.598255043577108D1*t8+0.2225569421150687D1*t5
     #+0.8004286349993634D0*t11+0.1897004325747559D0*t13
      t18 = 1.D0+0.1608197949869254D2/t15
      t19 = dlog(t18)
      t21 = 0.621814D-1*t7*t19
      t23 = 1.D0+0.6901399211255825D-1*t5
      t28 = 0.8157414703487641D1*t8+0.2247591863577616D1*t5
     #+0.4300972471276643D0*t11+0.1911512595127338D0*t13
      t31 = 1.D0+0.2960874997779344D2/t28
      t32 = dlog(t31)
      t33 = t23*t32
      t35 = rhoa-1.D0*rhob
      t36 = t35*t4
      t37 = 1.D0+t36
      t38 = t37**(1.D0/3.D0)
      t39 = t38*t37
      t41 = 1.D0-1.D0*t36
      t42 = t41**(1.D0/3.D0)
      t43 = t42*t41
      t44 = t39+t43-2.D0
      t45 = t35**2
      t46 = t45**2
      t47 = rho**2
      t48 = t47**2
      t49 = 1/t48
      t50 = t46*t49
      t52 = 1.D0-1.D0*t50
      t53 = t44*t52
      t55 = 0.37995525D-1*t33*t53
      t57 = 1.D0+0.1274696188700087D0*t5
      t62 = 0.1112037486309468D2*t8+0.3844746237447211D1*t5
     #+0.1644733775567609D1*t11+0.2405871291288192D0*t13
      t65 = 1.D0+0.3216395899738507D2/t62
      t66 = dlog(t65)
      t69 = -0.310907D-1*t57*t66+t21
      t70 = t69*t44
      t72 = 0.1923661050931536D1*t70*t50
      t73 = t38**2
      t75 = t42**2
      t77 = 0.5D0*t73+0.5D0*t75
      t78 = t77**2
      t79 = t78*t77
      t80 = 1/t78
      t81 = sigma*t80
      t82 = rho**(1.D0/3.D0)
      t84 = 1/t82/t47
      t85 = -t21+t55+t72
      t86 = 1/t79
      t89 = dexp(-0.3216396844291482D2*t85*t86)
      t90 = t89-1.D0
      t91 = 1/t90
      t92 = t91*sigma
      t93 = t80*t84
      t95 = 0.1362107888567592D0*t92*t93
      t96 = 1.D0+t95
      t98 = t90**2
      t99 = 1/t98
      t100 = sigma**2
      t101 = t99*t100
      t102 = t78**2
      t103 = 1/t102
      t104 = t82**2
      t106 = 1/t104/t48
      t107 = t103*t106
      t110 = 1.D0+t95+0.1855337900098064D-1*t101*t107
      t111 = 1/t110
      t115 = 1.D0+0.1362107888567592D0*t81*t84*t96*t111
      t116 = dlog(t115)
      t118 = 0.310906908696549D-1*t79*t116
      zk(i) = rho*(-t21+t55+t72+t118)
      t124 = 0.1333333333333333D1*t38*t4-0.1333333333333333D1*t42*t4
      t126 = t33*t124*t52
      t127 = 0.37995525D-1*t126
      t128 = t45*t35
      t129 = t44*t128
      t131 = t33*t129*t49
      t132 = 0.1519821D0*t131
      t133 = t69*t124
      t134 = t133*t50
      t135 = 0.1923661050931536D1*t134
      t136 = t128*t49
      t137 = t70*t136
      t138 = 0.7694644203726145D1*t137
      t139 = t78*t116
      t140 = 1/t38
      t143 = 1/t42
      t146 = 0.3333333333333333D0*t140*t4-0.3333333333333333D0*t143*t4
      t147 = t139*t146
      t148 = 0.9327207260896469D-1*t147
      t149 = sigma*t86
      t150 = t149*t84
      t151 = t96*t111
      t152 = t151*t146
      t155 = t99*sigma
      t156 = t155*t80
      t157 = t127-t132+t135+t138
      t160 = t85*t103
      t163 = -0.3216396844291482D2*t157*t86+0.9649190532874446D2*t160
     #*t146
      t164 = t84*t163
      t165 = t164*t89
      t167 = 0.1362107888567592D0*t156*t165
      t168 = t86*t84
      t171 = 0.2724215777135184D0*t92*t168*t146
      t172 = -t167-t171
      t177 = t81*t84
      t178 = t110**2
      t179 = 1/t178
      t180 = t96*t179
      t182 = 1/t98/t90
      t183 = t182*t100
      t184 = t183*t103
      t185 = t106*t163
      t186 = t185*t89
      t190 = 1/t102/t77
      t191 = t190*t106
      t192 = t191*t146
      t195 = -t167-t171-0.3710675800196129D-1*t184*t186
     #-0.7421351600392257D-1*t101*t192
      t196 = t180*t195
      t199 = -0.2724215777135184D0*t150*t152+0.1362107888567592D0*t81
     #*t84*t172*t111-0.1362107888567592D0*t177*t196
      t200 = t79*t199
      t201 = 1/t115
      t202 = t200*t201
      t203 = 0.310906908696549D-1*t202
      t206 = 1/t13
      t207 = 1/t47
      t208 = t206*t207
      t209 = t208*t19
      t210 = 0.2747773264188438D-2*t209
      t211 = t15**2
      t212 = 1/t211
      t213 = t7*t212
      t214 = t8**2
      t215 = t214**2
      t216 = t215*t8
      t217 = 1/t216
      t218 = t217*t207
      t221 = 1/t11
      t222 = t221*t207
      t224 = 1/t5
      t225 = t224*t207
      t227 = -0.99709173929518D0*t218-0.7418564737168958D0*t208
     #-0.4002143174996817D0*t222-0.1264669550498372D0*t225
      t228 = 1/t18
      t230 = t213*t227*t228
      t231 = 0.1D1*t230
      t232 = t32*t44
      t233 = t232*t52
      t234 = t208*t233
      t235 = 0.87407428755417D-3*t234
      t236 = t28**2
      t237 = 1/t236
      t238 = t23*t237
      t243 = -0.135956911724794D1*t218-0.7491972878592054D0*t208
     #-0.2150486235638321D0*t222-0.1274341730084892D0*t225
      t244 = t238*t243
      t245 = 1/t31
      t246 = t245*t44
      t247 = t246*t52
      t248 = t244*t247
      t249 = 0.1125D1*t248
      t250 = t38*t35
      t253 = t42*t35
      t256 = -0.1333333333333333D1*t250*t207+0.1333333333333333D1*t253
     #*t207
      t258 = t33*t256*t52
      t259 = 0.37995525D-1*t258
      t260 = t44*t46
      t261 = t48*rho
      t262 = 1/t261
      t264 = t33*t260*t262
      t265 = 0.1519821D0*t264
      t268 = t62**2
      t269 = 1/t268
      t270 = t57*t269
      t275 = -0.1853395810515781D1*t218-0.128158207914907D1*t208
     #-0.8223668877838045D0*t222-0.1603914194192128D0*t225
      t276 = 1/t65
      t282 = 0.1321039893133927D-2*t208*t66+0.1D1*t270*t275*t276
     #-0.2747773264188438D-2*t209-0.1D1*t230
      t283 = t282*t44
      t284 = t283*t50
      t285 = 0.1923661050931536D1*t284
      t286 = t69*t256
      t287 = t286*t50
      t288 = 0.1923661050931536D1*t287
      t289 = t46*t262
      t290 = t70*t289
      t291 = 0.7694644203726145D1*t290
      t292 = t140*t35
      t295 = t143*t35
      t298 = -0.3333333333333333D0*t292*t207+0.3333333333333333D0*t295
     #*t207
      t299 = t139*t298
      t301 = t151*t298
      t304 = t47*rho
      t306 = 1/t82/t304
      t311 = t210+t231-t235-t249+t259+t265+t285+t288-t291
      t316 = -0.3216396844291482D2*t311*t86+0.9649190532874446D2*t160
     #*t298
      t317 = t84*t316
      t318 = t317*t89
      t320 = 0.1362107888567592D0*t156*t318
      t323 = 0.2724215777135184D0*t92*t168*t298
      t324 = t80*t306
      t326 = 0.3178251739991049D0*t92*t324
      t327 = -t320-t323-t326
      t332 = t106*t316
      t333 = t332*t89
      t336 = t191*t298
      t340 = 1/t104/t261
      t341 = t103*t340
      t344 = -t320-t323-t326-0.3710675800196129D-1*t184*t333
     #-0.7421351600392257D-1*t101*t336-0.8658243533790967D-1*t101*t341
      t345 = t180*t344
      t348 = -0.2724215777135184D0*t150*t301-0.3178251739991049D0*t81
     #*t306*t96*t111+0.1362107888567592D0*t81*t84*t327*t111
     #-0.1362107888567592D0*t177*t345
      t349 = t79*t348
      t350 = t349*t201
      t352 = t210+t231-t235-t249+t259+t265+t285+t288-t291
     #+0.9327207260896469D-1*t299+0.310906908696549D-1*t350
      t353 = rho*t352
      vrhoa(i) = rho*(t127-t132+t135+t138+t148+t203)-t21+t55+t72+t118
     #+t353
      t354 = -t124
      t356 = t33*t354*t52
      t357 = 0.37995525D-1*t356
      t358 = 0.1519821D0*t131
      t359 = t69*t354
      t360 = t359*t50
      t361 = 0.1923661050931536D1*t360
      t362 = 0.7694644203726145D1*t137
      t363 = -t146
      t364 = t139*t363
      t365 = 0.9327207260896469D-1*t364
      t366 = t151*t363
      t369 = t357+t358+t361-t362
      t374 = -0.3216396844291482D2*t369*t86+0.9649190532874446D2*t160
     #*t363
      t375 = t84*t374
      t376 = t375*t89
      t378 = 0.1362107888567592D0*t156*t376
      t381 = 0.2724215777135184D0*t92*t168*t363
      t382 = -t378-t381
      t387 = t106*t374
      t388 = t387*t89
      t391 = t191*t363
      t394 = -t378-t381-0.3710675800196129D-1*t184*t388
     #-0.7421351600392257D-1*t101*t391
      t395 = t180*t394
      t398 = -0.2724215777135184D0*t150*t366+0.1362107888567592D0*t81
     #*t84*t382*t111-0.1362107888567592D0*t177*t395
      t399 = t79*t398
      t400 = t399*t201
      t401 = 0.310906908696549D-1*t400
      vrhob(i) = rho*(t357+t358+t361-t362+t365+t401)-t21+t55+t72+t118
     #+t353
      t404 = rho*t79
      t407 = sigma*t103
      t412 = t91*t80
      t417 = 0.1362107888567592D0*t412*t84+0.3710675800196129D-1*t155
     #*t107
      t418 = t180*t417
      t421 = 0.1362107888567592D0*t93*t151+0.1855337900098064D-1*t407
     #*t106*t91*t111-0.1362107888567592D0*t177*t418
      t423 = t404*t421*t201
      vsigmaaa(i) = 0.310906908696549D-1*t423
      vsigmaab(i) = 0.6218138173930979D-1*t423
      vsigmabb(i) = vsigmaaa(i)
      t427 = 0.87407428755417D-3*t208*t32*t124*t52
      t428 = t48*t47
      t429 = 1/t428
      t432 = t206*t429*t232*t128
      t433 = 0.349629715021668D-2*t432
      t437 = 0.1125D1*t244*t245*t124*t52
      t439 = t244*t246*t136
      t440 = 0.45D1*t439
      t441 = 1/t73
      t443 = 1/t304
      t448 = 1/t75
      t454 = -0.4444444444444444D0*t441*t35*t443-0.1333333333333333D1
     #*t38*t207-0.4444444444444444D0*t448*t35*t443
     #+0.1333333333333333D1*t42*t207
      t457 = 0.37995525D-1*t33*t454*t52
      t460 = t33*t256*t128*t49
      t461 = 0.1519821D0*t460
      t465 = 0.1519821D0*t33*t124*t46*t262
      t467 = t33*t129*t262
      t468 = 0.6079284D0*t467
      t471 = 0.1923661050931536D1*t282*t124*t50
      t472 = t283*t136
      t473 = 0.7694644203726145D1*t472
      t476 = 0.1923661050931536D1*t69*t454*t50
      t477 = t286*t136
      t478 = 0.7694644203726145D1*t477
      t480 = 0.7694644203726145D1*t133*t289
      t482 = t70*t128*t262
      t483 = 0.3077857681490458D2*t482
      t484 = t77*t116
      t485 = t298*t146
      t487 = 0.1865441452179294D0*t484*t485
      t488 = t78*t199
      t489 = t201*t298
      t491 = 0.9327207260896469D-1*t488*t489
      t492 = 1/t39
      t498 = 1/t43
      t504 = 0.1111111111111111D0*t492*t35*t443-0.3333333333333333D0
     #*t140*t207+0.1111111111111111D0*t498*t35*t443
     #+0.3333333333333333D0*t143*t207
      t506 = 0.9327207260896469D-1*t139*t504
      t507 = t78*t348
      t508 = t201*t146
      t510 = 0.9327207260896469D-1*t507*t508
      t511 = t407*t84
      t514 = 0.8172647331405553D0*t511*t151*t485
      t515 = t172*t111
      t518 = 0.2724215777135184D0*t150*t515*t298
      t522 = 0.2724215777135184D0*t150*t180*t298*t195
      t525 = 0.2724215777135184D0*t150*t151*t504
      t526 = t149*t306
      t528 = 0.6356503479982097D0*t526*t152
      t532 = 0.3178251739991049D0*t81*t306*t172*t111
      t533 = t81*t306
      t535 = 0.3178251739991049D0*t533*t196
      t536 = t327*t111
      t539 = 0.2724215777135184D0*t150*t536*t146
      t540 = t182*sigma
      t541 = t540*t80
      t542 = t89**2
      t543 = t542*t163
      t546 = 0.2724215777135184D0*t541*t317*t543
      t547 = t155*t86
      t548 = t89*t146
      t551 = 0.2724215777135184D0*t547*t317*t548
      t552 = -t427+t433-t437+t440+t457-t461+t465+t468+t471+t473+t476
     #+t478-t480-t483
      t555 = t311*t103
      t557 = 0.9649190532874446D2*t555*t146
      t558 = t157*t103
      t560 = 0.9649190532874446D2*t558*t298
      t561 = t85*t190
      t563 = 0.3859676213149779D3*t561*t485
      t565 = 0.9649190532874446D2*t160*t504
      t566 = -0.3216396844291482D2*t552*t86+t557+t560-t563+t565
      t570 = 0.1362107888567592D0*t156*t84*t566*t89
      t571 = t163*t89
      t574 = 0.1362107888567592D0*t156*t317*t571
      t575 = t84*t298
      t578 = 0.2724215777135184D0*t547*t575*t571
      t579 = t92*t103
      t582 = 0.8172647331405553D0*t579*t575*t146
      t585 = 0.2724215777135184D0*t92*t168*t504
      t589 = 0.3178251739991049D0*t156*t306*t163*t89
      t590 = t86*t306
      t593 = 0.6356503479982097D0*t92*t590*t146
      t599 = t327*t179
      t602 = 0.1362107888567592D0*t177*t599*t195
      t606 = 0.2724215777135184D0*t150*t180*t344*t146
      t607 = t172*t179
      t610 = 0.1362107888567592D0*t177*t607*t344
      t613 = t96/t178/t110
      t617 = 0.2724215777135184D0*t177*t613*t344*t195
      t618 = t98**2
      t621 = 1/t618*t100*t103
      t624 = 0.1113202740058839D0*t621*t332*t543
      t625 = t183*t190
      t628 = 0.1484270320078451D0*t625*t332*t548
      t635 = 0.3710675800196129D-1*t184*t332*t571
      t636 = t106*t298
      t639 = 0.1484270320078451D0*t625*t636*t571
      t641 = 1/t102/t78
      t642 = t101*t641
      t645 = 0.3710675800196129D0*t642*t636*t146
      t648 = 0.7421351600392257D-1*t101*t191*t504
      t652 = 0.1731648706758193D0*t184*t340*t163*t89
      t653 = t190*t340
      t656 = 0.3463297413516387D0*t101*t653*t146
      t657 = t546+t551-t570-t574+t578+t582-t585+t589+t593+t624+t628
     #-0.3710675800196129D-1*t184*t106*t566*t89-t635+t639+t645-t648
     #+t652+t656
      t661 = t514-t518+t522-t525+t528-t532+t535-t539
     #+0.1362107888567592D0*t81*t84*(t546+t551-t570-t574+t578+t582
     #-t585+t589+t593)*t111-t602+t606-t610+t617-0.1362107888567592D0
     #*t177*t180*t657
      t665 = t115**2
      t666 = 1/t665
      t669 = 0.310906908696549D-1*t349*t666*t199
      t670 = -t427+t433-t437+t440+t457-t461+t465+t468+t471+t473+t476
     #+t478-t480-t483+t487+t491+t506+t510+0.310906908696549D-1*t79
     #*t661*t201-t669
      t671 = rho*t670
      t677 = 0.4444444444444444D0*t441*t207+0.4444444444444444D0*t448
     #*t207
      t680 = 0.37995525D-1*t33*t677*t52
      t683 = t33*t124*t128*t49
      t684 = 0.3039642D0*t683
      t687 = t33*t44*t45*t49
      t688 = 0.4559463D0*t687
      t691 = 0.1923661050931536D1*t69*t677*t50
      t692 = t133*t136
      t693 = 0.1538928840745229D2*t692
      t695 = t70*t45*t49
      t696 = 0.2308393261117844D2*t695
      t697 = t146**2
      t706 = -0.1111111111111111D0*t492*t207-0.1111111111111111D0*t498
     #*t207
      t708 = 0.9327207260896469D-1*t139*t706
      t721 = 0.2724215777135184D0*t150*t151*t706
      t722 = t163**2
      t723 = t84*t722
      t726 = 0.2724215777135184D0*t541*t723*t542
      t729 = 0.5448431554270369D0*t547*t164*t548
      t738 = 0.9649190532874446D2*t160*t706
      t739 = -0.3216396844291482D2*(t680-t684-t688+t691+t693+t696)*t86
     #+0.1929838106574889D3*t558*t146-0.3859676213149779D3*t561*t697
     #+t738
      t743 = 0.1362107888567592D0*t156*t84*t739*t89
      t746 = 0.1362107888567592D0*t156*t723*t89
      t747 = t103*t84
      t750 = 0.8172647331405553D0*t92*t747*t697
      t753 = 0.2724215777135184D0*t92*t168*t706
      t762 = t195**2
      t766 = t106*t722
      t780 = t641*t106
      t786 = 0.7421351600392257D-1*t101*t191*t706
      t787 = t726+t729-t743-t746+t750-t753+0.1113202740058839D0*t621
     #*t766*t542+0.2968540640156903D0*t625*t185*t548
     #-0.3710675800196129D-1*t184*t106*t739*t89-0.3710675800196129D-1
     #*t184*t766*t89+0.3710675800196129D0*t101*t780*t697-t786
      t795 = t199**2
      t799 = t680-t684-t688+t691+t693+t696+0.1865441452179294D0*t484
     #*t697+0.1865441452179294D0*t488*t508+t708+0.310906908696549D-1
     #*t79*(0.8172647331405553D0*t511*t151*t697-0.5448431554270369D0
     #*t150*t515*t146+0.5448431554270369D0*t150*t180*t146*t195-t721
     #+0.1362107888567592D0*t81*t84*(t726+t729-t743-t746+t750-t753)
     #*t111-0.2724215777135184D0*t177*t607*t195+0.2724215777135184D0
     #*t177*t613*t762-0.1362107888567592D0*t177*t180*t787)*t201
     #-0.310906908696549D-1*t79*t795*t666
      t801 = 0.45D1*t439
      t802 = -t427-t437+t457+t465+t433+t801-t461+t468+t471+t476-t480
     #+t473+t478-t483
      t805 = -0.3216396844291482D2*t802*t86+t560+t557-t563+t565
      t809 = 0.1362107888567592D0*t156*t84*t805*t89
      t819 = t546+t578+t589-t809-t574+t551+t582+t593-t585+t624+t639
     #+t652-0.3710675800196129D-1*t184*t106*t805*t89-t635+t628+t645
     #+t656-t648
      t823 = t514+t528-t539+t606-t525-t518-t532+0.1362107888567592D0
     #*t81*t84*(t546+t578+t589-t809-t574+t551+t582+t593-t585)*t111
     #-t610+t522+t535-t602+t617-0.1362107888567592D0*t177*t180*t819
      t827 = -t427-t437+t457+t465+t433+t801-t461+t468+t471+t476-t480
     #+t473+t478-t483+t487+t510+t506+t491+0.310906908696549D-1*t79
     #*t823*t201-t669
      t828 = rho*t827
      t829 = 0.5495546528376876D-2*t209
      t831 = 0.2D1*t230
      t836 = 0.1538928840745229D2*t290
      t837 = 0.1865441452179294D0*t299
      t838 = 0.3847322101863073D1*t287
      t839 = 0.7599105D-1*t258
      t840 = 0.174814857510834D-2*t234
      t841 = 0.225D1*t248
      t842 = 0.3039642D0*t264
      t843 = 0.3847322101863073D1*t284
      t844 = 0.6218138173930979D-1*t350
      t848 = 1/t13/t4*t49
      t849 = t848*t19
      t850 = 0.1831848842792292D-2*t849
      t851 = t348**2
      t858 = t227**2
      t860 = t7/t211/t15*t858*t228
      t861 = 0.2D1*t860
      t864 = 1/t216/t4*t49
      t866 = t217*t443
      t869 = t206*t443
      t873 = 1/t11/t4*t49
      t875 = t221*t443
      t879 = 1/t5/t4*t49
      t881 = t224*t443
      t885 = t213*(-0.8309097827459833D0*t864+0.199418347859036D1*t866
     #-0.4945709824779306D0*t848+0.1483712947433792D1*t869
     #-0.2001071587498409D0*t873+0.8004286349993634D0*t875
     #-0.4215565168327908D-1*t879+0.2529339100996745D0*t881)*t228
      t886 = 0.1D1*t885
      t887 = t211**2
      t890 = t18**2
      t893 = t7/t887*t858/t890
      t894 = 0.1608197949869254D2*t893
      t898 = 0.225D1*t244*t245*t256*t52
      t900 = 0.174814857510834D-2*t869*t233
      t901 = t236**2
      t904 = t243**2
      t906 = t31**2
      t911 = 0.3330984372501762D2*t23/t901*t904/t906*t44*t52
      t914 = 0.9D1*t244*t246*t289
      t917 = t208*t212*t227*t228
      t918 = 0.8837926660346786D-1*t917
      t919 = t869*t19
      t920 = 0.5495546528376876D-2*t919
      t923 = 0.7599105D0*t33*t260*t429
      t927 = 0.3039642D0*t33*t256*t46*t262
      t939 = t275**2
      t955 = t268**2
      t958 = t65**2
      t969 = 0.8806932620892844D-3*t848*t66-0.2642079786267853D-2*t869
     #*t66-0.8497974591333914D-1*t208*t269*t275*t276-0.2D1*t57/t268
     #/t62*t939*t276+0.1D1*t270*(-0.1544496508763151D1*t864
     #+0.3706791621031562D1*t866-0.854388052766047D0*t848
     #+0.2563164158298141D1*t869-0.4111834438919023D0*t873
     #+0.1644733775567609D1*t875-0.5346380647307093D-1*t879
     #+0.3207828388384256D0*t881)*t276+0.3216395899738507D2*t57/t955
     #*t939/t958-0.1831848842792292D-2*t849+0.5495546528376876D-2*t919
     #+0.8837926660346786D-1*t917+0.2D1*t860-0.1D1*t885
     #-0.1608197949869254D2*t893
      t972 = 0.1923661050931536D1*t969*t44*t50
      t983 = 0.4444444444444444D0*t441*t45*t49+0.2666666666666667D1
     #*t250*t443+0.4444444444444444D0*t448*t45*t49
     #-0.2666666666666667D1*t253*t443
      t986 = 0.37995525D-1*t33*t983*t52
      t987 = t850-0.310906908696549D-1*t79*t851*t666-t861+t886+t894
     #-t898+t900-t911-t914-t918-t920-t923+t927+t972+t986
      t989 = 0.1538928840745229D2*t286*t289
      t992 = 0.1923661050931536D1*t69*t983*t50
      t994 = 0.58271619170278D-3*t848*t233
      t1006 = 0.1125D1*t238*(-0.1132974264373283D1*t864
     #+0.271913823449588D1*t866-0.4994648585728036D0*t848
     #+0.1498394575718411D1*t869-0.1075243117819161D0*t873
     #+0.4300972471276643D0*t875-0.4247805766949639D-1*t879
     #+0.2548683460169784D0*t881)*t247
      t1008 = 1/t48/t304
      t1012 = 0.699259430043336D-2*t206*t1008*t232*t46
      t1017 = 0.5176049408441869D-1*t208*t237*t243*t245*t53
      t1023 = 0.225D1*t23/t236/t28*t904*t247
      t1024 = t298**2
      t1029 = 0.3847322101863073D1*t282*t256*t50
      t1031 = 0.1538928840745229D2*t283*t289
      t1034 = 0.3847322101863073D2*t70*t46*t429
      t1038 = 0.174814857510834D-2*t208*t32*t256*t52
      t1054 = 1/t82/t48
      t1073 = -0.1111111111111111D0*t492*t45*t49+0.6666666666666667D0
     #*t292*t443-0.1111111111111111D0*t498*t45*t49
     #-0.6666666666666667D0*t295*t443
      t1077 = t316**2
      t1078 = t84*t1077
      t1081 = 0.2724215777135184D0*t541*t1078*t542
      t1082 = t89*t298
      t1085 = 0.5448431554270369D0*t547*t317*t1082
      t1089 = 0.6356503479982097D0*t156*t306*t316*t89
      t1090 = t850-t861+t886+t894-t898+t900-t911-t914-t918-t920-t923
     #+t927
      t1091 = t972+t986-t989+t992-t994-t1006-t1012+t1017+t1023+t1029
     #-t1031+t1034-t1038
      t1101 = -0.3216396844291482D2*(t1090+t1091)*t86
     #+0.1929838106574889D3*t555*t298-0.3859676213149779D3*t561*t1024
     #+0.9649190532874446D2*t160*t1073
      t1105 = 0.1362107888567592D0*t156*t84*t1101*t89
      t1108 = 0.1362107888567592D0*t156*t1078*t89
      t1111 = 0.8172647331405553D0*t92*t747*t1024
      t1114 = 0.1271300695996419D1*t92*t590*t298
      t1117 = 0.2724215777135184D0*t92*t168*t1073
      t1120 = 0.1059417246663683D1*t92*t80*t1054
      t1129 = t344**2
      t1150 = t106*t1077
      t1165 = t1111+t1114-t1117+t1120+t1089+t1081+t1085-t1108-t1105
     #+0.3710675800196129D0*t101*t780*t1024+0.6926594827032773D0*t101
     #*t653*t298+0.3463297413516387D0*t184*t340*t316*t89
     #-0.7421351600392257D-1*t101*t191*t1073-0.3710675800196129D-1
     #*t184*t106*t1101*t89+0.1113202740058839D0*t621*t1150*t542
     #+0.4906338002481548D0*t101*t103/t104/t428+0.2968540640156903D0
     #*t625*t332*t1082-0.3710675800196129D-1*t184*t1150*t89
      t1169 = 0.8172647331405553D0*t511*t151*t1024
     #+0.1271300695996419D1*t526*t301-0.5448431554270369D0*t150*t536
     #*t298+0.5448431554270369D0*t150*t180*t298*t344
     #+0.6356503479982097D0*t533*t345+0.1059417246663683D1*t81*t1054
     #*t96*t111-0.6356503479982097D0*t81*t306*t327*t111
     #-0.2724215777135184D0*t150*t151*t1073+0.1362107888567592D0*t81
     #*t84*(t1081+t1085+t1089-t1105-t1108+t1111+t1114-t1117+t1120)
     #*t111-0.2724215777135184D0*t177*t599*t344+0.2724215777135184D0
     #*t177*t613*t1129-0.1362107888567592D0*t177*t180*t1165
      t1177 = -t989+t992-t994-t1006-t1012+t1017+t1023
     #+0.1865441452179294D0*t484*t1024+t1029-t1031+t1034-t1038
     #+0.310906908696549D-1*t79*t1169*t201+0.1865441452179294D0*t507
     #*t489+0.9327207260896469D-1*t139*t1073
      t1179 = rho*(t987+t1177)
      t1180 = -t836+t837+t838+t839-t840-t841+t842+t843+t844
     #+0.3847322101863073D1*t134+t1179
      v2rhoa2(i) = t671-0.3039642D0*t131+rho*t799+t828+t829
     #+0.6218138173930979D-1*t202+t831+0.1538928840745229D2*t137
     #+0.1865441452179294D0*t147+0.7599105D-1*t126+t1180
      t1184 = t33*t354*t128*t49
      t1185 = 0.3039642D0*t1184
      t1186 = t359*t136
      t1187 = 0.1538928840745229D2*t1186
      t1188 = t363**2
      t1191 = t78*t398
      t1192 = t201*t363
      t1198 = t382*t111
      t1206 = t374**2
      t1207 = t84*t1206
      t1210 = 0.2724215777135184D0*t541*t1207*t542
      t1211 = t89*t363
      t1214 = 0.5448431554270369D0*t547*t375*t1211
      t1218 = t369*t103
      t1223 = -0.3216396844291482D2*(t680+t1185-t688+t691-t1187+t696)
     #*t86+0.1929838106574889D3*t1218*t363-0.3859676213149779D3*t561
     #*t1188+t738
      t1227 = 0.1362107888567592D0*t156*t84*t1223*t89
      t1230 = 0.1362107888567592D0*t156*t1207*t89
      t1233 = 0.8172647331405553D0*t92*t747*t1188
      t1239 = t382*t179
      t1243 = t394**2
      t1247 = t106*t1206
      t1264 = t1210+t1214-t1227-t1230+t1233-t753+0.1113202740058839D0
     #*t621*t1247*t542+0.2968540640156903D0*t625*t387*t1211
     #-0.3710675800196129D-1*t184*t106*t1223*t89-0.3710675800196129D-1
     #*t184*t1247*t89+0.3710675800196129D0*t101*t780*t1188-t786
      t1272 = t398**2
      t1276 = t680+t1185-t688+t691-t1187+t696+0.1865441452179294D0
     #*t484*t1188+0.1865441452179294D0*t1191*t1192+t708
     #+0.310906908696549D-1*t79*(0.8172647331405553D0*t511*t151*t1188
     #-0.5448431554270369D0*t150*t1198*t363+0.5448431554270369D0*t150
     #*t180*t363*t394-t721+0.1362107888567592D0*t81*t84*(t1210+t1214
     #-t1227-t1230+t1233-t753)*t111-0.2724215777135184D0*t177*t1239
     #*t394+0.2724215777135184D0*t177*t613*t1243-0.1362107888567592D0
     #*t177*t180*t1264)*t201-0.310906908696549D-1*t79*t1272*t666
      t1285 = 0.87407428755417D-3*t208*t32*t354*t52
      t1286 = 0.349629715021668D-2*t432
      t1290 = 0.1125D1*t244*t245*t354*t52
      t1291 = 0.45D1*t439
      t1292 = -t454
      t1295 = 0.37995525D-1*t33*t1292*t52
      t1296 = 0.1519821D0*t460
      t1300 = 0.1519821D0*t33*t354*t46*t262
      t1301 = 0.6079284D0*t467
      t1304 = 0.1923661050931536D1*t282*t354*t50
      t1305 = 0.7694644203726145D1*t472
      t1308 = 0.1923661050931536D1*t69*t1292*t50
      t1309 = 0.7694644203726145D1*t477
      t1311 = 0.7694644203726145D1*t359*t289
      t1312 = 0.3077857681490458D2*t482
      t1313 = t298*t363
      t1318 = -t504
      t1347 = t542*t374
      t1350 = 0.2724215777135184D0*t541*t317*t1347
      t1353 = 0.2724215777135184D0*t547*t317*t1211
      t1354 = -t1285-t1286-t1290-t1291+t1295+t1296+t1300-t1301+t1304
     #-t1305+t1308-t1309-t1311+t1312
      t1365 = -0.3216396844291482D2*t1354*t86+0.9649190532874446D2
     #*t555*t363+0.9649190532874446D2*t1218*t298-0.3859676213149779D3
     #*t561*t1313+0.9649190532874446D2*t160*t1318
      t1369 = 0.1362107888567592D0*t156*t84*t1365*t89
      t1370 = t374*t89
      t1373 = 0.1362107888567592D0*t156*t317*t1370
      t1376 = 0.2724215777135184D0*t547*t575*t1370
      t1379 = 0.8172647331405553D0*t579*t575*t363
      t1382 = 0.2724215777135184D0*t92*t168*t1318
      t1386 = 0.3178251739991049D0*t156*t306*t374*t89
      t1389 = 0.6356503479982097D0*t92*t590*t363
      t1438 = t1350+t1353-t1369-t1373+t1376+t1379-t1382+t1386+t1389
     #+0.1113202740058839D0*t621*t332*t1347+0.1484270320078451D0*t625
     #*t332*t1211-0.3710675800196129D-1*t184*t106*t1365*t89
     #-0.3710675800196129D-1*t184*t332*t1370+0.1484270320078451D0*t625
     #*t636*t1370+0.3710675800196129D0*t642*t636*t363
     #-0.7421351600392257D-1*t101*t191*t1318+0.1731648706758193D0*t184
     #*t340*t374*t89+0.3463297413516387D0*t101*t653*t363
      t1442 = 0.8172647331405553D0*t511*t151*t1313
     #-0.2724215777135184D0*t150*t1198*t298+0.2724215777135184D0*t150
     #*t180*t298*t394-0.2724215777135184D0*t150*t151*t1318
     #+0.6356503479982097D0*t526*t366-0.3178251739991049D0*t81*t306
     #*t382*t111+0.3178251739991049D0*t533*t395-0.2724215777135184D0
     #*t150*t536*t363+0.1362107888567592D0*t81*t84*(t1350+t1353-t1369
     #-t1373+t1376+t1379-t1382+t1386+t1389)*t111-0.1362107888567592D0
     #*t177*t599*t394+0.2724215777135184D0*t150*t180*t344*t363
     #-0.1362107888567592D0*t177*t1239*t344+0.2724215777135184D0*t177
     #*t613*t344*t394-0.1362107888567592D0*t177*t180*t1438
      t1446 = t666*t398
      t1449 = -t1285-t1286-t1290-t1291+t1295+t1296+t1300-t1301+t1304
     #-t1305+t1308-t1309-t1311+t1312+0.1865441452179294D0*t484*t1313
     #+0.9327207260896469D-1*t1191*t489+0.9327207260896469D-1*t139
     #*t1318+0.9327207260896469D-1*t507*t1192+0.310906908696549D-1*t79
     #*t1442*t201-0.310906908696549D-1*t349*t1446
      t1450 = rho*t1449
      v2rhob2(i) = 0.3039642D0*t131+t829+rho*t1276+t831
     #-0.1538928840745229D2*t137-t836+t837+t838+t839-t840-t841+t842
     #+t843+t844+0.3847322101863073D1*t360+0.7599105D-1*t356+2*t1450
     #+t1179+0.6218138173930979D-1*t400+0.1865441452179294D0*t364
      t1455 = -t677
      t1458 = 0.37995525D-1*t33*t1455*t52
      t1459 = 0.1519821D0*t683
      t1460 = 0.1519821D0*t1184
      t1461 = 0.4559463D0*t687
      t1464 = 0.1923661050931536D1*t69*t1455*t50
      t1465 = 0.7694644203726145D1*t692
      t1466 = 0.7694644203726145D1*t1186
      t1467 = 0.2308393261117844D2*t695
      t1468 = t146*t363
      t1473 = -t706
      t1496 = 0.2724215777135184D0*t541*t164*t1347
      t1499 = 0.2724215777135184D0*t547*t164*t1211
      t1511 = -0.3216396844291482D2*(t1458+t1459-t1460+t1461+t1464
     #-t1465+t1466-t1467)*t86+0.9649190532874446D2*t558*t363
     #+0.9649190532874446D2*t1218*t146-0.3859676213149779D3*t561*t1468
     #+0.9649190532874446D2*t160*t1473
      t1515 = 0.1362107888567592D0*t156*t84*t1511*t89
      t1518 = 0.1362107888567592D0*t156*t164*t1370
      t1519 = t84*t146
      t1522 = 0.2724215777135184D0*t547*t1519*t1370
      t1525 = 0.8172647331405553D0*t579*t1519*t363
      t1528 = 0.2724215777135184D0*t92*t168*t1473
      t1561 = t106*t146
      t1571 = t1496+t1499-t1515-t1518+t1522+t1525-t1528
     #+0.1113202740058839D0*t621*t185*t1347+0.1484270320078451D0*t625
     #*t185*t1211-0.3710675800196129D-1*t184*t106*t1511*t89
     #-0.3710675800196129D-1*t184*t185*t1370+0.1484270320078451D0*t625
     #*t1561*t1370+0.3710675800196129D0*t642*t1561*t363
     #-0.7421351600392257D-1*t101*t191*t1473
      t1575 = 0.8172647331405553D0*t511*t151*t1468
     #-0.2724215777135184D0*t150*t1198*t146+0.2724215777135184D0*t150
     #*t180*t146*t394-0.2724215777135184D0*t150*t151*t1473
     #-0.2724215777135184D0*t150*t515*t363+0.1362107888567592D0*t81
     #*t84*(t1496+t1499-t1515-t1518+t1522+t1525-t1528)*t111
     #-0.1362107888567592D0*t177*t607*t394+0.2724215777135184D0*t150
     #*t180*t195*t363-0.1362107888567592D0*t177*t1239*t195
     #+0.2724215777135184D0*t177*t613*t195*t394-0.1362107888567592D0
     #*t177*t180*t1571
      t1581 = t1458+t1459-t1460+t1461+t1464-t1465+t1466-t1467
     #+0.1865441452179294D0*t484*t1468+0.9327207260896469D-1*t1191
     #*t508+0.9327207260896469D-1*t139*t1473+0.9327207260896469D-1
     #*t488*t1192+0.310906908696549D-1*t79*t1575*t201
     #-0.310906908696549D-1*t200*t1446
      t1584 = t829+0.5D0*t671-t840+t148-t841+rho*t1581+t842+t1179+t831
     #+t839+t843+0.5D0*t828
      t1586 = -t836+t838+t203+t844+t357+t361+t365+t401+0.1D1*t1450
     #+t127+t135+t837
      v2rhoab(i) = t1584+t1586
      t1587 = t78*t421
      t1593 = sigma*t190*t106
      t1594 = t91*t111
      t1604 = t99*t80
      t1606 = 0.1362107888567592D0*t1604*t165
      t1607 = t91*t86
      t1609 = 0.2724215777135184D0*t1607*t1519
      t1620 = t407*t106
      t1621 = t91*t179
      t1629 = t540*t103
      t1642 = t666*t421
      t1646 = rho*(0.9327207260896469D-1*t1587*t508
     #+0.310906908696549D-1*t79*(-0.2724215777135184D0*t168*t152
     #-0.3710675800196129D-1*t1593*t1594*t146+0.2724215777135184D0
     #*t150*t180*t146*t417+0.1362107888567592D0*t93*t515
     #+0.1362107888567592D0*t81*t84*(-t1606-t1609)*t111
     #-0.1362107888567592D0*t177*t607*t417-0.1362107888567592D0*t93
     #*t196-0.1855337900098064D-1*t1620*t1621*t195
     #+0.2724215777135184D0*t177*t613*t195*t417-0.1362107888567592D0
     #*t177*t180*(-t1606-t1609-0.7421351600392257D-1*t1629*t186
     #-0.1484270320078451D0*t155*t192))*t201-0.310906908696549D-1*t200
     #*t1642)
      t1648 = t79*t421*t201
      t1649 = 0.310906908696549D-1*t1648
      t1672 = 0.1362107888567592D0*t1604*t318
      t1674 = 0.2724215777135184D0*t1607*t575
      t1676 = 0.3178251739991049D0*t412*t306
      t1704 = -0.2724215777135184D0*t168*t301-0.3710675800196129D-1
     #*t1593*t1594*t298+0.2724215777135184D0*t150*t180*t298*t417
     #-0.3178251739991049D0*t324*t151-0.4329121766895483D-1*t407*t340
     #*t91*t111+0.3178251739991049D0*t533*t418+0.1362107888567592D0
     #*t93*t536+0.1362107888567592D0*t81*t84*(-t1672-t1674-t1676)*t111
     #-0.1362107888567592D0*t177*t599*t417-0.1362107888567592D0*t93
     #*t345-0.1855337900098064D-1*t1620*t1621*t344
     #+0.2724215777135184D0*t177*t613*t344*t417-0.1362107888567592D0
     #*t177*t180*(-t1672-t1674-t1676-0.7421351600392257D-1*t1629*t333
     #-0.1484270320078451D0*t155*t336-0.1731648706758193D0*t155*t341)
      t1711 = rho*(0.9327207260896469D-1*t1587*t489
     #+0.310906908696549D-1*t79*t1704*t201-0.310906908696549D-1*t349
     #*t1642)
      v2rhoasigmaaa(i) = t1646+t1649+t1711
      t1713 = 0.6218138173930979D-1*t1648
      t1714 = 2.D0*t1711
      v2rhoasigmaab(i) = 2.D0*t1646+t1713+t1714
      v2rhoasigmabb(i) = v2rhoasigmaaa(i)
      t1729 = 0.1362107888567592D0*t1604*t376
      t1732 = 0.2724215777135184D0*t1607*t84*t363
      t1765 = rho*(0.9327207260896469D-1*t1587*t1192
     #+0.310906908696549D-1*t79*(-0.2724215777135184D0*t168*t366
     #-0.3710675800196129D-1*t1593*t1594*t363+0.2724215777135184D0
     #*t150*t180*t363*t417+0.1362107888567592D0*t93*t1198
     #+0.1362107888567592D0*t81*t84*(-t1729-t1732)*t111
     #-0.1362107888567592D0*t177*t1239*t417-0.1362107888567592D0*t93
     #*t395-0.1855337900098064D-1*t1620*t1621*t394
     #+0.2724215777135184D0*t177*t613*t394*t417-0.1362107888567592D0
     #*t177*t180*(-t1729-t1732-0.7421351600392257D-1*t1629*t388
     #-0.1484270320078451D0*t155*t391))*t201-0.310906908696549D-1*t399
     #*t1642)
      v2rhobsigmaaa(i) = t1765+t1649+t1711
      v2rhobsigmaab(i) = 2.D0*t1765+t1713+t1714
      v2rhobsigmabb(i) = v2rhobsigmaaa(i)
      t1774 = t417**2
      t1785 = t404*(0.3710675800196129D-1*t107*t1594
     #-0.2724215777135184D0*t93*t418-0.3710675800196129D-1*t1620*t1621
     #*t417+0.2724215777135184D0*t177*t613*t1774-0.5054340779364009D-2
     #*sigma*t641*t1008*t180*t99)*t201
      t1787 = t421**2
      t1789 = t404*t1787*t666
      v2sigmaaa2(i) = 0.310906908696549D-1*t1785-0.310906908696549D-1
     #*t1789
      v2sigmaaaab(i) = 0.6218138173930979D-1*t1785
     #-0.6218138173930979D-1*t1789
      v2sigmaaabb(i) = v2sigmaaa2(i)
      v2sigmaab2(i) = 0.1243627634786196D0*t1785-0.1243627634786196D0
     #*t1789
      v2sigmaabbb(i) = v2sigmaaaab(i)
      v2sigmabb2(i) = v2sigmaaabb(i)
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vrhob(i) = 0.0d0
      v2rhoa2(i) = 0.0d0
      v2rhob2(i) = 0.0d0
      v2rhoab(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      vsigmaab(i) = 0.0d0
      vsigmabb(i) = 0.0d0
      v2rhoasigmaaa(i) = 0.0d0
      v2rhoasigmaab(i) = 0.0d0
      v2rhoasigmabb(i) = 0.0d0
      v2rhobsigmaaa(i) = 0.0d0
      v2rhobsigmaab(i) = 0.0d0
      v2rhobsigmabb(i) = 0.0d0
      v2sigmaaa2(i) = 0.0d0
      v2sigmaab2(i) = 0.0d0
      v2sigmabb2(i) = 0.0d0
      v2sigmaaaab(i) = 0.0d0
      v2sigmaaabb(i) = 0.0d0
      v2sigmaabbb(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end
      
      
      subroutine rks_c_pbe
     & (ideriv,npt,rhoa1,sigmaaa1,
     &  zk,vrhoa,vsigmaaa,
     &  v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
c
c     J.P. Perdew, K. Burke, and M. Ernzerhof
c     Generalized gradient approximation made simple
c     Phys. Rev. Lett. 77 (1996) 3865-3868
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit real*8 (a-h,o-z)
      integer ideriv,npt
      real*8 rhoa1(npt)
      real*8 sigmaaa1(npt)
      real*8 zk(npt),vrhoa(npt),vsigmaaa(npt)
      real*8 v2rhoa2(npt),v2rhoasigmaaa(npt),v2sigmaaa2(npt)
      parameter(tol=1.0d-20)
      
      if(ideriv.eq.0) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = 1/rho
      t3 = t2**(1.D0/3.D0)
      t6 = t2**(1.D0/6.D0)
      t9 = dsqrt(t2)
      t11 = t3**2
      t17 = dlog(1.D0+0.1608197949869254D2/(0.598255043577108D1*t6
     #+0.2225569421150687D1*t3+0.8004286349993634D0*t9
     #+0.1897004325747559D0*t11))
      t18 = (1.D0+0.1325688999052018D0*t3)*t17
      t20 = rho**2
      t21 = rho**(1.D0/3.D0)
      t23 = 1/t21/t20
      t26 = dexp(0.2000000587336264D1*t18)
      t27 = t26-1.D0
      t31 = 0.1362107888567592D0/t27*sigma*t23
      t33 = t27**2
      t35 = sigma**2
      t37 = t20**2
      t38 = t21**2
      t49 = dlog(1.D0+0.1362107888567592D0*sigma*t23*(1.D0+t31)/(1.D0
     #+t31+0.1855337900098064D-1/t33*t35/t38/t37))
      zk(i) = rho*(-0.621814D-1*t18+0.310906908696549D-1*t49)
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = 1/rho
      t3 = t2**(1.D0/3.D0)
      t5 = 1.D0+0.1325688999052018D0*t3
      t6 = t2**(1.D0/6.D0)
      t9 = dsqrt(t2)
      t11 = t3**2
      t13 = 0.598255043577108D1*t6+0.2225569421150687D1*t3
     #+0.8004286349993634D0*t9+0.1897004325747559D0*t11
      t16 = 1.D0+0.1608197949869254D2/t13
      t17 = dlog(t16)
      t18 = t5*t17
      t19 = 0.621814D-1*t18
      t20 = rho**2
      t21 = rho**(1.D0/3.D0)
      t23 = 1/t21/t20
      t24 = sigma*t23
      t26 = dexp(0.2000000587336264D1*t18)
      t27 = t26-1.D0
      t28 = 1/t27
      t29 = t28*sigma
      t31 = 0.1362107888567592D0*t29*t23
      t32 = 1.D0+t31
      t33 = t27**2
      t34 = 1/t33
      t35 = sigma**2
      t36 = t34*t35
      t37 = t20**2
      t38 = t21**2
      t40 = 1/t38/t37
      t43 = 1.D0+t31+0.1855337900098064D-1*t36*t40
      t44 = 1/t43
      t45 = t32*t44
      t48 = 1.D0+0.1362107888567592D0*t24*t45
      t49 = dlog(t48)
      t50 = 0.310906908696549D-1*t49
      zk(i) = rho*(-t19+t50)
      t53 = 1/t20
      t54 = 1/t11*t53
      t56 = 0.2747773264188438D-2*t54*t17
      t57 = t13**2
      t60 = t6**2
      t61 = t60**2
      t77 = 0.1D1*t5/t57*(-0.99709173929518D0/t61/t6*t53
     #-0.7418564737168958D0*t54-0.4002143174996817D0/t9*t53
     #-0.1264669550498372D0/t3*t53)/t16
      t80 = 1/t21/t20/rho
      t84 = t34*sigma
      t85 = t56+t77
      t89 = 0.4381079514373337D1*t84*t23*t85*t26
      t91 = 0.3178251739991049D0*t29*t80
      t96 = t43**2
      t98 = t32/t96
      t116 = 1/t48
      vrhoa(i) = -t19+t50+rho*(t56+t77+0.310906908696549D-1*(
     #-0.3178251739991049D0*sigma*t80*t45+0.1362107888567592D0*t24*
     #(t89-t91)*t44-0.1362107888567592D0*t24*t98*(t89-t91
     #+0.119350059339396D1/t33/t27*t35*t40*t85*t26
     #-0.8658243533790967D-1*t36/t38/t37/rho))*t116)
      vsigmaaa(i) = 0.1243627634786196D0*rho*(0.1362107888567592D0*t23
     #*t32*t44+0.1855337900098064D-1*sigma*t40*t28*t44
     #-0.1362107888567592D0*t24*t98*(0.1362107888567592D0*t28*t23
     #+0.3710675800196129D-1*t84*t40))*t116
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = 1/rho
      t3 = t2**(1.D0/3.D0)
      t5 = 1.D0+0.1325688999052018D0*t3
      t6 = t2**(1.D0/6.D0)
      t9 = dsqrt(t2)
      t11 = t3**2
      t13 = 0.598255043577108D1*t6+0.2225569421150687D1*t3
     #+0.8004286349993634D0*t9+0.1897004325747559D0*t11
      t16 = 1.D0+0.1608197949869254D2/t13
      t17 = dlog(t16)
      t18 = t5*t17
      t19 = 0.621814D-1*t18
      t20 = rho**2
      t21 = rho**(1.D0/3.D0)
      t23 = 1/t21/t20
      t24 = sigma*t23
      t26 = dexp(0.2000000587336264D1*t18)
      t27 = t26-1.D0
      t28 = 1/t27
      t29 = t28*sigma
      t31 = 0.1362107888567592D0*t29*t23
      t32 = 1.D0+t31
      t33 = t27**2
      t34 = 1/t33
      t35 = sigma**2
      t36 = t34*t35
      t37 = t20**2
      t38 = t21**2
      t40 = 1/t38/t37
      t43 = 1.D0+t31+0.1855337900098064D-1*t36*t40
      t44 = 1/t43
      t45 = t32*t44
      t48 = 1.D0+0.1362107888567592D0*t24*t45
      t49 = dlog(t48)
      t50 = 0.310906908696549D-1*t49
      zk(i) = rho*(-t19+t50)
      t52 = 1/t11
      t53 = 1/t20
      t54 = t52*t53
      t55 = t54*t17
      t56 = 0.2747773264188438D-2*t55
      t57 = t13**2
      t58 = 1/t57
      t59 = t5*t58
      t60 = t6**2
      t61 = t60**2
      t62 = t61*t6
      t63 = 1/t62
      t67 = 1/t9
      t70 = 1/t3
      t73 = -0.99709173929518D0*t63*t53-0.7418564737168958D0*t54
     #-0.4002143174996817D0*t67*t53-0.1264669550498372D0*t70*t53
      t74 = 1/t16
      t76 = t59*t73*t74
      t77 = 0.1D1*t76
      t78 = t20*rho
      t80 = 1/t21/t78
      t81 = sigma*t80
      t84 = t34*sigma
      t85 = t56+t77
      t89 = 0.4381079514373337D1*t84*t23*t85*t26
      t91 = 0.3178251739991049D0*t29*t80
      t92 = t89-t91
      t93 = t92*t44
      t96 = t43**2
      t97 = 1/t96
      t98 = t32*t97
      t100 = 1/t33/t27
      t101 = t100*t35
      t103 = t40*t85*t26
      t108 = 1/t38/t37/rho
      t111 = t89-t91+0.119350059339396D1*t101*t103
     #-0.8658243533790967D-1*t36*t108
      t112 = t98*t111
      t115 = -0.3178251739991049D0*t81*t45+0.1362107888567592D0*t24
     #*t93-0.1362107888567592D0*t24*t112
      t116 = 1/t48
      t117 = t115*t116
      vrhoa(i) = -t19+t50+rho*(t56+t77+0.310906908696549D-1*t117)
      t121 = t23*t32
      t124 = sigma*t40
      t125 = t28*t44
      t132 = 0.1362107888567592D0*t28*t23+0.3710675800196129D-1*t84*t40
      t133 = t98*t132
      t136 = 0.1362107888567592D0*t121*t44+0.1855337900098064D-1*t124
     #*t125-0.1362107888567592D0*t24*t133
      vsigmaaa(i) = 0.1243627634786196D0*rho*t136*t116
      t151 = dlog(1.D0+0.2960874997779344D2/(0.8157414703487641D1*t6
     #+0.2247591863577616D1*t3+0.4300972471276643D0*t9
     #+0.1911512595127338D0*t11))
      t153 = (1.D0+0.6901399211255825D-1*t3)*t151*t53
      t155 = t49*t53
      t158 = 1/t21/t37
      t160 = sigma*t158*t45
      t165 = -0.1086299437397317D1*t153+0.1333333724890842D1*t18*t53
      t169 = 0.1362107888567592D0*t84*t23*t165*t26
      t170 = t29*t158
      t171 = 0.6053812838078188D-1*t170
      t183 = t36/t38/t37/t20
      t197 = -t165
      t201 = 0.1362107888567592D0*t84*t23*t197*t26
      t202 = 0.6053812838078188D-1*t170
      t223 = 1/t37
      t226 = 1/t78
      t231 = 1/t11/t2*t223
      t233 = t52*t226
      t250 = 0.1D1*t59*(-0.8309097827459833D0/t62/t2*t223
     #+0.199418347859036D1*t63*t226-0.4945709824779306D0*t231
     #+0.1483712947433792D1*t233-0.2001071587498409D0/t9/t2*t223
     #+0.8004286349993634D0*t67*t226-0.4215565168327908D-1/t3/t2*t223
     #+0.2529339100996745D0*t70*t226)*t74
      t251 = t57**2
      t254 = t73**2
      t255 = t16**2
      t259 = 0.1608197949869254D2*t5/t251*t254/t255
      t261 = 0.1831848842792292D-2*t231*t17
      t265 = 0.8837926660346786D-1*t54*t58*t73*t74
      t271 = 0.2D1*t5/t57/t13*t254*t74
      t273 = 0.5495546528376876D-2*t233*t17
      t274 = t115**2
      t275 = t48**2
      t276 = 1/t275
      t282 = t92*t97
      t287 = 1/t96/t43
      t288 = t32*t287
      t289 = t111**2
      t293 = t250+t259+t261-t265-t271-t273
      t297 = 0.4381079514373337D1*t84*t23*t293*t26
      t298 = t85**2
      t299 = t23*t298
      t302 = 0.1409129032462046D3*t84*t299*t26
      t311 = t33**2
      t314 = t40*t298
      t315 = t26**2
      t325 = 0.2044503773374224D2*t84*t80*t85*t26
      t326 = t100*sigma
      t329 = 0.2818258064924092D3*t326*t299*t315
      t330 = 0.1059417246663683D1*t170
      s1 = 0.1099109305675375D-1*t55+0.4D1*t76+rho*(0.337738D-1*t153
     #-0.207271272464366D-1*t155+0.310906908696549D-1*
     #(0.6053812838078188D-1*t160+0.1362107888567592D0*t24*(-t169+t171
     #)*t44-0.1362107888567592D0*t24*t98*(-t169+t171
     #-0.3710675800196129D-1*t101*t40*t165*t26+0.1649189244531613D-1
     #*t183))*t116)
      s2 = s1+rho*(-0.337738D-1*t153+0.207271272464366D-1*t155
     #+0.310906908696549D-1*(-0.6053812838078188D-1*t160
     #+0.1362107888567592D0*t24*(-t201-t202)*t44-0.1362107888567592D0
     #*t24*t98*(-t201-t202-0.3710675800196129D-1*t101*t40*t197*t26
     #-0.1649189244531613D-1*t183))*t116)
      v2rhoa2(i) = s2+2.D0*rho*(t250+t259+t261-t265-t271-t273
     #-0.310906908696549D-1*t274*t276+0.310906908696549D-1*
     #(0.1059417246663683D1*t160-0.6356503479982097D0*t81*t93
     #-0.2724215777135184D0*t24*t282*t111+0.2724215777135184D0*t24
     #*t288*t289-0.1362107888567592D0*t24*t98*(t297-t302
     #-0.1113933887167696D2*t101*t108*t85*t26+0.119350059339396D1*t101
     #*t40*t293*t26+0.1151631462675703D3/t311*t35*t314*t315
     #-0.3838771542252344D2*t101*t314*t26-t325+t329+t330
     #+0.4906338002481548D0*t183)+0.6356503479982097D0*t81*t112
     #+0.1362107888567592D0*t24*(t329-t325+t297-t302+t330)*t44)*t116)
     #+0.1243627634786196D0*t117
      t365 = 0.4381079514373337D1*t34*t23*t85*t26
      t367 = 0.3178251739991049D0*t28*t80
      t378 = t28*t97
      v2rhoasigmaaa(i) = 0.1243627634786196D0*t136*t116+4.D0*rho*
     #(0.310906908696549D-1*(-0.3178251739991049D0*t80*t32*t44
     #-0.4329121766895483D-1*sigma*t108*t125+0.3178251739991049D0*t81
     #*t133+0.1362107888567592D0*t23*t92*t44+0.1362107888567592D0*t24*
     #(t365-t367)*t44-0.1362107888567592D0*t24*t282*t132
     #-0.1362107888567592D0*t121*t97*t111-0.1855337900098064D-1*t124
     #*t378*t111+0.2724215777135184D0*t24*t32*t287*t111*t132
     #-0.1362107888567592D0*t24*t98*(t365-t367+0.238700118678792D1
     #*t326*t103-0.1731648706758193D0*t84*t108))*t116
     #-0.310906908696549D-1*t115*t276*t136)
      t413 = t132**2
      t427 = t136**2
      v2sigmaaa2(i) = 0.4974510539144783D0*rho*(0.3710675800196129D-1
     #*t40*t28*t44-0.2724215777135184D0*t121*t97*t132
     #-0.3710675800196129D-1*t124*t378*t132+0.2724215777135184D0*t24
     #*t288*t413-0.5054340779364009D-2*sigma/t37/t78*t98*t34)*t116
     #-0.4974510539144783D0*rho*t427*t276
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      v2rhoa2(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      v2rhoasigmaaa(i) = 0.0d0
      v2sigmaaa2(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end


      subroutine uks_xc_b97
     & (ideriv,npt,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,
     &  zk,vrhoa,vrhob,vsigmaaa,vsigmabb,vsigmaab,
     &  v2rhoa2,v2rhob2,v2rhoab,
     &  v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb,
     &  v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,
     &  v2sigmaaa2,v2sigmaaaab,v2sigmaaabb,
     &  v2sigmaab2,v2sigmaabbb,v2sigmabb2)
c
c     A.D. Becke
c     Density-functional thermochemistry. V. Systematic optimization of
c     exchange-correlation functionals
c     J. Chem. Phys. 107 (1997) 8554-8560
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit real*8 (a-h,o-z)
      integer ideriv,npt
      real*8 rhoa1(npt),rhob1(npt)
      real*8 sigmaaa1(npt),sigmabb1(npt),sigmaab1(npt)
      real*8 zk(npt),vrhoa(npt),vrhob(npt)
      real*8 vsigmaaa(npt),vsigmabb(npt),vsigmaab(npt)
      real*8 v2rhoa2(npt),v2rhob2(npt),v2rhoab(npt)
      real*8 v2rhoasigmaaa(npt),v2rhoasigmaab(npt)
      real*8 v2rhoasigmabb(npt),v2rhobsigmabb(npt)
      real*8 v2rhobsigmaab(npt),v2rhobsigmaaa(npt)
      real*8 v2sigmaaa2(npt),v2sigmaaaab(npt),v2sigmaaabb(npt)
      real*8 v2sigmaab2(npt),v2sigmaabbb(npt),v2sigmabb2(npt)
      parameter(tol=1.0d-20)
      
      if (ideriv.eq.0) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmabb
      t2 = rhob**(1.D0/3.D0)
      t4 = rhob**2
      t5 = t2**2
      t8 = sigmabb/t5/t4
      t10 = 1.D0+0.4D-2*t8
      t14 = sigmabb**2
      t15 = t4**2
      t19 = t14/t2/t15/rhob
      t20 = t10**2
      t27 = 1/rhob
      t28 = t27**(1.D0/3.D0)
      t32 = t27**(1.D0/6.D0)
      t35 = dsqrt(t27)
      t37 = t28**2
      t43 = dlog(1.D0+0.321646831778707D2/(0.1112037486309468D2*t32
     #+0.3844746237447211D1*t28+0.1644733775567609D1*t35
     #+0.2405871291288192D0*t37))
      t45 = 1.D0+0.2D0*t8
      t49 = t45**2
      zk(i) = -0.9305257363491D0*t2*rhob*(0.8094D0+0.20292D-2*t8/t10
     #+0.119696D-4*t19/t20)-0.3109D-1*rhob*(1.D0+0.1274696188700087D0
     #*t28)*t43*(0.1737D0+0.46974D0*t8/t45-0.99472D-1*t19/t49)
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = rhoa**(1.D0/3.D0)
      t4 = rhoa**2
      t5 = t2**2
      t8 = sigmaaa/t5/t4
      t10 = 1.D0+0.4D-2*t8
      t14 = sigmaaa**2
      t15 = t4**2
      t19 = t14/t2/t15/rhoa
      t20 = t10**2
      t27 = 1/rhoa
      t28 = t27**(1.D0/3.D0)
      t32 = t27**(1.D0/6.D0)
      t35 = dsqrt(t27)
      t37 = t28**2
      t43 = dlog(1.D0+0.321646831778707D2/(0.1112037486309468D2*t32
     #+0.3844746237447211D1*t28+0.1644733775567609D1*t35
     #+0.2405871291288192D0*t37))
      t45 = 1.D0+0.2D0*t8
      t49 = t45**2
      zk(i) = -0.9305257363491D0*t2*rhoa*(0.8094D0+0.20292D-2*t8/t10
     #+0.119696D-4*t19/t20)-0.3109D-1*rhoa*(1.D0+0.1274696188700087D0
     #*t28)*t43*(0.1737D0+0.46974D0*t8/t45-0.99472D-1*t19/t49)
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rhoa**(1.D0/3.D0)
      t6 = rhoa**2
      t7 = t4**2
      t10 = sigmaaa/t7/t6
      t12 = 1.D0+0.4D-2*t10
      t16 = sigmaaa**2
      t17 = t6**2
      t21 = t16/t4/t17/rhoa
      t22 = t12**2
      t29 = 1/rhoa
      t30 = t29**(1.D0/3.D0)
      t33 = rhoa*(1.D0+0.1274696188700087D0*t30)
      t34 = t29**(1.D0/6.D0)
      t37 = dsqrt(t29)
      t39 = t30**2
      t45 = dlog(1.D0+0.321646831778707D2/(0.1112037486309468D2*t34
     #+0.3844746237447211D1*t30+0.1644733775567609D1*t37
     #+0.2405871291288192D0*t39))
      t47 = 1.D0+0.2D0*t10
      t51 = t47**2
      t59 = rhob**(1.D0/3.D0)
      t61 = rhob**2
      t62 = t59**2
      t65 = sigmabb/t62/t61
      t67 = 1.D0+0.4D-2*t65
      t71 = sigmabb**2
      t72 = t61**2
      t76 = t71/t59/t72/rhob
      t77 = t67**2
      t84 = 1/rhob
      t85 = t84**(1.D0/3.D0)
      t88 = rhob*(1.D0+0.1274696188700087D0*t85)
      t89 = t84**(1.D0/6.D0)
      t92 = dsqrt(t84)
      t94 = t85**2
      t100 = dlog(1.D0+0.321646831778707D2/(0.1112037486309468D2*t89
     #+0.3844746237447211D1*t85+0.1644733775567609D1*t92
     #+0.2405871291288192D0*t94))
      t102 = 1.D0+0.2D0*t65
      t106 = t102**2
      t114 = rhoa+rhob
      t115 = 1/t114
      t116 = t115**(1.D0/3.D0)
      t119 = t115**(1.D0/6.D0)
      t122 = dsqrt(t115)
      t124 = t116**2
      t130 = dlog(1.D0+0.160818243221511D2/(0.598255043577108D1*t119
     #+0.2225569421150687D1*t116+0.8004286349993634D0*t122
     #+0.1897004325747559D0*t124))
      t132 = 0.62182D-1*(1.D0+0.1325688999052018D0*t116)*t130
      t143 = dlog(1.D0+0.2960857464321668D2/(0.8157414703487641D1*t119
     #+0.2247591863577616D1*t116+0.4300972471276643D0*t122
     #+0.1911512595127338D0*t124))
      t146 = rhoa-1.D0*rhob
      t147 = t146*t115
      t148 = 1.D0+t147
      t149 = t148**(1.D0/3.D0)
      t152 = 1.D0-1.D0*t147
      t153 = t152**(1.D0/3.D0)
      t155 = t149*t148+t153*t152-2.D0
      t156 = t146**2
      t157 = t156**2
      t158 = t114**2
      t159 = t158**2
      t161 = t157/t159
      t177 = dlog(1.D0+0.321646831778707D2/(0.1112037486309468D2*t119
     #+0.3844746237447211D1*t116+0.1644733775567609D1*t122
     #+0.2405871291288192D0*t124))
      t193 = 0.5D0*t10+0.5D0*t65
      t196 = 1.D0+0.3D-2*t10+0.3D-2*t65
      t200 = t193**2
      t201 = t196**2
      zk(i) = -0.9305257363491D0*t4*rhoa*(0.8094D0+0.20292D-2*t10/t12
     #+0.119696D-4*t21/t22)-0.3109D-1*t33*t45*(0.1737D0+0.46974D0*t10
     #/t47-0.99472D-1*t21/t51)-0.9305257363491D0*t59*rhob*(0.8094D0
     #+0.20292D-2*t65/t67+0.119696D-4*t76/t77)-0.3109D-1*t88*t100*
     #(0.1737D0+0.46974D0*t65/t102-0.99472D-1*t76/t106)+(t114*(-t132
     #+0.3799574853701528D-1*(1.D0+0.6901399211255825D-1*t116)*t143
     #*t155*(1.D0-1.D0*t161)+0.1923661050931536D1*(-0.3109D-1*(1.D0
     #+0.1274696188700087D0*t116)*t177+t132)*t155*t161)+0.3109D-1*t33
     #*t45+0.3109D-1*t88*t100)*(0.9454D0+0.44826D-2*t193/t196
     #-0.1654596D-3*t200/t201)
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmabb
      t2 = rhob**(1.D0/3.D0)
      t3 = t2*rhob
      t4 = rhob**2
      t5 = t2**2
      t7 = 1/t5/t4
      t8 = sigmabb*t7
      t10 = 1.D0+0.4D-2*t8
      t11 = 1/t10
      t14 = sigmabb**2
      t15 = t4**2
      t18 = 1/t2/t15/rhob
      t19 = t14*t18
      t20 = t10**2
      t21 = 1/t20
      t24 = 0.8094D0+0.20292D-2*t8*t11+0.119696D-4*t19*t21
      t27 = 1/rhob
      t28 = t27**(1.D0/3.D0)
      t30 = 1.D0+0.1274696188700087D0*t28
      t31 = rhob*t30
      t32 = t27**(1.D0/6.D0)
      t35 = dsqrt(t27)
      t37 = t28**2
      t39 = 0.1112037486309468D2*t32+0.3844746237447211D1*t28
     #+0.1644733775567609D1*t35+0.2405871291288192D0*t37
      t42 = 1.D0+0.321646831778707D2/t39
      t43 = dlog(t42)
      t45 = 1.D0+0.2D0*t8
      t46 = 1/t45
      t49 = t45**2
      t50 = 1/t49
      t53 = 0.1737D0+0.46974D0*t8*t46-0.99472D-1*t19*t50
      t54 = t43*t53
      zk(i) = -0.9305257363491D0*t3*t24-0.3109D-1*t31*t54
      vrhoa(i) = 0.D0
      t62 = sigmabb/t5/t4/rhob
      t68 = t14/t2/t15/t4
      t72 = t15**2
      t75 = t14*sigmabb/t72/rhob
      t77 = 1/t20/t10
      t86 = 1/t37
      t90 = t39**2
      t93 = t32**2
      t94 = t93**2
      t97 = 1/t4
      t119 = 1/t49/t45
      vrhob(i) = -0.12407009817988D1*t2*t24-0.9305257363491D0*t3*(
     #-0.54112D-2*t62*t11-0.4219306666666667D-4*t68*t21
     #+0.2553514666666667D-6*t75*t77)-0.3109D-1*t30*t43*t53
     #+0.1321010150222857D-2*t27*t86*t54+0.1D1*t31/t90*(
     #-0.1853395810515781D1/t94/t32*t97-0.128158207914907D1*t86*t97
     #-0.8223668877838045D0/t35*t97-0.1603914194192128D0/t28*t97)/t42
     #*t53-0.3109D-1*t31*t43*(-0.125264D1*t62*t46+0.7810453333333333D0
     #*t68*t50-0.1061034666666667D0*t75*t119)
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      t128 = sigmabb*t18
      t132 = t14/t72
      vsigmabb(i) = -0.9305257363491D0*t3*(0.20292D-2*t7*t11
     #+0.158224D-4*t128*t21-0.957568D-7*t132*t77)-0.3109D-1*t31*t43*
     #(0.46974D0*t7*t46-0.292892D0*t128*t50+0.397888D-1*t132*t119)
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = rhoa**(1.D0/3.D0)
      t3 = t2*rhoa
      t4 = rhoa**2
      t5 = t2**2
      t7 = 1/t5/t4
      t8 = sigmaaa*t7
      t10 = 1.D0+0.4D-2*t8
      t11 = 1/t10
      t14 = sigmaaa**2
      t15 = t4**2
      t18 = 1/t2/t15/rhoa
      t19 = t14*t18
      t20 = t10**2
      t21 = 1/t20
      t24 = 0.8094D0+0.20292D-2*t8*t11+0.119696D-4*t19*t21
      t27 = 1/rhoa
      t28 = t27**(1.D0/3.D0)
      t30 = 1.D0+0.1274696188700087D0*t28
      t31 = rhoa*t30
      t32 = t27**(1.D0/6.D0)
      t35 = dsqrt(t27)
      t37 = t28**2
      t39 = 0.1112037486309468D2*t32+0.3844746237447211D1*t28
     #+0.1644733775567609D1*t35+0.2405871291288192D0*t37
      t42 = 1.D0+0.321646831778707D2/t39
      t43 = dlog(t42)
      t45 = 1.D0+0.2D0*t8
      t46 = 1/t45
      t49 = t45**2
      t50 = 1/t49
      t53 = 0.1737D0+0.46974D0*t8*t46-0.99472D-1*t19*t50
      t54 = t43*t53
      zk(i) = -0.9305257363491D0*t3*t24-0.3109D-1*t31*t54
      t62 = sigmaaa/t5/t4/rhoa
      t68 = t14/t2/t15/t4
      t72 = t15**2
      t75 = t14*sigmaaa/t72/rhoa
      t77 = 1/t20/t10
      t86 = 1/t37
      t90 = t39**2
      t93 = t32**2
      t94 = t93**2
      t97 = 1/t4
      t119 = 1/t49/t45
      vrhoa(i) = -0.12407009817988D1*t2*t24-0.9305257363491D0*t3*(
     #-0.54112D-2*t62*t11-0.4219306666666667D-4*t68*t21
     #+0.2553514666666667D-6*t75*t77)-0.3109D-1*t30*t43*t53
     #+0.1321010150222857D-2*t27*t86*t54+0.1D1*t31/t90*(
     #-0.1853395810515781D1/t94/t32*t97-0.128158207914907D1*t86*t97
     #-0.8223668877838045D0/t35*t97-0.1603914194192128D0/t28*t97)/t42
     #*t53-0.3109D-1*t31*t43*(-0.125264D1*t62*t46+0.7810453333333333D0
     #*t68*t50-0.1061034666666667D0*t75*t119)
      vrhob(i) = 0.D0
      t128 = sigmaaa*t18
      t132 = t14/t72
      vsigmaaa(i) = -0.9305257363491D0*t3*(0.20292D-2*t7*t11
     #+0.158224D-4*t128*t21-0.957568D-7*t132*t77)-0.3109D-1*t31*t43*
     #(0.46974D0*t7*t46-0.292892D0*t128*t50+0.397888D-1*t132*t119)
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rhoa**(1.D0/3.D0)
      t5 = t4*rhoa
      t6 = rhoa**2
      t7 = t4**2
      t9 = 1/t7/t6
      t10 = sigmaaa*t9
      t12 = 1.D0+0.4D-2*t10
      t13 = 1/t12
      t16 = sigmaaa**2
      t17 = t6**2
      t20 = 1/t4/t17/rhoa
      t21 = t16*t20
      t22 = t12**2
      t23 = 1/t22
      t26 = 0.8094D0+0.20292D-2*t10*t13+0.119696D-4*t21*t23
      t29 = 1/rhoa
      t30 = t29**(1.D0/3.D0)
      t32 = 1.D0+0.1274696188700087D0*t30
      t33 = rhoa*t32
      t34 = t29**(1.D0/6.D0)
      t37 = dsqrt(t29)
      t39 = t30**2
      t41 = 0.1112037486309468D2*t34+0.3844746237447211D1*t30
     #+0.1644733775567609D1*t37+0.2405871291288192D0*t39
      t44 = 1.D0+0.321646831778707D2/t41
      t45 = dlog(t44)
      t47 = 1.D0+0.2D0*t10
      t48 = 1/t47
      t51 = t47**2
      t52 = 1/t51
      t55 = 0.1737D0+0.46974D0*t10*t48-0.99472D-1*t21*t52
      t56 = t45*t55
      t59 = rhob**(1.D0/3.D0)
      t60 = t59*rhob
      t61 = rhob**2
      t62 = t59**2
      t64 = 1/t62/t61
      t65 = sigmabb*t64
      t67 = 1.D0+0.4D-2*t65
      t68 = 1/t67
      t71 = sigmabb**2
      t72 = t61**2
      t75 = 1/t59/t72/rhob
      t76 = t71*t75
      t77 = t67**2
      t78 = 1/t77
      t81 = 0.8094D0+0.20292D-2*t65*t68+0.119696D-4*t76*t78
      t84 = 1/rhob
      t85 = t84**(1.D0/3.D0)
      t87 = 1.D0+0.1274696188700087D0*t85
      t88 = rhob*t87
      t89 = t84**(1.D0/6.D0)
      t92 = dsqrt(t84)
      t94 = t85**2
      t96 = 0.1112037486309468D2*t89+0.3844746237447211D1*t85
     #+0.1644733775567609D1*t92+0.2405871291288192D0*t94
      t99 = 1.D0+0.321646831778707D2/t96
      t100 = dlog(t99)
      t102 = 1.D0+0.2D0*t65
      t103 = 1/t102
      t106 = t102**2
      t107 = 1/t106
      t110 = 0.1737D0+0.46974D0*t65*t103-0.99472D-1*t76*t107
      t111 = t100*t110
      t114 = rhoa+rhob
      t115 = 1/t114
      t116 = t115**(1.D0/3.D0)
      t118 = 1.D0+0.1325688999052018D0*t116
      t119 = t115**(1.D0/6.D0)
      t122 = dsqrt(t115)
      t124 = t116**2
      t126 = 0.598255043577108D1*t119+0.2225569421150687D1*t116
     #+0.8004286349993634D0*t122+0.1897004325747559D0*t124
      t129 = 1.D0+0.160818243221511D2/t126
      t130 = dlog(t129)
      t132 = 0.62182D-1*t118*t130
      t134 = 1.D0+0.6901399211255825D-1*t116
      t139 = 0.8157414703487641D1*t119+0.2247591863577616D1*t116
     #+0.4300972471276643D0*t122+0.1911512595127338D0*t124
      t142 = 1.D0+0.2960857464321668D2/t139
      t143 = dlog(t142)
      t144 = t134*t143
      t146 = rhoa-1.D0*rhob
      t147 = t146*t115
      t148 = 1.D0+t147
      t149 = t148**(1.D0/3.D0)
      t152 = 1.D0-1.D0*t147
      t153 = t152**(1.D0/3.D0)
      t155 = t149*t148+t153*t152-2.D0
      t156 = t146**2
      t157 = t156**2
      t158 = t114**2
      t159 = t158**2
      t160 = 1/t159
      t161 = t157*t160
      t163 = 1.D0-1.D0*t161
      t166 = 0.3799574853701528D-1*t144*t155*t163
      t168 = 1.D0+0.1274696188700087D0*t116
      t173 = 0.1112037486309468D2*t119+0.3844746237447211D1*t116
     #+0.1644733775567609D1*t122+0.2405871291288192D0*t124
      t176 = 1.D0+0.321646831778707D2/t173
      t177 = dlog(t176)
      t180 = -0.3109D-1*t168*t177+t132
      t181 = t180*t155
      t183 = 0.1923661050931536D1*t181*t161
      t190 = t114*(-t132+t166+t183)+0.3109D-1*t33*t45+0.3109D-1*t88*t100
      t193 = 0.5D0*t10+0.5D0*t65
      t196 = 1.D0+0.3D-2*t10+0.3D-2*t65
      t197 = 1/t196
      t200 = t193**2
      t201 = t196**2
      t202 = 1/t201
      t205 = 0.9454D0+0.44826D-2*t193*t197-0.1654596D-3*t200*t202
      zk(i) = -0.9305257363491D0*t5*t26-0.3109D-1*t33*t56
     #-0.9305257363491D0*t60*t81-0.3109D-1*t88*t111+t190*t205
      t212 = sigmaaa/t7/t6/rhoa
      t218 = t16/t4/t17/t6
      t222 = t17**2
      t225 = t16*sigmaaa/t222/rhoa
      t227 = 1/t22/t12
      t233 = t32*t45
      t236 = 1/t39
      t237 = t29*t236
      t240 = t41**2
      t241 = 1/t240
      t243 = t34**2
      t244 = t243**2
      t247 = 1/t6
      t258 = -0.1853395810515781D1/t244/t34*t247-0.128158207914907D1
     #*t236*t247-0.8223668877838045D0/t37*t247-0.1603914194192128D0
     #/t30*t247
      t259 = 1/t44
      t269 = 1/t51/t47
      t277 = 1/t158
      t278 = 1/t124*t277
      t279 = t278*t130
      t280 = 0.2747799777968419D-2*t279
      t281 = t126**2
      t284 = t119**2
      t285 = t284**2
      t288 = 1/t285/t119*t277
      t292 = 1/t122*t277
      t295 = 1/t116*t277
      t300 = t118/t281*(-0.99709173929518D0*t288-0.7418564737168958D0
     #*t278-0.4002143174996817D0*t292-0.1264669550498372D0*t295)/t129
      t301 = 0.1D1*t300
      t305 = 0.8740794299481065D-3*t278*t143*t155*t163
      t306 = t139**2
      t319 = 0.1124999956683108D1*t134/t306*(-0.135956911724794D1*t288
     #-0.7491972878592054D0*t278-0.2150486235638321D0*t292
     #-0.1274341730084892D0*t295)/t142*t155*t163
      t320 = t146*t277
      t321 = 1.D0*t320
      t325 = 1.D0*t115
      t329 = 0.1333333333333333D1*t149*(t115-t321)
     #+0.1333333333333333D1*t153*(-t325+t320)
      t334 = t156*t146*t160
      t335 = 4.D0*t334
      t338 = t157/t159/t114
      t339 = 4.D0*t338
      t346 = t173**2
      t363 = 0.1923661050931536D1*(0.1321010150222857D-2*t278*t177
     #+0.1D1*t168/t346*(-0.1853395810515781D1*t288-0.128158207914907D1
     #*t278-0.8223668877838045D0*t292-0.1603914194192128D0*t295)/t176
     #-0.2747799777968419D-2*t279-0.1D1*t300)*t155*t161
      t367 = t181*t334
      t370 = 0.7694644203726145D1*t181*t338
      t384 = t193*t202
      t389 = t200/t201/t196
      vrhoa(i) = -0.12407009817988D1*t4*t26-0.9305257363491D0*t5*(
     #-0.54112D-2*t212*t13-0.4219306666666667D-4*t218*t23
     #+0.2553514666666667D-6*t225*t227)-0.3109D-1*t233*t55
     #+0.1321010150222857D-2*t237*t56+0.1D1*t33*t241*t258*t259*t55
     #-0.3109D-1*t33*t45*(-0.125264D1*t212*t48+0.7810453333333333D0
     #*t218*t52-0.1061034666666667D0*t225*t269)+(-t132+t166+t183+t114*
     #(t280+t301-t305-t319+0.3799574853701528D-1*t144*t329*t163
     #+0.3799574853701528D-1*t144*t155*(-t335+t339)+t363
     #+0.1923661050931536D1*t180*t329*t161+0.7694644203726145D1*t367
     #-t370)+0.3109D-1*t233-0.1321010150222857D-2*t237*t45-0.1D1*t33
     #*t241*t258*t259)*t205+t190*(-0.59768D-2*t212*t197+0.4770864D-3
     #*t384*t212-0.26473536D-5*t389*t212)
      t399 = sigmabb/t62/t61/rhob
      t405 = t71/t59/t72/t61
      t409 = t72**2
      t412 = t71*sigmabb/t409/rhob
      t414 = 1/t77/t67
      t420 = t87*t100
      t423 = 1/t94
      t424 = t84*t423
      t427 = t96**2
      t428 = 1/t427
      t430 = t89**2
      t431 = t430**2
      t434 = 1/t61
      t445 = -0.1853395810515781D1/t431/t89*t434-0.128158207914907D1
     #*t423*t434-0.8223668877838045D0/t92*t434-0.1603914194192128D0
     #/t85*t434
      t446 = 1/t99
      t456 = 1/t106/t102
      t469 = 0.1333333333333333D1*t149*(-t325-t321)
     #+0.1333333333333333D1*t153*(t115+t320)
      s1 = -0.12407009817988D1*t59*t81-0.9305257363491D0*t60*(
     #-0.54112D-2*t399*t68-0.4219306666666667D-4*t405*t78
     #+0.2553514666666667D-6*t412*t414)-0.3109D-1*t420*t110
     #+0.1321010150222857D-2*t424*t111
      vrhob(i) = s1+0.1D1*t88*t428*t445*t446*t110-0.3109D-1*t88*t100*(
     #-0.125264D1*t399*t103+0.7810453333333333D0*t405*t107
     #-0.1061034666666667D0*t412*t456)+(-t132+t166+t183+t114*(t280
     #+t301-t305-t319+0.3799574853701528D-1*t144*t469*t163
     #+0.3799574853701528D-1*t144*t155*(t335+t339)+t363
     #+0.1923661050931536D1*t180*t469*t161-0.7694644203726145D1*t367
     #-t370)+0.3109D-1*t420-0.1321010150222857D-2*t424*t100-0.1D1*t88
     #*t428*t445*t446)*t205+t190*(-0.59768D-2*t399*t197+0.4770864D-3
     #*t384*t399-0.26473536D-5*t389*t399)
      t502 = sigmaaa*t20
      t506 = t16/t222
      vsigmaaa(i) = -0.9305257363491D0*t5*(0.20292D-2*t9*t13
     #+0.158224D-4*t502*t23-0.957568D-7*t506*t227)-0.3109D-1*t33*t45*
     #(0.46974D0*t9*t48-0.292892D0*t502*t52+0.397888D-1*t506*t269)
     #+t190*(0.22413D-2*t9*t197-0.1789074D-3*t384*t9+0.9927576D-6*t389
     #*t9)
      vsigmaab(i) = 0.D0
      t532 = sigmabb*t75
      t536 = t71/t409
      vsigmabb(i) = -0.9305257363491D0*t60*(0.20292D-2*t64*t68
     #+0.158224D-4*t532*t78-0.957568D-7*t536*t414)-0.3109D-1*t88*t100*
     #(0.46974D0*t64*t103-0.292892D0*t532*t107+0.397888D-1*t536*t456)
     #+t190*(0.22413D-2*t64*t197-0.1789074D-3*t384*t64+0.9927576D-6
     #*t389*t64)
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vrhob(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      vsigmaab(i) = 0.0d0
      vsigmabb(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmabb
      t2 = rhob**(1.D0/3.D0)
      t3 = t2*rhob
      t4 = rhob**2
      t5 = t2**2
      t7 = 1/t5/t4
      t8 = sigmabb*t7
      t10 = 1.D0+0.4D-2*t8
      t11 = 1/t10
      t14 = sigmabb**2
      t15 = t4**2
      t18 = 1/t2/t15/rhob
      t19 = t14*t18
      t20 = t10**2
      t21 = 1/t20
      t24 = 0.8094D0+0.20292D-2*t8*t11+0.119696D-4*t19*t21
      t27 = 1/rhob
      t28 = t27**(1.D0/3.D0)
      t30 = 1.D0+0.1274696188700087D0*t28
      t31 = rhob*t30
      t32 = t27**(1.D0/6.D0)
      t35 = dsqrt(t27)
      t37 = t28**2
      t39 = 0.1112037486309468D2*t32+0.3844746237447211D1*t28
     #+0.1644733775567609D1*t35+0.2405871291288192D0*t37
      t42 = 1.D0+0.321646831778707D2/t39
      t43 = dlog(t42)
      t45 = 1.D0+0.2D0*t8
      t46 = 1/t45
      t49 = t45**2
      t50 = 1/t49
      t53 = 0.1737D0+0.46974D0*t8*t46-0.99472D-1*t19*t50
      t54 = t43*t53
      zk(i) = -0.9305257363491D0*t3*t24-0.3109D-1*t31*t54
      vrhoa(i) = 0.D0
      t59 = t4*rhob
      t62 = sigmabb/t5/t59
      t68 = t14/t2/t15/t4
      t71 = t14*sigmabb
      t72 = t15**2
      t75 = t71/t72/rhob
      t77 = 1/t20/t10
      t80 = -0.54112D-2*t62*t11-0.4219306666666667D-4*t68*t21
     #+0.2553514666666667D-6*t75*t77
      t83 = t30*t43
      t86 = 1/t37
      t87 = t27*t86
      t90 = t39**2
      t91 = 1/t90
      t92 = t31*t91
      t93 = t32**2
      t94 = t93**2
      t95 = t94*t32
      t96 = 1/t95
      t97 = 1/t4
      t102 = 1/t35
      t105 = 1/t28
      t108 = -0.1853395810515781D1*t96*t97-0.128158207914907D1*t86*t97
     #-0.8223668877838045D0*t102*t97-0.1603914194192128D0*t105*t97
      t109 = 1/t42
      t110 = t108*t109
      t111 = t110*t53
      t119 = 1/t49/t45
      t122 = -0.125264D1*t62*t46+0.7810453333333333D0*t68*t50
     #-0.1061034666666667D0*t75*t119
      t123 = t43*t122
      vrhob(i) = -0.12407009817988D1*t2*t24-0.9305257363491D0*t3*t80
     #-0.3109D-1*t83*t53+0.1321010150222857D-2*t87*t54+0.1D1*t92*t111
     #-0.3109D-1*t31*t123
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      t128 = sigmabb*t18
      t131 = 1/t72
      t132 = t14*t131
      vsigmabb(i) = -0.9305257363491D0*t3*(0.20292D-2*t7*t11
     #+0.158224D-4*t128*t21-0.957568D-7*t132*t77)-0.3109D-1*t31*t43*
     #(0.46974D0*t7*t46-0.292892D0*t128*t50+0.397888D-1*t132*t119)
      v2rhoa2(i) = 0.D0
      v2rhoab(i) = 0.D0
      t155 = sigmabb/t5/t15
      t161 = t14/t2/t15/t59
      t164 = t72*t4
      t166 = t71/t164
      t169 = t14**2
      t173 = t169/t5/t72/t15
      t174 = t20**2
      t175 = 1/t174
      t186 = 1/t59
      t188 = 1/t37/t27
      t200 = t108**2
      t207 = 1/t15
      t233 = t90**2
      t236 = t42**2
      t251 = t49**2
      t252 = 1/t251
      s1 = -0.4135669939329333D0/t5*t24-0.24814019635976D1*t2*t80
     #-0.9305257363491D0*t3*(0.1984106666666667D-1*t155*t11
     #+0.2095032888888889D-3*t161*t21-0.3198281955555556D-5*t166*t77
     #+0.8171246933333333D-8*t173*t175)+0.2D1*t30*t91*t111-0.6218D-1
     #*t83*t122+0.8806734334819047D-3*t186*t188*t54
      v2rhob2(i) = s1-0.8497974591333914D-1*t87*t91*t111
     #+0.2642020300445714D-2*t87*t123-0.2D1*t31/t90/t39*t200*t109*t53
     #+0.1D1*t92*(-0.1544496508763151D1/t95/t27*t207
     #+0.3706791621031562D1*t96*t186-0.854388052766047D0*t188*t207
     #+0.2563164158298141D1*t86*t186-0.4111834438919023D0/t35/t27*t207
     #+0.1644733775567609D1*t102*t186-0.5346380647307093D-1/t28/t27
     #*t207+0.3207828388384256D0*t105*t186)*t109*t53
     #+0.321646831778707D2*t31/t233*t200/t236*t53+0.2D1*t92*t110*t122
     #-0.3109D-1*t31*t43*(0.4593013333333333D1*t155*t46
     #-0.5614695111111111D1*t161*t50+0.1788046222222222D1*t166*t119
     #-0.1697655466666667D0*t173*t252)
      v2sigmaaa2(i) = 0.D0
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      t261 = sigmabb*t131
      t266 = t14/t5/t164
      v2sigmabb2(i) = -0.9305257363491D0*t3*(0.77056D-5*t18*t21
     #-0.3180928D-6*t261*t77+0.11490816D-8*t266*t175)-0.3109D-1*t31
     #*t43*(-0.38684D0*t18*t50+0.1967344D0*t261*t119-0.2387328D-1*t266
     #*t252)
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = rhoa**(1.D0/3.D0)
      t3 = t2*rhoa
      t4 = rhoa**2
      t5 = t2**2
      t7 = 1/t5/t4
      t8 = sigmaaa*t7
      t10 = 1.D0+0.4D-2*t8
      t11 = 1/t10
      t14 = sigmaaa**2
      t15 = t4**2
      t18 = 1/t2/t15/rhoa
      t19 = t14*t18
      t20 = t10**2
      t21 = 1/t20
      t24 = 0.8094D0+0.20292D-2*t8*t11+0.119696D-4*t19*t21
      t27 = 1/rhoa
      t28 = t27**(1.D0/3.D0)
      t30 = 1.D0+0.1274696188700087D0*t28
      t31 = rhoa*t30
      t32 = t27**(1.D0/6.D0)
      t35 = dsqrt(t27)
      t37 = t28**2
      t39 = 0.1112037486309468D2*t32+0.3844746237447211D1*t28
     #+0.1644733775567609D1*t35+0.2405871291288192D0*t37
      t42 = 1.D0+0.321646831778707D2/t39
      t43 = dlog(t42)
      t45 = 1.D0+0.2D0*t8
      t46 = 1/t45
      t49 = t45**2
      t50 = 1/t49
      t53 = 0.1737D0+0.46974D0*t8*t46-0.99472D-1*t19*t50
      t54 = t43*t53
      zk(i) = -0.9305257363491D0*t3*t24-0.3109D-1*t31*t54
      t59 = t4*rhoa
      t62 = sigmaaa/t5/t59
      t68 = t14/t2/t15/t4
      t71 = t14*sigmaaa
      t72 = t15**2
      t75 = t71/t72/rhoa
      t77 = 1/t20/t10
      t80 = -0.54112D-2*t62*t11-0.4219306666666667D-4*t68*t21
     #+0.2553514666666667D-6*t75*t77
      t83 = t30*t43
      t86 = 1/t37
      t87 = t27*t86
      t90 = t39**2
      t91 = 1/t90
      t92 = t31*t91
      t93 = t32**2
      t94 = t93**2
      t95 = t94*t32
      t96 = 1/t95
      t97 = 1/t4
      t102 = 1/t35
      t105 = 1/t28
      t108 = -0.1853395810515781D1*t96*t97-0.128158207914907D1*t86*t97
     #-0.8223668877838045D0*t102*t97-0.1603914194192128D0*t105*t97
      t109 = 1/t42
      t110 = t108*t109
      t111 = t110*t53
      t119 = 1/t49/t45
      t122 = -0.125264D1*t62*t46+0.7810453333333333D0*t68*t50
     #-0.1061034666666667D0*t75*t119
      t123 = t43*t122
      vrhoa(i) = -0.12407009817988D1*t2*t24-0.9305257363491D0*t3*t80
     #-0.3109D-1*t83*t53+0.1321010150222857D-2*t87*t54+0.1D1*t92*t111
     #-0.3109D-1*t31*t123
      vrhob(i) = 0.D0
      t128 = sigmaaa*t18
      t131 = 1/t72
      t132 = t14*t131
      vsigmaaa(i) = -0.9305257363491D0*t3*(0.20292D-2*t7*t11
     #+0.158224D-4*t128*t21-0.957568D-7*t132*t77)-0.3109D-1*t31*t43*
     #(0.46974D0*t7*t46-0.292892D0*t128*t50+0.397888D-1*t132*t119)
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      t155 = sigmaaa/t5/t15
      t161 = t14/t2/t15/t59
      t164 = t72*t4
      t166 = t71/t164
      t169 = t14**2
      t173 = t169/t5/t72/t15
      t174 = t20**2
      t175 = 1/t174
      t186 = 1/t59
      t188 = 1/t37/t27
      t200 = t108**2
      t207 = 1/t15
      t233 = t90**2
      t236 = t42**2
      t251 = t49**2
      t252 = 1/t251
      s1 = -0.4135669939329333D0/t5*t24-0.24814019635976D1*t2*t80
     #-0.9305257363491D0*t3*(0.1984106666666667D-1*t155*t11
     #+0.2095032888888889D-3*t161*t21-0.3198281955555556D-5*t166*t77
     #+0.8171246933333333D-8*t173*t175)+0.2D1*t30*t91*t111-0.6218D-1
     #*t83*t122+0.8806734334819047D-3*t186*t188*t54
      v2rhoa2(i) = s1-0.8497974591333914D-1*t87*t91*t111
     #+0.2642020300445714D-2*t87*t123-0.2D1*t31/t90/t39*t200*t109*t53
     #+0.1D1*t92*(-0.1544496508763151D1/t95/t27*t207
     #+0.3706791621031562D1*t96*t186-0.854388052766047D0*t188*t207
     #+0.2563164158298141D1*t86*t186-0.4111834438919023D0/t35/t27*t207
     #+0.1644733775567609D1*t102*t186-0.5346380647307093D-1/t28/t27
     #*t207+0.3207828388384256D0*t105*t186)*t109*t53
     #+0.321646831778707D2*t31/t233*t200/t236*t53+0.2D1*t92*t110*t122
     #-0.3109D-1*t31*t43*(0.4593013333333333D1*t155*t46
     #-0.5614695111111111D1*t161*t50+0.1788046222222222D1*t166*t119
     #-0.1697655466666667D0*t173*t252)
      v2rhob2(i) = 0.D0
      v2rhoab(i) = 0.D0
      t261 = sigmaaa*t131
      t266 = t14/t5/t164
      v2sigmaaa2(i) = -0.9305257363491D0*t3*(0.77056D-5*t18*t21
     #-0.3180928D-6*t261*t77+0.11490816D-8*t266*t175)-0.3109D-1*t31
     #*t43*(-0.38684D0*t18*t50+0.1967344D0*t261*t119-0.2387328D-1*t266
     #*t252)
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      v2sigmabb2(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rhoa**(1.D0/3.D0)
      t5 = t4*rhoa
      t6 = rhoa**2
      t7 = t4**2
      t9 = 1/t7/t6
      t10 = sigmaaa*t9
      t12 = 1.D0+0.4D-2*t10
      t13 = 1/t12
      t16 = sigmaaa**2
      t17 = t6**2
      t20 = 1/t4/t17/rhoa
      t21 = t16*t20
      t22 = t12**2
      t23 = 1/t22
      t26 = 0.8094D0+0.20292D-2*t10*t13+0.119696D-4*t21*t23
      t29 = 1/rhoa
      t30 = t29**(1.D0/3.D0)
      t32 = 1.D0+0.1274696188700087D0*t30
      t33 = rhoa*t32
      t34 = t29**(1.D0/6.D0)
      t37 = dsqrt(t29)
      t39 = t30**2
      t41 = 0.1112037486309468D2*t34+0.3844746237447211D1*t30
     #+0.1644733775567609D1*t37+0.2405871291288192D0*t39
      t44 = 1.D0+0.321646831778707D2/t41
      t45 = dlog(t44)
      t47 = 1.D0+0.2D0*t10
      t48 = 1/t47
      t51 = t47**2
      t52 = 1/t51
      t55 = 0.1737D0+0.46974D0*t10*t48-0.99472D-1*t21*t52
      t56 = t45*t55
      t59 = rhob**(1.D0/3.D0)
      t60 = t59*rhob
      t61 = rhob**2
      t62 = t59**2
      t64 = 1/t62/t61
      t65 = sigmabb*t64
      t67 = 1.D0+0.4D-2*t65
      t68 = 1/t67
      t71 = sigmabb**2
      t72 = t61**2
      t75 = 1/t59/t72/rhob
      t76 = t71*t75
      t77 = t67**2
      t78 = 1/t77
      t81 = 0.8094D0+0.20292D-2*t65*t68+0.119696D-4*t76*t78
      t84 = 1/rhob
      t85 = t84**(1.D0/3.D0)
      t87 = 1.D0+0.1274696188700087D0*t85
      t88 = rhob*t87
      t89 = t84**(1.D0/6.D0)
      t92 = dsqrt(t84)
      t94 = t85**2
      t96 = 0.1112037486309468D2*t89+0.3844746237447211D1*t85
     #+0.1644733775567609D1*t92+0.2405871291288192D0*t94
      t99 = 1.D0+0.321646831778707D2/t96
      t100 = dlog(t99)
      t102 = 1.D0+0.2D0*t65
      t103 = 1/t102
      t106 = t102**2
      t107 = 1/t106
      t110 = 0.1737D0+0.46974D0*t65*t103-0.99472D-1*t76*t107
      t111 = t100*t110
      t114 = rhoa+rhob
      t115 = 1/t114
      t116 = t115**(1.D0/3.D0)
      t118 = 1.D0+0.1325688999052018D0*t116
      t119 = t115**(1.D0/6.D0)
      t122 = dsqrt(t115)
      t124 = t116**2
      t126 = 0.598255043577108D1*t119+0.2225569421150687D1*t116
     #+0.8004286349993634D0*t122+0.1897004325747559D0*t124
      t129 = 1.D0+0.160818243221511D2/t126
      t130 = dlog(t129)
      t132 = 0.62182D-1*t118*t130
      t134 = 1.D0+0.6901399211255825D-1*t116
      t139 = 0.8157414703487641D1*t119+0.2247591863577616D1*t116
     #+0.4300972471276643D0*t122+0.1911512595127338D0*t124
      t142 = 1.D0+0.2960857464321668D2/t139
      t143 = dlog(t142)
      t144 = t134*t143
      t146 = rhoa-1.D0*rhob
      t147 = t146*t115
      t148 = 1.D0+t147
      t149 = t148**(1.D0/3.D0)
      t152 = 1.D0-1.D0*t147
      t153 = t152**(1.D0/3.D0)
      t155 = t149*t148+t153*t152-2.D0
      t156 = t146**2
      t157 = t156**2
      t158 = t114**2
      t159 = t158**2
      t160 = 1/t159
      t161 = t157*t160
      t163 = 1.D0-1.D0*t161
      t164 = t155*t163
      t166 = 0.3799574853701528D-1*t144*t164
      t168 = 1.D0+0.1274696188700087D0*t116
      t173 = 0.1112037486309468D2*t119+0.3844746237447211D1*t116
     #+0.1644733775567609D1*t122+0.2405871291288192D0*t124
      t176 = 1.D0+0.321646831778707D2/t173
      t177 = dlog(t176)
      t180 = -0.3109D-1*t168*t177+t132
      t181 = t180*t155
      t183 = 0.1923661050931536D1*t181*t161
      t190 = t114*(-t132+t166+t183)+0.3109D-1*t33*t45+0.3109D-1*t88*t100
      t193 = 0.5D0*t10+0.5D0*t65
      t196 = 1.D0+0.3D-2*t10+0.3D-2*t65
      t197 = 1/t196
      t200 = t193**2
      t201 = t196**2
      t202 = 1/t201
      t205 = 0.9454D0+0.44826D-2*t193*t197-0.1654596D-3*t200*t202
      zk(i) = -0.9305257363491D0*t5*t26-0.3109D-1*t33*t56
     #-0.9305257363491D0*t60*t81-0.3109D-1*t88*t111+t190*t205
      t209 = t6*rhoa
      t211 = 1/t7/t209
      t212 = sigmaaa*t211
      t217 = 1/t4/t17/t6
      t218 = t16*t217
      t221 = t16*sigmaaa
      t222 = t17**2
      t224 = 1/t222/rhoa
      t225 = t221*t224
      t227 = 1/t22/t12
      t230 = -0.54112D-2*t212*t13-0.4219306666666667D-4*t218*t23
     #+0.2553514666666667D-6*t225*t227
      t233 = t32*t45
      t236 = 1/t39
      t237 = t29*t236
      t240 = t41**2
      t241 = 1/t240
      t242 = t33*t241
      t243 = t34**2
      t244 = t243**2
      t245 = t244*t34
      t246 = 1/t245
      t247 = 1/t6
      t252 = 1/t37
      t255 = 1/t30
      t258 = -0.1853395810515781D1*t246*t247-0.128158207914907D1*t236
     #*t247-0.8223668877838045D0*t252*t247-0.1603914194192128D0*t255
     #*t247
      t259 = 1/t44
      t260 = t258*t259
      t261 = t260*t55
      t269 = 1/t51/t47
      t272 = -0.125264D1*t212*t48+0.7810453333333333D0*t218*t52
     #-0.1061034666666667D0*t225*t269
      t273 = t45*t272
      t276 = 1/t124
      t277 = 1/t158
      t278 = t276*t277
      t279 = t278*t130
      t280 = 0.2747799777968419D-2*t279
      t281 = t126**2
      t282 = 1/t281
      t283 = t118*t282
      t284 = t119**2
      t285 = t284**2
      t286 = t285*t119
      t287 = 1/t286
      t288 = t287*t277
      t291 = 1/t122
      t292 = t291*t277
      t294 = 1/t116
      t295 = t294*t277
      t297 = -0.99709173929518D0*t288-0.7418564737168958D0*t278
     #-0.4002143174996817D0*t292-0.1264669550498372D0*t295
      t298 = 1/t129
      t300 = t283*t297*t298
      t301 = 0.1D1*t300
      t302 = t143*t155
      t303 = t302*t163
      t304 = t278*t303
      t305 = 0.8740794299481065D-3*t304
      t306 = t139**2
      t307 = 1/t306
      t308 = t134*t307
      t313 = -0.135956911724794D1*t288-0.7491972878592054D0*t278
     #-0.2150486235638321D0*t292-0.1274341730084892D0*t295
      t314 = t308*t313
      t315 = 1/t142
      t316 = t315*t155
      t317 = t316*t163
      t318 = t314*t317
      t319 = 0.1124999956683108D1*t318
      t320 = t146*t277
      t321 = 1.D0*t320
      t322 = t115-t321
      t325 = 1.D0*t115
      t326 = -t325+t320
      t329 = 0.1333333333333333D1*t149*t322+0.1333333333333333D1*t153
     #*t326
      t331 = t144*t329*t163
      t332 = 0.3799574853701528D-1*t331
      t333 = t156*t146
      t334 = t333*t160
      t335 = 4.D0*t334
      t337 = 1/t159/t114
      t338 = t157*t337
      t339 = 4.D0*t338
      t340 = -t335+t339
      t342 = t144*t155*t340
      t343 = 0.3799574853701528D-1*t342
      t346 = t173**2
      t347 = 1/t346
      t348 = t168*t347
      t353 = -0.1853395810515781D1*t288-0.128158207914907D1*t278
     #-0.8223668877838045D0*t292-0.1603914194192128D0*t295
      t354 = 1/t176
      t360 = 0.1321010150222857D-2*t278*t177+0.1D1*t348*t353*t354
     #-0.2747799777968419D-2*t279-0.1D1*t300
      t361 = t360*t155
      t362 = t361*t161
      t363 = 0.1923661050931536D1*t362
      t364 = t180*t329
      t365 = t364*t161
      t366 = 0.1923661050931536D1*t365
      t367 = t181*t334
      t369 = t181*t338
      t370 = 0.7694644203726145D1*t369
      t377 = t241*t258*t259
      t380 = -t132+t166+t183+t114*(t280+t301-t305-t319+t332+t343+t363
     #+t366+0.7694644203726145D1*t367-t370)+0.3109D-1*t233
     #-0.1321010150222857D-2*t237*t45-0.1D1*t33*t377
      t384 = t193*t202
      t388 = 1/t201/t196
      t389 = t200*t388
      t392 = -0.59768D-2*t212*t197+0.4770864D-3*t384*t212
     #-0.26473536D-5*t389*t212
      vrhoa(i) = -0.12407009817988D1*t4*t26-0.9305257363491D0*t5*t230
     #-0.3109D-1*t233*t55+0.1321010150222857D-2*t237*t56+0.1D1*t242
     #*t261-0.3109D-1*t33*t273+t380*t205+t190*t392
      t396 = t61*rhob
      t398 = 1/t62/t396
      t399 = sigmabb*t398
      t404 = 1/t59/t72/t61
      t405 = t71*t404
      t408 = t71*sigmabb
      t409 = t72**2
      t411 = 1/t409/rhob
      t412 = t408*t411
      t414 = 1/t77/t67
      t417 = -0.54112D-2*t399*t68-0.4219306666666667D-4*t405*t78
     #+0.2553514666666667D-6*t412*t414
      t420 = t87*t100
      t423 = 1/t94
      t424 = t84*t423
      t427 = t96**2
      t428 = 1/t427
      t429 = t88*t428
      t430 = t89**2
      t431 = t430**2
      t432 = t431*t89
      t433 = 1/t432
      t434 = 1/t61
      t439 = 1/t92
      t442 = 1/t85
      t445 = -0.1853395810515781D1*t433*t434-0.128158207914907D1*t423
     #*t434-0.8223668877838045D0*t439*t434-0.1603914194192128D0*t442
     #*t434
      t446 = 1/t99
      t447 = t445*t446
      t448 = t447*t110
      t456 = 1/t106/t102
      t459 = -0.125264D1*t399*t103+0.7810453333333333D0*t405*t107
     #-0.1061034666666667D0*t412*t456
      t460 = t100*t459
      t463 = -t325-t321
      t466 = t115+t320
      t469 = 0.1333333333333333D1*t149*t463+0.1333333333333333D1*t153
     #*t466
      t471 = t144*t469*t163
      t472 = 0.3799574853701528D-1*t471
      t473 = t335+t339
      t475 = t144*t155*t473
      t476 = 0.3799574853701528D-1*t475
      t477 = t180*t469
      t478 = t477*t161
      t479 = 0.1923661050931536D1*t478
      t487 = t428*t445*t446
      t490 = -t132+t166+t183+t114*(t280+t301-t305-t319+t472+t476+t363
     #+t479-0.7694644203726145D1*t367-t370)+0.3109D-1*t420
     #-0.1321010150222857D-2*t424*t100-0.1D1*t88*t487
      t498 = -0.59768D-2*t399*t197+0.4770864D-3*t384*t399
     #-0.26473536D-5*t389*t399
      vrhob(i) = -0.12407009817988D1*t59*t81-0.9305257363491D0*t60
     #*t417-0.3109D-1*t420*t110+0.1321010150222857D-2*t424*t111+0.1D1
     #*t429*t448-0.3109D-1*t88*t460+t490*t205+t190*t498
      t502 = sigmaaa*t20
      t505 = 1/t222
      t506 = t16*t505
      t509 = 0.20292D-2*t9*t13+0.158224D-4*t502*t23-0.957568D-7*t506
     #*t227
      t518 = 0.46974D0*t9*t48-0.292892D0*t502*t52+0.397888D-1*t506*t269
      t519 = t45*t518
      t528 = 0.22413D-2*t9*t197-0.1789074D-3*t384*t9+0.9927576D-6*t389
     #*t9
      vsigmaaa(i) = -0.9305257363491D0*t5*t509-0.3109D-1*t33*t519+t190
     #*t528
      vsigmaab(i) = 0.D0
      t532 = sigmabb*t75
      t535 = 1/t409
      t536 = t71*t535
      t539 = 0.20292D-2*t64*t68+0.158224D-4*t532*t78-0.957568D-7*t536
     #*t414
      t548 = 0.46974D0*t64*t103-0.292892D0*t532*t107+0.397888D-1*t536
     #*t456
      t549 = t100*t548
      t558 = 0.22413D-2*t64*t197-0.1789074D-3*t384*t64+0.9927576D-6
     #*t389*t64
      vsigmabb(i) = -0.9305257363491D0*t60*t539-0.3109D-1*t88*t549
     #+t190*t558
      t562 = t32*t241
      t566 = 1/t240/t41
      t567 = t258**2
      t574 = 1/t209
      t576 = 1/t39/t29
      t577 = t574*t576
      t582 = 1/t17
      t603 = -0.1544496508763151D1/t245/t29*t582+0.3706791621031562D1
     #*t246*t574-0.854388052766047D0*t576*t582+0.2563164158298141D1
     #*t236*t574-0.4111834438919023D0/t37/t29*t582
     #+0.1644733775567609D1*t252*t574-0.5346380647307093D-1/t30/t29
     #*t582+0.3207828388384256D0*t255*t574
      t609 = t240**2
      t610 = 1/t609
      t612 = t44**2
      t613 = 1/t612
      t617 = 0.5495599555936838D-2*t279
      t618 = 0.2D1*t300
      t621 = 0.1538928840745229D2*t369
      t623 = 0.1748158859896213D-2*t304
      t624 = 0.2249999913366216D1*t318
      t625 = 0.3847322101863073D1*t362
      t628 = 1/t124/t115*t160
      t629 = t628*t130
      t630 = 0.1831866518645613D-2*t629
      t632 = 1/t158/t114
      t633 = t276*t632
      t634 = t633*t130
      t635 = 0.5495599555936838D-2*t634
      t638 = t278*t282*t297*t298
      t639 = 0.8837926660346786D-1*t638
      t641 = t278*t302*t340
      t643 = t281**2
      t646 = t297**2
      t647 = t129**2
      t650 = t118/t643*t646/t647
      t651 = 0.160818243221511D2*t650
      t652 = t156*t160
      t653 = 12.D0*t652
      t654 = t333*t337
      t655 = 32.D0*t654
      t658 = t157/t159/t158
      t659 = 20.D0*t658
      t664 = t149**2
      t665 = 1/t664
      t666 = t322**2
      t669 = 2.D0*t277
      t671 = 2.D0*t146*t632
      t672 = -t669+t671
      t675 = t153**2
      t676 = 1/t675
      t677 = t326**2
      t683 = 0.4444444444444444D0*t665*t666+0.1333333333333333D1*t149
     #*t672+0.4444444444444444D0*t676*t677-0.1333333333333333D1*t153
     #*t672
      t689 = t278*t143*t329*t163
      t692 = 0.5827196199654043D-3*t628*t303
      t694 = 0.1748158859896213D-2*t633*t303
      t699 = 0.5176049209143758D-1*t278*t307*t313*t315*t164
      t703 = t313**2
      t706 = 0.2249999913366216D1*t134/t306/t139*t703*t317
      t709 = 1/t286/t115*t160
      t711 = t287*t632
      t717 = 1/t122/t115*t160
      t719 = t291*t632
      t723 = 1/t116/t115*t160
      t725 = t294*t632
      t730 = 0.1124999956683108D1*t308*(-0.1132974264373283D1*t709
     #+0.271913823449588D1*t711-0.4994648585728036D0*t628
     #+0.1498394575718411D1*t633-0.1075243117819161D0*t717
     #+0.4300972471276643D0*t719-0.4247805766949639D-1*t723
     #+0.2548683460169784D0*t725)*t317
      t733 = t314*t315*t329*t163
      t735 = t630-t635-t639-0.1748158859896213D-2*t641+t651
     #+0.3799574853701528D-1*t144*t155*(-t653+t655-t659)
     #+0.3799574853701528D-1*t144*t683*t163-0.1748158859896213D-2*t689
     #-t692+t694+t699+t706-t730-0.2249999913366216D1*t733
      t737 = t314*t316*t340
      t742 = t364*t338
      t744 = t364*t334
      t746 = t181*t654
      t748 = t181*t652
      t749 = 0.2308393261117844D2*t748
      t761 = t353**2
      t777 = t346**2
      t780 = t176**2
      t792 = t118/t281/t126*t646*t298
      t804 = t283*(-0.8309097827459833D0*t709+0.199418347859036D1*t711
     #-0.4945709824779306D0*t628+0.1483712947433792D1*t633
     #-0.2001071587498409D0*t717+0.8004286349993634D0*t719
     #-0.4215565168327908D-1*t723+0.2529339100996745D0*t725)*t298
      t807 = 0.8806734334819047D-3*t628*t177-0.2642020300445714D-2
     #*t633*t177-0.8497974591333914D-1*t278*t347*t353*t354-0.2D1*t168
     #/t346/t173*t761*t354+0.1D1*t348*(-0.1544496508763151D1*t709
     #+0.3706791621031562D1*t711-0.854388052766047D0*t628
     #+0.2563164158298141D1*t633-0.4111834438919023D0*t717
     #+0.1644733775567609D1*t719-0.5346380647307093D-1*t723
     #+0.3207828388384256D0*t725)*t354+0.321646831778707D2*t168/t777
     #*t761/t780-0.1831866518645613D-2*t629+0.5495599555936838D-2*t634
     #+0.8837926660346786D-1*t638+0.2D1*t792-0.1D1*t804
     #-0.160818243221511D2*t650
      t810 = 0.1923661050931536D1*t807*t155*t161
      t812 = 0.3847322101863073D2*t181*t658
      t817 = t360*t329*t161
      t819 = 0.1D1*t804
      t820 = t306**2
      t824 = t142**2
      t829 = 0.3330964519106732D2*t134/t820*t703/t824*t155*t163
      t830 = 0.2D1*t792
      t831 = t361*t334
      t834 = 0.1538928840745229D2*t361*t338
      t835 = -0.2249999913366216D1*t737+0.1923661050931536D1*t180*t683
     #*t161-0.1538928840745229D2*t742+0.1538928840745229D2*t744
     #-0.6155715362980916D2*t746+t749+t810+t812+0.7599149707403056D-1
     #*t144*t329*t340+0.3847322101863073D1*t817+t819-t829-t830
     #+0.1538928840745229D2*t831-t834
      t838 = -0.2D1*t562*t260+0.2D1*t33*t566*t567*t259
     #+0.8497974591333914D-1*t237*t377-0.8806734334819047D-3*t577*t45
     #-0.1D1*t33*t241*t603*t259+0.7599149707403056D-1*t331
     #-0.321646831778707D2*t33*t610*t567*t613+t617+t618
     #+0.3847322101863073D1*t365+0.1538928840745229D2*t367-t621
     #+0.7599149707403056D-1*t342-t623-t624+t625+t114*(t735+t835)
      t842 = sigmaaa/t7/t17
      t848 = t16/t4/t17/t209
      t851 = t193*t388
      t856 = t201**2
      t858 = t200/t856
      t876 = t222*t6
      t878 = t221/t876
      t881 = t16**2
      t885 = t881/t7/t222/t17
      t886 = t22**2
      t887 = 1/t886
      t899 = t51**2
      t900 = 1/t899
      s1 = 2.D0*t380*t392+t838*t205+t190*(0.2191493333333333D-1*t842
     #*t197-0.6839296D-3*t848*t202+0.14692992D-4*t851*t848
     #-0.17493168D-2*t384*t842-0.635364864D-7*t858*t848+0.97069632D-5
     #*t389*t842)-0.6218D-1*t233*t272-0.4135669939329333D0/t7*t26
     #-0.24814019635976D1*t4*t230-0.9305257363491D0*t5*
     #(0.1984106666666667D-1*t842*t13+0.2095032888888889D-3*t848*t23
     #-0.3198281955555556D-5*t878*t227+0.8171246933333333D-8*t885*t887
     #)-0.3109D-1*t33*t45*(0.4593013333333333D1*t842*t48
     #-0.5614695111111111D1*t848*t52+0.1788046222222222D1*t878*t269
     #-0.1697655466666667D0*t885*t900)
      v2rhoa2(i) = s1+0.2D1*t562*t261+0.2642020300445714D-2*t237*t273
     #+0.8806734334819047D-3*t577*t56-0.8497974591333914D-1*t237*t241
     #*t261-0.2D1*t33*t566*t567*t259*t55+0.1D1*t242*t603*t259*t55
     #+0.2D1*t242*t260*t272+0.321646831778707D2*t33*t610*t567*t613*t55
      t937 = sigmabb/t62/t72
      t943 = t71/t59/t72/t396
      t962 = t409*t61
      t964 = t408/t962
      t967 = t71**2
      t971 = t967/t62/t409/t72
      t972 = t77**2
      t973 = 1/t972
      t987 = 1/t396
      t989 = 1/t94/t84
      t990 = t987*t989
      t995 = t87*t428
      t1000 = 1/t72
      t1021 = -0.1544496508763151D1/t432/t84*t1000
     #+0.3706791621031562D1*t433*t987-0.854388052766047D0*t989*t1000
     #+0.2563164158298141D1*t423*t987-0.4111834438919023D0/t92/t84
     #*t1000+0.1644733775567609D1*t439*t987-0.5346380647307093D-1/t85
     #/t84*t1000+0.3207828388384256D0*t442*t987
      t1027 = 1/t427/t96
      t1029 = t445**2
      t1034 = t427**2
      t1035 = 1/t1034
      t1037 = t99**2
      t1038 = 1/t1037
      t1068 = t477*t338
      t1070 = t477*t334
      t1073 = -0.1538928840745229D2*t1068-0.1538928840745229D2*t1070
     #+t630-t635-t639+t651-t692+t694+t699+t706-t730
     #+0.6155715362980916D2*t746+t749+t810
      t1078 = t463**2
      t1081 = t669+t671
      t1084 = t466**2
      t1090 = 0.4444444444444444D0*t665*t1078+0.1333333333333333D1
     #*t149*t1081+0.4444444444444444D0*t676*t1084-0.1333333333333333D1
     #*t153*t1081
      t1099 = t360*t469*t161
      t1105 = t278*t302*t473
      t1109 = t278*t143*t469*t163
      t1112 = t314*t316*t473
      t1116 = t314*t315*t469*t163
      t1118 = t812+t819-t829-t830-0.1538928840745229D2*t831-t834
     #+0.7599149707403056D-1*t144*t469*t473+0.3799574853701528D-1*t144
     #*t1090*t163+0.3799574853701528D-1*t144*t155*(-t653-t655-t659)
     #+0.3847322101863073D1*t1099+0.1923661050931536D1*t180*t1090*t161
     #-0.1748158859896213D-2*t1105-0.1748158859896213D-2*t1109
     #-0.2249999913366216D1*t1112-0.2249999913366216D1*t1116
      t1121 = -0.1D1*t88*t428*t1021*t446+0.2D1*t88*t1027*t1029*t446
     #-0.8806734334819047D-3*t990*t100+0.8497974591333914D-1*t424*t487
     #-0.2D1*t995*t447-0.321646831778707D2*t88*t1035*t1029*t1038+t617
     #+t618-0.1538928840745229D2*t367-t621-t623-t624+t625
     #+0.3847322101863073D1*t478+0.7599149707403056D-1*t471
     #+0.7599149707403056D-1*t475+t114*(t1073+t1118)
      t1129 = t106**2
      t1130 = 1/t1129
      s1 = 2.D0*t490*t498+t190*(0.2191493333333333D-1*t937*t197
     #-0.6839296D-3*t943*t202+0.14692992D-4*t851*t943-0.17493168D-2
     #*t384*t937-0.635364864D-7*t858*t943+0.97069632D-5*t389*t937)
     #-0.6218D-1*t420*t459-0.9305257363491D0*t60*
     #(0.1984106666666667D-1*t937*t68+0.2095032888888889D-3*t943*t78
     #-0.3198281955555556D-5*t964*t414+0.8171246933333333D-8*t971*t973
     #)-0.24814019635976D1*t59*t417-0.4135669939329333D0/t62*t81
     #-0.8497974591333914D-1*t424*t428*t448+0.8806734334819047D-3*t990
     #*t111
      v2rhob2(i) = s1+0.2642020300445714D-2*t424*t460+0.2D1*t995*t448
     #+0.1D1*t429*t1021*t446*t110-0.2D1*t88*t1027*t1029*t446*t110
     #+0.321646831778707D2*t88*t1035*t1029*t1038*t110+0.2D1*t429*t447
     #*t459+t1121*t205-0.3109D-1*t88*t100*(0.4593013333333333D1*t937
     #*t103-0.5614695111111111D1*t943*t107+0.1788046222222222D1*t964
     #*t456-0.1697655466666667D0*t971*t1130)
      t1160 = 0.4444444444444444D0*t665*t322*t463+0.2666666666666667D1
     #*t149*t146*t632+0.4444444444444444D0*t676*t326*t466
     #-0.2666666666666667D1*t153*t146*t632
      t1170 = -0.7694644203726145D1*t1068+0.7694644203726145D1*t1070
     #+t630-t635-t639-0.8740794299481065D-3*t641+t651
     #-0.8740794299481065D-3*t689+0.3799574853701528D-1*t144*t469*t340
     #+0.3799574853701528D-1*t144*t155*(t653-t659)
     #+0.3799574853701528D-1*t144*t1160*t163+0.3799574853701528D-1
     #*t144*t329*t473+0.1923661050931536D1*t180*t1160*t161-t692+t694
     #+t699+t706
      t1182 = -t730-0.1124999956683108D1*t733-0.1124999956683108D1
     #*t737-0.7694644203726145D1*t742-0.7694644203726145D1*t744
     #-0.2308393261117844D2*t748+t810+t812+0.1923661050931536D1*t817
     #+t819-t829-t830-t834+0.1923661050931536D1*t1099
     #-0.8740794299481065D-3*t1105-0.8740794299481065D-3*t1109
     #-0.1124999956683108D1*t1112-0.1124999956683108D1*t1116
      t1185 = t617+t618-t623-t624+t472+t476+t625+t479-t621+t332+t343
     #+t366+t114*(t1170+t1182)
      t1195 = t211*sigmabb*t398
      v2rhoab(i) = t1185*t205+t380*t498+t490*t392+t190*(-0.6839296D-3
     #*t212*t202*sigmabb*t398+0.14692992D-4*t851*sigmaaa*t1195
     #-0.635364864D-7*t858*sigmaaa*t1195)
      t1207 = sigmaaa*t217
      t1210 = t16*t224
      t1216 = t221/t7/t222/t209
      v2rhoasigmaaa(i) = -0.12407009817988D1*t4*t509-0.9305257363491D0
     #*t5*(-0.54112D-2*t211*t13-0.6274133333333333D-4*t1207*t23
     #+0.1103598933333333D-5*t1210*t227-0.30642176D-8*t1216*t887)
     #-0.3109D-1*t233*t518+0.1321010150222857D-2*t237*t519+0.1D1*t242
     #*t260*t518-0.3109D-1*t33*t45*(-0.125264D1*t211*t48
     #+0.1812618666666667D1*t1207*t52-0.6307285333333333D0*t1210*t269
     #+0.6366208D-1*t1216*t900)+t380*t528+t190*(-0.59768D-2*t211*t197
     #+0.2564736D-3*t1207*t202-0.5509872D-5*t851*t1207+0.4770864D-3
     #*t384*t211+0.238261824D-7*t858*t1207-0.26473536D-5*t389*t211)
      v2rhoasigmaab(i) = 0.D0
      t1260 = t212*t64
      v2rhoasigmabb(i) = t380*t558+t190*(0.2564736D-3*t212*t202*t64
     #-0.5509872D-5*t851*t1260+0.238261824D-7*t858*t1260)
      t1268 = t202*t9
      t1271 = t399*t9
      v2rhobsigmaaa(i) = t490*t528+t190*(0.2564736D-3*t399*t1268
     #-0.5509872D-5*t851*t1271+0.238261824D-7*t858*t1271)
      v2rhobsigmaab(i) = 0.D0
      t1282 = sigmabb*t404
      t1285 = t71*t411
      t1291 = t408/t62/t409/t396
      v2rhobsigmabb(i) = -0.12407009817988D1*t59*t539
     #-0.9305257363491D0*t60*(-0.54112D-2*t398*t68
     #-0.6274133333333333D-4*t1282*t78+0.1103598933333333D-5*t1285
     #*t414-0.30642176D-8*t1291*t973)-0.3109D-1*t420*t548
     #+0.1321010150222857D-2*t424*t549+0.1D1*t429*t447*t548-0.3109D-1
     #*t88*t100*(-0.125264D1*t398*t103+0.1812618666666667D1*t1282*t107
     #-0.6307285333333333D0*t1285*t456+0.6366208D-1*t1291*t1130)+t490
     #*t558+t190*(-0.59768D-2*t398*t197+0.2564736D-3*t1282*t202
     #-0.5509872D-5*t851*t1282+0.4770864D-3*t384*t398+0.238261824D-7
     #*t858*t1282-0.26473536D-5*t389*t398)
      t1333 = sigmaaa*t505
      t1338 = t16/t7/t876
      v2sigmaaa2(i) = -0.9305257363491D0*t5*(0.77056D-5*t20*t23
     #-0.3180928D-6*t1333*t227+0.11490816D-8*t1338*t887)-0.3109D-1*t33
     #*t45*(-0.38684D0*t20*t52+0.1967344D0*t1333*t269-0.2387328D-1
     #*t1338*t900)+t190*(-0.961776D-4*t20*t202+0.2066202D-5*t851*t20
     #-0.89348184D-8*t858*t20)
      v2sigmaaaab(i) = 0.D0
      t1364 = t9*t64
      v2sigmaaabb(i) = t190*(-0.961776D-4*t1268*t64+0.2066202D-5*t851
     #*t1364-0.89348184D-8*t858*t1364)
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      t1372 = sigmabb*t535
      t1377 = t71/t62/t962
      v2sigmabb2(i) = -0.9305257363491D0*t60*(0.77056D-5*t75*t78
     #-0.3180928D-6*t1372*t414+0.11490816D-8*t1377*t973)-0.3109D-1*t88
     #*t100*(-0.38684D0*t75*t107+0.1967344D0*t1372*t456-0.2387328D-1
     #*t1377*t1130)+t190*(-0.961776D-4*t75*t202+0.2066202D-5*t851*t75
     #-0.89348184D-8*t858*t75)
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vrhob(i) = 0.0d0
      v2rhoa2(i) = 0.0d0
      v2rhob2(i) = 0.0d0
      v2rhoab(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      vsigmaab(i) = 0.0d0
      vsigmabb(i) = 0.0d0
      v2rhoasigmaaa(i) = 0.0d0
      v2rhoasigmaab(i) = 0.0d0
      v2rhoasigmabb(i) = 0.0d0
      v2rhobsigmaaa(i) = 0.0d0
      v2rhobsigmaab(i) = 0.0d0
      v2rhobsigmabb(i) = 0.0d0
      v2sigmaaa2(i) = 0.0d0
      v2sigmaab2(i) = 0.0d0
      v2sigmabb2(i) = 0.0d0
      v2sigmaaaab(i) = 0.0d0
      v2sigmaaabb(i) = 0.0d0
      v2sigmaabbb(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end
      
      
      subroutine rks_xc_b97
     & (ideriv,npt,rhoa1,sigmaaa1,
     &  zk,vrhoa,vsigmaaa,
     &  v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
c
c     A.D. Becke
c     Density-functional thermochemistry. V. Systematic optimization of
c     exchange-correlation functionals
c     J. Chem. Phys. 107 (1997) 8554-8560
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit real*8 (a-h,o-z)
      integer ideriv,npt
      real*8 rhoa1(npt)
      real*8 sigmaaa1(npt)
      real*8 zk(npt),vrhoa(npt),vsigmaaa(npt)
      real*8 v2rhoa2(npt),v2rhoasigmaaa(npt),v2sigmaaa2(npt)
      parameter(tol=1.0d-20)
      
      if(ideriv.eq.0) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = rho**(1.D0/3.D0)
      t4 = rho**2
      t5 = t2**2
      t8 = sigma/t5/t4
      t10 = 1.D0+0.6349604207872798D-2*t8
      t14 = sigma**2
      t15 = t4**2
      t19 = t14/t2/t15/rho
      t20 = t10**2
      t27 = 1/rho
      t28 = t27**(1.D0/3.D0)
      t31 = rho*(1.D0+0.1606016560364007D0*t28)
      t32 = t27**(1.D0/6.D0)
      t35 = dsqrt(t27)
      t37 = t28**2
      t43 = dlog(1.D0+0.321646831778707D2/(0.1248219874679732D2*t32
     #+0.4844076716063854D1*t28+0.2326004811900819D1*t35
     #+0.3819082618690966D0*t37))
      t45 = 1.D0+0.3174802103936399D0*t8
      t49 = t45**2
      t68 = dlog(1.D0+0.160818243221511D2/(0.598255043577108D1*t32
     #+0.2225569421150687D1*t28+0.8004286349993634D0*t35
     #+0.1897004325747559D0*t37))
      t75 = 1.D0+0.9524406311809197D-2*t8
      t79 = t75**2
      zk(i) = -0.7385587663820224D0*t2*rho*(0.8094D0
     #+0.322115421465387D-2*t8/t10+0.3016150199764335D-4*t19/t20)
     #-0.3109D-1*t31*t43*(0.1737D0+0.745665770151542D0*t8/t45
     #-0.2506537333502856D0*t19/t49)+(-0.62182D-1*rho*(1.D0
     #+0.1325688999052018D0*t28)*t68+0.3109D-1*t31*t43)*(0.9454D0
     #+0.7115683955552651D-2*t8/t75-0.4169320658943715D-3*t19/t79)
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = rho**(1.D0/3.D0)
      t3 = t2*rho
      t4 = rho**2
      t5 = t2**2
      t7 = 1/t5/t4
      t8 = sigma*t7
      t10 = 1.D0+0.6349604207872798D-2*t8
      t11 = 1/t10
      t14 = sigma**2
      t15 = t4**2
      t18 = 1/t2/t15/rho
      t19 = t14*t18
      t20 = t10**2
      t21 = 1/t20
      t24 = 0.8094D0+0.322115421465387D-2*t8*t11+0.3016150199764335D-4
     #*t19*t21
      t27 = 1/rho
      t28 = t27**(1.D0/3.D0)
      t30 = 1.D0+0.1606016560364007D0*t28
      t31 = rho*t30
      t32 = t27**(1.D0/6.D0)
      t35 = dsqrt(t27)
      t37 = t28**2
      t39 = 0.1248219874679732D2*t32+0.4844076716063854D1*t28
     #+0.2326004811900819D1*t35+0.3819082618690966D0*t37
      t42 = 1.D0+0.321646831778707D2/t39
      t43 = dlog(t42)
      t45 = 1.D0+0.3174802103936399D0*t8
      t46 = 1/t45
      t49 = t45**2
      t50 = 1/t49
      t53 = 0.1737D0+0.745665770151542D0*t8*t46-0.2506537333502856D0
     #*t19*t50
      t54 = t43*t53
      t58 = 1.D0+0.1325688999052018D0*t28
      t64 = 0.598255043577108D1*t32+0.2225569421150687D1*t28
     #+0.8004286349993634D0*t35+0.1897004325747559D0*t37
      t67 = 1.D0+0.160818243221511D2/t64
      t68 = dlog(t67)
      t73 = -0.62182D-1*rho*t58*t68+0.3109D-1*t31*t43
      t75 = 1.D0+0.9524406311809197D-2*t8
      t76 = 1/t75
      t79 = t75**2
      t80 = 1/t79
      t83 = 0.9454D0+0.7115683955552651D-2*t8*t76
     #-0.4169320658943715D-3*t19*t80
      zk(i) = -0.7385587663820224D0*t3*t24-0.3109D-1*t31*t54+t73*t83
      t90 = sigma/t5/t4/rho
      t96 = t14/t2/t15/t4
      t100 = t15**2
      t103 = t14*sigma/t100/rho
      t105 = 1/t20/t10
      t111 = t30*t43
      t114 = 1/t37
      t115 = t27*t114
      t118 = t39**2
      t119 = 1/t118
      t121 = t32**2
      t122 = t121**2
      t125 = 1/t4
      t126 = 1/t122/t32*t125
      t128 = t114*t125
      t131 = 1/t35*t125
      t134 = 1/t28*t125
      t136 = -0.4160732915599108D1*t126-0.3229384477375903D1*t128
     #-0.2326004811900819D1*t131-0.5092110158254621D0*t134
      t137 = 1/t42
      t147 = 1/t49/t45
      t158 = t64**2
      t186 = 1/t79/t75
      vrhoa(i) = -0.9847450218426965D0*t2*t24-0.3692793831910112D0*t3*
     #(-0.1717948914482064D-1*t90*t11-0.2126397314118042D-3*t96*t21
     #+0.2042811733333333D-5*t103*t105)-0.3109D-1*t111*t53
     #+0.1664368495390566D-2*t115*t54+0.5D0*t31*t119*t136*t137*t53
     #-0.15545D-1*t31*t43*(-0.3976884107474891D1*t90*t46
     #+0.3936221825555298D1*t96*t50-0.8488277333333333D0*t103*t147)+(
     #-0.62182D-1*t58*t68+rho*(0.2747799777968419D-2*t128*t68+0.1D1
     #*t58/t158*(-0.99709173929518D0*t126-0.7418564737168958D0*t128
     #-0.4002143174996817D0*t131-0.1264669550498372D0*t134)/t67)
     #+0.3109D-1*t111-0.1664368495390566D-2*t115*t43-0.5D0*t31*t119
     #*t136*t137)*t83+t73*(-0.1897515721480707D-1*t90*t76
     #+0.2404364791914262D-2*t96*t80-0.211788288D-4*t103*t186)
      t193 = sigma*t18
      t197 = t14/t100
      vsigmaaa(i) = -0.7385587663820224D0*t3*(0.1288461685861548D-1*t7
     #*t11+0.1594797985588531D-3*t193*t21-0.15321088D-5*t197*t105)
     #-0.3109D-1*t31*t43*(0.2982663080606168D1*t7*t46
     #-0.2952166369166474D1*t193*t50+0.6366208D0*t197*t147)+2.D0*t73*
     #(0.142313679111053D-1*t7*t76-0.1803273593935696D-2*t193*t80
     #+0.158841216D-4*t197*t186)
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = rho**(1.D0/3.D0)
      t3 = t2*rho
      t4 = rho**2
      t5 = t2**2
      t7 = 1/t5/t4
      t8 = sigma*t7
      t10 = 1.D0+0.6349604207872798D-2*t8
      t11 = 1/t10
      t14 = sigma**2
      t15 = t4**2
      t18 = 1/t2/t15/rho
      t19 = t14*t18
      t20 = t10**2
      t21 = 1/t20
      t24 = 0.8094D0+0.322115421465387D-2*t8*t11+0.3016150199764335D-4
     #*t19*t21
      t27 = 1/rho
      t28 = t27**(1.D0/3.D0)
      t30 = 1.D0+0.1606016560364007D0*t28
      t31 = rho*t30
      t32 = t27**(1.D0/6.D0)
      t35 = dsqrt(t27)
      t37 = t28**2
      t39 = 0.1248219874679732D2*t32+0.4844076716063854D1*t28
     #+0.2326004811900819D1*t35+0.3819082618690966D0*t37
      t42 = 1.D0+0.321646831778707D2/t39
      t43 = dlog(t42)
      t45 = 1.D0+0.3174802103936399D0*t8
      t46 = 1/t45
      t49 = t45**2
      t50 = 1/t49
      t53 = 0.1737D0+0.745665770151542D0*t8*t46-0.2506537333502856D0
     #*t19*t50
      t54 = t43*t53
      t58 = 1.D0+0.1325688999052018D0*t28
      t64 = 0.598255043577108D1*t32+0.2225569421150687D1*t28
     #+0.8004286349993634D0*t35+0.1897004325747559D0*t37
      t67 = 1.D0+0.160818243221511D2/t64
      t68 = dlog(t67)
      t73 = -0.62182D-1*rho*t58*t68+0.3109D-1*t31*t43
      t75 = 1.D0+0.9524406311809197D-2*t8
      t76 = 1/t75
      t79 = t75**2
      t80 = 1/t79
      t83 = 0.9454D0+0.7115683955552651D-2*t8*t76
     #-0.4169320658943715D-3*t19*t80
      zk(i) = -0.7385587663820224D0*t3*t24-0.3109D-1*t31*t54+t73*t83
      t87 = t4*rho
      t89 = 1/t5/t87
      t90 = sigma*t89
      t95 = 1/t2/t15/t4
      t96 = t14*t95
      t99 = t14*sigma
      t100 = t15**2
      t102 = 1/t100/rho
      t103 = t99*t102
      t105 = 1/t20/t10
      t108 = -0.1717948914482064D-1*t90*t11-0.2126397314118042D-3*t96
     #*t21+0.2042811733333333D-5*t103*t105
      t111 = t30*t43
      t114 = 1/t37
      t115 = t27*t114
      t118 = t39**2
      t119 = 1/t118
      t120 = t31*t119
      t121 = t32**2
      t122 = t121**2
      t123 = t122*t32
      t124 = 1/t123
      t125 = 1/t4
      t126 = t124*t125
      t128 = t114*t125
      t130 = 1/t35
      t131 = t130*t125
      t133 = 1/t28
      t134 = t133*t125
      t136 = -0.4160732915599108D1*t126-0.3229384477375903D1*t128
     #-0.2326004811900819D1*t131-0.5092110158254621D0*t134
      t137 = 1/t42
      t138 = t136*t137
      t139 = t138*t53
      t147 = 1/t49/t45
      t150 = -0.3976884107474891D1*t90*t46+0.3936221825555298D1*t96
     #*t50-0.8488277333333333D0*t103*t147
      t151 = t43*t150
      t156 = t128*t68
      t158 = t64**2
      t159 = 1/t158
      t160 = t58*t159
      t165 = -0.99709173929518D0*t126-0.7418564737168958D0*t128
     #-0.4002143174996817D0*t131-0.1264669550498372D0*t134
      t166 = 1/t67
      t168 = t160*t165*t166
      t176 = t119*t136*t137
      t179 = -0.62182D-1*t58*t68+rho*(0.2747799777968419D-2*t156+0.1D1
     #*t168)+0.3109D-1*t111-0.1664368495390566D-2*t115*t43-0.5D0*t31
     #*t176
      t186 = 1/t79/t75
      t189 = -0.1897515721480707D-1*t90*t76+0.2404364791914262D-2*t96
     #*t80-0.211788288D-4*t103*t186
      vrhoa(i) = -0.9847450218426965D0*t2*t24-0.3692793831910112D0*t3
     #*t108-0.3109D-1*t111*t53+0.1664368495390566D-2*t115*t54+0.5D0
     #*t120*t139-0.15545D-1*t31*t151+t179*t83+t73*t189
      t193 = sigma*t18
      t196 = 1/t100
      t197 = t14*t196
      t200 = 0.1288461685861548D-1*t7*t11+0.1594797985588531D-3*t193
     #*t21-0.15321088D-5*t197*t105
      t209 = 0.2982663080606168D1*t7*t46-0.2952166369166474D1*t193*t50
     #+0.6366208D0*t197*t147
      t210 = t43*t209
      t219 = 0.142313679111053D-1*t7*t76-0.1803273593935696D-2*t193
     #*t80+0.158841216D-4*t197*t186
      vsigmaaa(i) = -0.7385587663820224D0*t3*t200-0.3109D-1*t31*t210
     #+2.D0*t73*t219
      t224 = sigma/t5/t15
      t230 = t14/t2/t15/t87
      t233 = t100*t4
      t235 = t99/t233
      t238 = t14**2
      t242 = t238/t5/t100/t15
      t243 = t20**2
      t244 = 1/t243
      t252 = 1/t87
      t254 = 1/t37/t27
      t255 = t252*t254
      t262 = 1/t118/t39
      t264 = t136**2
      t271 = t30*t119
      t280 = t49**2
      t281 = 1/t280
      t290 = 1/t15
      t291 = 1/t123/t27*t290
      t293 = t124*t252
      t295 = t254*t290
      t297 = t114*t252
      t301 = 1/t35/t27*t290
      t303 = t130*t252
      t307 = 1/t28/t27*t290
      t309 = t133*t252
      t311 = -0.6934554859331846D1*t291+0.1664293166239643D2*t293
     #-0.4305845969834537D1*t295+0.1291753790950361D2*t297
     #-0.2326004811900819D1*t301+0.9304019247603276D1*t303
     #-0.3394740105503081D0*t307+0.2036844063301848D1*t309
      t321 = 0.5495599555936838D-2*t156
      t322 = 0.2D1*t168
      t333 = dlog(1.D0+0.2960857464321668D2/(0.8157414703487641D1*t32
     #+0.2247591863577616D1*t28+0.4300972471276643D0*t35
     #+0.1911512595127338D0*t37))
      t335 = (1.D0+0.6901399211255825D-1*t28)*t333*t125
      t340 = 0.8837926660346786D-1*t128*t159*t165*t166
      t342 = 0.5495599555936838D-2*t297*t68
      t344 = 0.1831866518645613D-2*t295*t68
      t345 = t158**2
      t348 = t165**2
      t349 = t67**2
      t353 = 0.160818243221511D2*t58/t345*t348/t349
      t365 = 0.1D1*t160*(-0.8309097827459833D0*t291
     #+0.199418347859036D1*t293-0.4945709824779306D0*t295
     #+0.1483712947433792D1*t297-0.2001071587498409D0*t301
     #+0.8004286349993634D0*t303-0.4215565168327908D-1*t307
     #+0.2529339100996745D0*t309)*t166
      t371 = 0.2D1*t58/t158/t64*t348*t166
      t380 = t118**2
      t381 = 1/t380
      t383 = t42**2
      t384 = 1/t383
      t410 = t230*t80
      t412 = t235*t186
      t414 = t79**2
      t415 = 1/t414
      t417 = 0.1613726165595571D-5*t242*t415
      s1 = -0.3692793831910112D0*t3*(0.1259829203953514D0*t224*t11
     #+0.2111660829546542D-2*t230*t21-0.5117251128888889D-4*t235*t105
     #+0.2075367356458441D-6*t242*t244)-0.6218D-1*t111*t150
     #+0.2219157993854088D-2*t255*t54-0.1070677706909338D0*t115*t119
     #*t139-0.1D1*t31*t262*t264*t137*t53+0.3328736990781132D-2*t115
     #*t151+0.2D1*t271*t139-0.15545D-1*t31*t43*(0.291638167881492D2
     #*t224*t46-0.5659258047384578D2*t230*t50+0.2860873955555556D2
     #*t235*t147-0.431177611786597D1*t242*t281)+0.5D0*t120*t311*t137*t53
      v2rhoa2(i) = s1-0.6564966812284644D0/t5*t24-0.1969490043685393D1
     #*t2*t108+(t321+t322+rho*(0.3377399869956914D-1*t335-t340-t342
     #+t344+t353+t365-t371)-0.2219157993854088D-2*t255*t43
     #+0.1070677706909338D0*t115*t176-0.2D1*t271*t138
     #-0.1608234158893535D2*t31*t381*t264*t384-0.5D0*t31*t119*t311
     #*t137+0.1D1*t31*t262*t264*t137)*t83+0.1608234158893535D2*t31
     #*t381*t264*t384*t53+0.1D1*t120*t138*t150+4.D0*t179*t189+t73*
     #(0.1391511529085852D0*t224*t76-0.2452558687152736D-1*t410
     #+0.3903992832D-3*t412-t417)+(t321+t322+rho*(
     #-0.3377399869956914D-1*t335-t340-t342+t344+t353+t365-t371))*t83
     #+t73*(-0.6893578397489445D-2*t410+0.235087872D-3*t412-t417)
      t433 = sigma*t95
      t436 = t14*t102
      t442 = t99/t5/t100/t87
      t471 = t433*t80
      t473 = t436*t186
      t476 = 0.1210294624196678D-5*t442*t415
      v2rhoasigmaaa(i) = -0.9847450218426965D0*t2*t200
     #-0.3692793831910112D0*t3*(-0.6871795657928257D-1*t89*t11
     #-0.1264786025042201D-2*t433*t21+0.3531516586666667D-4*t436*t105
     #-0.1556525517343831D-6*t442*t244)-0.3109D-1*t111*t209
     #+0.1664368495390566D-2*t115*t210+0.5D0*t120*t138*t209-0.15545D-1
     #*t31*t43*(-0.1590753642989956D2*t89*t46+0.3654010261705139D2
     #*t433*t50-0.2018331306666667D2*t436*t147+0.3233832088399478D1
     #*t442*t281)+2.D0*t179*t219+t73*(-0.7590062885922828D-1*t89*t76
     #+0.1478764296577413D-1*t471-0.2610312192D-3*t473+t476)+t73*
     #(0.5170183798117084D-2*t471-0.176315904D-3*t473+t476)
      t485 = sigma*t196
      t490 = t14/t5/t233
      v2sigmaaa2(i) = -0.7385587663820224D0*t3*(0.3106703245462379D-3
     #*t18*t21-0.203579392D-4*t485*t105+0.1167394138007873D-6*t490
     #*t244)-0.3109D-1*t31*t43*(-0.1559641148612265D2*t18*t50
     #+0.125910016D2*t485*t147-0.2425374066299608D1*t490*t281)+4.D0
     #*t73*(-0.3877637848587813D-2*t18*t80+0.132236928D-3*t485*t186
     #-0.9077209681475088D-6*t490*t415)
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      v2rhoa2(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      v2rhoasigmaaa(i) = 0.0d0
      v2sigmaaa2(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end


      subroutine uks_xc_hcth407
     & (ideriv,npt,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,
     &  zk,vrhoa,vrhob,vsigmaaa,vsigmabb,vsigmaab,
     &  v2rhoa2,v2rhob2,v2rhoab,
     &  v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb,
     &  v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,
     &  v2sigmaaa2,v2sigmaaaab,v2sigmaaabb,
     &  v2sigmaab2,v2sigmaabbb,v2sigmabb2)
c
c     A.D. Boese, N.C. Handy
c     A new parametriztion of exchange-correlation generalized gradient
c     approximation functionals
c     J. Chem. Phys. 114 (2001) 5497-5503.
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit real*8 (a-h,o-z)
      integer ideriv,npt
      real*8 rhoa1(npt),rhob1(npt)
      real*8 sigmaaa1(npt),sigmabb1(npt),sigmaab1(npt)
      real*8 zk(npt),vrhoa(npt),vrhob(npt)
      real*8 vsigmaaa(npt),vsigmabb(npt),vsigmaab(npt)
      real*8 v2rhoa2(npt),v2rhob2(npt),v2rhoab(npt)
      real*8 v2rhoasigmaaa(npt),v2rhoasigmaab(npt)
      real*8 v2rhoasigmabb(npt),v2rhobsigmabb(npt)
      real*8 v2rhobsigmaab(npt),v2rhobsigmaaa(npt)
      real*8 v2sigmaaa2(npt),v2sigmaaaab(npt),v2sigmaaabb(npt)
      real*8 v2sigmaab2(npt),v2sigmaabbb(npt),v2sigmabb2(npt)
      parameter(tol=1.0d-20)
      
      if (ideriv.eq.0) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmabb
      t2 = rhob**(1.D0/3.D0)
      t4 = rhob**2
      t5 = t2**2
      t8 = sigmabb/t5/t4
      t10 = 1.D0+0.4D-2*t8
      t14 = sigmabb**2
      t15 = t4**2
      t19 = t14/t2/t15/rhob
      t20 = t10**2
      t25 = t15**2
      t27 = t14*sigmabb/t25
      t32 = t14**2
      t36 = t32/t5/t25/t4
      t37 = t20**2
      t44 = 1/rhob
      t45 = t44**(1.D0/3.D0)
      t49 = t44**(1.D0/6.D0)
      t52 = dsqrt(t44)
      t54 = t45**2
      t60 = dlog(1.D0+0.321646831778707D2/(0.1112037486309468D2*t49
     #+0.3844746237447211D1*t45+0.1644733775567609D1*t52
     #+0.2405871291288192D0*t54))
      t62 = 1.D0+0.2D0*t8
      t66 = t62**2
      t74 = t66**2
      zk(i) = -0.9305257363491D0*t2*rhob*(0.108184D1-0.2073356D-2*t8
     #/t10+0.5480992D-4*t19/t20-0.16825664D-6*t27/t20/t10+0.5858688D-9
     #*t36/t37)-0.3109D-1*rhob*(1.D0+0.1274696188700087D0*t45)*t60*
     #(0.118777D1-0.480584D0*t8/t62+0.2246964D0*t19/t66-0.7343384D-1
     #*t27/t66/t62+0.9996768D-2*t36/t74)
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = rhoa**(1.D0/3.D0)
      t4 = rhoa**2
      t5 = t2**2
      t8 = sigmaaa/t5/t4
      t10 = 1.D0+0.4D-2*t8
      t14 = sigmaaa**2
      t15 = t4**2
      t19 = t14/t2/t15/rhoa
      t20 = t10**2
      t25 = t15**2
      t27 = t14*sigmaaa/t25
      t32 = t14**2
      t36 = t32/t5/t25/t4
      t37 = t20**2
      t44 = 1/rhoa
      t45 = t44**(1.D0/3.D0)
      t49 = t44**(1.D0/6.D0)
      t52 = dsqrt(t44)
      t54 = t45**2
      t60 = dlog(1.D0+0.321646831778707D2/(0.1112037486309468D2*t49
     #+0.3844746237447211D1*t45+0.1644733775567609D1*t52
     #+0.2405871291288192D0*t54))
      t62 = 1.D0+0.2D0*t8
      t66 = t62**2
      t74 = t66**2
      zk(i) = -0.9305257363491D0*t2*rhoa*(0.108184D1-0.2073356D-2*t8
     #/t10+0.5480992D-4*t19/t20-0.16825664D-6*t27/t20/t10+0.5858688D-9
     #*t36/t37)-0.3109D-1*rhoa*(1.D0+0.1274696188700087D0*t45)*t60*
     #(0.118777D1-0.480584D0*t8/t62+0.2246964D0*t19/t66-0.7343384D-1
     #*t27/t66/t62+0.9996768D-2*t36/t74)
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rhoa**(1.D0/3.D0)
      t6 = rhoa**2
      t7 = t4**2
      t10 = sigmaaa/t7/t6
      t12 = 1.D0+0.4D-2*t10
      t16 = sigmaaa**2
      t17 = t6**2
      t21 = t16/t4/t17/rhoa
      t22 = t12**2
      t27 = t17**2
      t29 = t16*sigmaaa/t27
      t34 = t16**2
      t38 = t34/t7/t27/t6
      t39 = t22**2
      t46 = 1/rhoa
      t47 = t46**(1.D0/3.D0)
      t50 = rhoa*(1.D0+0.1274696188700087D0*t47)
      t51 = t46**(1.D0/6.D0)
      t54 = dsqrt(t46)
      t56 = t47**2
      t62 = dlog(1.D0+0.321646831778707D2/(0.1112037486309468D2*t51
     #+0.3844746237447211D1*t47+0.1644733775567609D1*t54
     #+0.2405871291288192D0*t56))
      t64 = 1.D0+0.2D0*t10
      t68 = t64**2
      t76 = t68**2
      t84 = rhob**(1.D0/3.D0)
      t86 = rhob**2
      t87 = t84**2
      t90 = sigmabb/t87/t86
      t92 = 1.D0+0.4D-2*t90
      t96 = sigmabb**2
      t97 = t86**2
      t101 = t96/t84/t97/rhob
      t102 = t92**2
      t107 = t97**2
      t109 = t96*sigmabb/t107
      t114 = t96**2
      t118 = t114/t87/t107/t86
      t119 = t102**2
      t126 = 1/rhob
      t127 = t126**(1.D0/3.D0)
      t130 = rhob*(1.D0+0.1274696188700087D0*t127)
      t131 = t126**(1.D0/6.D0)
      t134 = dsqrt(t126)
      t136 = t127**2
      t142 = dlog(1.D0+0.321646831778707D2/(0.1112037486309468D2*t131
     #+0.3844746237447211D1*t127+0.1644733775567609D1*t134
     #+0.2405871291288192D0*t136))
      t144 = 1.D0+0.2D0*t90
      t148 = t144**2
      t156 = t148**2
      t164 = rhoa+rhob
      t165 = 1/t164
      t166 = t165**(1.D0/3.D0)
      t169 = t165**(1.D0/6.D0)
      t172 = dsqrt(t165)
      t174 = t166**2
      t180 = dlog(1.D0+0.160818243221511D2/(0.598255043577108D1*t169
     #+0.2225569421150687D1*t166+0.8004286349993634D0*t172
     #+0.1897004325747559D0*t174))
      t182 = 0.62182D-1*(1.D0+0.1325688999052018D0*t166)*t180
      t193 = dlog(1.D0+0.2960857464321668D2/(0.8157414703487641D1*t169
     #+0.2247591863577616D1*t166+0.4300972471276643D0*t172
     #+0.1911512595127338D0*t174))
      t196 = rhoa-1.D0*rhob
      t197 = t196*t165
      t198 = 1.D0+t197
      t199 = t198**(1.D0/3.D0)
      t202 = 1.D0-1.D0*t197
      t203 = t202**(1.D0/3.D0)
      t205 = t199*t198+t203*t202-2.D0
      t206 = t196**2
      t207 = t206**2
      t208 = t164**2
      t209 = t208**2
      t211 = t207/t209
      t227 = dlog(1.D0+0.321646831778707D2/(0.1112037486309468D2*t169
     #+0.3844746237447211D1*t166+0.1644733775567609D1*t172
     #+0.2405871291288192D0*t174))
      t243 = 0.5D0*t10+0.5D0*t90
      t246 = 1.D0+0.3D-2*t10+0.3D-2*t90
      t250 = t243**2
      t251 = t246**2
      t260 = t250**2
      t261 = t251**2
      s1 = -0.9305257363491D0*t4*rhoa*(0.108184D1-0.2073356D-2*t10/t12
     #+0.5480992D-4*t21/t22-0.16825664D-6*t29/t22/t12+0.5858688D-9*t38
     #/t39)-0.3109D-1*t50*t62*(0.118777D1-0.480584D0*t10/t64
     #+0.2246964D0*t21/t68-0.7343384D-1*t29/t68/t64+0.9996768D-2*t38
     #/t76)
      s2 = s1-0.9305257363491D0*t84*rhob*(0.108184D1-0.2073356D-2*t90
     #/t92+0.5480992D-4*t101/t102-0.16825664D-6*t109/t102/t92
     #+0.5858688D-9*t118/t119)
      zk(i) = s2-0.3109D-1*t130*t142*(0.118777D1-0.480584D0*t90/t144
     #+0.2246964D0*t101/t148-0.7343384D-1*t109/t148/t144+0.9996768D-2
     #*t118/t156)+(t164*(-t182+0.3799574853701528D-1*(1.D0
     #+0.6901399211255825D-1*t166)*t193*t205*(1.D0-1.D0*t211)
     #+0.1923661050931536D1*(-0.3109D-1*(1.D0+0.1274696188700087D0
     #*t166)*t227+t182)*t205*t211)+0.3109D-1*t50*t62+0.3109D-1*t130
     #*t142)*(0.589076D0+0.2654244D-1*t243/t246-0.6919848D-3*t250/t251
     #+0.91955736D-5*t250*t243/t251/t246-0.544387392D-7*t260/t261)
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmabb
      t2 = rhob**(1.D0/3.D0)
      t3 = t2*rhob
      t4 = rhob**2
      t5 = t2**2
      t7 = 1/t5/t4
      t8 = sigmabb*t7
      t10 = 1.D0+0.4D-2*t8
      t11 = 1/t10
      t14 = sigmabb**2
      t15 = t4**2
      t16 = t15*rhob
      t18 = 1/t2/t16
      t19 = t14*t18
      t20 = t10**2
      t21 = 1/t20
      t24 = t14*sigmabb
      t25 = t15**2
      t26 = 1/t25
      t27 = t24*t26
      t29 = 1/t20/t10
      t32 = t14**2
      t35 = 1/t5/t25/t4
      t36 = t32*t35
      t37 = t20**2
      t38 = 1/t37
      t41 = 0.108184D1-0.2073356D-2*t8*t11+0.5480992D-4*t19*t21
     #-0.16825664D-6*t27*t29+0.5858688D-9*t36*t38
      t44 = 1/rhob
      t45 = t44**(1.D0/3.D0)
      t47 = 1.D0+0.1274696188700087D0*t45
      t48 = rhob*t47
      t49 = t44**(1.D0/6.D0)
      t52 = dsqrt(t44)
      t54 = t45**2
      t56 = 0.1112037486309468D2*t49+0.3844746237447211D1*t45
     #+0.1644733775567609D1*t52+0.2405871291288192D0*t54
      t59 = 1.D0+0.321646831778707D2/t56
      t60 = dlog(t59)
      t62 = 1.D0+0.2D0*t8
      t63 = 1/t62
      t66 = t62**2
      t67 = 1/t66
      t71 = 1/t66/t62
      t74 = t66**2
      t75 = 1/t74
      t78 = 0.118777D1-0.480584D0*t8*t63+0.2246964D0*t19*t67
     #-0.7343384D-1*t27*t71+0.9996768D-2*t36*t75
      t79 = t60*t78
      zk(i) = -0.9305257363491D0*t3*t41-0.3109D-1*t48*t79
      vrhoa(i) = 0.D0
      t84 = t4*rhob
      t87 = sigmabb/t5/t84
      t90 = t15*t4
      t93 = t14/t2/t90
      t98 = t24/t25/rhob
      t104 = t32/t5/t25/t84
      t111 = t32*sigmabb/t2/t25/t90
      t113 = 1/t37/t10
      t122 = 1/t54
      t126 = t56**2
      t129 = t49**2
      t130 = t129**2
      t133 = 1/t4
      t159 = 1/t74/t62
      vrhob(i) = -0.12407009817988D1*t2*t41-0.9305257363491D0*t3*
     #(0.5528949333333333D-2*t87*t11-0.3144353706666667D-3*t93*t21
     #+0.2515331413333333D-5*t98*t29-0.1163347968D-7*t104*t38
     #+0.249970688D-10*t111*t113)-0.3109D-1*t47*t60*t78
     #+0.1321010150222857D-2*t44*t122*t79+0.1D1*t48/t126*(
     #-0.1853395810515781D1/t130/t49*t133-0.128158207914907D1*t122
     #*t133-0.8223668877838045D0/t52*t133-0.1603914194192128D0/t45
     #*t133)/t59*t78-0.3109D-1*t48*t60*(0.1281557333333333D1*t87*t63
     #-0.1454692266666667D1*t93*t67+0.82714688D0*t98*t71-0.224126336D0
     #*t104*t75+0.213264384D-1*t111*t159)
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      t168 = sigmabb*t18
      t171 = t14*t26
      t174 = t24*t35
      t180 = t32/t2/t25/t16
      vsigmabb(i) = -0.9305257363491D0*t3*(-0.2073356D-2*t7*t11
     #+0.117913264D-3*t168*t21-0.94324928D-6*t171*t29+0.436255488D-8
     #*t174*t38-0.93739008D-11*t180*t113)-0.3109D-1*t48*t60*(
     #-0.480584D0*t7*t63+0.5455096D0*t168*t67-0.31018008D0*t171*t71
     #+0.84047376D-1*t174*t75-0.79974144D-2*t180*t159)
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = rhoa**(1.D0/3.D0)
      t3 = t2*rhoa
      t4 = rhoa**2
      t5 = t2**2
      t7 = 1/t5/t4
      t8 = sigmaaa*t7
      t10 = 1.D0+0.4D-2*t8
      t11 = 1/t10
      t14 = sigmaaa**2
      t15 = t4**2
      t16 = t15*rhoa
      t18 = 1/t2/t16
      t19 = t14*t18
      t20 = t10**2
      t21 = 1/t20
      t24 = t14*sigmaaa
      t25 = t15**2
      t26 = 1/t25
      t27 = t24*t26
      t29 = 1/t20/t10
      t32 = t14**2
      t35 = 1/t5/t25/t4
      t36 = t32*t35
      t37 = t20**2
      t38 = 1/t37
      t41 = 0.108184D1-0.2073356D-2*t8*t11+0.5480992D-4*t19*t21
     #-0.16825664D-6*t27*t29+0.5858688D-9*t36*t38
      t44 = 1/rhoa
      t45 = t44**(1.D0/3.D0)
      t47 = 1.D0+0.1274696188700087D0*t45
      t48 = rhoa*t47
      t49 = t44**(1.D0/6.D0)
      t52 = dsqrt(t44)
      t54 = t45**2
      t56 = 0.1112037486309468D2*t49+0.3844746237447211D1*t45
     #+0.1644733775567609D1*t52+0.2405871291288192D0*t54
      t59 = 1.D0+0.321646831778707D2/t56
      t60 = dlog(t59)
      t62 = 1.D0+0.2D0*t8
      t63 = 1/t62
      t66 = t62**2
      t67 = 1/t66
      t71 = 1/t66/t62
      t74 = t66**2
      t75 = 1/t74
      t78 = 0.118777D1-0.480584D0*t8*t63+0.2246964D0*t19*t67
     #-0.7343384D-1*t27*t71+0.9996768D-2*t36*t75
      t79 = t60*t78
      zk(i) = -0.9305257363491D0*t3*t41-0.3109D-1*t48*t79
      t84 = t4*rhoa
      t87 = sigmaaa/t5/t84
      t90 = t15*t4
      t93 = t14/t2/t90
      t98 = t24/t25/rhoa
      t104 = t32/t5/t25/t84
      t111 = t32*sigmaaa/t2/t25/t90
      t113 = 1/t37/t10
      t122 = 1/t54
      t126 = t56**2
      t129 = t49**2
      t130 = t129**2
      t133 = 1/t4
      t159 = 1/t74/t62
      vrhoa(i) = -0.12407009817988D1*t2*t41-0.9305257363491D0*t3*
     #(0.5528949333333333D-2*t87*t11-0.3144353706666667D-3*t93*t21
     #+0.2515331413333333D-5*t98*t29-0.1163347968D-7*t104*t38
     #+0.249970688D-10*t111*t113)-0.3109D-1*t47*t60*t78
     #+0.1321010150222857D-2*t44*t122*t79+0.1D1*t48/t126*(
     #-0.1853395810515781D1/t130/t49*t133-0.128158207914907D1*t122
     #*t133-0.8223668877838045D0/t52*t133-0.1603914194192128D0/t45
     #*t133)/t59*t78-0.3109D-1*t48*t60*(0.1281557333333333D1*t87*t63
     #-0.1454692266666667D1*t93*t67+0.82714688D0*t98*t71-0.224126336D0
     #*t104*t75+0.213264384D-1*t111*t159)
      vrhob(i) = 0.D0
      t168 = sigmaaa*t18
      t171 = t14*t26
      t174 = t24*t35
      t180 = t32/t2/t25/t16
      vsigmaaa(i) = -0.9305257363491D0*t3*(-0.2073356D-2*t7*t11
     #+0.117913264D-3*t168*t21-0.94324928D-6*t171*t29+0.436255488D-8
     #*t174*t38-0.93739008D-11*t180*t113)-0.3109D-1*t48*t60*(
     #-0.480584D0*t7*t63+0.5455096D0*t168*t67-0.31018008D0*t171*t71
     #+0.84047376D-1*t174*t75-0.79974144D-2*t180*t159)
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rhoa**(1.D0/3.D0)
      t5 = t4*rhoa
      t6 = rhoa**2
      t7 = t4**2
      t9 = 1/t7/t6
      t10 = sigmaaa*t9
      t12 = 1.D0+0.4D-2*t10
      t13 = 1/t12
      t16 = sigmaaa**2
      t17 = t6**2
      t18 = t17*rhoa
      t20 = 1/t4/t18
      t21 = t16*t20
      t22 = t12**2
      t23 = 1/t22
      t26 = t16*sigmaaa
      t27 = t17**2
      t28 = 1/t27
      t29 = t26*t28
      t31 = 1/t22/t12
      t34 = t16**2
      t37 = 1/t7/t27/t6
      t38 = t34*t37
      t39 = t22**2
      t40 = 1/t39
      t43 = 0.108184D1-0.2073356D-2*t10*t13+0.5480992D-4*t21*t23
     #-0.16825664D-6*t29*t31+0.5858688D-9*t38*t40
      t46 = 1/rhoa
      t47 = t46**(1.D0/3.D0)
      t49 = 1.D0+0.1274696188700087D0*t47
      t50 = rhoa*t49
      t51 = t46**(1.D0/6.D0)
      t54 = dsqrt(t46)
      t56 = t47**2
      t58 = 0.1112037486309468D2*t51+0.3844746237447211D1*t47
     #+0.1644733775567609D1*t54+0.2405871291288192D0*t56
      t61 = 1.D0+0.321646831778707D2/t58
      t62 = dlog(t61)
      t64 = 1.D0+0.2D0*t10
      t65 = 1/t64
      t68 = t64**2
      t69 = 1/t68
      t73 = 1/t68/t64
      t76 = t68**2
      t77 = 1/t76
      t80 = 0.118777D1-0.480584D0*t10*t65+0.2246964D0*t21*t69
     #-0.7343384D-1*t29*t73+0.9996768D-2*t38*t77
      t81 = t62*t80
      t84 = rhob**(1.D0/3.D0)
      t85 = t84*rhob
      t86 = rhob**2
      t87 = t84**2
      t89 = 1/t87/t86
      t90 = sigmabb*t89
      t92 = 1.D0+0.4D-2*t90
      t93 = 1/t92
      t96 = sigmabb**2
      t97 = t86**2
      t98 = t97*rhob
      t100 = 1/t84/t98
      t101 = t96*t100
      t102 = t92**2
      t103 = 1/t102
      t106 = t96*sigmabb
      t107 = t97**2
      t108 = 1/t107
      t109 = t106*t108
      t111 = 1/t102/t92
      t114 = t96**2
      t117 = 1/t87/t107/t86
      t118 = t114*t117
      t119 = t102**2
      t120 = 1/t119
      t123 = 0.108184D1-0.2073356D-2*t90*t93+0.5480992D-4*t101*t103
     #-0.16825664D-6*t109*t111+0.5858688D-9*t118*t120
      t126 = 1/rhob
      t127 = t126**(1.D0/3.D0)
      t129 = 1.D0+0.1274696188700087D0*t127
      t130 = rhob*t129
      t131 = t126**(1.D0/6.D0)
      t134 = dsqrt(t126)
      t136 = t127**2
      t138 = 0.1112037486309468D2*t131+0.3844746237447211D1*t127
     #+0.1644733775567609D1*t134+0.2405871291288192D0*t136
      t141 = 1.D0+0.321646831778707D2/t138
      t142 = dlog(t141)
      t144 = 1.D0+0.2D0*t90
      t145 = 1/t144
      t148 = t144**2
      t149 = 1/t148
      t153 = 1/t148/t144
      t156 = t148**2
      t157 = 1/t156
      t160 = 0.118777D1-0.480584D0*t90*t145+0.2246964D0*t101*t149
     #-0.7343384D-1*t109*t153+0.9996768D-2*t118*t157
      t161 = t142*t160
      t164 = rhoa+rhob
      t165 = 1/t164
      t166 = t165**(1.D0/3.D0)
      t168 = 1.D0+0.1325688999052018D0*t166
      t169 = t165**(1.D0/6.D0)
      t172 = dsqrt(t165)
      t174 = t166**2
      t176 = 0.598255043577108D1*t169+0.2225569421150687D1*t166
     #+0.8004286349993634D0*t172+0.1897004325747559D0*t174
      t179 = 1.D0+0.160818243221511D2/t176
      t180 = dlog(t179)
      t182 = 0.62182D-1*t168*t180
      t184 = 1.D0+0.6901399211255825D-1*t166
      t189 = 0.8157414703487641D1*t169+0.2247591863577616D1*t166
     #+0.4300972471276643D0*t172+0.1911512595127338D0*t174
      t192 = 1.D0+0.2960857464321668D2/t189
      t193 = dlog(t192)
      t194 = t184*t193
      t196 = rhoa-1.D0*rhob
      t197 = t196*t165
      t198 = 1.D0+t197
      t199 = t198**(1.D0/3.D0)
      t202 = 1.D0-1.D0*t197
      t203 = t202**(1.D0/3.D0)
      t205 = t199*t198+t203*t202-2.D0
      t206 = t196**2
      t207 = t206**2
      t208 = t164**2
      t209 = t208**2
      t210 = 1/t209
      t211 = t207*t210
      t213 = 1.D0-1.D0*t211
      t216 = 0.3799574853701528D-1*t194*t205*t213
      t218 = 1.D0+0.1274696188700087D0*t166
      t223 = 0.1112037486309468D2*t169+0.3844746237447211D1*t166
     #+0.1644733775567609D1*t172+0.2405871291288192D0*t174
      t226 = 1.D0+0.321646831778707D2/t223
      t227 = dlog(t226)
      t230 = -0.3109D-1*t218*t227+t182
      t231 = t230*t205
      t233 = 0.1923661050931536D1*t231*t211
      t240 = t164*(-t182+t216+t233)+0.3109D-1*t50*t62+0.3109D-1*t130
     #*t142
      t243 = 0.5D0*t10+0.5D0*t90
      t246 = 1.D0+0.3D-2*t10+0.3D-2*t90
      t247 = 1/t246
      t250 = t243**2
      t251 = t246**2
      t252 = 1/t251
      t255 = t250*t243
      t257 = 1/t251/t246
      t260 = t250**2
      t261 = t251**2
      t262 = 1/t261
      t265 = 0.589076D0+0.2654244D-1*t243*t247-0.6919848D-3*t250*t252
     #+0.91955736D-5*t255*t257-0.544387392D-7*t260*t262
      zk(i) = -0.9305257363491D0*t5*t43-0.3109D-1*t50*t81
     #-0.9305257363491D0*t85*t123-0.3109D-1*t130*t161+t240*t265
      t269 = t6*rhoa
      t272 = sigmaaa/t7/t269
      t275 = t17*t6
      t278 = t16/t4/t275
      t283 = t26/t27/rhoa
      t289 = t34/t7/t27/t269
      t296 = t34*sigmaaa/t4/t27/t275
      t298 = 1/t39/t12
      t304 = t49*t62
      t307 = 1/t56
      t308 = t46*t307
      t311 = t58**2
      t312 = 1/t311
      t314 = t51**2
      t315 = t314**2
      t318 = 1/t6
      t329 = -0.1853395810515781D1/t315/t51*t318-0.128158207914907D1
     #*t307*t318-0.8223668877838045D0/t54*t318-0.1603914194192128D0
     #/t47*t318
      t330 = 1/t61
      t344 = 1/t76/t64
      t352 = 1/t208
      t353 = 1/t174*t352
      t354 = t353*t180
      t355 = 0.2747799777968419D-2*t354
      t356 = t176**2
      t359 = t169**2
      t360 = t359**2
      t363 = 1/t360/t169*t352
      t367 = 1/t172*t352
      t370 = 1/t166*t352
      t375 = t168/t356*(-0.99709173929518D0*t363-0.7418564737168958D0
     #*t353-0.4002143174996817D0*t367-0.1264669550498372D0*t370)/t179
      t376 = 0.1D1*t375
      t380 = 0.8740794299481065D-3*t353*t193*t205*t213
      t381 = t189**2
      t394 = 0.1124999956683108D1*t184/t381*(-0.135956911724794D1*t363
     #-0.7491972878592054D0*t353-0.2150486235638321D0*t367
     #-0.1274341730084892D0*t370)/t192*t205*t213
      t395 = t196*t352
      t396 = 1.D0*t395
      t400 = 1.D0*t165
      t404 = 0.1333333333333333D1*t199*(t165-t396)
     #+0.1333333333333333D1*t203*(-t400+t395)
      t409 = t206*t196*t210
      t410 = 4.D0*t409
      t413 = t207/t209/t164
      t414 = 4.D0*t413
      t421 = t223**2
      t438 = 0.1923661050931536D1*(0.1321010150222857D-2*t353*t227
     #+0.1D1*t218/t421*(-0.1853395810515781D1*t363-0.128158207914907D1
     #*t353-0.8223668877838045D0*t367-0.1603914194192128D0*t370)/t226
     #-0.2747799777968419D-2*t354-0.1D1*t375)*t205*t211
      t442 = t231*t409
      t445 = 0.7694644203726145D1*t231*t413
      t459 = t243*t252
      t462 = t250*t257
      t465 = t255*t262
      t470 = t260/t261/t246
      s1 = -0.12407009817988D1*t4*t43-0.9305257363491D0*t5*
     #(0.5528949333333333D-2*t272*t13-0.3144353706666667D-3*t278*t23
     #+0.2515331413333333D-5*t283*t31-0.1163347968D-7*t289*t40
     #+0.249970688D-10*t296*t298)-0.3109D-1*t304*t80
     #+0.1321010150222857D-2*t308*t81
      vrhoa(i) = s1+0.1D1*t50*t312*t329*t330*t80-0.3109D-1*t50*t62*
     #(0.1281557333333333D1*t272*t65-0.1454692266666667D1*t278*t69
     #+0.82714688D0*t283*t73-0.224126336D0*t289*t77+0.213264384D-1
     #*t296*t344)+(-t182+t216+t233+t164*(t355+t376-t380-t394
     #+0.3799574853701528D-1*t194*t404*t213+0.3799574853701528D-1*t194
     #*t205*(-t410+t414)+t438+0.1923661050931536D1*t230*t404*t211
     #+0.7694644203726145D1*t442-t445)+0.3109D-1*t304
     #-0.1321010150222857D-2*t308*t62-0.1D1*t50*t312*t329*t330)*t265
     #+t240*(-0.3538992D-1*t272*t247+0.205763232D-2*t459*t272
     #-0.478540512D-4*t462*t272+0.5110337088D-6*t465*t272
     #-0.17420396544D-8*t470*t272)
      t477 = t86*rhob
      t480 = sigmabb/t87/t477
      t483 = t97*t86
      t486 = t96/t84/t483
      t491 = t106/t107/rhob
      t497 = t114/t87/t107/t477
      t504 = t114*sigmabb/t84/t107/t483
      t506 = 1/t119/t92
      t512 = t129*t142
      t515 = 1/t136
      t516 = t126*t515
      t519 = t138**2
      t520 = 1/t519
      t522 = t131**2
      t523 = t522**2
      t526 = 1/t86
      t537 = -0.1853395810515781D1/t523/t131*t526-0.128158207914907D1
     #*t515*t526-0.8223668877838045D0/t134*t526-0.1603914194192128D0
     #/t127*t526
      t538 = 1/t141
      t552 = 1/t156/t144
      t565 = 0.1333333333333333D1*t199*(-t400-t396)
     #+0.1333333333333333D1*t203*(t165+t395)
      s1 = -0.12407009817988D1*t84*t123-0.9305257363491D0*t85*
     #(0.5528949333333333D-2*t480*t93-0.3144353706666667D-3*t486*t103
     #+0.2515331413333333D-5*t491*t111-0.1163347968D-7*t497*t120
     #+0.249970688D-10*t504*t506)-0.3109D-1*t512*t160
     #+0.1321010150222857D-2*t516*t161
      vrhob(i) = s1+0.1D1*t130*t520*t537*t538*t160-0.3109D-1*t130*t142
     #*(0.1281557333333333D1*t480*t145-0.1454692266666667D1*t486*t149
     #+0.82714688D0*t491*t153-0.224126336D0*t497*t157+0.213264384D-1
     #*t504*t552)+(-t182+t216+t233+t164*(t355+t376-t380-t394
     #+0.3799574853701528D-1*t194*t565*t213+0.3799574853701528D-1*t194
     #*t205*(t410+t414)+t438+0.1923661050931536D1*t230*t565*t211
     #-0.7694644203726145D1*t442-t445)+0.3109D-1*t512
     #-0.1321010150222857D-2*t516*t142-0.1D1*t130*t520*t537*t538)*t265
     #+t240*(-0.3538992D-1*t480*t247+0.205763232D-2*t459*t480
     #-0.478540512D-4*t462*t480+0.5110337088D-6*t465*t480
     #-0.17420396544D-8*t470*t480)
      t602 = sigmaaa*t20
      t605 = t16*t28
      t608 = t26*t37
      t614 = t34/t4/t27/t18
      vsigmaaa(i) = -0.9305257363491D0*t5*(-0.2073356D-2*t9*t13
     #+0.117913264D-3*t602*t23-0.94324928D-6*t605*t31+0.436255488D-8
     #*t608*t40-0.93739008D-11*t614*t298)-0.3109D-1*t50*t62*(
     #-0.480584D0*t9*t65+0.5455096D0*t602*t69-0.31018008D0*t605*t73
     #+0.84047376D-1*t608*t77-0.79974144D-2*t614*t344)+t240*
     #(0.1327122D-1*t9*t247-0.77161212D-3*t459*t9+0.179452692D-4*t462
     #*t9-0.1916376408D-6*t465*t9+0.6532648704D-9*t470*t9)
      vsigmaab(i) = 0.D0
      t648 = sigmabb*t100
      t651 = t96*t108
      t654 = t106*t117
      t660 = t114/t84/t107/t98
      vsigmabb(i) = -0.9305257363491D0*t85*(-0.2073356D-2*t89*t93
     #+0.117913264D-3*t648*t103-0.94324928D-6*t651*t111+0.436255488D-8
     #*t654*t120-0.93739008D-11*t660*t506)-0.3109D-1*t130*t142*(
     #-0.480584D0*t89*t145+0.5455096D0*t648*t149-0.31018008D0*t651
     #*t153+0.84047376D-1*t654*t157-0.79974144D-2*t660*t552)+t240*
     #(0.1327122D-1*t89*t247-0.77161212D-3*t459*t89+0.179452692D-4
     #*t462*t89-0.1916376408D-6*t465*t89+0.6532648704D-9*t470*t89)
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vrhob(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      vsigmaab(i) = 0.0d0
      vsigmabb(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmabb
      t2 = rhob**(1.D0/3.D0)
      t3 = t2*rhob
      t4 = rhob**2
      t5 = t2**2
      t7 = 1/t5/t4
      t8 = sigmabb*t7
      t10 = 1.D0+0.4D-2*t8
      t11 = 1/t10
      t14 = sigmabb**2
      t15 = t4**2
      t16 = t15*rhob
      t18 = 1/t2/t16
      t19 = t14*t18
      t20 = t10**2
      t21 = 1/t20
      t24 = t14*sigmabb
      t25 = t15**2
      t26 = 1/t25
      t27 = t24*t26
      t29 = 1/t20/t10
      t32 = t14**2
      t33 = t25*t4
      t35 = 1/t5/t33
      t36 = t32*t35
      t37 = t20**2
      t38 = 1/t37
      t41 = 0.108184D1-0.2073356D-2*t8*t11+0.5480992D-4*t19*t21
     #-0.16825664D-6*t27*t29+0.5858688D-9*t36*t38
      t44 = 1/rhob
      t45 = t44**(1.D0/3.D0)
      t47 = 1.D0+0.1274696188700087D0*t45
      t48 = rhob*t47
      t49 = t44**(1.D0/6.D0)
      t52 = dsqrt(t44)
      t54 = t45**2
      t56 = 0.1112037486309468D2*t49+0.3844746237447211D1*t45
     #+0.1644733775567609D1*t52+0.2405871291288192D0*t54
      t59 = 1.D0+0.321646831778707D2/t56
      t60 = dlog(t59)
      t62 = 1.D0+0.2D0*t8
      t63 = 1/t62
      t66 = t62**2
      t67 = 1/t66
      t71 = 1/t66/t62
      t74 = t66**2
      t75 = 1/t74
      t78 = 0.118777D1-0.480584D0*t8*t63+0.2246964D0*t19*t67
     #-0.7343384D-1*t27*t71+0.9996768D-2*t36*t75
      t79 = t60*t78
      zk(i) = -0.9305257363491D0*t3*t41-0.3109D-1*t48*t79
      vrhoa(i) = 0.D0
      t84 = t4*rhob
      t87 = sigmabb/t5/t84
      t90 = t15*t4
      t93 = t14/t2/t90
      t98 = t24/t25/rhob
      t104 = t32/t5/t25/t84
      t107 = t32*sigmabb
      t111 = t107/t2/t25/t90
      t113 = 1/t37/t10
      t116 = 0.5528949333333333D-2*t87*t11-0.3144353706666667D-3*t93
     #*t21+0.2515331413333333D-5*t98*t29-0.1163347968D-7*t104*t38
     #+0.249970688D-10*t111*t113
      t119 = t47*t60
      t122 = 1/t54
      t123 = t44*t122
      t126 = t56**2
      t127 = 1/t126
      t128 = t48*t127
      t129 = t49**2
      t130 = t129**2
      t131 = t130*t49
      t132 = 1/t131
      t133 = 1/t4
      t138 = 1/t52
      t141 = 1/t45
      t144 = -0.1853395810515781D1*t132*t133-0.128158207914907D1*t122
     #*t133-0.8223668877838045D0*t138*t133-0.1603914194192128D0*t141
     #*t133
      t145 = 1/t59
      t146 = t144*t145
      t147 = t146*t78
      t159 = 1/t74/t62
      t162 = 0.1281557333333333D1*t87*t63-0.1454692266666667D1*t93*t67
     #+0.82714688D0*t98*t71-0.224126336D0*t104*t75+0.213264384D-1*t111
     #*t159
      t163 = t60*t162
      vrhob(i) = -0.12407009817988D1*t2*t41-0.9305257363491D0*t3*t116
     #-0.3109D-1*t119*t78+0.1321010150222857D-2*t123*t79+0.1D1*t128
     #*t147-0.3109D-1*t48*t163
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      t168 = sigmabb*t18
      t171 = t14*t26
      t174 = t24*t35
      t179 = 1/t2/t25/t16
      t180 = t32*t179
      vsigmabb(i) = -0.9305257363491D0*t3*(-0.2073356D-2*t7*t11
     #+0.117913264D-3*t168*t21-0.94324928D-6*t171*t29+0.436255488D-8
     #*t174*t38-0.93739008D-11*t180*t113)-0.3109D-1*t48*t60*(
     #-0.480584D0*t7*t63+0.5455096D0*t168*t67-0.31018008D0*t171*t71
     #+0.84047376D-1*t174*t75-0.79974144D-2*t180*t159)
      v2rhoa2(i) = 0.D0
      v2rhoab(i) = 0.D0
      t207 = sigmabb/t5/t15
      t210 = t15*t84
      t213 = t14/t2/t210
      t217 = t24/t33
      t223 = t32/t5/t25/t15
      t229 = t107/t2/t25/t210
      t233 = t25**2
      t236 = t32*t14/t233/t4
      t238 = 1/t37/t20
      t249 = 1/t84
      t251 = 1/t54/t44
      t263 = t144**2
      t270 = 1/t15
      t296 = t126**2
      t299 = t59**2
      t319 = 1/t74/t66
      s1 = -0.4135669939329333D0/t5*t41-0.24814019635976D1*t2*t116
     #-0.9305257363491D0*t3*(-0.2027281422222222D-1*t207*t11
     #+0.2050399473777778D-2*t213*t21-0.2934593729422222D-4*t217*t29
     #+0.2162145348266667D-6*t223*t38-0.8546531191466667D-9*t229*t113
     #+0.1333177002666667D-11*t236*t238)+0.2D1*t47*t127*t147-0.6218D-1
     #*t119*t162+0.8806734334819047D-3*t249*t251*t79
      v2rhob2(i) = s1-0.8497974591333914D-1*t123*t127*t147
     #+0.2642020300445714D-2*t123*t163-0.2D1*t48/t126/t56*t263*t145
     #*t78+0.1D1*t128*(-0.1544496508763151D1/t131/t44*t270
     #+0.3706791621031562D1*t132*t249-0.854388052766047D0*t251*t270
     #+0.2563164158298141D1*t122*t249-0.4111834438919023D0/t52/t44
     #*t270+0.1644733775567609D1*t138*t249-0.5346380647307093D-1/t45
     #/t44*t270+0.3207828388384256D0*t141*t249)*t145*t78
     #+0.321646831778707D2*t48/t296*t263/t299*t78+0.2D1*t128*t146*t162
     #-0.3109D-1*t48*t60*(-0.4699043555555556D1*t207*t63
     #+0.9896548266666667D1*t213*t67-0.8995993671111111D1*t217*t71
     #+0.3938242261333333D1*t223*t75-0.7838151338666667D0*t229*t159
     #+0.568705024D-1*t236*t319)
      v2sigmaaa2(i) = 0.D0
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      t328 = sigmabb*t26
      t331 = t14*t35
      t334 = t24*t179
      t338 = t32/t233
      v2sigmabb2(i) = -0.9305257363491D0*t3*(0.126206688D-3*t18*t21
     #-0.2829804672D-5*t328*t29+0.24406656D-7*t331*t38
     #-0.10729648128D-9*t334*t113+0.187478016D-12*t338*t238)-0.3109D-1
     #*t48*t60*(0.6416264D0*t18*t67-0.838564D0*t328*t71+0.438250176D0
     #*t331*t75-0.992275584D-1*t334*t159+0.79974144D-2*t338*t319)
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = rhoa**(1.D0/3.D0)
      t3 = t2*rhoa
      t4 = rhoa**2
      t5 = t2**2
      t7 = 1/t5/t4
      t8 = sigmaaa*t7
      t10 = 1.D0+0.4D-2*t8
      t11 = 1/t10
      t14 = sigmaaa**2
      t15 = t4**2
      t16 = t15*rhoa
      t18 = 1/t2/t16
      t19 = t14*t18
      t20 = t10**2
      t21 = 1/t20
      t24 = t14*sigmaaa
      t25 = t15**2
      t26 = 1/t25
      t27 = t24*t26
      t29 = 1/t20/t10
      t32 = t14**2
      t33 = t25*t4
      t35 = 1/t5/t33
      t36 = t32*t35
      t37 = t20**2
      t38 = 1/t37
      t41 = 0.108184D1-0.2073356D-2*t8*t11+0.5480992D-4*t19*t21
     #-0.16825664D-6*t27*t29+0.5858688D-9*t36*t38
      t44 = 1/rhoa
      t45 = t44**(1.D0/3.D0)
      t47 = 1.D0+0.1274696188700087D0*t45
      t48 = rhoa*t47
      t49 = t44**(1.D0/6.D0)
      t52 = dsqrt(t44)
      t54 = t45**2
      t56 = 0.1112037486309468D2*t49+0.3844746237447211D1*t45
     #+0.1644733775567609D1*t52+0.2405871291288192D0*t54
      t59 = 1.D0+0.321646831778707D2/t56
      t60 = dlog(t59)
      t62 = 1.D0+0.2D0*t8
      t63 = 1/t62
      t66 = t62**2
      t67 = 1/t66
      t71 = 1/t66/t62
      t74 = t66**2
      t75 = 1/t74
      t78 = 0.118777D1-0.480584D0*t8*t63+0.2246964D0*t19*t67
     #-0.7343384D-1*t27*t71+0.9996768D-2*t36*t75
      t79 = t60*t78
      zk(i) = -0.9305257363491D0*t3*t41-0.3109D-1*t48*t79
      t84 = t4*rhoa
      t87 = sigmaaa/t5/t84
      t90 = t15*t4
      t93 = t14/t2/t90
      t98 = t24/t25/rhoa
      t104 = t32/t5/t25/t84
      t107 = t32*sigmaaa
      t111 = t107/t2/t25/t90
      t113 = 1/t37/t10
      t116 = 0.5528949333333333D-2*t87*t11-0.3144353706666667D-3*t93
     #*t21+0.2515331413333333D-5*t98*t29-0.1163347968D-7*t104*t38
     #+0.249970688D-10*t111*t113
      t119 = t47*t60
      t122 = 1/t54
      t123 = t44*t122
      t126 = t56**2
      t127 = 1/t126
      t128 = t48*t127
      t129 = t49**2
      t130 = t129**2
      t131 = t130*t49
      t132 = 1/t131
      t133 = 1/t4
      t138 = 1/t52
      t141 = 1/t45
      t144 = -0.1853395810515781D1*t132*t133-0.128158207914907D1*t122
     #*t133-0.8223668877838045D0*t138*t133-0.1603914194192128D0*t141
     #*t133
      t145 = 1/t59
      t146 = t144*t145
      t147 = t146*t78
      t159 = 1/t74/t62
      t162 = 0.1281557333333333D1*t87*t63-0.1454692266666667D1*t93*t67
     #+0.82714688D0*t98*t71-0.224126336D0*t104*t75+0.213264384D-1*t111
     #*t159
      t163 = t60*t162
      vrhoa(i) = -0.12407009817988D1*t2*t41-0.9305257363491D0*t3*t116
     #-0.3109D-1*t119*t78+0.1321010150222857D-2*t123*t79+0.1D1*t128
     #*t147-0.3109D-1*t48*t163
      vrhob(i) = 0.D0
      t168 = sigmaaa*t18
      t171 = t14*t26
      t174 = t24*t35
      t179 = 1/t2/t25/t16
      t180 = t32*t179
      vsigmaaa(i) = -0.9305257363491D0*t3*(-0.2073356D-2*t7*t11
     #+0.117913264D-3*t168*t21-0.94324928D-6*t171*t29+0.436255488D-8
     #*t174*t38-0.93739008D-11*t180*t113)-0.3109D-1*t48*t60*(
     #-0.480584D0*t7*t63+0.5455096D0*t168*t67-0.31018008D0*t171*t71
     #+0.84047376D-1*t174*t75-0.79974144D-2*t180*t159)
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      t207 = sigmaaa/t5/t15
      t210 = t15*t84
      t213 = t14/t2/t210
      t217 = t24/t33
      t223 = t32/t5/t25/t15
      t229 = t107/t2/t25/t210
      t233 = t25**2
      t236 = t32*t14/t233/t4
      t238 = 1/t37/t20
      t249 = 1/t84
      t251 = 1/t54/t44
      t263 = t144**2
      t270 = 1/t15
      t296 = t126**2
      t299 = t59**2
      t319 = 1/t74/t66
      s1 = -0.4135669939329333D0/t5*t41-0.24814019635976D1*t2*t116
     #-0.9305257363491D0*t3*(-0.2027281422222222D-1*t207*t11
     #+0.2050399473777778D-2*t213*t21-0.2934593729422222D-4*t217*t29
     #+0.2162145348266667D-6*t223*t38-0.8546531191466667D-9*t229*t113
     #+0.1333177002666667D-11*t236*t238)+0.2D1*t47*t127*t147-0.6218D-1
     #*t119*t162+0.8806734334819047D-3*t249*t251*t79
      v2rhoa2(i) = s1-0.8497974591333914D-1*t123*t127*t147
     #+0.2642020300445714D-2*t123*t163-0.2D1*t48/t126/t56*t263*t145
     #*t78+0.1D1*t128*(-0.1544496508763151D1/t131/t44*t270
     #+0.3706791621031562D1*t132*t249-0.854388052766047D0*t251*t270
     #+0.2563164158298141D1*t122*t249-0.4111834438919023D0/t52/t44
     #*t270+0.1644733775567609D1*t138*t249-0.5346380647307093D-1/t45
     #/t44*t270+0.3207828388384256D0*t141*t249)*t145*t78
     #+0.321646831778707D2*t48/t296*t263/t299*t78+0.2D1*t128*t146*t162
     #-0.3109D-1*t48*t60*(-0.4699043555555556D1*t207*t63
     #+0.9896548266666667D1*t213*t67-0.8995993671111111D1*t217*t71
     #+0.3938242261333333D1*t223*t75-0.7838151338666667D0*t229*t159
     #+0.568705024D-1*t236*t319)
      v2rhob2(i) = 0.D0
      v2rhoab(i) = 0.D0
      t328 = sigmaaa*t26
      t331 = t14*t35
      t334 = t24*t179
      t338 = t32/t233
      v2sigmaaa2(i) = -0.9305257363491D0*t3*(0.126206688D-3*t18*t21
     #-0.2829804672D-5*t328*t29+0.24406656D-7*t331*t38
     #-0.10729648128D-9*t334*t113+0.187478016D-12*t338*t238)-0.3109D-1
     #*t48*t60*(0.6416264D0*t18*t67-0.838564D0*t328*t71+0.438250176D0
     #*t331*t75-0.992275584D-1*t334*t159+0.79974144D-2*t338*t319)
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      v2sigmabb2(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rhoa**(1.D0/3.D0)
      t5 = t4*rhoa
      t6 = rhoa**2
      t7 = t4**2
      t9 = 1/t7/t6
      t10 = sigmaaa*t9
      t12 = 1.D0+0.4D-2*t10
      t13 = 1/t12
      t16 = sigmaaa**2
      t17 = t6**2
      t18 = t17*rhoa
      t20 = 1/t4/t18
      t21 = t16*t20
      t22 = t12**2
      t23 = 1/t22
      t26 = t16*sigmaaa
      t27 = t17**2
      t28 = 1/t27
      t29 = t26*t28
      t31 = 1/t22/t12
      t34 = t16**2
      t35 = t27*t6
      t37 = 1/t7/t35
      t38 = t34*t37
      t39 = t22**2
      t40 = 1/t39
      t43 = 0.108184D1-0.2073356D-2*t10*t13+0.5480992D-4*t21*t23
     #-0.16825664D-6*t29*t31+0.5858688D-9*t38*t40
      t46 = 1/rhoa
      t47 = t46**(1.D0/3.D0)
      t49 = 1.D0+0.1274696188700087D0*t47
      t50 = rhoa*t49
      t51 = t46**(1.D0/6.D0)
      t54 = dsqrt(t46)
      t56 = t47**2
      t58 = 0.1112037486309468D2*t51+0.3844746237447211D1*t47
     #+0.1644733775567609D1*t54+0.2405871291288192D0*t56
      t61 = 1.D0+0.321646831778707D2/t58
      t62 = dlog(t61)
      t64 = 1.D0+0.2D0*t10
      t65 = 1/t64
      t68 = t64**2
      t69 = 1/t68
      t73 = 1/t68/t64
      t76 = t68**2
      t77 = 1/t76
      t80 = 0.118777D1-0.480584D0*t10*t65+0.2246964D0*t21*t69
     #-0.7343384D-1*t29*t73+0.9996768D-2*t38*t77
      t81 = t62*t80
      t84 = rhob**(1.D0/3.D0)
      t85 = t84*rhob
      t86 = rhob**2
      t87 = t84**2
      t89 = 1/t87/t86
      t90 = sigmabb*t89
      t92 = 1.D0+0.4D-2*t90
      t93 = 1/t92
      t96 = sigmabb**2
      t97 = t86**2
      t98 = t97*rhob
      t100 = 1/t84/t98
      t101 = t96*t100
      t102 = t92**2
      t103 = 1/t102
      t106 = t96*sigmabb
      t107 = t97**2
      t108 = 1/t107
      t109 = t106*t108
      t111 = 1/t102/t92
      t114 = t96**2
      t115 = t107*t86
      t117 = 1/t87/t115
      t118 = t114*t117
      t119 = t102**2
      t120 = 1/t119
      t123 = 0.108184D1-0.2073356D-2*t90*t93+0.5480992D-4*t101*t103
     #-0.16825664D-6*t109*t111+0.5858688D-9*t118*t120
      t126 = 1/rhob
      t127 = t126**(1.D0/3.D0)
      t129 = 1.D0+0.1274696188700087D0*t127
      t130 = rhob*t129
      t131 = t126**(1.D0/6.D0)
      t134 = dsqrt(t126)
      t136 = t127**2
      t138 = 0.1112037486309468D2*t131+0.3844746237447211D1*t127
     #+0.1644733775567609D1*t134+0.2405871291288192D0*t136
      t141 = 1.D0+0.321646831778707D2/t138
      t142 = dlog(t141)
      t144 = 1.D0+0.2D0*t90
      t145 = 1/t144
      t148 = t144**2
      t149 = 1/t148
      t153 = 1/t148/t144
      t156 = t148**2
      t157 = 1/t156
      t160 = 0.118777D1-0.480584D0*t90*t145+0.2246964D0*t101*t149
     #-0.7343384D-1*t109*t153+0.9996768D-2*t118*t157
      t161 = t142*t160
      t164 = rhoa+rhob
      t165 = 1/t164
      t166 = t165**(1.D0/3.D0)
      t168 = 1.D0+0.1325688999052018D0*t166
      t169 = t165**(1.D0/6.D0)
      t172 = dsqrt(t165)
      t174 = t166**2
      t176 = 0.598255043577108D1*t169+0.2225569421150687D1*t166
     #+0.8004286349993634D0*t172+0.1897004325747559D0*t174
      t179 = 1.D0+0.160818243221511D2/t176
      t180 = dlog(t179)
      t182 = 0.62182D-1*t168*t180
      t184 = 1.D0+0.6901399211255825D-1*t166
      t189 = 0.8157414703487641D1*t169+0.2247591863577616D1*t166
     #+0.4300972471276643D0*t172+0.1911512595127338D0*t174
      t192 = 1.D0+0.2960857464321668D2/t189
      t193 = dlog(t192)
      t194 = t184*t193
      t196 = rhoa-1.D0*rhob
      t197 = t196*t165
      t198 = 1.D0+t197
      t199 = t198**(1.D0/3.D0)
      t202 = 1.D0-1.D0*t197
      t203 = t202**(1.D0/3.D0)
      t205 = t199*t198+t203*t202-2.D0
      t206 = t196**2
      t207 = t206**2
      t208 = t164**2
      t209 = t208**2
      t210 = 1/t209
      t211 = t207*t210
      t213 = 1.D0-1.D0*t211
      t214 = t205*t213
      t216 = 0.3799574853701528D-1*t194*t214
      t218 = 1.D0+0.1274696188700087D0*t166
      t223 = 0.1112037486309468D2*t169+0.3844746237447211D1*t166
     #+0.1644733775567609D1*t172+0.2405871291288192D0*t174
      t226 = 1.D0+0.321646831778707D2/t223
      t227 = dlog(t226)
      t230 = -0.3109D-1*t218*t227+t182
      t231 = t230*t205
      t233 = 0.1923661050931536D1*t231*t211
      t240 = t164*(-t182+t216+t233)+0.3109D-1*t50*t62+0.3109D-1*t130
     #*t142
      t243 = 0.5D0*t10+0.5D0*t90
      t246 = 1.D0+0.3D-2*t10+0.3D-2*t90
      t247 = 1/t246
      t250 = t243**2
      t251 = t246**2
      t252 = 1/t251
      t255 = t250*t243
      t257 = 1/t251/t246
      t260 = t250**2
      t261 = t251**2
      t262 = 1/t261
      t265 = 0.589076D0+0.2654244D-1*t243*t247-0.6919848D-3*t250*t252
     #+0.91955736D-5*t255*t257-0.544387392D-7*t260*t262
      zk(i) = -0.9305257363491D0*t5*t43-0.3109D-1*t50*t81
     #-0.9305257363491D0*t85*t123-0.3109D-1*t130*t161+t240*t265
      t269 = t6*rhoa
      t271 = 1/t7/t269
      t272 = sigmaaa*t271
      t275 = t17*t6
      t277 = 1/t4/t275
      t278 = t16*t277
      t282 = 1/t27/rhoa
      t283 = t26*t282
      t288 = 1/t7/t27/t269
      t289 = t34*t288
      t292 = t34*sigmaaa
      t295 = 1/t4/t27/t275
      t296 = t292*t295
      t298 = 1/t39/t12
      t301 = 0.5528949333333333D-2*t272*t13-0.3144353706666667D-3*t278
     #*t23+0.2515331413333333D-5*t283*t31-0.1163347968D-7*t289*t40
     #+0.249970688D-10*t296*t298
      t304 = t49*t62
      t307 = 1/t56
      t308 = t46*t307
      t311 = t58**2
      t312 = 1/t311
      t313 = t50*t312
      t314 = t51**2
      t315 = t314**2
      t316 = t315*t51
      t317 = 1/t316
      t318 = 1/t6
      t323 = 1/t54
      t326 = 1/t47
      t329 = -0.1853395810515781D1*t317*t318-0.128158207914907D1*t307
     #*t318-0.8223668877838045D0*t323*t318-0.1603914194192128D0*t326
     #*t318
      t330 = 1/t61
      t331 = t329*t330
      t332 = t331*t80
      t344 = 1/t76/t64
      t347 = 0.1281557333333333D1*t272*t65-0.1454692266666667D1*t278
     #*t69+0.82714688D0*t283*t73-0.224126336D0*t289*t77+0.213264384D-1
     #*t296*t344
      t348 = t62*t347
      t351 = 1/t174
      t352 = 1/t208
      t353 = t351*t352
      t354 = t353*t180
      t355 = 0.2747799777968419D-2*t354
      t356 = t176**2
      t357 = 1/t356
      t358 = t168*t357
      t359 = t169**2
      t360 = t359**2
      t361 = t360*t169
      t362 = 1/t361
      t363 = t362*t352
      t366 = 1/t172
      t367 = t366*t352
      t369 = 1/t166
      t370 = t369*t352
      t372 = -0.99709173929518D0*t363-0.7418564737168958D0*t353
     #-0.4002143174996817D0*t367-0.1264669550498372D0*t370
      t373 = 1/t179
      t375 = t358*t372*t373
      t376 = 0.1D1*t375
      t377 = t193*t205
      t378 = t377*t213
      t379 = t353*t378
      t380 = 0.8740794299481065D-3*t379
      t381 = t189**2
      t382 = 1/t381
      t383 = t184*t382
      t388 = -0.135956911724794D1*t363-0.7491972878592054D0*t353
     #-0.2150486235638321D0*t367-0.1274341730084892D0*t370
      t389 = t383*t388
      t390 = 1/t192
      t391 = t390*t205
      t392 = t391*t213
      t393 = t389*t392
      t394 = 0.1124999956683108D1*t393
      t395 = t196*t352
      t396 = 1.D0*t395
      t397 = t165-t396
      t400 = 1.D0*t165
      t401 = -t400+t395
      t404 = 0.1333333333333333D1*t199*t397+0.1333333333333333D1*t203
     #*t401
      t406 = t194*t404*t213
      t407 = 0.3799574853701528D-1*t406
      t408 = t206*t196
      t409 = t408*t210
      t410 = 4.D0*t409
      t412 = 1/t209/t164
      t413 = t207*t412
      t414 = 4.D0*t413
      t415 = -t410+t414
      t417 = t194*t205*t415
      t418 = 0.3799574853701528D-1*t417
      t421 = t223**2
      t422 = 1/t421
      t423 = t218*t422
      t428 = -0.1853395810515781D1*t363-0.128158207914907D1*t353
     #-0.8223668877838045D0*t367-0.1603914194192128D0*t370
      t429 = 1/t226
      t435 = 0.1321010150222857D-2*t353*t227+0.1D1*t423*t428*t429
     #-0.2747799777968419D-2*t354-0.1D1*t375
      t436 = t435*t205
      t437 = t436*t211
      t438 = 0.1923661050931536D1*t437
      t439 = t230*t404
      t440 = t439*t211
      t441 = 0.1923661050931536D1*t440
      t442 = t231*t409
      t444 = t231*t413
      t445 = 0.7694644203726145D1*t444
      t452 = t312*t329*t330
      t455 = -t182+t216+t233+t164*(t355+t376-t380-t394+t407+t418+t438
     #+t441+0.7694644203726145D1*t442-t445)+0.3109D-1*t304
     #-0.1321010150222857D-2*t308*t62-0.1D1*t50*t452
      t459 = t243*t252
      t462 = t250*t257
      t465 = t255*t262
      t469 = 1/t261/t246
      t470 = t260*t469
      t473 = -0.3538992D-1*t272*t247+0.205763232D-2*t459*t272
     #-0.478540512D-4*t462*t272+0.5110337088D-6*t465*t272
     #-0.17420396544D-8*t470*t272
      vrhoa(i) = -0.12407009817988D1*t4*t43-0.9305257363491D0*t5*t301
     #-0.3109D-1*t304*t80+0.1321010150222857D-2*t308*t81+0.1D1*t313
     #*t332-0.3109D-1*t50*t348+t455*t265+t240*t473
      t477 = t86*rhob
      t479 = 1/t87/t477
      t480 = sigmabb*t479
      t483 = t97*t86
      t485 = 1/t84/t483
      t486 = t96*t485
      t490 = 1/t107/rhob
      t491 = t106*t490
      t496 = 1/t87/t107/t477
      t497 = t114*t496
      t500 = t114*sigmabb
      t503 = 1/t84/t107/t483
      t504 = t500*t503
      t506 = 1/t119/t92
      t509 = 0.5528949333333333D-2*t480*t93-0.3144353706666667D-3*t486
     #*t103+0.2515331413333333D-5*t491*t111-0.1163347968D-7*t497*t120
     #+0.249970688D-10*t504*t506
      t512 = t129*t142
      t515 = 1/t136
      t516 = t126*t515
      t519 = t138**2
      t520 = 1/t519
      t521 = t130*t520
      t522 = t131**2
      t523 = t522**2
      t524 = t523*t131
      t525 = 1/t524
      t526 = 1/t86
      t531 = 1/t134
      t534 = 1/t127
      t537 = -0.1853395810515781D1*t525*t526-0.128158207914907D1*t515
     #*t526-0.8223668877838045D0*t531*t526-0.1603914194192128D0*t534
     #*t526
      t538 = 1/t141
      t539 = t537*t538
      t540 = t539*t160
      t552 = 1/t156/t144
      t555 = 0.1281557333333333D1*t480*t145-0.1454692266666667D1*t486
     #*t149+0.82714688D0*t491*t153-0.224126336D0*t497*t157
     #+0.213264384D-1*t504*t552
      t556 = t142*t555
      t559 = -t400-t396
      t562 = t165+t395
      t565 = 0.1333333333333333D1*t199*t559+0.1333333333333333D1*t203
     #*t562
      t567 = t194*t565*t213
      t568 = 0.3799574853701528D-1*t567
      t569 = t410+t414
      t571 = t194*t205*t569
      t572 = 0.3799574853701528D-1*t571
      t573 = t230*t565
      t574 = t573*t211
      t575 = 0.1923661050931536D1*t574
      t583 = t520*t537*t538
      t586 = -t182+t216+t233+t164*(t355+t376-t380-t394+t568+t572+t438
     #+t575-0.7694644203726145D1*t442-t445)+0.3109D-1*t512
     #-0.1321010150222857D-2*t516*t142-0.1D1*t130*t583
      t598 = -0.3538992D-1*t480*t247+0.205763232D-2*t459*t480
     #-0.478540512D-4*t462*t480+0.5110337088D-6*t465*t480
     #-0.17420396544D-8*t470*t480
      vrhob(i) = -0.12407009817988D1*t84*t123-0.9305257363491D0*t85
     #*t509-0.3109D-1*t512*t160+0.1321010150222857D-2*t516*t161+0.1D1
     #*t521*t540-0.3109D-1*t130*t556+t586*t265+t240*t598
      t602 = sigmaaa*t20
      t605 = t16*t28
      t608 = t26*t37
      t613 = 1/t4/t27/t18
      t614 = t34*t613
      t617 = -0.2073356D-2*t9*t13+0.117913264D-3*t602*t23
     #-0.94324928D-6*t605*t31+0.436255488D-8*t608*t40-0.93739008D-11
     #*t614*t298
      t630 = -0.480584D0*t9*t65+0.5455096D0*t602*t69-0.31018008D0*t605
     #*t73+0.84047376D-1*t608*t77-0.79974144D-2*t614*t344
      t631 = t62*t630
      t644 = 0.1327122D-1*t9*t247-0.77161212D-3*t459*t9+0.179452692D-4
     #*t462*t9-0.1916376408D-6*t465*t9+0.6532648704D-9*t470*t9
      vsigmaaa(i) = -0.9305257363491D0*t5*t617-0.3109D-1*t50*t631+t240
     #*t644
      vsigmaab(i) = 0.D0
      t648 = sigmabb*t100
      t651 = t96*t108
      t654 = t106*t117
      t659 = 1/t84/t107/t98
      t660 = t114*t659
      t663 = -0.2073356D-2*t89*t93+0.117913264D-3*t648*t103
     #-0.94324928D-6*t651*t111+0.436255488D-8*t654*t120-0.93739008D-11
     #*t660*t506
      t676 = -0.480584D0*t89*t145+0.5455096D0*t648*t149-0.31018008D0
     #*t651*t153+0.84047376D-1*t654*t157-0.79974144D-2*t660*t552
      t677 = t142*t676
      t690 = 0.1327122D-1*t89*t247-0.77161212D-3*t459*t89
     #+0.179452692D-4*t462*t89-0.1916376408D-6*t465*t89
     #+0.6532648704D-9*t470*t89
      vsigmabb(i) = -0.9305257363491D0*t85*t663-0.3109D-1*t130*t677
     #+t240*t690
      t696 = sigmaaa/t7/t17
      t699 = t17*t269
      t702 = t16/t4/t699
      t705 = t243*t257
      t710 = t250*t262
      t715 = t255*t469
      t722 = t260/t261/t251
      t739 = t26/t35
      t745 = t34/t7/t27/t17
      t751 = t292/t4/t27/t699
      t755 = t27**2
      t758 = t34*t16/t755/t6
      t760 = 1/t39/t22
      t768 = t49*t312
      t774 = 1/t311/t58
      t776 = t329**2
      t781 = 1/t269
      t783 = 1/t56/t46
      t784 = t781*t783
      t792 = 1/t17
      t813 = -0.1544496508763151D1/t316/t46*t792+0.3706791621031562D1
     #*t317*t781-0.854388052766047D0*t783*t792+0.2563164158298141D1
     #*t307*t781-0.4111834438919023D0/t54/t46*t792
     #+0.1644733775567609D1*t323*t781-0.5346380647307093D-1/t47/t46
     #*t792+0.3207828388384256D0*t326*t781
      t821 = t311**2
      t822 = 1/t821
      t824 = t61**2
      t825 = 1/t824
      t841 = 1/t76/t68
      t848 = 0.2D1*t375
      t849 = 0.5495599555936838D-2*t354
      t851 = 0.1748158859896213D-2*t379
      t852 = 0.2249999913366216D1*t393
      t853 = 0.3847322101863073D1*t437
      t856 = 0.1538928840745229D2*t444
      t877 = 1/t208/t164
      t878 = t351*t877
      t879 = t878*t180
      t880 = 0.5495599555936838D-2*t879
      t883 = 1/t174/t165*t210
      t885 = 0.5827196199654043D-3*t883*t378
      t886 = t883*t180
      t887 = 0.1831866518645613D-2*t886
      t888 = t408*t412
      t889 = t231*t888
      t892 = t435*t404*t211
      t895 = t353*t377*t415
      t898 = 0.1748158859896213D-2*t878*t378
      t901 = t353*t193*t404*t213
      t903 = t199**2
      t904 = 1/t903
      t905 = t397**2
      t908 = 2.D0*t352
      t910 = 2.D0*t196*t877
      t911 = -t908+t910
      t914 = t203**2
      t915 = 1/t914
      t916 = t401**2
      t922 = 0.4444444444444444D0*t904*t905+0.1333333333333333D1*t199
     #*t911+0.4444444444444444D0*t915*t916-0.1333333333333333D1*t203
     #*t911
      t928 = t389*t390*t404*t213
      t932 = t353*t357*t372*t373
      t933 = 0.8837926660346786D-1*t932
      t935 = t389*t391*t415
      t937 = t381**2
      t940 = t388**2
      t942 = t192**2
      t947 = 0.3330964519106732D2*t184/t937*t940/t942*t205*t213
      t950 = 1/t361/t165*t210
      t952 = t362*t877
      t958 = 1/t172/t165*t210
      t960 = t366*t877
      t964 = 1/t166/t165*t210
      t966 = t369*t877
      t971 = 0.1124999956683108D1*t383*(-0.1132974264373283D1*t950
     #+0.271913823449588D1*t952-0.4994648585728036D0*t883
     #+0.1498394575718411D1*t878-0.1075243117819161D0*t958
     #+0.4300972471276643D0*t960-0.4247805766949639D-1*t964
     #+0.2548683460169784D0*t966)*t392
      t972 = -t880-t885+t887-0.6155715362980916D2*t889
     #+0.3847322101863073D1*t892-0.1748158859896213D-2*t895+t898
     #-0.1748158859896213D-2*t901+0.3799574853701528D-1*t194*t922*t213
     #-0.2249999913366216D1*t928-t933-0.2249999913366216D1*t935-t947
     #-t971
      t977 = 0.5176049209143758D-1*t353*t382*t388*t390*t214
      t983 = 0.2249999913366216D1*t184/t381/t189*t940*t392
      t987 = t372**2
      t989 = t168/t356/t176*t987*t373
      t990 = 0.2D1*t989
      t1001 = t358*(-0.8309097827459833D0*t950+0.199418347859036D1
     #*t952-0.4945709824779306D0*t883+0.1483712947433792D1*t878
     #-0.2001071587498409D0*t958+0.8004286349993634D0*t960
     #-0.4215565168327908D-1*t964+0.2529339100996745D0*t966)*t373
      t1002 = 0.1D1*t1001
      t1003 = t356**2
      t1006 = t179**2
      t1009 = t168/t1003*t987/t1006
      t1010 = 0.160818243221511D2*t1009
      t1022 = t428**2
      t1038 = t421**2
      t1041 = t226**2
      t1052 = 0.8806734334819047D-3*t883*t227-0.2642020300445714D-2
     #*t878*t227-0.8497974591333914D-1*t353*t422*t428*t429-0.2D1*t218
     #/t421/t223*t1022*t429+0.1D1*t423*(-0.1544496508763151D1*t950
     #+0.3706791621031562D1*t952-0.854388052766047D0*t883
     #+0.2563164158298141D1*t878-0.4111834438919023D0*t958
     #+0.1644733775567609D1*t960-0.5346380647307093D-1*t964
     #+0.3207828388384256D0*t966)*t429+0.321646831778707D2*t218/t1038
     #*t1022/t1041-0.1831866518645613D-2*t886+0.5495599555936838D-2
     #*t879+0.8837926660346786D-1*t932+0.2D1*t989-0.1D1*t1001
     #-0.160818243221511D2*t1009
      t1055 = 0.1923661050931536D1*t1052*t205*t211
      t1056 = t436*t409
      t1061 = t439*t409
      t1063 = t206*t210
      t1064 = 12.D0*t1063
      t1065 = 32.D0*t888
      t1068 = t207/t209/t208
      t1069 = 20.D0*t1068
      t1075 = 0.3847322101863073D2*t231*t1068
      t1076 = t439*t413
      t1078 = t231*t1063
      t1079 = 0.2308393261117844D2*t1078
      t1084 = 0.1538928840745229D2*t436*t413
      t1085 = t977+t983-t990+t1002+t1010+t1055+0.1538928840745229D2
     #*t1056+0.1923661050931536D1*t230*t922*t211+0.1538928840745229D2
     #*t1061+0.3799574853701528D-1*t194*t205*(-t1064+t1065-t1069)
     #+t1075-0.1538928840745229D2*t1076+t1079+0.7599149707403056D-1
     #*t194*t404*t415-t1084
      t1088 = t848+t849+0.7599149707403056D-1*t417-t851-t852+t853
     #+0.7599149707403056D-1*t406+0.1538928840745229D2*t442-t856
     #+0.3847322101863073D1*t440-0.2D1*t768*t331-0.321646831778707D2
     #*t50*t822*t776*t825-0.1D1*t50*t312*t813*t330
     #-0.8806734334819047D-3*t784*t62+0.8497974591333914D-1*t308*t452
     #+0.2D1*t50*t774*t776*t330+t164*(t972+t1085)
      s1 = 2.D0*t455*t473+t240*(0.12976304D0*t696*t247-0.302662912D-2
     #*t702*t252+0.16053292032D-3*t705*t702-0.754465184D-2*t459*t696
     #-0.3192632064D-5*t710*t702+0.1754648544D-3*t462*t696
     #+0.256439568384D-7*t715*t702-0.18737902656D-5*t465*t696
     #-0.69681586176D-10*t722*t702+0.63874787328D-8*t470*t696)
     #-0.4135669939329333D0/t7*t43-0.24814019635976D1*t4*t301
     #-0.9305257363491D0*t5*(-0.2027281422222222D-1*t696*t13
     #+0.2050399473777778D-2*t702*t23-0.2934593729422222D-4*t739*t31
     #+0.2162145348266667D-6*t745*t40-0.8546531191466667D-9*t751*t298
     #+0.1333177002666667D-11*t758*t760)-0.6218D-1*t304*t347+0.2D1
     #*t768*t332+0.2642020300445714D-2*t308*t348
      v2rhoa2(i) = s1-0.2D1*t50*t774*t776*t330*t80
     #+0.8806734334819047D-3*t784*t81-0.8497974591333914D-1*t308*t312
     #*t332+0.1D1*t313*t813*t330*t80+0.2D1*t313*t331*t347
     #+0.321646831778707D2*t50*t822*t776*t825*t80-0.3109D-1*t50*t62*(
     #-0.4699043555555556D1*t696*t65+0.9896548266666667D1*t702*t69
     #-0.8995993671111111D1*t739*t73+0.3938242261333333D1*t745*t77
     #-0.7838151338666667D0*t751*t344+0.568705024D-1*t758*t841)+t1088
     #*t265
      t1097 = sigmabb/t87/t97
      t1100 = t97*t477
      t1103 = t96/t84/t1100
      t1107 = t106/t115
      t1113 = t114/t87/t107/t97
      t1119 = t500/t84/t107/t1100
      t1123 = t107**2
      t1126 = t114*t96/t1123/t86
      t1128 = 1/t119/t102
      t1139 = t519**2
      t1140 = 1/t1139
      t1142 = t537**2
      t1143 = t141**2
      t1144 = 1/t1143
      t1154 = t559**2
      t1157 = t908+t910
      t1160 = t562**2
      t1166 = 0.4444444444444444D0*t904*t1154+0.1333333333333333D1
     #*t199*t1157+0.4444444444444444D0*t915*t1160-0.1333333333333333D1
     #*t203*t1157
      t1171 = t435*t565*t211
      t1173 = -t880-t885+t887+0.6155715362980916D2*t889+t898-t933-t947
     #-t971+t977+t983-t990+t1002+0.3799574853701528D-1*t194*t1166*t213
     #+0.3847322101863073D1*t1171
      t1184 = t573*t413
      t1186 = t573*t409
      t1189 = t389*t391*t569
      t1194 = t353*t193*t565*t213
      t1197 = t353*t377*t569
      t1201 = t389*t390*t565*t213
      t1203 = 0.7599149707403056D-1*t194*t565*t569
     #+0.1923661050931536D1*t230*t1166*t211+0.3799574853701528D-1*t194
     #*t205*(-t1064-t1065-t1069)-0.1538928840745229D2*t1184+t1010
     #+t1055-0.1538928840745229D2*t1186-0.2249999913366216D1*t1189
     #-0.1538928840745229D2*t1056-0.1748158859896213D-2*t1194
     #-0.1748158859896213D-2*t1197-0.2249999913366216D1*t1201+t1075
     #+t1079-t1084
      t1206 = t129*t520
      t1209 = 1/t477
      t1211 = 1/t136/t126
      t1212 = t1209*t1211
      t1218 = 1/t519/t138
      t1225 = 1/t97
      t1246 = -0.1544496508763151D1/t524/t126*t1225
     #+0.3706791621031562D1*t525*t1209-0.854388052766047D0*t1211*t1225
     #+0.2563164158298141D1*t515*t1209-0.4111834438919023D0/t134/t126
     #*t1225+0.1644733775567609D1*t531*t1209-0.5346380647307093D-1
     #/t127/t126*t1225+0.3207828388384256D0*t534*t1209
      t1255 = t848+t849-t851-t852+t853-0.1538928840745229D2*t442-t856
     #+0.7599149707403056D-1*t567+0.7599149707403056D-1*t571
     #+0.3847322101863073D1*t574+t164*(t1173+t1203)-0.2D1*t1206*t539
     #-0.8806734334819047D-3*t1212*t142+0.8497974591333914D-1*t516
     #*t583+0.2D1*t130*t1218*t1142*t538-0.1D1*t130*t520*t1246*t538
     #-0.321646831778707D2*t130*t1140*t1142*t1144
      t1272 = 1/t156/t148
      s1 = -0.4135669939329333D0/t87*t123-0.24814019635976D1*t84*t509
     #-0.9305257363491D0*t85*(-0.2027281422222222D-1*t1097*t93
     #+0.2050399473777778D-2*t1103*t103-0.2934593729422222D-4*t1107
     #*t111+0.2162145348266667D-6*t1113*t120-0.8546531191466667D-9
     #*t1119*t506+0.1333177002666667D-11*t1126*t1128)-0.6218D-1*t512
     #*t555+0.2D1*t521*t539*t555+0.321646831778707D2*t130*t1140*t1142
     #*t1144*t160+t1255*t265+0.2D1*t1206*t540
      s2 = s1+0.2642020300445714D-2*t516*t556-0.3109D-1*t130*t142*(
     #-0.4699043555555556D1*t1097*t145+0.9896548266666667D1*t1103*t149
     #-0.8995993671111111D1*t1107*t153+0.3938242261333333D1*t1113*t157
     #-0.7838151338666667D0*t1119*t552+0.568705024D-1*t1126*t1272)
     #+0.8806734334819047D-3*t1212*t161
      v2rhob2(i) = s2-0.8497974591333914D-1*t516*t520*t540-0.2D1*t130
     #*t1218*t1142*t538*t160+0.1D1*t521*t1246*t538*t160+2.D0*t586*t598
     #+t240*(0.12976304D0*t1097*t247-0.302662912D-2*t1103*t252
     #+0.16053292032D-3*t705*t1103-0.754465184D-2*t459*t1097
     #-0.3192632064D-5*t710*t1103+0.1754648544D-3*t462*t1097
     #+0.256439568384D-7*t715*t1103-0.18737902656D-5*t465*t1097
     #-0.69681586176D-10*t722*t1103+0.63874787328D-8*t470*t1097)
      t1323 = -t880-t885+t887+0.1923661050931536D1*t892
     #-0.8740794299481065D-3*t895+t898-0.8740794299481065D-3*t901
     #-0.1124999956683108D1*t928-t933-0.1124999956683108D1*t935-t947
     #-t971+t977+t983-t990+t1002+0.1923661050931536D1*t1171
      t1352 = 0.4444444444444444D0*t904*t397*t559+0.2666666666666667D1
     #*t199*t196*t877+0.4444444444444444D0*t915*t401*t562
     #-0.2666666666666667D1*t203*t196*t877
      t1362 = -0.7694644203726145D1*t1184+t1010+t1055
     #+0.7694644203726145D1*t1186-0.1124999956683108D1*t1189
     #+0.3799574853701528D-1*t194*t205*(t1064-t1069)
     #-0.7694644203726145D1*t1061-0.8740794299481065D-3*t1194
     #-0.8740794299481065D-3*t1197-0.1124999956683108D1*t1201
     #+0.3799574853701528D-1*t194*t404*t569+t1075-0.7694644203726145D1
     #*t1076-0.2308393261117844D2*t1078-t1084+0.3799574853701528D-1
     #*t194*t1352*t213+0.3799574853701528D-1*t194*t565*t415
     #+0.1923661050931536D1*t230*t1352*t211
      t1365 = t849+t848-t851-t852+t568+t572+t853+t575-t856+t407+t418
     #+t441+t164*(t1323+t1362)
      t1375 = t271*sigmabb*t479
      v2rhoab(i) = t1365*t265+t455*t598+t586*t473+t240*(
     #-0.302662912D-2*t272*t252*sigmabb*t479+0.16053292032D-3*t705
     #*sigmaaa*t1375-0.3192632064D-5*t710*sigmaaa*t1375
     #+0.256439568384D-7*t715*sigmaaa*t1375-0.69681586176D-10*t722
     #*sigmaaa*t1375)
      t1393 = sigmaaa*t277
      t1396 = t16*t282
      t1399 = t26*t288
      t1402 = t34*t295
      t1407 = t292/t755/rhoa
      s1 = -0.12407009817988D1*t4*t617-0.9305257363491D0*t5*
     #(0.5528949333333333D-2*t271*t13-0.6509865386666667D-3*t1393*t23
     #+0.1006147720533333D-4*t1396*t31-0.7671789568D-7*t1399*t40
     #+0.31112101888D-9*t1402*t298-0.499941376D-12*t1407*t760)
     #-0.3109D-1*t304*t630+0.1321010150222857D-2*t308*t631
      v2rhoasigmaaa(i) = s1+0.1D1*t313*t331*t630-0.3109D-1*t50*t62*
     #(0.1281557333333333D1*t271*t65-0.3165696D1*t1393*t69
     #+0.3063317546666667D1*t1396*t73-0.1392793472D1*t1399*t77
     #+0.2859332608D0*t1402*t344-0.213264384D-1*t1407*t841)+t455*t644
     #+t240*(-0.3538992D-1*t271*t247+0.113498592D-2*t1393*t252
     #-0.6019984512D-4*t705*t1393+0.205763232D-2*t459*t271
     #+0.1197237024D-5*t710*t1393-0.478540512D-4*t462*t271
     #-0.96164838144D-8*t715*t1393+0.5110337088D-6*t465*t271
     #+0.26130594816D-10*t722*t1393-0.17420396544D-8*t470*t271)
      v2rhoasigmaab(i) = 0.D0
      t1463 = t272*t89
      v2rhoasigmabb(i) = t455*t690+t240*(0.113498592D-2*t272*t252*t89
     #-0.6019984512D-4*t705*t1463+0.1197237024D-5*t710*t1463
     #-0.96164838144D-8*t715*t1463+0.26130594816D-10*t722*t1463)
      t1475 = t252*t9
      t1478 = t480*t9
      v2rhobsigmaaa(i) = t586*t644+t240*(0.113498592D-2*t480*t1475
     #-0.6019984512D-4*t705*t1478+0.1197237024D-5*t710*t1478
     #-0.96164838144D-8*t715*t1478+0.26130594816D-10*t722*t1478)
      v2rhobsigmaab(i) = 0.D0
      t1493 = sigmabb*t485
      t1496 = t96*t490
      t1499 = t106*t496
      t1502 = t114*t503
      t1507 = t500/t1123/rhob
      s1 = -0.12407009817988D1*t84*t663-0.9305257363491D0*t85*
     #(0.5528949333333333D-2*t479*t93-0.6509865386666667D-3*t1493*t103
     #+0.1006147720533333D-4*t1496*t111-0.7671789568D-7*t1499*t120
     #+0.31112101888D-9*t1502*t506-0.499941376D-12*t1507*t1128)
     #-0.3109D-1*t512*t676+0.1321010150222857D-2*t516*t677
      v2rhobsigmabb(i) = s1+0.1D1*t521*t539*t676-0.3109D-1*t130*t142*
     #(0.1281557333333333D1*t479*t145-0.3165696D1*t1493*t149
     #+0.3063317546666667D1*t1496*t153-0.1392793472D1*t1499*t157
     #+0.2859332608D0*t1502*t552-0.213264384D-1*t1507*t1272)+t586*t690
     #+t240*(-0.3538992D-1*t479*t247+0.113498592D-2*t1493*t252
     #-0.6019984512D-4*t705*t1493+0.205763232D-2*t459*t479
     #+0.1197237024D-5*t710*t1493-0.478540512D-4*t462*t479
     #-0.96164838144D-8*t715*t1493+0.5110337088D-6*t465*t479
     #+0.26130594816D-10*t722*t1493-0.17420396544D-8*t470*t479)
      t1561 = sigmaaa*t28
      t1564 = t16*t37
      t1567 = t26*t613
      t1571 = t34/t755
      v2sigmaaa2(i) = -0.9305257363491D0*t5*(0.126206688D-3*t20*t23
     #-0.2829804672D-5*t1561*t31+0.24406656D-7*t1564*t40
     #-0.10729648128D-9*t1567*t298+0.187478016D-12*t1571*t760)
     #-0.3109D-1*t50*t62*(0.6416264D0*t20*t69-0.838564D0*t1561*t73
     #+0.438250176D0*t1564*t77-0.992275584D-1*t1567*t344+0.79974144D-2
     #*t1571*t841)+t240*(-0.42561972D-3*t20*t252+0.2257494192D-4*t705
     #*t20-0.448963884D-6*t710*t20+0.36061814304D-8*t715*t20
     #-0.9798973056D-11*t722*t20)
      v2sigmaaaab(i) = 0.D0
      t1605 = t9*t89
      v2sigmaaabb(i) = t240*(-0.42561972D-3*t1475*t89+0.2257494192D-4
     #*t705*t1605-0.448963884D-6*t710*t1605+0.36061814304D-8*t715
     #*t1605-0.9798973056D-11*t722*t1605)
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      t1617 = sigmabb*t108
      t1620 = t96*t117
      t1623 = t106*t659
      t1627 = t114/t1123
      v2sigmabb2(i) = -0.9305257363491D0*t85*(0.126206688D-3*t100*t103
     #-0.2829804672D-5*t1617*t111+0.24406656D-7*t1620*t120
     #-0.10729648128D-9*t1623*t506+0.187478016D-12*t1627*t1128)
     #-0.3109D-1*t130*t142*(0.6416264D0*t100*t149-0.838564D0*t1617
     #*t153+0.438250176D0*t1620*t157-0.992275584D-1*t1623*t552
     #+0.79974144D-2*t1627*t1272)+t240*(-0.42561972D-3*t100*t252
     #+0.2257494192D-4*t705*t100-0.448963884D-6*t710*t100
     #+0.36061814304D-8*t715*t100-0.9798973056D-11*t722*t100)
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vrhob(i) = 0.0d0
      v2rhoa2(i) = 0.0d0
      v2rhob2(i) = 0.0d0
      v2rhoab(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      vsigmaab(i) = 0.0d0
      vsigmabb(i) = 0.0d0
      v2rhoasigmaaa(i) = 0.0d0
      v2rhoasigmaab(i) = 0.0d0
      v2rhoasigmabb(i) = 0.0d0
      v2rhobsigmaaa(i) = 0.0d0
      v2rhobsigmaab(i) = 0.0d0
      v2rhobsigmabb(i) = 0.0d0
      v2sigmaaa2(i) = 0.0d0
      v2sigmaab2(i) = 0.0d0
      v2sigmabb2(i) = 0.0d0
      v2sigmaaaab(i) = 0.0d0
      v2sigmaaabb(i) = 0.0d0
      v2sigmaabbb(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end
      
      
      subroutine rks_xc_hcth407
     & (ideriv,npt,rhoa1,sigmaaa1,
     &  zk,vrhoa,vsigmaaa,
     &  v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
c
c     A.D. Boese, N.C. Handy
c     A new parametriztion of exchange-correlation generalized gradient
c     approximation functionals
c     J. Chem. Phys. 114 (2001) 5497-5503.
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit real*8 (a-h,o-z)
      integer ideriv,npt
      real*8 rhoa1(npt)
      real*8 sigmaaa1(npt)
      real*8 zk(npt),vrhoa(npt),vsigmaaa(npt)
      real*8 v2rhoa2(npt),v2rhoasigmaaa(npt),v2sigmaaa2(npt)
      parameter(tol=1.0d-20)
      
      if(ideriv.eq.0) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = rho**(1.D0/3.D0)
      t4 = rho**2
      t5 = t2**2
      t8 = sigma/t5/t4
      t10 = 1.D0+0.6349604207872798D-2*t8
      t14 = sigma**2
      t15 = t4**2
      t19 = t14/t2/t15/rho
      t20 = t10**2
      t25 = t15**2
      t27 = t14*sigma/t25
      t32 = t14**2
      t36 = t32/t5/t25/t4
      t37 = t20**2
      t44 = 1/rho
      t45 = t44**(1.D0/3.D0)
      t48 = rho*(1.D0+0.1606016560364007D0*t45)
      t49 = t44**(1.D0/6.D0)
      t52 = dsqrt(t44)
      t54 = t45**2
      t60 = dlog(1.D0+0.321646831778707D2/(0.1248219874679732D2*t49
     #+0.4844076716063854D1*t45+0.2326004811900819D1*t52
     #+0.3819082618690966D0*t54))
      t62 = 1.D0+0.3174802103936399D0*t8
      t66 = t62**2
      t74 = t66**2
      t93 = dlog(1.D0+0.160818243221511D2/(0.598255043577108D1*t49
     #+0.2225569421150687D1*t45+0.8004286349993634D0*t52
     #+0.1897004325747559D0*t54))
      t100 = 1.D0+0.9524406311809197D-2*t8
      t104 = t100**2
      t112 = t104**2
      zk(i) = -0.7385587663820224D0*t2*rho*(0.108184D1
     #-0.3291247495504578D-2*t8/t10+0.138112343902108D-3*t19/t20
     #-0.67302656D-6*t27/t20/t10+0.3720034997741387D-8*t36/t37)
     #-0.3109D-1*t48*t60*(0.118777D1-0.7628795471590852D0*t8/t62
     #+0.5661994483911968D0*t19/t66-0.29373536D0*t27/t66/t62
     #+0.6347552015792813D-1*t36/t74)+(-0.62182D-1*rho*(1.D0
     #+0.1325688999052018D0*t45)*t93+0.3109D-1*t48*t60)*(0.589076D0
     #+0.4213349717780282D-1*t8/t100-0.1743692431454588D-2*t19/t104
     #+0.367822944D-4*t27/t104/t100-0.3456644474956098D-6*t36/t112)
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = rho**(1.D0/3.D0)
      t3 = t2*rho
      t4 = rho**2
      t5 = t2**2
      t7 = 1/t5/t4
      t8 = sigma*t7
      t10 = 1.D0+0.6349604207872798D-2*t8
      t11 = 1/t10
      t14 = sigma**2
      t15 = t4**2
      t16 = t15*rho
      t18 = 1/t2/t16
      t19 = t14*t18
      t20 = t10**2
      t21 = 1/t20
      t24 = t14*sigma
      t25 = t15**2
      t26 = 1/t25
      t27 = t24*t26
      t29 = 1/t20/t10
      t32 = t14**2
      t35 = 1/t5/t25/t4
      t36 = t32*t35
      t37 = t20**2
      t38 = 1/t37
      t41 = 0.108184D1-0.3291247495504578D-2*t8*t11
     #+0.138112343902108D-3*t19*t21-0.67302656D-6*t27*t29
     #+0.3720034997741387D-8*t36*t38
      t44 = 1/rho
      t45 = t44**(1.D0/3.D0)
      t47 = 1.D0+0.1606016560364007D0*t45
      t48 = rho*t47
      t49 = t44**(1.D0/6.D0)
      t52 = dsqrt(t44)
      t54 = t45**2
      t56 = 0.1248219874679732D2*t49+0.4844076716063854D1*t45
     #+0.2326004811900819D1*t52+0.3819082618690966D0*t54
      t59 = 1.D0+0.321646831778707D2/t56
      t60 = dlog(t59)
      t62 = 1.D0+0.3174802103936399D0*t8
      t63 = 1/t62
      t66 = t62**2
      t67 = 1/t66
      t71 = 1/t66/t62
      t74 = t66**2
      t75 = 1/t74
      t78 = 0.118777D1-0.7628795471590852D0*t8*t63
     #+0.5661994483911968D0*t19*t67-0.29373536D0*t27*t71
     #+0.6347552015792813D-1*t36*t75
      t79 = t60*t78
      t83 = 1.D0+0.1325688999052018D0*t45
      t89 = 0.598255043577108D1*t49+0.2225569421150687D1*t45
     #+0.8004286349993634D0*t52+0.1897004325747559D0*t54
      t92 = 1.D0+0.160818243221511D2/t89
      t93 = dlog(t92)
      t98 = -0.62182D-1*rho*t83*t93+0.3109D-1*t48*t60
      t100 = 1.D0+0.9524406311809197D-2*t8
      t101 = 1/t100
      t104 = t100**2
      t105 = 1/t104
      t109 = 1/t104/t100
      t112 = t104**2
      t113 = 1/t112
      t116 = 0.589076D0+0.4213349717780282D-1*t8*t101
     #-0.1743692431454588D-2*t19*t105+0.367822944D-4*t27*t109
     #-0.3456644474956098D-6*t36*t113
      zk(i) = -0.7385587663820224D0*t3*t41-0.3109D-1*t48*t79+t98*t116
      t120 = t4*rho
      t123 = sigma/t5/t120
      t126 = t15*t4
      t129 = t14/t2/t126
      t134 = t24/t25/rho
      t140 = t32/t5/t25/t120
      t147 = t32*sigma/t2/t25/t126
      t149 = 1/t37/t10
      t155 = t47*t60
      t158 = 1/t54
      t159 = t44*t158
      t162 = t56**2
      t163 = 1/t162
      t165 = t49**2
      t166 = t165**2
      t169 = 1/t4
      t170 = 1/t166/t49*t169
      t172 = t158*t169
      t175 = 1/t52*t169
      t178 = 1/t45*t169
      t180 = -0.4160732915599108D1*t170-0.3229384477375903D1*t172
     #-0.2326004811900819D1*t175-0.5092110158254621D0*t178
      t181 = 1/t59
      t195 = 1/t74/t62
      t206 = t89**2
      t238 = 1/t112/t100
      s1 = -0.9847450218426965D0*t2*t41-0.3692793831910112D0*t3*
     #(0.1755331997602442D-1*t123*t11-0.1584654969337721D-2*t129*t21
     #+0.2012265130666667D-4*t134*t29-0.1477359830566614D-6*t140*t38
     #+0.503909330668646D-9*t147*t149)-0.3109D-1*t155*t78
     #+0.1664368495390566D-2*t159*t79
      vrhoa(i) = s1+0.5D0*t48*t163*t180*t181*t78-0.15545D-1*t48*t60*
     #(0.4068690918181788D1*t123*t63-0.7331189631570478D1*t129*t67
     #+0.661717504D1*t134*t71-0.2846227052321425D1*t140*t75
     #+0.4299140585511414D0*t147*t195)+(-0.62182D-1*t83*t93+rho*
     #(0.2747799777968419D-2*t172*t93+0.1D1*t83/t206*(
     #-0.99709173929518D0*t170-0.7418564737168958D0*t172
     #-0.4002143174996817D0*t175-0.1264669550498372D0*t178)/t92)
     #+0.3109D-1*t155-0.1664368495390566D-2*t159*t60-0.5D0*t48*t163
     #*t180*t181)*t116+t98*(-0.1123559924741408D0*t123*t101
     #+0.1036981709164809D-1*t129*t105-0.3828324096D-3*t134*t109
     #+0.6489723575522644D-5*t140*t113-0.351173188852824D-7*t147*t238)
      t245 = sigma*t18
      t248 = t14*t26
      t251 = t24*t35
      t257 = t32/t2/t25/t16
      vsigmaaa(i) = -0.7385587663820224D0*t3*(-0.1316498998201831D-1
     #*t7*t11+0.1188491227003291D-2*t245*t21-0.1509198848D-4*t248*t29
     #+0.110801987292496D-6*t251*t38-0.3779319980014845D-9*t257*t149)
     #-0.3109D-1*t48*t60*(-0.3051518188636341D1*t7*t63
     #+0.5498392223677858D1*t245*t67-0.496288128D1*t248*t71
     #+0.2134670289241069D1*t251*t75-0.3224355439133561D0*t257*t195)
     #+2.D0*t98*(0.8426699435560563D-1*t7*t101-0.7777362818736071D-2
     #*t245*t105+0.2871243072D-3*t248*t109-0.4867292681641983D-5*t251
     #*t113+0.263379891639618D-7*t257*t238)
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = rho**(1.D0/3.D0)
      t3 = t2*rho
      t4 = rho**2
      t5 = t2**2
      t7 = 1/t5/t4
      t8 = sigma*t7
      t10 = 1.D0+0.6349604207872798D-2*t8
      t11 = 1/t10
      t14 = sigma**2
      t15 = t4**2
      t16 = t15*rho
      t18 = 1/t2/t16
      t19 = t14*t18
      t20 = t10**2
      t21 = 1/t20
      t24 = t14*sigma
      t25 = t15**2
      t26 = 1/t25
      t27 = t24*t26
      t29 = 1/t20/t10
      t32 = t14**2
      t33 = t25*t4
      t35 = 1/t5/t33
      t36 = t32*t35
      t37 = t20**2
      t38 = 1/t37
      t41 = 0.108184D1-0.3291247495504578D-2*t8*t11
     #+0.138112343902108D-3*t19*t21-0.67302656D-6*t27*t29
     #+0.3720034997741387D-8*t36*t38
      t44 = 1/rho
      t45 = t44**(1.D0/3.D0)
      t47 = 1.D0+0.1606016560364007D0*t45
      t48 = rho*t47
      t49 = t44**(1.D0/6.D0)
      t52 = dsqrt(t44)
      t54 = t45**2
      t56 = 0.1248219874679732D2*t49+0.4844076716063854D1*t45
     #+0.2326004811900819D1*t52+0.3819082618690966D0*t54
      t59 = 1.D0+0.321646831778707D2/t56
      t60 = dlog(t59)
      t62 = 1.D0+0.3174802103936399D0*t8
      t63 = 1/t62
      t66 = t62**2
      t67 = 1/t66
      t71 = 1/t66/t62
      t74 = t66**2
      t75 = 1/t74
      t78 = 0.118777D1-0.7628795471590852D0*t8*t63
     #+0.5661994483911968D0*t19*t67-0.29373536D0*t27*t71
     #+0.6347552015792813D-1*t36*t75
      t79 = t60*t78
      t83 = 1.D0+0.1325688999052018D0*t45
      t89 = 0.598255043577108D1*t49+0.2225569421150687D1*t45
     #+0.8004286349993634D0*t52+0.1897004325747559D0*t54
      t92 = 1.D0+0.160818243221511D2/t89
      t93 = dlog(t92)
      t98 = -0.62182D-1*rho*t83*t93+0.3109D-1*t48*t60
      t100 = 1.D0+0.9524406311809197D-2*t8
      t101 = 1/t100
      t104 = t100**2
      t105 = 1/t104
      t109 = 1/t104/t100
      t112 = t104**2
      t113 = 1/t112
      t116 = 0.589076D0+0.4213349717780282D-1*t8*t101
     #-0.1743692431454588D-2*t19*t105+0.367822944D-4*t27*t109
     #-0.3456644474956098D-6*t36*t113
      zk(i) = -0.7385587663820224D0*t3*t41-0.3109D-1*t48*t79+t98*t116
      t120 = t4*rho
      t122 = 1/t5/t120
      t123 = sigma*t122
      t126 = t15*t4
      t128 = 1/t2/t126
      t129 = t14*t128
      t133 = 1/t25/rho
      t134 = t24*t133
      t139 = 1/t5/t25/t120
      t140 = t32*t139
      t143 = t32*sigma
      t146 = 1/t2/t25/t126
      t147 = t143*t146
      t149 = 1/t37/t10
      t152 = 0.1755331997602442D-1*t123*t11-0.1584654969337721D-2*t129
     #*t21+0.2012265130666667D-4*t134*t29-0.1477359830566614D-6*t140
     #*t38+0.503909330668646D-9*t147*t149
      t155 = t47*t60
      t158 = 1/t54
      t159 = t44*t158
      t162 = t56**2
      t163 = 1/t162
      t164 = t48*t163
      t165 = t49**2
      t166 = t165**2
      t167 = t166*t49
      t168 = 1/t167
      t169 = 1/t4
      t170 = t168*t169
      t172 = t158*t169
      t174 = 1/t52
      t175 = t174*t169
      t177 = 1/t45
      t178 = t177*t169
      t180 = -0.4160732915599108D1*t170-0.3229384477375903D1*t172
     #-0.2326004811900819D1*t175-0.5092110158254621D0*t178
      t181 = 1/t59
      t182 = t180*t181
      t183 = t182*t78
      t195 = 1/t74/t62
      t198 = 0.4068690918181788D1*t123*t63-0.7331189631570478D1*t129
     #*t67+0.661717504D1*t134*t71-0.2846227052321425D1*t140*t75
     #+0.4299140585511414D0*t147*t195
      t199 = t60*t198
      t204 = t172*t93
      t206 = t89**2
      t207 = 1/t206
      t208 = t83*t207
      t213 = -0.99709173929518D0*t170-0.7418564737168958D0*t172
     #-0.4002143174996817D0*t175-0.1264669550498372D0*t178
      t214 = 1/t92
      t216 = t208*t213*t214
      t224 = t163*t180*t181
      t227 = -0.62182D-1*t83*t93+rho*(0.2747799777968419D-2*t204+0.1D1
     #*t216)+0.3109D-1*t155-0.1664368495390566D-2*t159*t60-0.5D0*t48
     #*t224
      t230 = 0.1123559924741408D0*t123*t101
      t232 = 0.1036981709164809D-1*t129*t105
      t234 = 0.3828324096D-3*t134*t109
      t235 = t140*t113
      t238 = 1/t112/t100
      t240 = 0.351173188852824D-7*t147*t238
      t241 = -t230+t232-t234+0.6489723575522644D-5*t235-t240
      vrhoa(i) = -0.9847450218426965D0*t2*t41-0.3692793831910112D0*t3
     #*t152-0.3109D-1*t155*t78+0.1664368495390566D-2*t159*t79+0.5D0
     #*t164*t183-0.15545D-1*t48*t199+t227*t116+t98*t241
      t245 = sigma*t18
      t248 = t14*t26
      t251 = t24*t35
      t256 = 1/t2/t25/t16
      t257 = t32*t256
      t260 = -0.1316498998201831D-1*t7*t11+0.1188491227003291D-2*t245
     #*t21-0.1509198848D-4*t248*t29+0.110801987292496D-6*t251*t38
     #-0.3779319980014845D-9*t257*t149
      t273 = -0.3051518188636341D1*t7*t63+0.5498392223677858D1*t245
     #*t67-0.496288128D1*t248*t71+0.2134670289241069D1*t251*t75
     #-0.3224355439133561D0*t257*t195
      t274 = t60*t273
      t287 = 0.8426699435560563D-1*t7*t101-0.7777362818736071D-2*t245
     #*t105+0.2871243072D-3*t248*t109-0.4867292681641983D-5*t251*t113
     #+0.263379891639618D-7*t257*t238
      vsigmaaa(i) = -0.7385587663820224D0*t3*t260-0.3109D-1*t48*t274
     #+2.D0*t98*t287
      t294 = sigma/t5/t15
      t297 = t15*t120
      t300 = t14/t2/t297
      t301 = t300*t105
      t304 = t24/t33
      t305 = t304*t109
      t310 = t32/t5/t25/t15
      t311 = t310*t113
      t316 = t143/t2/t25/t297
      t317 = t316*t238
      t320 = t25**2
      t323 = t32*t14/t320/t4
      t325 = 1/t112/t104
      t327 = 0.4459621515264D-8*t323*t325
      t330 = t47*t163
      t346 = 1/t37/t20
      t357 = t162**2
      t358 = 1/t357
      t360 = t180**2
      t361 = t59**2
      t362 = 1/t361
      t381 = 1/t74/t66
      t390 = 1/t15
      t391 = 1/t167/t44*t390
      t393 = 1/t120
      t394 = t168*t393
      t397 = 1/t54/t44
      t398 = t397*t390
      t400 = t158*t393
      t404 = 1/t52/t44*t390
      t406 = t174*t393
      t410 = 1/t45/t44*t390
      t412 = t177*t393
      t414 = -0.6934554859331846D1*t391+0.1664293166239643D2*t394
     #-0.4305845969834537D1*t398+0.1291753790950361D2*t400
     #-0.2326004811900819D1*t404+0.9304019247603276D1*t406
     #-0.3394740105503081D0*t410+0.2036844063301848D1*t412
      t419 = t393*t397
      t426 = 1/t162/t56
      t452 = 0.5495599555936838D-2*t204
      t453 = 0.2D1*t216
      t457 = t213**2
      t460 = 0.2D1*t83/t206/t89*t457*t214
      t472 = 0.1D1*t208*(-0.8309097827459833D0*t391
     #+0.199418347859036D1*t394-0.4945709824779306D0*t398
     #+0.1483712947433792D1*t400-0.2001071587498409D0*t404
     #+0.8004286349993634D0*t406-0.4215565168327908D-1*t410
     #+0.2529339100996745D0*t412)*t214
      t473 = t206**2
      t476 = t92**2
      t480 = 0.160818243221511D2*t83/t473*t457/t476
      t482 = 0.1831866518645613D-2*t398*t93
      t484 = 0.5495599555936838D-2*t400*t93
      t495 = dlog(1.D0+0.2960857464321668D2/(0.8157414703487641D1*t49
     #+0.2247591863577616D1*t45+0.4300972471276643D0*t52
     #+0.1911512595127338D0*t54))
      t497 = (1.D0+0.6901399211255825D-1*t45)*t495*t169
      t502 = 0.8837926660346786D-1*t172*t207*t213*t214
      s1 = 3.D0*t227*t241+t98*(0.8239439448103662D0*t294*t101
     #-0.1065518352468551D0*t301+0.537596439552D-2*t305
     #-0.1286791061715555D-3*t311+0.1291426557904201D-5*t317-t327)
     #+0.2D1*t330*t183-0.6218D-1*t155*t198-0.3692793831910112D0*t3*(
     #-0.1287243464908457D0*t294*t11+0.2066673166164795D-1*t300*t21
     #-0.4695349967075556D-3*t304*t29+0.5491506880554649D-5*t310*t38
     #-0.3445745456547828D-7*t316*t149+0.8532332817066667D-10*t323
     #*t346)-0.6564966812284644D0/t5*t41-0.1969490043685393D1*t2*t152
     #+0.1608234158893535D2*t48*t358*t360*t362*t78+0.1D1*t164*t182*t198
      s2 = s1-0.15545D-1*t48*t60*(-0.2983706673333311D2*t294*t63
     #+0.9975095585979163D2*t300*t67-0.1439358987377778D3*t304*t71
     #+0.1000251185367385D3*t310*t75-0.31601445964313D2*t316*t195
     #+0.36397121536D1*t323*t381)+0.5D0*t164*t414*t181*t78
     #+0.2219157993854088D-2*t419*t79-0.1070677706909338D0*t159*t163
     #*t183
      v2rhoa2(i) = s2-0.1D1*t48*t426*t360*t181*t78
     #+0.3328736990781132D-2*t159*t199+(-0.2D1*t330*t182+0.1D1*t48
     #*t426*t360*t181-0.2219157993854088D-2*t419*t60
     #+0.1070677706909338D0*t159*t224-0.1608234158893535D2*t48*t358
     #*t360*t362-0.5D0*t48*t163*t414*t181+t452+t453+rho*(-t460+t472
     #+t480+t482-t484+0.3377399869956914D-1*t497-t502))*t116+(t452
     #+t453+rho*(-t460+t472+t480+t482-t484-0.3377399869956914D-1*t497
     #-t502))*t116+t227*(-t230+t232-t234+0.6489723575522644D-5*t235
     #-t240)+t98*(-0.3050650990810237D-1*t301+0.256852672512D-2*t305
     #-0.8108779995105606D-4*t311+0.1033899552745464D-5*t317-t327)
      t525 = sigma*t128
      t528 = t14*t133
      t531 = t24*t139
      t534 = t32*t146
      t539 = t143/t320/rho
      t572 = t525*t105
      t574 = t528*t109
      t576 = t531*t113
      t578 = t534*t238
      t581 = 0.3344716136448D-8*t539*t325
      s1 = -0.9847450218426965D0*t2*t260-0.3692793831910112D0*t3*
     #(0.7021327990409767D-1*t122*t11-0.1312306629222938D-1*t525*t21
     #+0.3219672705706667D-3*t528*t29-0.3897026185830995D-5*t531*t38
     #+0.2508722692810574D-7*t534*t149-0.63992496128D-10*t539*t346)
     #-0.3109D-1*t155*t273+0.1664368495390566D-2*t159*t274
      v2rhoasigmaaa(i) = s1+0.5D0*t164*t182*t273-0.15545D-1*t48*t60*
     #(0.1627476367272715D2*t122*t63-0.6381643244748801D2*t525*t67
     #+0.9802616149333333D2*t528*t71-0.7074949832407171D2*t531*t75
     #+0.2305621338540804D2*t534*t195-0.27297841152D1*t539*t381)+2.D0
     #*t227*t287+t98*(-0.4494239698965634D0*t122*t101
     #+0.6435915079766915D-1*t572-0.345772468224D-2*t574
     #+0.8677474426538262D-4*t576-0.9158939401002274D-6*t578+t581)+t98
     #*(0.2287988243107678D-1*t572-0.192639504384D-2*t574
     #+0.6081584996329205D-4*t576-0.7754246645590978D-6*t578+t581)
      t592 = sigma*t26
      t595 = t14*t35
      t598 = t24*t256
      t602 = t32/t320
      v2sigmaaa2(i) = -0.7385587663820224D0*t3*(0.508833481115887D-2
     #*t18*t21-0.181107499008D-3*t592*t29+0.2479561690203262D-5*t595
     #*t38-0.1730369220407337D-7*t598*t149+0.47994372096D-10*t602*t346
     #)-0.3109D-1*t48*t60*(0.2586875544090457D2*t18*t67-0.53668096D2
     #*t592*t71+0.4452344258608951D2*t595*t75-0.160024178634026D2*t598
     #*t195+0.20473380864D1*t602*t381)+4.D0*t98*(
     #-0.1715991182330758D-1*t18*t105+0.144479628288D-2*t592*t109
     #-0.4561188747246904D-4*t595*t113+0.5815684984193233D-6*t598*t238
     #-0.2508537102336D-8*t602*t325)
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      v2rhoa2(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      v2rhoasigmaaa(i) = 0.0d0
      v2sigmaaa2(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end