!/****************************************************************/
!/**                       DISLIN.F90                           **/
!/**                                                            **/
!/** Empty implementation of DISLIN Fortran 90                  **/
!/**                            for avoiding source processing. **/
!/**                                                            **/
!/** Date     :  05.05.2022                                     **/
!/** Routines :  795                                            **/
!/** Version  :  11.1 / explicit-shape / double precision       **/
!/****************************************************************/

  subroutine abs3pt(x,y,z,xp,yp)
    implicit none
    double precision, intent (in)  :: x,y,z
    double precision, intent (out) :: xp,yp
    call doesnotexist("abs3pt")
  end subroutine abs3pt

  subroutine addlab(cstr,v,itic,cax)
    implicit none
    character (len = *), intent (in) :: cstr,cax
    double precision, intent (in) :: v
    integer, intent (in) :: itic
    call doesnotexist("addlab")
  end subroutine addlab

  subroutine angle(i)
    implicit none
    integer, intent (in) :: i
    call doesnotexist("angle")
  end subroutine angle

  subroutine arcell(nx,ny,na,nb,alpha,beta,theta)
    implicit none
    integer, intent (in) :: nx,ny,na,nb
    double precision, intent (in)   :: alpha,beta,theta
    call doesnotexist("arcell")
  end subroutine arcell

  subroutine areaf(ix,iy,n)
    implicit none
    integer, intent (in) :: n
    integer, dimension (n), intent (in) :: ix,iy
    call doesnotexist("areaf")
  end subroutine areaf

  subroutine autres(i,j)
    implicit none
    integer, intent (in) :: i,j
    call doesnotexist("autres")
  end subroutine autres

  subroutine ax2grf()
    call doesnotexist("ax2grf")
  end subroutine ax2grf

  subroutine ax3len(i,j,k)
    implicit none
    integer, intent (in) :: i,j,k
    call doesnotexist("ax3len")
  end subroutine ax3len

  subroutine axclrs(n,copt,cax)
    implicit none
    integer, intent(in) :: n
    character (len = *) , intent (in) :: copt, cax
    call doesnotexist("axclrs")
  end subroutine axclrs

  subroutine axends(copt,cax)
    implicit none
    character (len = *) , intent (in) :: copt, cax
    call doesnotexist("axends")
  end subroutine axends

  subroutine axgit()
    call doesnotexist("axgit")
  end subroutine axgit

  subroutine axis3d(x,y,z)
    implicit none
    double precision, intent (in) :: x,y,z
    call doesnotexist("axis3d")
  end subroutine axis3d

  subroutine axsbgd(n)
    implicit none
    integer, intent (in) :: n
    call doesnotexist("axsbgd")
  end subroutine axsbgd

  subroutine axsers()
    call doesnotexist("axsers")
  end subroutine axsers

  subroutine axslen(i,j)
    implicit none
    integer, intent (in) :: i,j
    call doesnotexist("axslen")
  end subroutine axslen

  subroutine axsorg(i,j)
    implicit none
    integer, intent (in) :: i,j
    call doesnotexist("axsorg")
  end subroutine axsorg

  subroutine axspos(i,j)
    implicit none
    integer, intent (in) :: i,j
    call doesnotexist("axspos")
  end subroutine axspos

  subroutine axsscl(copt,cax)
    implicit none
    character (len = *), intent (in) :: copt,cax
    call doesnotexist("axsscl")
  end subroutine axsscl

  subroutine axstyp(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("axstyp")
  end subroutine axstyp

  subroutine barbor(iclr)
    implicit none
    integer, intent (in) :: iclr
    call doesnotexist("barbor")
  end subroutine barbor

  subroutine barclr(ic1,ic2,ic3)
    implicit none
    integer, intent (in) :: ic1,ic2,ic3
    call doesnotexist("barclr")
  end subroutine barclr

  subroutine bargrp(n,xgap)
    implicit none
    integer, intent (in) :: n
    double precision, intent (in) :: xgap
    call doesnotexist("bargrp")
  end subroutine bargrp

  subroutine barmod(cmode,copt)
    implicit none
    character (len = *), intent (in) :: cmode,copt
    call doesnotexist("barmod")
  end subroutine barmod

  subroutine baropt(x1,x2)
    implicit none
    double precision, intent (in) :: x1,x2
    call doesnotexist("baropt")
  end subroutine baropt

  subroutine barpos(cpos)
    implicit none
    character (len = *), intent (in) :: cpos
    call doesnotexist("barpos")
  end subroutine barpos

  subroutine bars(xray,y1ray,y2ray,n)
    implicit none
    integer, intent (in) :: n
    double precision, intent (in out), dimension (n) :: xray,y1ray,y2ray
    call doesnotexist("bars")
  end subroutine bars

  subroutine bars3d(xray,yray,z1ray,z2ray,xwray,ywray,icray,n)
    implicit none
    integer, intent (in) :: n
    double precision, intent (in), dimension (n) :: xray,yray,z1ray,z2ray,xwray,ywray
    integer, intent (in), dimension (n) :: icray
    call doesnotexist("bars3d")
  end subroutine bars3d

  subroutine bartyp(ctyp)
    implicit none
    character (len = *), intent (in) :: ctyp
    call doesnotexist("bartyp")
  end subroutine bartyp

  subroutine barwth(fact)
    implicit none
    double precision, intent (in) :: fact
    call doesnotexist("barwth")
  end subroutine barwth

  subroutine basalf(calph)
    implicit none
    character (len = *), intent (in) :: calph
    call doesnotexist("basalf")
  end subroutine basalf

  subroutine basdat(id,im,iy)
    implicit none
    integer, intent (in) :: id, im, iy
    call doesnotexist("basdat")
  end subroutine basdat

  subroutine bezier(xray,yray,nray,x,y,n)
    implicit none
    integer, intent (in)  :: nray,n
    double precision, dimension (nray), intent (in)  :: xray, yray
    double precision, dimension (n), intent (out) :: x, y
    call doesnotexist("bezier")
  end subroutine bezier

  subroutine bfcclr(ic)
    implicit none
    integer, intent (in) :: ic
    call doesnotexist("bfcclr")
  end subroutine bfcclr

  subroutine bfcmsh(ic)
    implicit none
    integer, intent (in) :: ic
    call doesnotexist("bfcmsh")
  end subroutine bfcmsh

  subroutine bitsi2(nbits,mher,iher,mhin,ihin,lob)
    implicit none
    integer, intent (in) :: nbits,iher,ihin,lob
    integer (kind=selected_int_kind(4)), intent (in) :: mher
    integer (kind=selected_int_kind(4)), intent (in out) :: mhin
    call doesnotexist("bitsi2")
  end subroutine bitsi2

  subroutine bitsi4(nbits,mher,iher,mhin,ihin,lob)
    implicit none
    integer, intent (in) :: nbits,mher,iher,ihin,lob
    integer, intent (in out) :: mhin
    call doesnotexist("bitsi4")
  end subroutine bitsi4

  subroutine bmpfnt(cfnt)
    implicit none
    character (len = *), intent (in) :: cfnt
    call doesnotexist("bmpfnt")
  end subroutine bmpfnt

  subroutine bmpmod(n,cval,copt)
    implicit none
    character (len = *), intent (in) :: cval,copt
    integer, intent (in) :: n
    call doesnotexist("bmpmod")
  end subroutine bmpmod

  subroutine box2d()
    call doesnotexist("box2d")
  end subroutine box2d

  subroutine box3d()
    call doesnotexist("box3d")
  end subroutine box3d

  subroutine bufmod(cmod,ckey)
    implicit none
    character (len = *), intent (in) :: cmod,ckey
    call doesnotexist("bufmod")
  end subroutine bufmod

  subroutine center()
    call doesnotexist("center")
  end subroutine center

  subroutine cgmbgd(xr,xg,xb)
    implicit none
    double precision, intent (in) :: xr,xg,xb
    call doesnotexist("cgmbgd")
  end subroutine cgmbgd

  subroutine cgmpic(ct)
    implicit none
    character (len = *), intent (in) :: ct
    call doesnotexist("cgmpic")
  end subroutine cgmpic

  subroutine cgmver(n)
    implicit none
    integer, intent (in) :: n
    call doesnotexist("cgmver")
  end subroutine cgmver

  subroutine chaang(x)
    implicit none
    double precision, intent (in) :: x
    call doesnotexist("chaang")
  end subroutine chaang

  subroutine chacod(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("chacod")
  end subroutine chacod

  subroutine chaspc(x)
    implicit none
    double precision, intent (in) :: x
    call doesnotexist("chaspc")
  end subroutine chaspc

  subroutine chawth(x)
    implicit none
    double precision, intent (in) :: x
    call doesnotexist("chawth")
  end subroutine chawth

  subroutine chnatt()
    call doesnotexist("chnatt")
  end subroutine chnatt

  subroutine chncrv(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("chncrv")
  end subroutine chncrv

  subroutine chndot()
    call doesnotexist("chndot")
  end subroutine chndot

  subroutine chndsh()
    call doesnotexist("chndsh")
  end subroutine chndsh

  subroutine chnbar(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("chnbar")
  end subroutine chnbar

  subroutine chnpie(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("chnpie")
  end subroutine chnpie

  subroutine circ3p(x1,y1,x2,y2,x3,y3,xm,ym,r)
    implicit none
    double precision, intent (in) :: x1,y1,x2,y2,x3,y3
    double precision, intent (out) :: xm,ym,r
    call doesnotexist("circ3p")
  end subroutine circ3p

  subroutine circle(nx,ny,nr)
    implicit none
    integer, intent (in) :: nx,ny,nr
    call doesnotexist("circle")
  end subroutine circle

  subroutine circsp(n)
    implicit none
    integer, intent (in) :: n
    call doesnotexist("circsp")
  end subroutine circsp

  subroutine clip3d(ctyp)
    implicit none
    character (len = *), intent (in) :: ctyp
    call doesnotexist("clip3d")
  end subroutine clip3d

  subroutine closfl(nlu)
    implicit none
    integer, intent (in) :: nlu
    call doesnotexist("closfl")
  end subroutine closfl

  subroutine clpbor(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("clpbor")
  end subroutine clpbor

  subroutine clpmod(cmod)
    implicit none
    character (len = *), intent (in) :: cmod
    call doesnotexist("clpmod")
  end subroutine clpmod

  subroutine clpwin(nx,ny,nw,nh)
    implicit none
    integer, intent (in) :: nx,ny,nw,nh
    call doesnotexist("clpwin")
  end subroutine clpwin

  subroutine clrcyc(i,iclr)
    implicit none
    integer, intent (in) :: i,iclr
    call doesnotexist("clrcyc")
  end subroutine clrcyc

  subroutine clrmod(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("clrmod")
  end subroutine clrmod

  subroutine clswin(id)
    implicit none
    integer, intent (in) :: id
    call doesnotexist("clswin")
  end subroutine clswin

  subroutine color(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("color")
  end subroutine color

  subroutine colran(i,j)
    implicit none
    integer, intent (in) :: i,j
    call doesnotexist("colran")
  end subroutine colran

  subroutine colray(z,ncol,n)
    implicit none
    integer, intent (in) :: n
    double precision, dimension (n), intent (in) :: z
    integer, dimension (n), intent (out) :: ncol
    call doesnotexist("colray")
  end subroutine colray

  subroutine complx()
    call doesnotexist("complx")
  end subroutine complx

  subroutine conclr(iray,n)
    implicit none
    integer, intent (in) :: n
    integer, dimension (n), intent (in) :: iray
    call doesnotexist("conclr")
  end subroutine conclr

  subroutine concrv(x,y,n,zlev)
    implicit none
    integer, intent (in) :: n
    double precision, dimension (n), intent (in) :: x,y
    double precision, intent (in) :: zlev
    call doesnotexist("concrv")
  end subroutine concrv

  subroutine cone3d(x,y,z,r,h1,h2,nsk1,nsk2)
    implicit none
    double precision, intent (in) :: x,y,z,r,h1,h2
    integer, intent (in) :: nsk1,nsk2
    call doesnotexist("cone3d")
  end subroutine cone3d

  subroutine confll(xray,yray,zray,n,i1ray,i2ray,i3ray,ntri,zlev,nlev)
    implicit none
    integer, intent (in) :: n,ntri,nlev
    double precision, dimension (n), intent (in) :: xray,yray,zray
    double precision, dimension (nlev), intent (in) :: zlev
    integer, dimension (ntri), intent (in) :: i1ray,i2ray,i3ray
    call doesnotexist("confll")
  end subroutine confll

  subroutine congap(xgap)
    implicit none
    double precision, intent (in) :: xgap
    call doesnotexist("congap")
  end subroutine congap

  subroutine conlab(cstr)
    implicit none
    character (len = *), intent (in) :: cstr
    call doesnotexist("conlab")
  end subroutine conlab

  subroutine conmat(zmat,n,m,zlev)
    implicit none
    integer, intent (in) :: n,m
    double precision, dimension (n,m), intent (in) :: zmat
    double precision, intent (in) :: zlev
    call doesnotexist("conmat")
  end subroutine conmat

  subroutine conmod (xf1,xf2)
    implicit none
    double precision, intent (in) :: xf1,xf2
    call doesnotexist("conmod")
  end subroutine conmod

  subroutine conn3d(x2,y2,z2)
    implicit none
    double precision, intent (in) :: x2,y2,z2
    call doesnotexist("conn3d")
  end subroutine conn3d

  subroutine connpt(x,y)
    implicit none
    double precision, intent (in) :: x,y
    call doesnotexist("connpt")
  end subroutine connpt

  subroutine conpts(x,n,y,m,z,zlev,xpts,ypts,maxpts,nray,maxray,nlins)
    implicit none
    integer, intent (in) :: n,m,maxpts,maxray
    double precision, dimension (n), intent (in) :: x
    double precision, dimension (m), intent (in) :: y
    double precision, dimension (n,m), intent (in) :: z
    double precision, intent (in) :: zlev
    integer, dimension (maxray), intent (out) :: nray
    integer, intent (out) :: nlins
    double precision, dimension (maxpts), intent (out) :: xpts,ypts
    call doesnotexist("conpts")
  end subroutine conpts

  subroutine conshd(xray,n,yray,m,zmat,zlev,nlray)
    implicit none
    integer, intent (in) :: n,m,nlray
    double precision, dimension (n), intent (in) :: xray
    double precision, dimension (m), intent (in) :: yray
    double precision, dimension (nlray), intent (in) :: zlev
    double precision, dimension (n,m), intent (in) :: zmat
    call doesnotexist("conshd")
  end subroutine conshd

  subroutine conshd2(xmat,ymat,zmat,n,m,zlev,nlray)
    implicit none
    integer, intent (in) :: n,m,nlray
    double precision, dimension (n,m), intent (in) :: xmat,ymat,zmat
    double precision, dimension (nlray), intent (in) :: zlev
    call doesnotexist("conshd2")
  end subroutine conshd2

  subroutine conshd3d(xray,n,yray,m,zmat,zlev,nlray)
    implicit none
    integer, intent (in) :: n,m,nlray
    double precision, dimension (n), intent (in) :: xray
    double precision, dimension (m), intent (in) :: yray
    double precision, dimension (nlray), intent (in) :: zlev
    double precision, dimension (n,m), intent (in) :: zmat
    call doesnotexist("conshd3d")
  end subroutine conshd3d

  subroutine contri(xray,yray,zray,n,i1ray,i2ray,i3ray,ntri,zlev)
    implicit none
    integer, intent (in) :: n,ntri
    double precision, dimension (n), intent (in) :: xray,yray,zray
    integer, dimension (ntri), intent (in) :: i1ray,i2ray,i3ray
    double precision, intent (in) :: zlev
    call doesnotexist("contri")
  end subroutine contri

  subroutine contur(x,n,y,m,z,zlev)
    implicit none
    integer, intent (in) :: n,m
    double precision, dimension (n), intent (in) :: x
    double precision, dimension (m), intent (in) :: y
    double precision, dimension (n,m), intent (in) :: z
    double precision, intent (in) :: zlev
    call doesnotexist("contur")
  end subroutine contur

  subroutine contur2(x,y,z,n,m,zlev)
    implicit none
    integer, intent (in) :: n,m
    double precision, dimension (n,m), intent (in) :: x,y,z
    double precision, intent (in) :: zlev
    call doesnotexist("contur2")
  end subroutine contur2

  subroutine cross()
    call doesnotexist("cross")
  end subroutine cross

  subroutine crvmat(zmat,ixdim,iydim,ixpts,iypts)
    implicit none
    integer, intent (in) :: ixdim,iydim,ixpts,iypts
    double precision, dimension (ixdim,iydim), intent (in) :: zmat
    call doesnotexist("crvmat")
  end subroutine crvmat

  subroutine crvqdr(xray,yray,zray,n)
    implicit none
    integer, intent (in) :: n
    double precision, dimension (n), intent (in) :: xray,yray,zray
    call doesnotexist("crvqdr")
  end subroutine crvqdr

  subroutine crvt3d(x,y,z,r,ic,n)
    implicit none
    integer, intent (in) :: n
    double precision, dimension (n), intent (in) :: x,y,z,r
    integer, dimension (n), intent (in) :: ic
    call doesnotexist("crvt3d")
  end subroutine crvt3d

  subroutine crvtri(xray,yray,zray,n,i1ray,i2ray,i3ray,ntri)
    implicit none
    integer, intent (in) :: n,ntri
    double precision, dimension (n), intent (in) :: xray,yray,zray
    integer, dimension (ntri), intent (in) :: i1ray,i2ray,i3ray
    call doesnotexist("crvtri")
  end subroutine crvtri

  subroutine curv3d(x,y,z,n)
    implicit none
    integer, intent (in) :: n
    double precision, dimension (n), intent (in) :: x,y,z
    call doesnotexist("curv3d")
  end subroutine curv3d

  subroutine curv4d(x,y,z,w,n)
    implicit none
    integer, intent (in) :: n
    double precision, dimension (n), intent (in) :: x,y,z,w
    call doesnotexist("curv4d")
  end subroutine curv4d

  subroutine csrkey(ik)
    implicit none
    integer, intent (out) :: ik
    call doesnotexist("csrkey")
  end subroutine csrkey

  subroutine csrlin(ix1,iy1,ix2,iy2)
    implicit none
    integer, intent (out) :: ix1,iy1,ix2,iy2
    call doesnotexist("csrlin")
  end subroutine csrlin

  subroutine csrmod(cmod,ckey)
    implicit none
    character (len = *), intent (in) :: cmod,ckey
    call doesnotexist("csrmod")
  end subroutine csrmod

  subroutine csrpol(ixray,iyray,nmax,n,iret)
    implicit none
    integer, intent (in) :: nmax
    integer, dimension (nmax), intent (out) :: ixray,iyray
    integer, intent (out) :: n, iret
    call doesnotexist("csrpol")
  end subroutine csrpol

  subroutine csrpos(ix,iy,ik)
    implicit none
    integer, intent (in out) :: ix,iy
    integer, intent (out) :: ik
    call doesnotexist("csrpos")
  end subroutine csrpos

  subroutine csrpt1(ix,iy)
    implicit none
    integer, intent (out) :: ix,iy
    call doesnotexist("csrpt1")
  end subroutine csrpt1

  subroutine csrmov(ixray,iyray,nmax,n,iret)
    implicit none
    integer, intent (in) :: nmax
    integer, dimension (nmax), intent (out) :: ixray,iyray
    integer, intent (out) :: n, iret
    call doesnotexist("csrmov")
  end subroutine csrmov

  subroutine csrpts(ixray,iyray,nmax,n,iret)
    implicit none
    integer, intent (in) :: nmax
    integer, dimension (nmax), intent (out) :: ixray,iyray
    integer, intent (out) :: n, iret
    call doesnotexist("csrpts")
  end subroutine csrpts

  subroutine csrrec(ix1,iy1,ix2,iy2)
    implicit none
    integer, intent (out) :: ix1,iy1,ix2,iy2
    call doesnotexist("csrrec")
  end subroutine csrrec

  subroutine csrtyp(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("csrtyp")
  end subroutine csrtyp

  subroutine csruni(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("csruni")
  end subroutine csruni

  subroutine curve(x,y,n)
    implicit none
    integer, intent (in) :: n
    double precision, dimension (n), intent (in) :: x,y
    call doesnotexist("curve")
  end subroutine curve

  subroutine curve3(x,y,z,n)
    implicit none
    integer, intent (in) :: n
    double precision, dimension (n), intent (in) :: x,y,z
    call doesnotexist("curve3")
  end subroutine curve3

  subroutine curvmp(x,y,n)
    implicit none
    integer, intent (in) :: n
    double precision, dimension (n), intent (in) :: x,y
    call doesnotexist("curvmp")
  end subroutine curvmp

  subroutine curvx3(x,y,z,n)
    implicit none
    integer, intent (in) :: n
    double precision, dimension (n), intent (in) :: x,z
    double precision, intent (in) :: y
    call doesnotexist("curvx3")
  end subroutine curvx3

  subroutine curvy3(x,y,z,n)
    implicit none
    integer, intent (in) :: n
    double precision, dimension (n), intent (in) :: y,z
    double precision, intent (in) :: x
    call doesnotexist("curvy3")
  end subroutine curvy3

  subroutine cyli3d(x,y,z,r,h,nsk1,nsk2)
    implicit none
    double precision, intent (in) :: x,y,z,r,h
    integer, intent (in) :: nsk1,nsk2
    call doesnotexist("cyli3d")
  end subroutine cyli3d

  subroutine dash()
    call doesnotexist("dash")
  end subroutine dash

  subroutine dashl()
    call doesnotexist("dashl")
  end subroutine dashl

  subroutine dashm()
    call doesnotexist("dashm")
  end subroutine dashm

  subroutine dattim(cdat,ctim)
    implicit none
    character (len = *), intent (out) :: cdat,ctim
    call doesnotexist("dattim")
  end subroutine dattim

  subroutine dbffin()
    call doesnotexist("dbffin")
  end subroutine dbffin

  subroutine dbfini(iret)
    implicit none
    integer, intent (out) :: iret
    call doesnotexist("dbfini")
  end subroutine dbfini

  subroutine dbfmod(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("dbfmod")
  end subroutine dbfmod

  subroutine delglb()
    implicit none
    call doesnotexist("delglb")
  end subroutine delglb

  subroutine digits(i,cax)
    implicit none
    integer, intent (in) :: i
    character (len = *), intent (in) :: cax
    call doesnotexist("digits")
  end subroutine digits

  subroutine disalf()
    call doesnotexist("disalf")
  end subroutine disalf

  subroutine disenv(cstr)
    implicit none
    character (len = *), intent (in) :: cstr
    call doesnotexist("disenv")
  end subroutine disenv

  subroutine disfin()
    call doesnotexist("disfin")
  end subroutine disfin

  subroutine disini()
    call doesnotexist("disini")
  end subroutine disini

  subroutine disk3d(x,y,z,r1,r2,nsk1,nsk2)
    implicit none
    double precision, intent (in) :: x,y,z,r1,r2
    integer, intent (in) :: nsk1,nsk2
    call doesnotexist("disk3d")
  end subroutine disk3d

  subroutine doevnt()
    call doesnotexist("doevnt")
  end subroutine doevnt

  subroutine dot()
    call doesnotexist("dot")
  end subroutine dot

  subroutine dotl()
    call doesnotexist("dotl")
  end subroutine dotl

  subroutine duplx()
    call doesnotexist("duplx")
  end subroutine duplx

  subroutine dwgbut(cstr,ival)
    implicit none
    character (len=*), intent (in) :: cstr
    integer, intent (in out) :: ival
    call doesnotexist("dwgbut")
  end subroutine dwgbut

  subroutine dwgerr(ival)
    implicit none
    integer, intent (in out) :: ival
    call doesnotexist("dwgerr")
  end subroutine dwgerr

  subroutine dwgfil(clab,cstr,cmask)
    implicit none
    character (len=*), intent (in) :: clab,cmask
    character (len=*), intent (in out) :: cstr
    call doesnotexist("dwgfil")
  end subroutine dwgfil

  subroutine dwglis(clab,clis,ilis)
    implicit none
    character (len=*), intent(in) :: clab,clis
    integer, intent (in out) :: ilis
    call doesnotexist("dwglis")
  end subroutine dwglis

  subroutine dwgmsg(cstr)
    implicit none
    character (len=*), intent (in) :: cstr
    call doesnotexist("dwgmsg")
  end subroutine dwgmsg

  subroutine dwgtxt(clab,cstr)
    implicit none
    character (len=*), intent(in) :: clab
    character (len=*), intent(in out) :: cstr
    call doesnotexist("dwgtxt")
  end subroutine dwgtxt

  subroutine ellips(nx,ny,na,nb)
    implicit none
    integer, intent (in) :: nx,ny,na,nb
    call doesnotexist("ellips")
  end subroutine ellips

  subroutine endgrf()
    call doesnotexist("endgrf")
  end subroutine endgrf

  subroutine erase()
    call doesnotexist("erase")
  end subroutine erase

  subroutine errbar(x,y,err1,err2,n)
    implicit none
    integer, intent (in) :: n
    double precision, dimension (n), intent (in) :: x,y,err1,err2
    call doesnotexist("errbar")
  end subroutine errbar

  subroutine errdev(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("errdev")
  end subroutine errdev

  subroutine errfil(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("errfil")
  end subroutine errfil

  subroutine errmod(cstr,cmode)
    implicit none
    character (len = *), intent (in) :: cstr,cmode
    call doesnotexist("errmod")
  end subroutine errmod

  subroutine eushft(calph,csft)
    implicit none
    character (len = *), intent (in) :: calph,csft
    call doesnotexist("eushft")
  end subroutine eushft

  subroutine expimg(cfl,copt)
    implicit none
    character (len = *), intent (in) :: cfl,copt
    call doesnotexist("expimg")
  end subroutine expimg

  subroutine expzlb(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("expzlb")
  end subroutine expzlb

  subroutine fbars(x,y1,y2,y3,y4,n)
    implicit none
    integer, intent (in) :: n
    double precision, dimension (n), intent (in) :: x,y1,y2,y3,y4
    call doesnotexist("fbars")
  end subroutine fbars

  subroutine fcha(x,ndez,nl,cstr)
    implicit none
    double precision, intent (in) :: x
    integer, intent (in) :: ndez
    integer, intent (out) :: nl
    character (len = *), intent (out) :: cstr
    call doesnotexist("fcha")
  end subroutine fcha

  subroutine field(xray,yray,uray,vray,n,ivec)
    implicit none
    integer, intent (in) :: n,ivec
    double precision, dimension (n), intent (in) :: xray,yray,uray,vray
    call doesnotexist("field")
  end subroutine field

  subroutine field3d(x1ray,y1ray,z1ray,x2ray,y2ray,z2ray,n,ivec)
    implicit none
    integer, intent (in) :: n,ivec
    double precision, dimension (n), intent (in) :: x1ray,y1ray,z1ray, &
           x2ray,y2ray,z2ray
    call doesnotexist("field3d")
  end subroutine field3d

  subroutine filbox(nx,ny,nw,nh)
    implicit none
    integer, intent (in) :: nx,ny,nw,nh
    call doesnotexist("filbox")
  end subroutine filbox

  subroutine filclr(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("filclr")
  end subroutine filclr

  subroutine filmod(cmod)
    implicit none
    character (len = *), intent (in) :: cmod
    call doesnotexist("filmod")
  end subroutine filmod

  subroutine filopt(copt,ckey)
    implicit none
    character (len = *), intent (in) :: copt,ckey
    call doesnotexist("filopt")
  end subroutine filopt

  subroutine filsiz(cfl,nw,nh,iret)
    implicit none
    character (len = *), intent (in) :: cfl
    integer, intent (out) :: nw,nh,iret
    call doesnotexist("filsiz")
  end subroutine filsiz

  subroutine filtyp(cfl,iret)
    implicit none
    character (len = *), intent (in) :: cfl
    integer, intent (out) :: iret
    call doesnotexist("filtyp")
  end subroutine filtyp

  subroutine filwin(nx,ny,nw,nh)
    implicit none
    integer, intent (in) :: nx,ny,nw,nh
    call doesnotexist("filwin")
  end subroutine filwin

  subroutine fitscls()
    call doesnotexist("fitscls")
  end subroutine fitscls

  subroutine fitsflt(ckey,xv)
    implicit none
    character (len = *), intent (in) :: ckey
    double precision, intent (out) :: xv
    call doesnotexist("fitsflt")
  end subroutine fitsflt

  subroutine fitshdu(nhdu,n)
    implicit none
    integer, intent (in) :: nhdu
    integer, intent (out) :: n
    call doesnotexist("fitshdu")
  end subroutine fitshdu

  subroutine fitsimg(iray,nmax,n)
    implicit none
    character (len=1), intent (in), dimension (*) :: iray
    integer, intent (in) :: nmax
    integer, intent (out) :: n
    call doesnotexist("fitsimg")
  end subroutine fitsimg

  subroutine fitsopn(cfl,iret)
    implicit none
    character (len = *), intent (in) :: cfl
    integer, intent (out) :: iret
    call doesnotexist("fitsopn")
  end subroutine fitsopn

  subroutine fitsstr(ckey,cval,nmax)
    implicit none
    character (len = *), intent (in) :: ckey
    character (len = *), intent (out) :: cval
    integer, intent (in) :: nmax
    call doesnotexist("fitsstr")
  end subroutine fitsstr

  subroutine fitstyp(ckey,iv)
    implicit none
    character (len = *), intent (in) :: ckey
    integer, intent (out) :: iv
    call doesnotexist("fitstyp")
  end subroutine fitstyp

  subroutine fitsval(ckey,iv)
    implicit none
    character (len = *), intent (in) :: ckey
    integer, intent (out) :: iv
    call doesnotexist("fitsval")
  end subroutine fitsval

  subroutine fixspc(x)
    implicit none
    double precision, intent (in) :: x
    call doesnotexist("fixspc")
  end subroutine fixspc

  subroutine flab3d()
    call doesnotexist("flab3d")
  end subroutine flab3d

  subroutine flen(x,ndez,nx)
    implicit none
    double precision, intent (in) :: x
    integer, intent (in) :: ndez
    integer, intent (out) :: nx
    call doesnotexist("flen")
  end subroutine flen

  subroutine frame(i)
    implicit none
    integer, intent (in):: i
    call doesnotexist("frame")
  end subroutine frame

  subroutine frmbar(i)
    implicit none
    integer, intent (in):: i
    call doesnotexist("frmbar")
  end subroutine frmbar

  subroutine frmclr(i)
    implicit none
    integer, intent (in):: i
    call doesnotexist("frmclr")
  end subroutine frmclr

  subroutine frmess(i)
    implicit none
    integer, intent (in):: i
    call doesnotexist("frmess")
  end subroutine frmess

  subroutine gapcrv(x)
    implicit none
    double precision, intent (in) :: x
    call doesnotexist("gapcrv")
  end subroutine gapcrv

  subroutine gapsiz(x,cax)
    implicit none
    double precision, intent (in) :: x
    character (len = *), intent (in) :: cax
    call doesnotexist("gapsiz")
  end subroutine gapsiz

  subroutine gaxpar(a1,a2,copt,cax,a,b,or,stp,ndig)
    implicit none
    double precision, intent (in) :: a1,a2
    double precision, intent (out) :: a,b,or,stp
    integer, intent (out) :: ndig
    character (len=*), intent (in) :: copt, cax
    call doesnotexist("gaxpar")
  end subroutine gaxpar

  subroutine getalf(cstr)
    implicit none
    character (len = *), intent (out) :: cstr
    call doesnotexist("getalf")
  end subroutine getalf

  subroutine getang(i)
    implicit none
    integer, intent (out) :: i
    call doesnotexist("getang")
  end subroutine getang

  subroutine getbpp(i)
    implicit none
    integer, intent (out) :: i
    call doesnotexist("getbpp")
  end subroutine getbpp

  subroutine getclp(nx,ny,nw,nh)
    implicit none
    integer, intent (out) :: nx,ny,nw,nh
    call doesnotexist("getclp")
  end subroutine getclp

  subroutine getclr(i)
    implicit none
    integer, intent (out) :: i
    call doesnotexist("getclr")
  end subroutine getclr

  subroutine getdig(i,j,k)
    implicit none
    integer, intent (out) :: i,j,k
    call doesnotexist("getdig")
  end subroutine getdig

  subroutine getdsp(cdsp)
    implicit none
    character (len = *), intent (out) :: cdsp
    call doesnotexist("getdsp")
  end subroutine getdsp

  subroutine getfil(cstr)
    implicit none
    character (len = *), intent (out) :: cstr
    call doesnotexist("getfil")
  end subroutine getfil

  subroutine getgrf(a,e,or,step,copt)
    implicit none
    double precision, intent (out) :: a,e,or,step
    character (len = *), intent (in) :: copt
    call doesnotexist("getgrf")
  end subroutine getgrf

  subroutine gethgt(i)
    implicit none
    integer, intent (out) :: i
    call doesnotexist("gethgt")
  end subroutine gethgt

  subroutine gethnm(i)
    implicit none
    integer, intent (out) :: i
    call doesnotexist("gethnm")
  end subroutine gethnm

  subroutine getind(i,xr,xg,xb)
    implicit none
    integer, intent (in) :: i
    double precision, intent (out) :: xr,xg,xb
    call doesnotexist("getind")
  end subroutine getind

  subroutine getico(x,y,xp,yp)
    implicit none
    double precision, intent (in) :: x,y
    double precision, intent (out) :: xp,yp
    call doesnotexist("getico")
  end subroutine getico

  subroutine getlab(c1,c2,c3)
    implicit none
    character (len = *), intent (out) :: c1,c2,c3
    call doesnotexist("getlab")
  end subroutine getlab

  subroutine getlen(i,j,k)
    implicit none
    integer, intent (out) :: i,j,k
    call doesnotexist("getlen")
  end subroutine getlen

  subroutine getlev(i)
    implicit none
    integer, intent (out) :: i
    call doesnotexist("getlev")
  end subroutine getlev

  subroutine getlin(i)
    implicit none
    integer, intent (out) :: i
    call doesnotexist("getlin")
  end subroutine getlin

  subroutine getlit(xp,yp,zp,xn,yn,zn,i)
    implicit none
    integer, intent (in) :: xp,yp,zp,xn,yn,zn
    integer, intent (out) :: i
    call doesnotexist("getlit")
  end subroutine getlit

  subroutine getmat(x,y,z,n,zmat,nx,ny,zval,imat,wmat)
    implicit none
    integer, intent (in) :: n,nx,ny
    double precision, dimension (n), intent (in) :: x,y,z
    double precision, dimension (nx,ny), intent (out) :: zmat
    double precision, intent (in) :: zval
    integer, dimension (nx,ny), intent (in out) :: imat
    double precision, dimension (nx,ny), intent (in out) :: wmat
    call doesnotexist("getmat")
  end subroutine getmat

  subroutine getmfl(cstr)
    implicit none
    character (len = *), intent (out) :: cstr
    call doesnotexist("getmfl")
  end subroutine getmfl

  subroutine getmix(c,cstr)
    implicit none
    character (len = *), intent (out) :: c
    character (len = *), intent (in) :: cstr
    call doesnotexist("getmix")
  end subroutine getmix

  subroutine getor(i,j)
    implicit none
    integer, intent (out) :: i,j
    call doesnotexist("getor")
  end subroutine getor

  subroutine getpag(i,j)
    implicit none
    integer, intent (out) :: i,j
    call doesnotexist("getpag")
  end subroutine getpag

  subroutine getpat(i)
    implicit none
    integer, intent (out) :: i
    call doesnotexist("getpat")
  end subroutine getpat

  subroutine getplv(i)
    implicit none
    integer, intent (out) :: i
    call doesnotexist("getplv")
  end subroutine getplv

  subroutine getpos(i,j)
    implicit none
    integer, intent (out) :: i,j
    call doesnotexist("getpos")
  end subroutine getpos

  subroutine getran(i,j)
    implicit none
    integer, intent (out) :: i,j
    call doesnotexist("getran")
  end subroutine getran

  subroutine getrco(x,y,xp,yp)
    implicit none
    double precision, intent (in) :: x,y
    double precision, intent (out) :: xp,yp
    call doesnotexist("getrco")
  end subroutine getrco

  subroutine getres(i,j)
    implicit none
    integer, intent (out) :: i,j
    call doesnotexist("getres")
  end subroutine getres

  subroutine getrgb(xr,xg,xb)
    implicit none
    double precision, intent (out) :: xr,xg,xb
    call doesnotexist("getrgb")
  end subroutine getrgb

  subroutine getscl(i,j,k)
    implicit none
    integer, intent (out) :: i,j,k
    call doesnotexist("getscl")
  end subroutine getscl

  subroutine getscm(i,j,k)
    implicit none
    integer, intent (out) :: i,j,k
    call doesnotexist("getscm")
  end subroutine getscm

  subroutine getscr(i,j)
    implicit none
    integer, intent (out) :: i,j
    call doesnotexist("getscr")
  end subroutine getscr

  subroutine getshf(cstr,c)
    implicit none
    character (len = *), intent (out) :: c
    character (len = *), intent (in) :: cstr
    call doesnotexist("getshf")
  end subroutine getshf

  subroutine getsp1(i,j,k)
    implicit none
    integer, intent (out) :: i,j,k
    call doesnotexist("getsp1")
  end subroutine getsp1

  subroutine getsp2(i,j,k)
    implicit none
    integer, intent (out) :: i,j,k
    call doesnotexist("getsp2")
  end subroutine getsp2

  subroutine getsym(i,j)
    implicit none
    integer, intent (out) :: i,j
    call doesnotexist("getsym")
  end subroutine getsym

  subroutine gettcl(i,j)
    implicit none
    integer, intent (out) :: i,j
    call doesnotexist("gettcl")
  end subroutine gettcl

  subroutine gettic(i,j,k)
    implicit none
    integer, intent (out) :: i,j,k
    call doesnotexist("gettic")
  end subroutine gettic

  subroutine gettyp(i)
    implicit none
    integer, intent (out) :: i
    call doesnotexist("gettyp")
  end subroutine gettyp

  subroutine getuni(i)
    implicit none
    integer, intent (out) :: i
    call doesnotexist("getuni")
  end subroutine getuni

  subroutine getver(xver)
    implicit none
    double precision, intent (out) :: xver
    call doesnotexist("getver")
  end subroutine getver

  subroutine getvk(i,j,k)
    implicit none
    integer, intent (out) :: i,j,k
    call doesnotexist("getvk")
  end subroutine getvk

  subroutine getvlt(ctab)
    implicit none
    character (len = *), intent (out) :: ctab
    call doesnotexist("getvlt")
  end subroutine getvlt

  subroutine getwid(i)
    implicit none
    integer, intent (out) :: i
    call doesnotexist("getwid")
  end subroutine getwid

  subroutine getwin(ix,iy,nw,nh)
    implicit none
    integer, intent (out) :: ix,iy,nw,nh
    call doesnotexist("getwin")
  end subroutine getwin

  subroutine getxid (ival, copt)
    implicit none
    integer, intent (out) :: ival
    character (len = *), intent (in) :: copt
    call doesnotexist("getxid")
  end subroutine getxid

  subroutine gifmod(cmod,ckey)
    implicit none
    character (len = *), intent (in) :: cmod,ckey
    call doesnotexist("gifmod")
  end subroutine gifmod

  subroutine gmxalf(calph,ca,cb,n)
    implicit none
    character (len = *), intent (in) :: calph
    character (len = *), intent (out) :: ca,cb
    integer, intent (out) :: n
    call doesnotexist("gmxalf")
  end subroutine gmxalf

  subroutine gothic()
    call doesnotexist("gothic")
  end subroutine gothic

  subroutine grace(i)
    implicit none
    integer, intent (in) :: i
    call doesnotexist("grace")
  end subroutine grace

  subroutine graf(ax,ex,orx,stepx,ay,ey,ory,stepy)
    implicit none
    double precision, intent (in) :: ax,ex,orx,stepx,ay,ey,ory,stepy
    call doesnotexist("graf")
  end subroutine graf

  subroutine graf3(ax,ex,orx,stepx,ay,ey,ory,stepy,az,ez,orz,stepz)
    implicit none
    double precision, intent (in) :: ax,ex,orx,stepx,ay,ey,ory,stepy, &
                                     az,ez,orz,stepz
    call doesnotexist("graf3")
  end subroutine graf3

  subroutine graf3d(ax,ex,orx,stepx,ay,ey,ory,stepy,az,ez,orz,stepz)
    implicit none
    double precision, intent (in) :: ax,ex,orx,stepx,ay,ey,ory,stepy, &
                                     az,ez,orz,stepz
    call doesnotexist("graf3d")
  end subroutine graf3d

  subroutine grafmp(ax,ex,orx,stepx,ay,ey,ory,stepy)
    implicit none
    double precision, intent (in) :: ax,ex,orx,stepx,ay,ey,ory,stepy
    call doesnotexist("grafmp")
  end subroutine grafmp

  subroutine grafp(ex,orx,stepx,ory,stepy)
    implicit none
    double precision, intent (in) :: ex,orx,stepx,ory,stepy
    call doesnotexist("grafp")
  end subroutine grafp

  subroutine grafr(zre,nre,zimg,nimg)
    implicit none
    integer, intent (in) :: nre,nimg
    double precision, dimension (nre), intent (in) :: zre
    double precision, dimension (nimg), intent (in) :: zimg
    call doesnotexist("grafr")
  end subroutine grafr

  subroutine grdpol(igrd,jgrd)
    implicit none
    integer, intent (in) :: igrd,jgrd
    call doesnotexist("grdpol")
  end subroutine grdpol

  subroutine grffin()
    implicit none
    call doesnotexist("grffin")
  end subroutine grffin

  subroutine grfimg(cstr)
    implicit none
    character (len = *), intent (in) :: cstr
    call doesnotexist("grfimg")
  end subroutine grfimg

  subroutine grfini(x1,y1,z1,x2,y2,z2,x3,y3,z3)
    implicit none
    double precision, intent (in) :: x1,y1,z1,x2,y2,z2,x3,y3,z3
    call doesnotexist("grfini")
  end subroutine grfini

  subroutine grid(i,j)
    implicit none
    integer, intent (in) :: i,j
    call doesnotexist("grid")
  end subroutine grid

  subroutine grid3d(igrid,jgrid,copt)
    implicit none
    integer, intent (in) :: igrid,jgrid
    character (len = *), intent (in) :: copt
    call doesnotexist("grid3d")
  end subroutine grid3d

  subroutine gridim(zim,zre1,zre2,n)
    implicit none
    double precision, intent (in) :: zim,zre1,zre2
    integer, intent (in) :: n
    call doesnotexist("gridim")
  end subroutine gridim

  subroutine gridmp(i,j)
    implicit none
    integer, intent (in) :: i,j
    call doesnotexist("gridmp")
  end subroutine gridmp

  subroutine gridre(zre,zim1,zim2,n)
    implicit none
    double precision, intent (in) :: zre,zim1,zim2
    integer, intent (in) :: n
    call doesnotexist("gridre")
  end subroutine gridre

  subroutine gwgatt (id,ival,copt)
    implicit none
    integer, intent (in) :: id
    integer, intent (out) :: ival
    character (len=*), intent (in) :: copt
    call doesnotexist("gwgatt")
  end subroutine gwgatt

  subroutine gwgbox(id,ival)
    implicit none
    integer, intent (in) :: id
    integer, intent (out) :: ival
    call doesnotexist("gwgbox")
  end subroutine gwgbox

  subroutine gwgbut(id,ival)
    implicit none
    integer, intent (in) :: id
    integer, intent (out) :: ival
    call doesnotexist("gwgbut")
  end subroutine gwgbut

  subroutine gwgfil(id,cstr)
    implicit none
    integer, intent (in) :: id
    character (len=*), intent (out) :: cstr
    call doesnotexist("gwgfil")
  end subroutine gwgfil

  subroutine gwgflt(id,xv)
    implicit none
    integer, intent (in) :: id
    double precision, intent (out) :: xv
    call doesnotexist("gwgflt")
  end subroutine gwgflt

  subroutine gwggui(ival)
    implicit none
    integer, intent (out) :: ival
    call doesnotexist("gwggui")
  end subroutine gwggui

  subroutine gwgint(id,iv)
    implicit none
    integer, intent (in) :: id
    integer, intent (out) :: iv
    call doesnotexist("gwgint")
  end subroutine gwgint

  subroutine gwglis(id,ival)
    implicit none
    integer, intent (in) :: id
    integer, intent (out) :: ival
    call doesnotexist("gwglis")
  end subroutine gwglis

  subroutine gwgscl(id,xval)
    implicit none
    integer, intent (in) :: id
    double precision, intent (out) :: xval
    call doesnotexist("gwgscl")
  end subroutine gwgscl

  subroutine gwgsiz (id,nw,nh)
    implicit none
    integer, intent (in) :: id
    integer, intent (out) :: nw,nh
    call doesnotexist("gwgsiz")
  end subroutine gwgsiz

  subroutine gwgtbf(id,i,j,xv)
    implicit none
    integer, intent (in) :: id,i,j
    double precision, intent (out) :: xv
    call doesnotexist("gwgtbf")
  end subroutine gwgtbf

  subroutine gwgtbi(id,i,j,iv)
    implicit none
    integer, intent (in) :: id,i,j
    integer, intent (out) :: iv
    call doesnotexist("gwgtbi")
  end subroutine gwgtbi

  subroutine gwgtbl(id,xray,n,idx,copt)
    implicit none
    integer, intent (in) :: id,n,idx
    double precision, dimension (n), intent (out) :: xray
    character (len=*), intent (in) :: copt
    call doesnotexist("gwgtbl")
  end subroutine gwgtbl

  subroutine gwgtbs(id,i,j,cstr)
    implicit none
    integer, intent (in) :: id,i,j
    character (len=*), intent (out) :: cstr
    call doesnotexist("gwgtbs")
  end subroutine gwgtbs

  subroutine gwgtxt(id,cstr)
    implicit none
    integer, intent (in) :: id
    character (len=*), intent (out) :: cstr
    call doesnotexist("gwgtxt")
  end subroutine gwgtxt

  subroutine gwgxid(id,ival)
    implicit none
    integer, intent (in) :: id
    integer, intent (out) :: ival
    call doesnotexist("gwgxid")
  end subroutine gwgxid

  subroutine height(i)
    implicit none
    integer, intent (in) :: i
    call doesnotexist("height")
  end subroutine height

  subroutine helve()
    call doesnotexist("helve")
  end subroutine helve

  subroutine helves()
    call doesnotexist("helves")
  end subroutine helves

  subroutine helvet()
    call doesnotexist("helvet")
  end subroutine helvet

  subroutine hidwin(i,copt)
    implicit none
    integer, intent (in) :: i
    character (len = *), intent (in) :: copt
    call doesnotexist("hidwin")
  end subroutine hidwin

  subroutine histog(xray,n,x,y,m)
    implicit none
    integer, intent (in) :: n
    double precision, dimension (n), intent (in) :: xray
    integer, intent (out) :: m
    double precision, dimension (n), intent (out) :: x,y
    call doesnotexist("histog")
  end subroutine histog

  subroutine hname(i)
    implicit none
    integer, intent (in) :: i
    call doesnotexist("hname")
  end subroutine hname

  subroutine hpgmod(cmod,ckey)
    implicit none
    character (len = *), intent (in) :: cmod,ckey
    call doesnotexist("hpgmod")
  end subroutine hpgmod

  subroutine hsvrgb(xh,xs,xv,r,g,b)
    implicit none
    double precision, intent (in)  :: xh,xs,xv
    double precision, intent (out) :: r,g,b
    call doesnotexist("hsvrgb")
  end subroutine hsvrgb

  subroutine hsym3d(xh)
    implicit none
    double precision, intent (in) :: xh
    call doesnotexist("hsym3d")
  end subroutine hsym3d

  subroutine hsymbl(i)
    implicit none
    integer, intent (in) :: i
    call doesnotexist("hsymbl")
  end subroutine hsymbl

  subroutine htitle(i)
    implicit none
    integer, intent (in) :: i
    call doesnotexist("htitle")
  end subroutine htitle

  subroutine hwfont()
    call doesnotexist("hwfont")
  end subroutine hwfont

  subroutine hwmode(copt,ckey)
    implicit none
    character (len = *), intent (in) :: copt,ckey
    call doesnotexist("hwmode")
  end subroutine hwmode

  subroutine hworig(nx,ny)
    implicit none
    integer, intent (in) :: nx,ny
    call doesnotexist("hworig")
  end subroutine hworig

  subroutine hwpage(nxp,nyp)
    implicit none
    integer, intent (in) :: nxp,nyp
    call doesnotexist("hwpage")
  end subroutine hwpage

  subroutine hwscal(x)
    implicit none
    double precision, intent (in) :: x
    call doesnotexist("hwscal")
  end subroutine hwscal

  subroutine imgbox(nx,ny,nw,nh)
    implicit none
    integer, intent (in) :: nx,ny,nw,nh
    call doesnotexist("imgbox")
  end subroutine imgbox

  subroutine imgclp(nx,ny,nw,nh)
    implicit none
    integer, intent (in) :: nx,ny,nw,nh
    call doesnotexist("imgclp")
  end subroutine imgclp

  subroutine imgfin()
    call doesnotexist("imgfin")
  end subroutine imgfin

  subroutine imgfmt(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("imgfmt")
  end subroutine imgfmt

  subroutine imgini()
    call doesnotexist("imgini")
  end subroutine imgini

  subroutine imgmod(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("imgmod")
  end subroutine imgmod

  subroutine imgsiz(nw,nh)
    implicit none
    integer, intent (in) :: nw,nh
    call doesnotexist("imgsiz")
  end subroutine imgsiz

  subroutine imgtpr(n)
    implicit none
    integer, intent (in) :: n
    call doesnotexist("imgtpr")
  end subroutine imgtpr

  subroutine inccrv(i)
    implicit none
    integer, intent (in) :: i
    call doesnotexist("inccrv")
  end subroutine inccrv

  function incdat(id,im,iy)
    implicit none
    integer, intent (in) :: id, im, iy
    integer :: incdat
    call doesnotexist("incdat")
  end function incdat

  subroutine incfil(cstr)
    implicit none
    character (len = *), intent (in) :: cstr
    call doesnotexist("incfil")
  end subroutine incfil

  subroutine incmrk(i)
    implicit none
    integer, intent (in) :: i
    call doesnotexist("incmrk")
  end subroutine incmrk

  function indrgb(xr,xg,xb)
    implicit none
    double precision, intent (in) :: xr,xg,xb
    integer :: indrgb
    call doesnotexist("indrgb")
  end function indrgb

  subroutine intax()
    call doesnotexist("intax")
  end subroutine intax

  subroutine intcha(num,n,cnum)
    implicit none
    integer, intent (in)  :: num
    integer, intent (out) :: n
    character (len = *), intent (out) :: cnum
    call doesnotexist("intcha")
  end subroutine intcha

  subroutine intlen(nm,nlaen)
    implicit none
    integer, intent (in)  :: nm
    integer, intent (out) :: nlaen
    call doesnotexist("intlen")
  end subroutine intlen

  function intrgb(xr,xg,xb)
    implicit none
    double precision, intent (in) :: xr,xg,xb
    integer :: intrgb
    call doesnotexist("intrgb")
  end function intrgb

  subroutine intutf(iray,nray,cstr,nmax,nl)
    implicit none
    integer, intent (in) :: nray,nmax
    integer, dimension (nray), intent (in) :: iray
    character (len=*), intent (out) :: cstr
    integer, intent (out) :: nl
    call doesnotexist("intutf")
  end subroutine intutf

  subroutine isopts(xray,nx,yray,ny,zray,nz,wmat,wlev, &
                    xtri,ytri,ztri,nmax,ntri)
    implicit none
    integer, intent (in) :: nx,ny,nz,nmax
    integer, intent (out) :: ntri
    double precision, dimension (nx), intent (in) :: xray
    double precision, dimension (ny), intent (in) :: yray
    double precision, dimension (nz), intent (in) :: zray
    double precision, dimension (nx,ny,nz), intent (in) :: wmat
    double precision, intent (in) :: wlev
    double precision, dimension (nmax), intent (out) :: xtri,ytri,ztri
    call doesnotexist("isopts")
  end subroutine isopts

  subroutine itmcat(clis,cstr)
    implicit none
    character (len=*), intent (in out) :: clis
    character (len=*), intent (in) :: cstr
    call doesnotexist("itmcat")
  end subroutine itmcat

  subroutine itmncat(clis,nmx,cstr)
    implicit none
    character (len=*), intent (in out) :: clis
    character (len=*), intent (in) :: cstr
    integer, intent (in) :: nmx
    call doesnotexist("itmncat")
  end subroutine itmncat

  function itmcnt(clis)
    implicit none
    character (len=*), intent (in) :: clis
    integer :: itmcnt
    call doesnotexist("itmcnt")
  end function itmcnt

  subroutine itmstr(clis,nlis,cstr)
    implicit none
    character (len=*), intent (in) :: clis
    character (len=*), intent (out) :: cstr
    integer, intent (in) :: nlis
    call doesnotexist("itmstr")
  end subroutine itmstr

  subroutine jusbar(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("jusbar")
  end subroutine jusbar

  subroutine labclr(iclr,copt)
    implicit none
    integer, intent (in) :: iclr
    character (len = *), intent (in) :: copt
    call doesnotexist("labclr")
  end subroutine labclr

  subroutine labdig(i,cax)
    implicit none
    integer, intent (in) :: i
    character (len = *), intent (in) :: cax
    call doesnotexist("labdig")
  end subroutine labdig

  subroutine labdis(i,cax)
    implicit none
    integer, intent (in) :: i
    character (len = *), intent (in) :: cax
    call doesnotexist("labdis")
  end subroutine labdis

  subroutine labels(copt,cax)
    implicit none
    character (len = *), intent (in) :: copt,cax
    call doesnotexist("labels")
  end subroutine labels

  subroutine labjus(copt,cax)
    implicit none
    character (len = *), intent (in) :: copt,cax
    call doesnotexist("labjus")
  end subroutine labjus

  subroutine labl3d(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("labl3d")
  end subroutine labl3d

  subroutine labmod(ckey,copt,cax)
    implicit none
    character (len = *), intent (in) :: ckey,copt,cax
    call doesnotexist("labmod")
  end subroutine labmod

  subroutine labpos(copt,cax)
    implicit none
    character (len = *), intent (in) :: copt,cax
    call doesnotexist("labpos")
  end subroutine labpos

  subroutine labtyp(copt,cax)
    implicit none
    character (len = *), intent (in) :: copt,cax
    call doesnotexist("labtyp")
  end subroutine labtyp

  subroutine ldimg(cstr,iray,nmax,nc,iret)
    implicit none
    character (len = *), intent (in) :: cstr
    integer, intent (in) :: nmax,nc
    integer, intent (out) :: iret
    integer (kind=selected_int_kind(4)), dimension (*), &
            intent (out) :: iray
    call doesnotexist("ldimg")
  end subroutine ldimg

  subroutine legbgd(n)
    implicit none
    integer, intent (in) :: n
    call doesnotexist("legbgd")
  end subroutine legbgd

  subroutine legclr()
    call doesnotexist("legclr")
  end subroutine legclr

  subroutine legend(cbf,ncor)
    implicit none
    character (len = *), intent (in) :: cbf
    integer, intent (in) :: ncor
    call doesnotexist("legend")
  end subroutine legend

  subroutine legini(cbf,nlin,nmax)
    implicit none
    character (len = *), intent (in out) :: cbf
    integer, intent (in) :: nlin, nmax
    call doesnotexist("legini")
  end subroutine legini

  subroutine leglin(cbf,cstr,n)
    implicit none
    character (len = *), intent (in out) :: cbf
    character (len = *), intent (in) :: cstr
    integer, intent (in) :: n
    call doesnotexist("leglin")
  end subroutine leglin

  subroutine legopt(x1,x2,x3)
    implicit none
    double precision, intent (in) :: x1,x2,x3
    call doesnotexist("legopt")
  end subroutine legopt

  subroutine legpat(ilin,ithk,isym,iclr,ipat,i)
    implicit none
    integer, intent (in) :: ilin,ithk,isym,iclr,ipat,i
    call doesnotexist("legpat")
  end subroutine legpat

  subroutine legpos(nx,ny)
    implicit none
    integer, intent (in) :: nx,ny
    call doesnotexist("legpos")
  end subroutine legpos

  subroutine legsel(nray,n)
    implicit none
    integer, intent (in) :: n
    integer, dimension (n), intent (in) :: nray
    call doesnotexist("legsel")
  end subroutine legsel

  subroutine legtit(cstr)
    implicit none
    character (len = *), intent (in) :: cstr
    call doesnotexist("legtit")
  end subroutine legtit

  subroutine legtyp(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("legtyp")
  end subroutine legtyp

  subroutine legval(x,copt)
    implicit none
    double precision, intent (in) :: x
    character (len = *), intent (in) :: copt
    call doesnotexist("legval")
  end subroutine legval

  subroutine lfttit()
    call doesnotexist("lfttit")
  end subroutine lfttit

  subroutine licmod(cmod,ckey)
    implicit none
    character (len = *), intent (in) :: cmod,ckey
    call doesnotexist("licmod")
  end subroutine licmod

  subroutine licpts(xmat,ymat,nx,ny,itmat,iwmat,wmat)
    implicit none
    integer, intent (in) :: nx,ny
    double precision, dimension (nx,ny), intent (in) :: xmat,ymat
    integer, dimension (nx,ny), intent (in) :: itmat
    integer, dimension (nx,ny), intent (out) :: iwmat
    double precision, dimension (nx,ny), intent (out) :: wmat
    call doesnotexist("licpts")
  end subroutine licpts

  subroutine linclr(nray,n)
    implicit none
    integer, intent (in) :: n
    integer, dimension (n), intent (in) :: nray
    call doesnotexist("linclr")
  end subroutine linclr

  subroutine lincyc(i,ilin)
    implicit none
    integer, intent (in) :: i,ilin
    call doesnotexist("lincyc")
  end subroutine lincyc

  subroutine line(nx,ny,nu,nv)
    implicit none
    integer, intent (in) :: nx,ny,nu,nv
    call doesnotexist("line")
  end subroutine line

  subroutine linesp(x)
    implicit none
    double precision, intent (in) :: x
    call doesnotexist("linesp")
  end subroutine linesp

  subroutine linfit(x,y,n,a,b,r,copt)
    implicit none
    integer, intent (in) :: n
    double precision, dimension (n), intent (in) :: x,y
    double precision, intent (out) :: a,b,r
    character (len = *), intent (in) :: copt
    call doesnotexist("linfit")
  end subroutine linfit

  subroutine linmod(cmod,ckey)
    implicit none
    character (len = *), intent (in) :: cmod,ckey
    call doesnotexist("linmod")
  end subroutine linmod

  subroutine lintyp(i)
    implicit none
    integer, intent (in) :: i
    call doesnotexist("lintyp")
  end subroutine lintyp

  subroutine linwid(i)
    implicit none
    integer, intent (in) :: i
    call doesnotexist("linwid")
  end subroutine linwid

  subroutine light(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("light")
  end subroutine light

  subroutine litmod(id,copt)
    implicit none
    integer, intent (in) :: id
    character (len = *), intent (in) :: copt
    call doesnotexist("litmod")
  end subroutine litmod

  subroutine litopt(id,x,copt)
    implicit none
    integer, intent (in) :: id
    double precision, intent (in) :: x
    character (len = *), intent (in) :: copt
    call doesnotexist("litopt")
  end subroutine litopt

  subroutine litop3(id,xr,xg,xb,copt)
    implicit none
    integer, intent (in) :: id
    double precision, intent (in) :: xr,xg,xb
    character (len = *), intent (in) :: copt
    call doesnotexist("litop3")
  end subroutine litop3

  subroutine litpos(id,xp,yp,zp,copt)
    implicit none
    integer, intent (in) :: id
    double precision, intent (in) :: xp,yp,zp
    character (len = *), intent (in) :: copt
    call doesnotexist("litpos")
  end subroutine litpos

  subroutine lncap(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("lncap")
  end subroutine lncap

  subroutine lnjoin(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("lnjoin")
  end subroutine lnjoin

  subroutine lnmlt(x)
    implicit none
    double precision, intent (in) :: x
    call doesnotexist("lnmlt")
  end subroutine lnmlt

  subroutine logtic(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("logtic")
  end subroutine logtic

  subroutine lsechk(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("lsechk")
  end subroutine lsechk

  subroutine mapbas(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("mapbas")
  end subroutine mapbas

  subroutine mapfil(cfl, copt)
    implicit none
    character (len = *), intent (in) :: cfl, copt
    call doesnotexist("mapfil")
  end subroutine mapfil

  subroutine mapimg(cfl,x1,x2,x3,x4,x5,x6)
    implicit none
    character (len = *), intent (in) :: cfl
    double precision, intent (in) :: x1,x2,x3,x4,x5,x6
    call doesnotexist("mapimg")
  end subroutine mapimg

  subroutine maplab(copt,ckey)
    implicit none
    character (len = *), intent (in) :: copt,ckey
    call doesnotexist("maplab")
  end subroutine maplab

  subroutine maplev(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("maplev")
  end subroutine maplev

  subroutine mapmod(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("mapmod")
  end subroutine mapmod

  subroutine mappol(x,y)
    implicit none
    double precision, intent (in) :: x,y
    call doesnotexist("mappol")
  end subroutine mappol

  subroutine mapopt(copt,ckey)
    implicit none
    character (len = *), intent (in) :: copt,ckey
    call doesnotexist("mapopt")
  end subroutine mapopt

  subroutine mapref(ylw,yup)
    implicit none
    double precision, intent (in) :: ylw,yup
    call doesnotexist("mapref")
  end subroutine mapref

  subroutine mapsph(xrad)
    implicit none
    double precision, intent (in) :: xrad
    call doesnotexist("mapsph")
  end subroutine mapsph

  subroutine marker(i)
    implicit none
    integer, intent (in) :: i
    call doesnotexist("marker")
  end subroutine marker

  subroutine matopt(x,copt)
    implicit none
    double precision, intent (in) :: x
    character (len = *), intent (in) :: copt
    call doesnotexist("matopt")
  end subroutine matopt

  subroutine matop3(xr,xg,xb,copt)
    implicit none
    double precision, intent (in) :: xr,xg,xb
    character (len = *), intent (in) :: copt
    call doesnotexist("matop3")
  end subroutine matop3

  subroutine mdfmat(i,j,x)
    implicit none
    integer, intent (in) :: i,j
    double precision, intent (in) :: x
    call doesnotexist("mdfmat")
  end subroutine mdfmat

  subroutine messag(cstr,nx,ny)
    implicit none
    character (len = *), intent (in) :: cstr
    integer, intent (in) :: nx,ny
    call doesnotexist("messag")
  end subroutine messag

  subroutine metafl(ct)
    implicit none
    character (len = *), intent (in) :: ct
    call doesnotexist("metafl")
  end subroutine metafl

  subroutine mixalf()
    call doesnotexist("mixalf")
  end subroutine mixalf

  subroutine mixleg()
    call doesnotexist("mixleg")
  end subroutine mixleg

  subroutine moment(xray,n,copt,xv)
    implicit none
    integer, intent (in) :: n
    double precision, dimension (n), intent (in) :: xray
    character (len = *), intent (in) :: copt
    double precision, intent (out) :: xv
    call doesnotexist("moment")
  end subroutine moment

  subroutine mpaepl(i)
    implicit none
    integer, intent (in) :: i
    call doesnotexist("mpaepl")
  end subroutine mpaepl

  subroutine mplang(x)
    implicit none
    double precision, intent (in) :: x
    call doesnotexist("mplang")
  end subroutine mplang

  subroutine mplclr(i,j)
    implicit none
    integer, intent (in) :: i,j
    call doesnotexist("mplclr")
  end subroutine mplclr

  subroutine mplpos(i,j)
    implicit none
    integer, intent (in) :: i,j
    call doesnotexist("mplpos")
  end subroutine mplpos

  subroutine mplsiz(i)
    implicit none
    integer, intent (in) :: i
    call doesnotexist("mplsiz")
  end subroutine mplsiz

  subroutine mpslogo(nx,ny,nsize,copt)
    implicit none
    integer, intent (in) :: nx,ny,nsize
    character (len = *), intent (in) :: copt
    call doesnotexist("mpslogo")
  end subroutine mpslogo

  subroutine mrkclr(i)
    implicit none
    integer, intent (in) :: i
    call doesnotexist("mrkclr")
  end subroutine mrkclr

  subroutine msgbox(cstr)
    implicit none
    character (len = *), intent (in) :: cstr
    call doesnotexist("msgbox")
  end subroutine msgbox

  subroutine mshclr(ic)
    implicit none
    integer, intent (in) :: ic
    call doesnotexist("mshclr")
  end subroutine mshclr

  subroutine mshcrv(n)
    implicit none
    integer, intent (in) :: n
    call doesnotexist("mshcrv")
  end subroutine mshcrv

  subroutine mylab(cstr,i,cax)
    implicit none
    character (len = *), intent (in) :: cstr,cax
    integer, intent (in) :: i
    call doesnotexist("mylab")
  end subroutine mylab

  subroutine myline(nray,n)
    implicit none
    integer, intent (in) :: n
    integer, dimension (n), intent (in) :: nray
    call doesnotexist("myline")
  end subroutine myline

  subroutine mypat(iang,itype,idens,icross)
    implicit none
    integer, intent (in) :: iang,itype,idens,icross
    call doesnotexist("mypat")
  end subroutine mypat

  subroutine mypie(iseg,xdat,xper,nrad,noff,ang,nvbox,idrw,iann)
    implicit none
    integer, intent (in out) :: iseg,nrad,noff,nvbox,idrw,iann
    double precision, intent (in out)    :: xdat,xper,ang
    call doesnotexist("mypie")
  end subroutine mypie

  subroutine mysymb(xray,yray,n,isym,iflag)
    implicit none
    integer, intent (in) :: n,isym,iflag
    double precision, dimension (n), intent (in) :: xray,yray
    call doesnotexist("mysymb")
  end subroutine mysymb

  subroutine myvlt(xr,xg,xb,n)
    implicit none
    integer, intent (in) :: n
    double precision, dimension (n), intent (in) :: xr,xg,xb
    call doesnotexist("myvlt")
  end subroutine myvlt

  subroutine namdis(i,cax)
    implicit none
    integer, intent (in) :: i
    character (len = *), intent (in) :: cax
    call doesnotexist("namdis")
  end subroutine namdis

  subroutine name(cnam,cax)
    implicit none
    character (len = *), intent (in) :: cnam,cax
    call doesnotexist("name")
  end subroutine name

  subroutine namjus(copt,cax)
    implicit none
    character (len = *), intent (in) :: copt,cax
    call doesnotexist("namjus")
  end subroutine namjus

  subroutine nancrv(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("nancrv")
  end subroutine nancrv

  subroutine neglog(e)
    implicit none
    double precision, intent (in) :: e
    call doesnotexist("neglog")
  end subroutine neglog

  subroutine newmix()
    call doesnotexist("newmix")
  end subroutine newmix

  subroutine newpag()
    call doesnotexist("newpag")
  end subroutine newpag

  function nlmess(cstr)
    implicit none
    character (len = *), intent (in) :: cstr
    integer :: nlmess
    call doesnotexist("nlmess")
  end function nlmess

  function nlnumb(x,ndez)
    implicit none
    double precision, intent (in) :: x
    integer, intent (in) :: ndez
    integer :: nlnumb
    call doesnotexist("nlnumb")
  end function nlnumb

  subroutine noarln()
    call doesnotexist("noarln")
  end subroutine noarln

  subroutine nobar()
    call doesnotexist("nobar")
  end subroutine nobar

  subroutine nobgd()
    call doesnotexist("nobgd")
  end subroutine nobgd

  subroutine nochek()
    call doesnotexist("nochek")
  end subroutine nochek

  subroutine noclip()
    call doesnotexist("noclip")
  end subroutine noclip

  subroutine nofill()
    call doesnotexist("nofill")
  end subroutine nofill

  subroutine nograf()
    call doesnotexist("nograf")
  end subroutine nograf

  subroutine nohide()
    call doesnotexist("nohide")
  end subroutine nohide

  subroutine noline(cax)
    implicit none
    character (len = *), intent (in) :: cax
    call doesnotexist("noline")
  end subroutine noline

  subroutine number(x,ndez,nx,ny)
    implicit none
    double precision, intent (in) :: x
    integer, intent (in) :: ndez,nx,ny
    call doesnotexist("number")
  end subroutine number

  subroutine numfmt(ct)
    implicit none
    character (len = *), intent (in) :: ct
    call doesnotexist("numfmt")
  end subroutine numfmt

  subroutine numode(cdec,cgrp,cpos,cspc)
    implicit none
    character (len = *), intent (in) :: cdec,cgrp,cpos,cspc
    call doesnotexist("numode")
  end subroutine numode

  function nwkday(id,im,iy)
    implicit none
    integer, intent (in) :: id,im,iy
    integer :: nwkday
    call doesnotexist("nwkday")
  end function nwkday

  function nxlegn(cbf)
    implicit none
    character (len = *), intent (in) :: cbf
    integer :: nxlegn
    call doesnotexist("nxlegn")
  end function nxlegn

  function nxpixl(ix,iy)
    implicit none
    integer, intent (in) :: ix,iy
    integer :: nxpixl
    call doesnotexist("nxpixl")
  end function nxpixl

  function nxposn(x)
    implicit none
    double precision, intent (in) :: x
    integer :: nxposn
    call doesnotexist("nxposn")
  end function nxposn

  function nylegn(cbf)
    implicit none
    character (len = *), intent (in) :: cbf
    integer :: nylegn
    call doesnotexist("nylegn")
  end function nylegn

  function nypixl(ix,iy)
    implicit none
    integer, intent (in) :: ix,iy
    integer :: nypixl
    call doesnotexist("nypixl")
  end function nypixl

  function nyposn(y)
    implicit none
    double precision, intent (in) :: y
    integer :: nyposn
    call doesnotexist("nyposn")
  end function nyposn

  function nzposn(z)
    implicit none
    double precision, intent (in) :: z
    integer :: nzposn
    call doesnotexist("nzposn")
  end function nzposn

  subroutine openfl(cstr,nlu,irw,istat)
    implicit none
    character (len = *), intent (in) :: cstr
    integer, intent (in) :: nlu,irw
    integer, intent (out) :: istat
    call doesnotexist("openfl")
  end subroutine openfl

  subroutine opnwin(id)
    implicit none
    integer, intent (in) :: id
    call doesnotexist("opnwin")
  end subroutine opnwin

  subroutine origin(i,j)
    implicit none
    integer, intent (in) :: i,j
    call doesnotexist("origin")
  end subroutine origin

  subroutine page(nxp,nyp)
    implicit none
    integer, intent (in) :: nxp,nyp
    call doesnotexist("page")
  end subroutine page

  subroutine pagera()
    call doesnotexist("pagera")
  end subroutine pagera

  subroutine pagfll(n)
    implicit none
    integer, intent (in) :: n
    call doesnotexist("pagfll")
  end subroutine pagfll

  subroutine paghdr(c1,c2,iopt,idir)
    implicit none
    character (len = *), intent (in) :: c1,c2
    integer, intent (in) :: iopt,idir
    call doesnotexist("paghdr")
  end subroutine paghdr

  subroutine pagmod(cmod)
    implicit none
    character (len = *), intent (in) :: cmod
    call doesnotexist("pagmod")
  end subroutine pagmod

  subroutine pagorg(cmod)
    implicit none
    character (len = *), intent (in) :: cmod
    call doesnotexist("pagorg")
  end subroutine pagorg

  subroutine pagwin(nxp,nyp)
    implicit none
    integer, intent (in) :: nxp,nyp
    call doesnotexist("pagwin")
  end subroutine pagwin

  subroutine patcyc(i,ipat)
    implicit none
    integer, intent (in) :: i,ipat
    call doesnotexist("patcyc")
  end subroutine patcyc

  subroutine pdfbuf(iray,nmax,nn)
    implicit none
    integer, intent (in) :: nmax
    character (len=1), intent (out), dimension (*) :: iray
    integer, intent (out) :: nn
    call doesnotexist("pdfbuf")
  end subroutine pdfbuf

  subroutine pdfmod(cmod,ckey)
    implicit none
    character (len = *), intent (in) :: cmod,ckey
    call doesnotexist("pdfmod")
  end subroutine pdfmod

  subroutine pdfmrk(cstr,copt)
    implicit none
    character (len = *), intent (in) :: cstr,copt
    call doesnotexist("pdfmrk")
  end subroutine pdfmrk

  subroutine penwid(x)
    implicit none
    double precision, intent (in) :: x
    call doesnotexist("penwid")
  end subroutine penwid

  subroutine pie(nx,ny,nr,a,b)
    implicit none
    integer, intent (in) :: nx,ny,nr
    double precision, intent (in) :: a,b
    call doesnotexist("pie")
  end subroutine pie

  subroutine piebor(iclr)
    implicit none
    integer, intent (in) :: iclr
    call doesnotexist("piebor")
  end subroutine piebor

  subroutine piecbk (callbk)
    implicit none

    interface
       subroutine callbk(iseg,xdat,xper,nrad,noff,a,nvx,nvy,idrw,iann)
         implicit none
         integer, intent (in) :: iseg
         double precision, intent (in) :: xdat,xper
         integer, intent (in out) :: nrad,noff,nvx,nvy,idrw,iann
         double precision, intent (in out) :: a
       end subroutine callbk
    end interface
    call doesnotexist("piecbk")
  end subroutine piecbk

  subroutine pieclr(ic1,ic2,n)
    implicit none
    integer, intent (in) :: n
    integer, dimension (n), intent (in) :: ic1,ic2
    call doesnotexist("pieclr")
  end subroutine pieclr

  subroutine pieexp()
    call doesnotexist("pieexp")
  end subroutine pieexp

  subroutine piegrf(cstr,nlin,xray,n)
    implicit none
    character (len = *), intent (in) :: cstr
    integer, intent (in) :: nlin,n
    double precision, dimension (n), intent (in) :: xray
    call doesnotexist("piegrf")
  end subroutine piegrf

  subroutine pielab(cdat,cstr)
    implicit none
    character(len = *), intent (in) :: cdat,cstr
    call doesnotexist("pielab")
  end subroutine pielab

  subroutine pieopt(x1,x2)
    implicit none
    double precision, intent (in) :: x1,x2
    call doesnotexist("pieopt")
  end subroutine pieopt

  subroutine pierot(xrot)
    implicit none
    double precision, intent (in) :: xrot
    call doesnotexist("pierot")
  end subroutine pierot

  subroutine pietyp(ctyp)
    implicit none
    character (len = *), intent (in) :: ctyp
    call doesnotexist("pietyp")
  end subroutine pietyp

  subroutine pieval(x,copt)
    implicit none
    double precision, intent (in) :: x
    character (len = *), intent (in) :: copt
    call doesnotexist("pieval")
  end subroutine pieval

  subroutine pievec(ivec,copt)
    implicit none
    integer, intent (in) :: ivec
    character (len = *), intent (in) :: copt
    call doesnotexist("pievec")
  end subroutine pievec

  subroutine pike3d(x1,y1,z1,x2,y2,z2,r,nsk1,nsk2)
    implicit none
    double precision, intent (in) :: x1,y1,z1,x2,y2,z2,r
    integer, intent (in) :: nsk1,nsk2
    call doesnotexist("pike3d")
  end subroutine pike3d

  subroutine plat3d(x,y,z,xl,copt)
    implicit none
    double precision, intent (in) :: x,y,z,xl
    character (len = *), intent (in) :: copt
    call doesnotexist("plat3d")
  end subroutine plat3d

  subroutine plyfin(cfl,cobj)
    implicit none
    character (len = *), intent (in) :: cfl,cobj
    call doesnotexist("plyfin")
  end subroutine plyfin

  subroutine plyini(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("plyini")
  end subroutine plyini

  subroutine pngmod(cmod,ckey)
    implicit none
    character (len = *), intent (in) :: cmod,ckey
    call doesnotexist("pngmod")
  end subroutine pngmod

  subroutine point(nx,ny,nb,nh,ncol)
    implicit none
    integer, intent (in) :: nx,ny,nb,nh,ncol
    call doesnotexist("point")
  end subroutine point

  subroutine polar(ex,orx,stepx,ory,stepy)
    implicit none
    double precision, intent (in) :: ex,orx,stepx,ory,stepy
    call doesnotexist("polar")
  end subroutine polar

  subroutine polclp(xray,yray,n,xout,yout,nmax,nout,xv,cedge)
    implicit none
    integer, intent (in) :: n,nmax
    integer, intent (out) :: nout
    double precision, dimension (n), intent (in) :: xray,yray
    double precision, dimension (nmax), intent (out) :: xout,yout
    double precision, intent (in) :: xv
    character (len = *), intent (in) :: cedge
    call doesnotexist("polclp")
  end subroutine polclp

  subroutine polcrv(cstr)
    implicit none
    character (len = *), intent (in) :: cstr
    call doesnotexist("polcrv")
  end subroutine polcrv

  subroutine polmod(cpos,cdir)
    implicit none
    character (len = *), intent (in) :: cpos,cdir
    call doesnotexist("polmod")
  end subroutine polmod

  subroutine pos2pt(x,y,xp,yp)
    implicit none
    double precision, intent (in) :: x,y
    double precision, intent (out) :: xp,yp
    call doesnotexist("pos2pt")
  end subroutine pos2pt

  subroutine pos3pt(x,y,z,xp,yp,zp)
    implicit none
    double precision, intent (in) :: x,y,z
    double precision, intent (out) :: xp,yp,zp
    call doesnotexist("pos3pt")
  end subroutine pos3pt

  subroutine posbar(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("posbar")
  end subroutine posbar

  subroutine posifl(nlu,nbyt,istat)
    implicit none
    integer, intent (in) :: nlu,nbyt
    integer, intent (out) :: istat
    call doesnotexist("posifl")
  end subroutine posifl

  subroutine proj3d(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("proj3d")
  end subroutine proj3d

  subroutine projct(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("projct")
  end subroutine projct

  subroutine psfont(cstr)
    implicit none
    character (len = *), intent (in) :: cstr
    call doesnotexist("psfont")
  end subroutine psfont

  subroutine psmode(cstr)
    implicit none
    character (len = *), intent (in) :: cstr
    call doesnotexist("psmode")
  end subroutine psmode

  subroutine pt2pos(x,y,xp,yp)
    implicit none
    double precision, intent (in) :: x,y
    double precision, intent (out) :: xp,yp
    call doesnotexist("pt2pos")
  end subroutine pt2pos

  subroutine pyra3d(x,y,z,xl,h1,h2,n)
    implicit none
    double precision, intent (in) :: x,y,z,xl,h1,h2
    integer, intent (in) :: n
    call doesnotexist("pyra3d")
  end subroutine pyra3d

  subroutine qplbar(x,n)
    implicit none
    integer, intent (in) :: n
    double precision, dimension (n), intent (in) :: x
    call doesnotexist("qplbar")
  end subroutine qplbar

  subroutine qplclr(x,n,m)
    implicit none
    integer, intent (in) :: n,m
    double precision, dimension (n,m), intent (in) :: x
    call doesnotexist("qplclr")
  end subroutine qplclr

  subroutine qplcon(x,n,m,nlev)
    implicit none
    integer, intent (in) :: n,m,nlev
    double precision, dimension (n,m), intent (in) :: x
    call doesnotexist("qplcon")
  end subroutine qplcon

  subroutine qplcrv(x,y,n,copt)
    implicit none
    integer, intent (in) :: n
    double precision, dimension (n), intent (in) :: x,y
    character (len = *), intent (in) :: copt
    call doesnotexist("qplcrv")
  end subroutine qplcrv

  subroutine qplot(x,y,n)
    implicit none
    integer, intent (in) :: n
    double precision, dimension (n), intent (in) :: x,y
    call doesnotexist("qplot")
  end subroutine qplot

  subroutine qplpie(x,n)
    implicit none
    integer, intent (in) :: n
    double precision, dimension (n), intent (in) :: x
    call doesnotexist("qplpie")
  end subroutine qplpie

  subroutine qplsca(x,y,n)
    implicit none
    integer, intent (in) :: n
    double precision, dimension (n), intent (in) :: x,y
    call doesnotexist("qplsca")
  end subroutine qplsca

  subroutine qplscl(a,e,or,step,copt)
    implicit none
    double precision, intent (in) :: a,e,or,step
    character (len = *), intent (in) :: copt
    call doesnotexist("qplscl")
  end subroutine qplscl

  subroutine qplsur(x,n,m)
    implicit none
    integer, intent (in) :: n,m
    double precision, dimension (n,m), intent (in) :: x
    call doesnotexist("qplsur")
  end subroutine qplsur

  subroutine quad3d(x,y,z,xl,yl,zl)
    implicit none
    double precision, intent (in) :: x,y,z,xl,yl,zl
    call doesnotexist("quad3d")
  end subroutine quad3d

  subroutine rbfpng(iray,nmax,nn)
    implicit none
    integer, intent (in) :: nmax
    character (len=1), intent (out), dimension (*) :: iray
    integer, intent (out) :: nn
    call doesnotexist("rbfpng")
  end subroutine rbfpng

  subroutine rbmp(cfl)
    implicit none
    character (len=*), intent (in) :: cfl
    call doesnotexist("rbmp")
  end subroutine rbmp

  subroutine readfl(nlu,iray,nbyt,istat)
    implicit none
    integer, intent (in) :: nlu,nbyt
    character (len=1), intent (out), dimension (nbyt) :: iray
    integer, intent (out) :: istat
    call doesnotexist("readfl")
  end subroutine readfl

  subroutine reawgt()
    call doesnotexist("reawgt")
  end subroutine reawgt

  subroutine recfll(nx,ny,nb,nh,ncol)
    implicit none
    integer, intent (in) :: nx,ny,nb,nh,ncol
    call doesnotexist("recfll")
  end subroutine recfll

  subroutine rectan(nx,ny,nb,nh)
    implicit none
    integer, intent (in) :: nx,ny,nb,nh
    call doesnotexist("rectan")
  end subroutine rectan

  subroutine rel3pt(x,y,z,xp,yp)
    implicit none
    double precision, intent (in) :: x,y,z
    double precision, intent (out) :: xp,yp
    call doesnotexist("rel3pt")
  end subroutine rel3pt

  subroutine resatt()
    call doesnotexist("resatt")
  end subroutine resatt

  subroutine reset(cw)
    implicit none
    character (len = *), intent (in) :: cw
    call doesnotexist("reset")
  end subroutine reset

  subroutine revscr()
    call doesnotexist("revscr")
  end subroutine revscr

  subroutine rgbhsv(r,g,b,xh,xs,xv)
    implicit none
    double precision, intent (in) :: r,g,b
    double precision, intent (out) :: xh,xs,xv
    call doesnotexist("rgbhsv")
  end subroutine rgbhsv

  subroutine rgif(cfl)
    implicit none
    character (len=*), intent (in) :: cfl
    call doesnotexist("rgif")
  end subroutine rgif

  subroutine rgtlab()
    call doesnotexist("rgtlab")
  end subroutine rgtlab

  subroutine rimage(cfl)
    implicit none
    character (len = *), intent (in) :: cfl
    call doesnotexist("rimage")
  end subroutine rimage

  subroutine rlarc(xm,ym,a,b,alpha,beta,theta)
    implicit none
    double precision, intent (in) :: xm,ym,a,b,alpha,beta,theta
    call doesnotexist("rlarc")
  end subroutine rlarc

  subroutine rlarea(x,y,n)
    implicit none
    integer, intent (in) :: n
    double precision, dimension (n), intent (in) :: x,y
    call doesnotexist("rlarea")
  end subroutine rlarea

  subroutine rlcirc(xm,ym,r)
    implicit none
    double precision, intent (in) :: xm,ym,r
    call doesnotexist("rlcirc")
  end subroutine rlcirc

  subroutine rlconn(x,y)
    implicit none
    double precision, intent (in) :: x,y
    call doesnotexist("rlconn")
  end subroutine rlconn

  subroutine rlell(xm,ym,a,b)
    implicit none
    double precision, intent (in) :: xm,ym,a,b
    call doesnotexist("rlell")
  end subroutine rlell

  subroutine rline(x,y,u,v)
    implicit none
    double precision, intent (in) :: x,y,u,v
    call doesnotexist("rline")
  end subroutine rline

  subroutine rlmess(cstr,x,y)
    implicit none
    character (len = *), intent (in) :: cstr
    double precision, intent (in) :: x,y
    call doesnotexist("rlmess")
  end subroutine rlmess

  subroutine rlnumb(z,ndez,x,y)
    implicit none
    double precision, intent (in) :: z,x,y
    integer, intent (in) :: ndez
    call doesnotexist("rlnumb")
  end subroutine rlnumb

  subroutine rlpie(xm,ym,r,alpha,beta)
    implicit none
    double precision, intent (in) :: xm,ym,r,alpha,beta
    call doesnotexist("rlpie")
  end subroutine rlpie

  subroutine rlpoin(x,y,nb,nh,ncol)
    implicit none
    double precision, intent (in) :: x,y
    integer, intent (in) :: nb,nh,ncol
    call doesnotexist("rlpoin")
  end subroutine rlpoin

  subroutine rlrec(x,y,xb,xh)
    implicit none
    double precision, intent (in) :: x,y,xb,xh
    call doesnotexist("rlrec")
  end subroutine rlrec

  subroutine rlrnd(x,y,xb,xh,irnd)
    implicit none
    double precision, intent (in) :: x,y,xb,xh
    integer, intent (in) :: irnd
    call doesnotexist("rlrnd")
  end subroutine rlrnd

  subroutine rlsec(xm,ym,r1,r,beta,alpha,ncol)
    implicit none
    double precision, intent (in) :: xm,ym,r1,r,beta,alpha
    integer, intent (in) :: ncol
    call doesnotexist("rlsec")
  end subroutine rlsec

  subroutine rlstrt(x,y)
    implicit none
    double precision, intent (in) :: x,y
    call doesnotexist("rlstrt")
  end subroutine rlstrt

  subroutine rlsymb(i,x,y)
    implicit none
    integer, intent (in) :: i
    double precision, intent (in) :: x,y
    call doesnotexist("rlsymb")
  end subroutine rlsymb

  subroutine rlvec(x,y,u,v,ivec)
    implicit none
    double precision, intent (in) :: x,y,u,v
    integer, intent (in) :: ivec
    call doesnotexist("rlvec")
  end subroutine rlvec

  subroutine rlwind(xk,x,y,nw,a)
    implicit none
    double precision, intent (in) :: xk,x,y,a
    integer, intent (in) :: nw
    call doesnotexist("rlwind")
  end subroutine rlwind

  subroutine rndrec(nx,ny,nb,nh,irnd)
    implicit none
    integer, intent (in) :: nx,ny,nb,nh,irnd
    call doesnotexist("rndrec")
  end subroutine rndrec

  subroutine rot3d(xa,ya,za)
    implicit none
    double precision, intent (in) :: xa,ya,za
    call doesnotexist("rot3d")
  end subroutine rot3d

  subroutine rpixel(ix,iy,n)
    implicit none
    integer, intent (in)  :: ix,iy
    integer, intent (out) :: n
    call doesnotexist("rpixel")
  end subroutine rpixel

  subroutine rpixls(iray,ix,iy,nw,nh)
    implicit none
    integer, intent (in) :: ix,iy,nw,nh
    character (len=1), intent (out), dimension (*) :: iray
    call doesnotexist("rpixls")
  end subroutine rpixls

  subroutine rpng(cfl)
    implicit none
    character (len=*), intent (in) :: cfl
    call doesnotexist("rpng")
  end subroutine rpng

  subroutine rppm(cfl)
    implicit none
    character (len=*), intent (in) :: cfl
    call doesnotexist("rppm")
  end subroutine rppm

  subroutine rpxrow(iray,ix,iy,n)
    implicit none
    integer, intent (in) :: ix,iy,n
    character (len=1), intent (out), dimension (n) :: iray
    call doesnotexist("rpxrow")
  end subroutine rpxrow

  subroutine rtiff(cfl)
    implicit none
    character (len=*), intent (in) :: cfl
    call doesnotexist("rtiff")
  end subroutine rtiff

  subroutine rvynam()
    call doesnotexist("rvynam")
  end subroutine rvynam

  subroutine scale(copt,cax)
    implicit none
    character (len = *), intent (in) :: copt,cax
    call doesnotexist("scale")
  end subroutine scale

  subroutine sclfac(x)
    implicit none
    double precision, intent (in) :: x
    call doesnotexist("sclfac")
  end subroutine sclfac

  subroutine sclmod(cmod)
    implicit none
    character (len = *), intent (in) :: cmod
    call doesnotexist("sclmod")
  end subroutine sclmod

  subroutine scmplx()
    call doesnotexist("scmplx")
  end subroutine scmplx

  subroutine scrmod(cmod)
    implicit none
    character (len = *), intent (in) :: cmod
    call doesnotexist("scrmod")
  end subroutine scrmod

  subroutine sector(nx,ny,nr1,nr2,alpha,beta,ncol)
    implicit none
    integer, intent (in) :: nx,ny,nr1,nr2,ncol
    double precision, intent (in) :: alpha,beta
    call doesnotexist("sector")
  end subroutine sector

  subroutine selwin(id)
    implicit none
    integer, intent (in) :: id
    call doesnotexist("selwin")
  end subroutine selwin

  subroutine sendbf()
    call doesnotexist("sendbf")
  end subroutine sendbf

  subroutine sendmb()
    call doesnotexist("sendmb")
  end subroutine sendmb

  subroutine sendok()
    call doesnotexist("sendok")
  end subroutine sendok

  subroutine serif()
    call doesnotexist("serif")
  end subroutine serif

  subroutine setbas(f)
    implicit none
    double precision, intent (in) :: f
    call doesnotexist("setbas")
  end subroutine setbas

  subroutine setcbk (callbk, copt)
    implicit none
    character (len = *), intent (in) :: copt

    interface
       subroutine callbk (xp, yp)
         implicit none
         double precision, intent (in out) :: xp,yp
       end subroutine callbk
    end interface
    call doesnotexist("setcbk")
  end subroutine setcbk

  subroutine setclr(n)
    implicit none
    integer, intent (in) :: n
    call doesnotexist("setclr")
  end subroutine setclr

  subroutine setcsr(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("setcsr")
  end subroutine setcsr

  subroutine setexp(f)
    implicit none
    double precision, intent (in) :: f
    call doesnotexist("setexp")
  end subroutine setexp

  subroutine setfce(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("setfce")
  end subroutine setfce

  subroutine setfil(ct)
    implicit none
    character (len = *), intent (in) :: ct
    call doesnotexist("setfil")
  end subroutine setfil

  subroutine setgrf(c1,c2,c3,c4)
    implicit none
    character (len = *), intent (in) :: c1,c2,c3,c4
    call doesnotexist("setgrf")
  end subroutine setgrf

  subroutine setind(i,xr,xg,xb)
    implicit none
    integer, intent (in) :: i
    double precision, intent (in) :: xr,xg,xb
    call doesnotexist("setind")
  end subroutine setind

  subroutine setmix(c,cstr)
    implicit none
    character (len = *), intent (in) :: c,cstr
    call doesnotexist("setmix")
  end subroutine setmix

  subroutine setpag(ct)
    implicit none
    character (len = *), intent (in) :: ct
    call doesnotexist("setpag")
  end subroutine setpag

  subroutine setres(i,j)
    implicit none
    integer, intent (in) :: i,j
    call doesnotexist("setres")
  end subroutine setres

  subroutine setrgb(xr,xg,xb)
    implicit none
    double precision, intent (in) :: xr,xg,xb
    call doesnotexist("setrgb")
  end subroutine setrgb

  subroutine setscl(xray,n,cstr)
    implicit none
    integer, intent (in) :: n
    double precision, dimension (n), intent (in) :: xray
    character (len = *), intent (in) :: cstr
    call doesnotexist("setscl")
  end subroutine setscl

  subroutine setvlt(ctab)
    implicit none
    character (len = *), intent (in) :: ctab
    call doesnotexist("setvlt")
  end subroutine setvlt

  subroutine setxid(i,copt)
    implicit none
    integer, intent (in) :: i
    character (len = *), intent (in) :: copt
    call doesnotexist("setxid")
  end subroutine setxid

  subroutine shdafr(inat,ishd,iclr,n)
    implicit none
    integer, intent (in) :: n
    integer, dimension (n), intent (in) :: inat,ishd,iclr
    call doesnotexist("shdafr")
  end subroutine shdafr

  subroutine shdasi(inat,ishd,iclr,n)
    implicit none
    integer, intent (in) :: n
    integer, dimension (n), intent (in) :: inat,ishd,iclr
    call doesnotexist("shdasi")
  end subroutine shdasi

  subroutine shdaus(inat,ishd,iclr,n)
    implicit none
    integer, intent (in) :: n
    integer, dimension (n), intent (in) :: inat,ishd,iclr
    call doesnotexist("shdaus")
  end subroutine shdaus

  subroutine shdcha()
    call doesnotexist("shdcha")
  end subroutine shdcha

  subroutine shdcrv(x1,y1,n1,x2,y2,n2)
    implicit none
    integer, intent (in) :: n1,n2
    double precision, dimension (n1), intent (in) :: x1,y1
    double precision, dimension (n2), intent (in) :: x2,y2
    call doesnotexist("shdcrv")
  end subroutine shdcrv

  subroutine shdeur(inat,ishd,iclr,n)
    implicit none
    integer, intent (in) :: n
    integer, dimension (n), intent (in) :: inat,ishd,iclr
    call doesnotexist("shdeur")
  end subroutine shdeur

  subroutine shdfac(x)
    implicit none
    double precision, intent (in) :: x
    call doesnotexist("shdfac")
  end subroutine shdfac

  subroutine shdmap(cmap)
    implicit none
    character (len = *), intent (in) :: cmap
    call doesnotexist("shdmap")
  end subroutine shdmap

  subroutine shdmod(copt,ctype)
    implicit none
    character (len = *), intent (in) :: copt,ctype
    call doesnotexist("shdmod")
  end subroutine shdmod

  subroutine shdnor(inat,ishd,iclr,n)
    implicit none
    integer, intent (in) :: n
    integer, dimension (n), intent (in) :: inat,ishd,iclr
    call doesnotexist("shdnor")
  end subroutine shdnor

  subroutine shdpat(i)
    implicit none
    integer, intent (in) :: i
    call doesnotexist("shdpat")
  end subroutine shdpat

  subroutine shdsou(inat,ishd,iclr,n)
    implicit none
    integer, intent (in) :: n
    integer, dimension (n), intent (in) :: inat,ishd,iclr
    call doesnotexist("shdsou")
  end subroutine shdsou

  subroutine shdusa(inat,ishd,iclr,n)
    implicit none
    integer, intent (in) :: n
    integer, dimension (n), intent (in) :: inat,ishd,iclr
    call doesnotexist("shdusa")
  end subroutine shdusa

  subroutine shield(cblnk,cmode)
    implicit none
    character (len = *), intent (in) :: cblnk,cmode
    call doesnotexist("shield")
  end subroutine shield

  subroutine shlcir(nx,ny,nr)
    implicit none
    integer, intent (in) :: nx,ny,nr
    call doesnotexist("shlcir")
  end subroutine shlcir

  subroutine shldel(id)
    implicit none
    integer, intent (in) :: id
    call doesnotexist("shldel")
  end subroutine shldel

  subroutine shlell(nx,ny,na,nb,ang)
    implicit none
    integer, intent (in) :: nx,ny,na,nb
    double precision, intent (in) :: ang
    call doesnotexist("shlell")
  end subroutine shlell

  subroutine shlind(id)
    implicit none
    integer, intent (out) :: id
    call doesnotexist("shlind")
  end subroutine shlind

  subroutine shlpie(nx,ny,nr,alph,beta)
    implicit none
    integer, intent (in) :: nx,ny,nr
    double precision, intent (in) :: alph,beta
    call doesnotexist("shlpie")
  end subroutine shlpie

  subroutine shlpol(nxray,nyray,n)
    implicit none
    integer, intent (in) :: n
    integer, dimension (n), intent (in) :: nxray,nyray
    call doesnotexist("shlpol")
  end subroutine shlpol

  subroutine shlrct(nx,ny,nw,nh,ang)
    implicit none
    integer, intent (in) :: nx,ny,nw,nh
    double precision,    intent (in) :: ang
    call doesnotexist("shlrct")
  end subroutine shlrct

  subroutine shlrec(nx,ny,nw,nh)
    implicit none
    integer, intent (in) :: nx,ny,nw,nh
    call doesnotexist("shlrec")
  end subroutine shlrec

  subroutine shlres(nn)
    implicit none
    integer, intent (in) :: nn
    call doesnotexist("shlres")
  end subroutine shlres

  subroutine shlsur()
    call doesnotexist("shlsur")
  end subroutine shlsur

  subroutine shlvis(id,cvis)
    implicit none
    integer, intent (in) :: id
    character (len = *), intent (in) :: cvis
    call doesnotexist("shlvis")
  end subroutine shlvis

  subroutine simplx()
    call doesnotexist("simplx")
  end subroutine simplx

  subroutine skipfl(nlu,nbyt,istat)
    implicit none
    integer, intent (in) :: nlu,nbyt
    integer, intent (out) :: istat
    call doesnotexist("skipfl")
  end subroutine skipfl

  subroutine smxalf(calph,ca,cb,n)
    implicit none
    character (len = *), intent (in) :: calph,ca,cb
    integer, intent (in) :: n
    call doesnotexist("smxalf")
  end subroutine smxalf

  subroutine solid()
    call doesnotexist("solid")
  end subroutine solid

  subroutine sortr1(x,n,copt)
    implicit none
    integer, intent (in) :: n
    double precision, dimension (n), intent (in out) :: x
    character (len = *), intent (in) :: copt
    call doesnotexist("sortr1")
  end subroutine sortr1

  subroutine sortr2(x,y,n,copt)
    implicit none
    integer, intent (in) :: n
    double precision, dimension (n), intent (in out) :: x,y
    character (len = *), intent (in) :: copt
    call doesnotexist("sortr2")
  end subroutine sortr2

  subroutine spcbar(i)
    implicit none
    integer, intent (in) :: i
    call doesnotexist("spcbar")
  end subroutine spcbar

  subroutine sphe3d(xm,ym,zm,r,nsk1,nsk2)
    implicit none
    double precision, intent (in) :: xm,ym,zm,r
    integer, intent (in) :: nsk1,nsk2
    call doesnotexist("sphe3d")
  end subroutine sphe3d

  subroutine spline(x,y,n,xray,yray,ndat)
    implicit none
    integer, intent (in) :: n
    double precision, dimension (n),intent (in) :: x,y
    double precision, dimension (n),intent (out) :: xray,yray
    integer, intent (out) :: ndat
    call doesnotexist("spline")
  end subroutine spline

  subroutine splmod(k,n)
    implicit none
    integer, intent (in) :: k,n
    call doesnotexist("splmod")
  end subroutine splmod

  subroutine stmmod(cmod,ckey)
    implicit none
    character (len = *), intent (in) :: cmod,ckey
    call doesnotexist("stmmod")
  end subroutine stmmod

  subroutine stmopt(n,copt)
    implicit none
    integer, intent (in) :: n
    character (len = *), intent (in) :: copt
    call doesnotexist("stmopt")
  end subroutine stmopt

  subroutine stmpts(xmat,ymat,nx,ny,xp,yp,x0,y0,xray,yray,nmax,nray)
    implicit none
    integer, intent (in) :: nx,ny,nmax
    integer, intent (out) :: nray
    double precision, dimension (nx,ny), intent (in) :: xmat,ymat
    double precision, dimension (nx), intent (in) :: xp
    double precision, dimension (ny), intent (in) :: yp
    double precision, intent (in) :: x0,y0
    double precision, dimension (nmax),  intent (out) :: xray,yray
    call doesnotexist("stmpts")
  end subroutine stmpts

  subroutine stmpts3d(xv,yv,zv,nx,ny,nz,xp,yp,zp,x0,y0,z0, &
                      xray,yray,zray,nmax,nray)
    implicit none
    integer, intent (in) :: nx,ny,nz,nmax
    integer, intent (out) :: nray
    double precision, dimension (nx,ny,nz), intent (in) :: xv,yv,zv
    double precision, dimension (nx), intent (in) :: xp
    double precision, dimension (ny), intent (in) :: yp
    double precision, dimension (nz), intent (in) :: zp
    double precision, intent (in) :: x0,y0,z0
    double precision, dimension (nmax),  intent (out) :: xray,yray,zray
    call doesnotexist("stmpts3d")
  end subroutine stmpts3d

  subroutine stmtri(xvray,yvray,xpray,ypray,n, &
                    i1ray,i2ray,i3ray,ntri,xs,ys,nray)
    implicit none
    integer, intent (in) :: n,ntri,nray
    double precision, dimension (n), intent (in) :: xvray,yvray,xpray,ypray
    integer, dimension (ntri), intent (in) :: i1ray,i2ray,i3ray
    double precision, dimension (nray),  intent (in) :: xs,ys
    call doesnotexist("stmtri")
  end subroutine stmtri

  subroutine stmval(x,copt)
    implicit none
    double precision, intent (in) :: x
    character (len = *), intent (in) :: copt
    call doesnotexist("stmval")
  end subroutine stmval

  subroutine stream(xmat,ymat,nx,ny,xp,yp,xs,ys,n)
    implicit none
    integer, intent (in) :: nx,ny,n
    double precision, dimension (nx,ny), intent (in) :: xmat,ymat
    double precision, dimension (nx), intent (in) :: xp
    double precision, dimension (ny), intent (in) :: yp
    double precision, dimension (n),  intent (in) :: xs,ys
    call doesnotexist("stream")
  end subroutine stream

  subroutine stream3d(xv,yv,zv,nx,ny,nz,xp,yp,zp,xs,ys,zs,n)
    implicit none
    integer, intent (in) :: nx,ny,nz,n
    double precision, dimension (nx,ny,nz), intent (in) :: xv,yv,zv
    double precision, dimension (nx), intent (in) :: xp
    double precision, dimension (ny), intent (in) :: yp
    double precision, dimension (nz), intent (in) :: zp
    double precision, dimension (n),  intent (in) :: xs,ys,zs
    call doesnotexist("stream3d")
  end subroutine stream3d

  subroutine strt3d(x,y,z)
    implicit none
    double precision, intent (in) :: x,y,z
    call doesnotexist("strt3d")
  end subroutine strt3d

  subroutine strtpt(x,y)
    implicit none
    double precision, intent (in) :: x,y
    call doesnotexist("strtpt")
  end subroutine strtpt

  subroutine surclr(itop,ibot)
    implicit none
    integer, intent (in) :: itop,ibot
    call doesnotexist("surclr")
  end subroutine surclr

  subroutine surfce(xray,ixdim,yray,iydim,zmat)
    implicit none
    integer, intent (in) :: ixdim,iydim
    double precision, dimension (ixdim), intent (in) :: xray
    double precision, dimension (iydim), intent (in) :: yray
    double precision, dimension (ixdim,iydim), intent (in) :: zmat
    call doesnotexist("surfce")
  end subroutine surfce

  subroutine surfcp(zfun,a1,a2,astp,b1,b2,bstp)
    implicit none
    interface
       function zfun(x,y,iopt)
         implicit none
         double precision, intent (in) :: x,y
         integer, intent (in) :: iopt
         double precision :: zfun
       end function zfun
    end interface
    double precision, intent (in) :: a1,a2,astp,b1,b2,bstp
    call doesnotexist("surfcp")
  end subroutine surfcp

  subroutine suriso(xray,nx,yray,ny,zray,nz,wmat,wlev)
    implicit none
    integer, intent (in) :: nx,ny,nz
    double precision, dimension (nx), intent (in) :: xray
    double precision, dimension (ny), intent (in) :: yray
    double precision, dimension (nz), intent (in) :: zray
    double precision, dimension (nx,ny,nz), intent (in) :: wmat
    double precision, intent (in) :: wlev
    call doesnotexist("suriso")
  end subroutine suriso

  subroutine surfun(zfun,ixpts,xdel,iypts,ydel)
    implicit none
    interface
       function zfun(x,y)
         implicit none
         double precision, intent (in) :: x,y
         double precision :: zfun
       end function zfun
    end interface
    integer, intent (in) :: ixpts,iypts
    double precision, intent (in) :: xdel,ydel
    call doesnotexist("surfun")
  end subroutine surfun

  subroutine surmat(zmat,ixdim,iydim,ixp,iyp)
    implicit none
    integer, intent (in) :: ixdim,iydim,ixp,iyp
    double precision, dimension (ixdim,iydim), intent (in) :: zmat
    call doesnotexist("surmat")
  end subroutine surmat

  subroutine surmsh(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("surmsh")
  end subroutine surmsh

  subroutine suropt(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("suropt")
  end subroutine suropt

  subroutine surshc(xray,ixdim,yray,iydim,zmat,wmat)
    implicit none
    integer, intent (in) :: ixdim,iydim
    double precision, dimension (ixdim), intent (in) :: xray
    double precision, dimension (iydim), intent (in) :: yray
    double precision, dimension (ixdim,iydim), intent (in) :: zmat,wmat
    call doesnotexist("surshc")
  end subroutine surshc

  subroutine surshd(xray,ixdim,yray,iydim,zmat)
    implicit none
    integer, intent (in) :: ixdim,iydim
    double precision, dimension (ixdim), intent (in) :: xray
    double precision, dimension (iydim), intent (in) :: yray
    double precision, dimension (ixdim,iydim), intent (in) :: zmat
    call doesnotexist("surshd")
  end subroutine surshd

  subroutine sursze(ax,ex,ay,ey)
    implicit none
    double precision, intent (in) :: ax,ex,ay,ey
    call doesnotexist("sursze")
  end subroutine sursze

  subroutine surtri(xray,yray,zray,n,i1ray,i2ray,i3ray,ntri)
    implicit none
    integer, intent (in) :: n,ntri
    double precision, dimension (n), intent (in) :: xray,yray,zray
    integer, dimension (ntri), intent (in) :: i1ray,i2ray,i3ray
    call doesnotexist("surtri")
  end subroutine surtri

  subroutine survis(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("survis")
  end subroutine survis

  subroutine swapi2(iray,n)
    implicit none
    integer, intent (in) :: n
    integer (kind=selected_int_kind(4)), dimension (n), &
            intent (in out) :: iray
    call doesnotexist("swapi2")
  end subroutine swapi2

  subroutine swapi4(iray,n)
    implicit none
    integer, intent (in) :: n
    integer (kind=selected_int_kind(9)), dimension (n), &
            intent (in out) :: iray
    call doesnotexist("swapi4")
  end subroutine swapi4

  subroutine swgatt (id,cval,copt)
    implicit none
    integer, intent (in) :: id
    character (len=*), intent (in) :: cval,copt
    call doesnotexist("swgatt")
  end subroutine swgatt

  subroutine swgbgd(id,xr,xg,xb)
    implicit none
    integer, intent (in) :: id
    double precision, intent (in) :: xr,xg,xb
    call doesnotexist("swgbgd")
  end subroutine swgbgd

  subroutine swgbox(id,ival)
    implicit none
    integer, intent (in) :: id,ival
    call doesnotexist("swgbox")
  end subroutine swgbox

  subroutine swgbut(id,ival)
    implicit none
    integer, intent (in) :: id,ival
    call doesnotexist("swgbut")
  end subroutine swgbut

  subroutine swgcb (id, callbk, iray)
    implicit none
    integer, intent (in) :: id
    integer, dimension (1), intent (in) :: iray

    interface
       subroutine callbk (id,iray)
         implicit none
         integer, intent (in) :: id
         integer, dimension (1), intent (in) :: iray
       end subroutine callbk
    end interface
    call doesnotexist("swgcb")
  end subroutine swgcb

  subroutine swgcb2 (id, callbk)
    implicit none
    integer, intent (in) :: id

    interface
       subroutine callbk (id,irow,icol)
         implicit none
         integer, intent (in) :: id,irow,icol
       end subroutine callbk
    end interface
    call doesnotexist("swgcb2")
  end subroutine swgcb2

  subroutine swgcb3 (id, callbk)
    implicit none
    integer, intent (in) :: id

    interface
       subroutine callbk (id,ival)
         implicit none
         integer, intent (in) :: id,ival
       end subroutine callbk
    end interface
    call doesnotexist("swgcb3")
  end subroutine swgcb3

  subroutine swgcbk (id, callbk)
    implicit none
    integer, intent (in) :: id

    interface
       subroutine callbk (id)
         implicit none
         integer, intent (in) :: id
       end subroutine callbk
    end interface
    call doesnotexist("swgcbk")
  end subroutine swgcbk

  subroutine swgclr(xr,xg,xb,copt)
    implicit none
    double precision, intent (in) :: xr,xg,xb
    character (len=*), intent (in) :: copt
    call doesnotexist("swgclr")
  end subroutine swgclr

  subroutine swgdrw(x)
    implicit none
    double precision, intent (in) :: x
    call doesnotexist("swgdrw")
  end subroutine swgdrw

  subroutine swgfgd(id,xr,xg,xb)
    implicit none
    integer, intent (in) :: id
    double precision, intent (in) :: xr,xg,xb
    call doesnotexist("swgfgd")
  end subroutine swgfgd

  subroutine swgfil(id,cstr)
    implicit none
    integer, intent (in) :: id
    character (len=*), intent (in) :: cstr
    call doesnotexist("swgfil")
  end subroutine swgfil

  subroutine swgflt(id,xval,ndig)
    implicit none
    integer, intent (in) :: id,ndig
    double precision, intent (in) :: xval
    call doesnotexist("swgflt")
  end subroutine swgflt

  subroutine swgfnt(cstr,n)
    implicit none
    character (len=*), intent (in) :: cstr
    integer, intent (in) :: n
    call doesnotexist("swgfnt")
  end subroutine swgfnt

  subroutine swgfoc(id)
    implicit none
    integer, intent (in) :: id
    call doesnotexist("swgfoc")
  end subroutine swgfoc

  subroutine swghlp(cstr)
    implicit none
    character (len=*), intent (in) :: cstr
    call doesnotexist("swghlp")
  end subroutine swghlp

  subroutine swgint(id,iv)
    implicit none
    integer, intent (in) :: id,iv
    call doesnotexist("swgint")
  end subroutine swgint

  subroutine swgiop (ival,copt)
    implicit none
    integer, intent (in) :: ival
    character (len=*), intent (in) :: copt
    call doesnotexist("swgiop")
  end subroutine swgiop

  subroutine swgjus (ctype,cwidg)
    implicit none
    character (len=*), intent (in) :: ctype,cwidg
    call doesnotexist("swgjus")
  end subroutine swgjus

  subroutine swglis(id,ival)
    implicit none
    integer, intent (in) :: id,ival
    call doesnotexist("swglis")
  end subroutine swglis

  subroutine swgmix(c,cstr)
    implicit none
    character (len=*), intent (in) :: c,cstr
    call doesnotexist("swgmix")
  end subroutine swgmix

  subroutine swgmod(cmod)
    implicit none
    character (len=*), intent (in) :: cmod
    call doesnotexist("swgmod")
  end subroutine swgmod

  subroutine swgmrg(ival,cstr)
    implicit none
    integer, intent (in) :: ival
    character (len=*), intent (in) :: cstr
    call doesnotexist("swgmrg")
  end subroutine swgmrg

  subroutine swgoff(nx,ny)
    implicit none
    integer, intent (in) :: nx,ny
    call doesnotexist("swgoff")
  end subroutine swgoff

  subroutine swgopt (cval,copt)
    implicit none
    character (len=*), intent (in) :: cval,copt
    call doesnotexist("swgopt")
  end subroutine swgopt

  subroutine swgpop (copt)
    implicit none
    character (len=*), intent (in) :: copt
    call doesnotexist("swgpop")
  end subroutine swgpop

  subroutine swgpos(nx,ny)
    implicit none
    integer, intent (in) :: nx,ny
    call doesnotexist("swgpos")
  end subroutine swgpos

  subroutine swgray(xray,n,copt)
    implicit none
    integer, intent (in) :: n
    double precision, dimension (n), intent (in) :: xray
    character (len=*), intent (in) :: copt
    call doesnotexist("swgray")
  end subroutine swgray

  subroutine swgscl(id,xval)
    implicit none
    integer, intent (in) :: id
    double precision, intent (in) :: xval
    call doesnotexist("swgscl")
  end subroutine swgscl

  subroutine swgsiz(nx,ny)
    implicit none
    integer, intent (in) :: nx,ny
    call doesnotexist("swgsiz")
  end subroutine swgsiz

  subroutine swgspc(x, y)
    implicit none
    double precision, intent (in) :: x,y
    call doesnotexist("swgspc")
  end subroutine swgspc

  subroutine swgstp(x)
    implicit none
    double precision, intent (in) :: x
    call doesnotexist("swgstp")
  end subroutine swgstp

  subroutine swgtbf(id,xval,ndig,irow,icol,copt)
    implicit none
    integer, intent (in) :: id,ndig,irow,icol
    double precision, intent (in) :: xval
    character (len=*), intent (in) :: copt
    call doesnotexist("swgtbf")
  end subroutine swgtbf

  subroutine swgtbi(id,ival,irow,icol,copt)
    implicit none
    integer, intent (in) :: id,ival,irow,icol
    character (len=*), intent (in) :: copt
    call doesnotexist("swgtbi")
  end subroutine swgtbi

  subroutine swgtbl(id,xray,n,ndig,idx,copt)
    implicit none
    integer, intent (in) :: id,n,ndig,idx
    double precision, dimension (n), intent (in) :: xray
    character (len=*), intent (in) :: copt
    call doesnotexist("swgtbl")
  end subroutine swgtbl

  subroutine swgtbs(id,cstr,irow,icol,copt)
    implicit none
    integer, intent (in) :: id,irow,icol
    character (len=*), intent (in) :: cstr,copt
    call doesnotexist("swgtbs")
  end subroutine swgtbs

  subroutine swgtit(cstr)
    implicit none
    character (len=*), intent (in) :: cstr
    call doesnotexist("swgtit")
  end subroutine swgtit

  subroutine wgtbl(ip,nrows,ncols,id)
    implicit none
    integer, intent (in)  :: ip,nrows,ncols
    integer, intent (out) :: id
    call doesnotexist("wgtbl")
  end subroutine wgtbl

  subroutine swgtxt(id,cstr)
    implicit none
    integer, intent (in) :: id
    character (len=*), intent (in) :: cstr
    call doesnotexist("swgtxt")
  end subroutine swgtxt

  subroutine swgtyp (ctype,cwidg)
    implicit none
    character (len=*), intent (in) :: ctype,cwidg
    call doesnotexist("swgtyp")
  end subroutine swgtyp

  subroutine swgval(id,xval)
    implicit none
    integer, intent (in) :: id
    double precision, intent (in) :: xval
    call doesnotexist("swgval")
  end subroutine swgval

  subroutine swgwin(nx,ny,nw,nh)
    implicit none
    integer, intent (in) :: nx,ny,nw,nh
    call doesnotexist("swgwin")
  end subroutine swgwin

  subroutine swgwth (nwth)
    implicit none
    integer, intent (in) :: nwth
    call doesnotexist("swgwth")
  end subroutine swgwth

  subroutine symb3d(i,x,y,z)
    implicit none
    integer, intent (in) :: i
    double precision, intent (in) :: x,y,z
    call doesnotexist("symb3d")
  end subroutine symb3d

  subroutine symbol(i,nx,ny)
    implicit none
    integer, intent (in) :: i,nx,ny
    call doesnotexist("symbol")
  end subroutine symbol

  subroutine symfil(cdv,cst)
    implicit none
    character (len=*), intent (in) :: cdv,cst
    call doesnotexist("symfil")
  end subroutine symfil

  subroutine symrot(xrot)
    implicit none
    double precision, intent (in) :: xrot
    call doesnotexist("symrot")
  end subroutine symrot

  subroutine tellfl(nlu,nbyt)
    implicit none
    integer, intent (in) :: nlu
    integer, intent (out) :: nbyt
    call doesnotexist("tellfl")
  end subroutine tellfl

  subroutine texmod(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("texmod")
  end subroutine texmod

  subroutine texopt(copt,ctype)
    implicit none
    character (len = *), intent (in) :: copt,ctype
    call doesnotexist("texopt")
  end subroutine texopt

  subroutine texval(x,copt)
    implicit none
    double precision, intent (in) :: x
    character (len = *), intent (in) :: copt
    call doesnotexist("texval")
  end subroutine texval

  subroutine thkc3d(x)
    implicit none
    double precision, intent (in) :: x
    call doesnotexist("thkc3d")
  end subroutine thkc3d

  subroutine thkcrv(i)
    implicit none
    integer, intent (in) :: i
    call doesnotexist("thkcrv")
  end subroutine thkcrv

  subroutine thrfin()
    implicit none
    call doesnotexist("thrfin")
  end subroutine thrfin

  subroutine thrini(i)
    implicit none
    integer, intent (in) :: i
    call doesnotexist("thrini")
  end subroutine thrini

  subroutine ticks(i,cax)
    implicit none
    integer, intent (in) :: i
    character (len = *), intent (in) :: cax
    call doesnotexist("ticks")
  end subroutine ticks

  subroutine ticlen(i1,i2)
    implicit none
    integer, intent (in) :: i1,i2
    call doesnotexist("ticlen")
  end subroutine ticlen

  subroutine ticmod(copt,cax)
    implicit none
    character (len = *), intent (in) :: copt,cax
    call doesnotexist("ticmod")
  end subroutine ticmod

  subroutine ticpos(copt,cax)
    implicit none
    character (len = *), intent (in) :: copt,cax
    call doesnotexist("ticpos")
  end subroutine ticpos

  subroutine tifmod(n,cval,copt)
    implicit none
    character (len = *), intent (in) :: cval,copt
    integer, intent (in) :: n
    call doesnotexist("tifmod")
  end subroutine tifmod

  subroutine tiforg(nx,ny)
    implicit none
    integer, intent (in) :: nx,ny
    call doesnotexist("tiforg")
  end subroutine tiforg

  subroutine tifwin(nx,ny,nw,nh)
    implicit none
    integer, intent (in) :: nx,ny,nw,nh
    call doesnotexist("tifwin")
  end subroutine tifwin

  subroutine timopt()
    call doesnotexist("timopt")
  end subroutine timopt

  subroutine titjus(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("titjus")
  end subroutine titjus

  subroutine title()
    call doesnotexist("title")
  end subroutine title

  subroutine titlin(cstr,j)
    implicit none
    character (len = *), intent (in) :: cstr
    integer, intent (in) :: j
    call doesnotexist("titlin")
  end subroutine titlin

  subroutine titpos(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("titpos")
  end subroutine titpos

  subroutine torus3d(xm,ym,zm,r1,r2,h,a1,a2,nsk1,nsk2)
    implicit none
    double precision, intent (in) :: xm,ym,zm,r1,r2,h,a1,a2
    integer, intent (in) :: nsk1,nsk2
    call doesnotexist("torus3d")
  end subroutine torus3d

  subroutine tprfin()
    call doesnotexist("tprfin")
  end subroutine tprfin

  subroutine tprini()
    call doesnotexist("tprini")
  end subroutine tprini

  subroutine tprmod(copt,ckey)
    implicit none
    character (len = *), intent (in) :: copt,ckey
    call doesnotexist("tprmod")
  end subroutine tprmod

  subroutine tprval(x)
    implicit none
    double precision, intent (in) :: x
    call doesnotexist("tprval")
  end subroutine tprval

  subroutine tr3axs(x,y,z,a)
    implicit none
    double precision, intent (in) :: x,y,z,a
    call doesnotexist("tr3axs")
  end subroutine tr3axs

  subroutine tr3res()
    call doesnotexist("tr3res")
  end subroutine tr3res

  subroutine tr3rot(a,b,c)
    implicit none
    double precision, intent (in) :: a,b,c
    call doesnotexist("tr3rot")
  end subroutine tr3rot

  subroutine tr3scl(xscl,yscl,zscl)
    implicit none
    double precision, intent (in) :: xscl,yscl,zscl
    call doesnotexist("tr3scl")
  end subroutine tr3scl

  subroutine tr3shf(xshf,yshf,zshf)
    implicit none
    double precision, intent (in) :: xshf,yshf,zshf
    call doesnotexist("tr3shf")
  end subroutine tr3shf

  subroutine trfco1(xray,n,cfrom,cto)
    implicit none
    integer, intent (in) :: n
    double precision, dimension (n), intent (in out) :: xray
    character (len = *), intent(in) :: cfrom, cto
    call doesnotexist("trfco1")
  end subroutine trfco1

  subroutine trfco2(xray,yray,n,cfrom,cto)
    implicit none
    integer, intent (in) :: n
    double precision, dimension (n), intent (in out) :: xray,yray
    character (len = *), intent(in) :: cfrom, cto
    call doesnotexist("trfco2")
  end subroutine trfco2

  subroutine trfco3(xray,yray,zray,n,cfrom,cto)
    implicit none
    integer, intent (in) :: n
    double precision, dimension (n), intent (in out) :: xray,yray,zray
    character (len = *), intent(in) :: cfrom, cto
    call doesnotexist("trfco3")
  end subroutine trfco3

  subroutine trfdat(ndays,id,im,iy)
    implicit none
    integer, intent (in) :: ndays
    integer, intent (out) :: id,im,iy
    call doesnotexist("trfdat")
  end subroutine trfdat

  subroutine trfmat(zmat,nx,ny,zmat2,nx2,ny2)
    implicit none
    integer, intent (in) :: nx,ny,nx2,ny2
    double precision, dimension (nx,ny), intent (in) :: zmat
    double precision, dimension (nx2,ny2), intent (out) :: zmat2
    call doesnotexist("trfmat")
  end subroutine trfmat

  subroutine trfrel(x,y,n)
    implicit none
    integer, intent (in) :: n
    double precision, dimension (n), intent (in out) :: x,y
    call doesnotexist("trfrel")
  end subroutine trfrel

  subroutine trfres()
    call doesnotexist("trfres")
  end subroutine trfres

  subroutine trfrot(xang,nx,ny)
    implicit none
    double precision, intent (in) :: xang
    integer, intent (in) :: nx,ny
    call doesnotexist("trfrot")
  end subroutine trfrot

  subroutine trfscl(xscl,yscl)
    implicit none
    double precision, intent (in) :: xscl,yscl
    call doesnotexist("trfscl")
  end subroutine trfscl

  subroutine trfshf(nxshf,nyshf)
    implicit none
    integer, intent (in) :: nxshf,nyshf
    call doesnotexist("trfshf")
  end subroutine trfshf

  subroutine tria3d(x,y,z)
    implicit none
    double precision, dimension (3), intent (in) :: x,y,z
    call doesnotexist("tria3d")
  end subroutine tria3d

  subroutine triang(xray,yray,n,i1ray,i2ray,i3ray,nmax,ntri)
    implicit none
    integer, intent (in) :: n,nmax
    double precision, dimension (n), intent (in) :: xray,yray
    integer, dimension (nmax), intent (out) :: i1ray,i2ray,i3ray
    integer, intent (out) :: ntri
    call doesnotexist("triang")
  end subroutine triang

  subroutine triflc(xray,yray,iray,n)
    implicit none
    integer, intent (in) :: n
    double precision, dimension (n), intent (in) :: xray,yray
    integer, dimension (n), intent (in) :: iray
    call doesnotexist("triflc")
  end subroutine triflc

  subroutine trifll(x,y)
    implicit none
    double precision, dimension (3), intent (in) :: x,y
    call doesnotexist("trifll")
  end subroutine trifll

  subroutine triplx()
    call doesnotexist("triplx")
  end subroutine triplx

  subroutine tripts(xray,yray,zray,n,i1ray,i2ray,i3ray,ntri,zlev, &
                    xpts, ypts, maxpts, iray, maxray, nlins)
    implicit none
    integer, intent (in) :: n,ntri,maxpts,maxray
    double precision, dimension (n), intent (in) :: xray,yray,zray
    double precision, dimension (maxpts), intent (out) :: xpts,ypts
    integer, dimension (ntri), intent (in) :: i1ray,i2ray,i3ray
    integer, dimension (maxray), intent (out) :: iray
    double precision, intent (in) :: zlev
    integer, intent (out) :: nlins
    call doesnotexist("tripts")
  end subroutine tripts

  function trmlen(cstr)
    implicit none
    character (len = *), intent (in) :: cstr
    double precision :: trmlen
    call doesnotexist("trmlen")
  end function trmlen

  subroutine ttfont(cfnt)
    implicit none
    character (len = *), intent (in) :: cfnt
    call doesnotexist("ttfont")
  end subroutine ttfont

  subroutine tube3d(x1,y1,z1,x2,y2,z2,r,nsk1,nsk2)
    implicit none
    double precision, intent (in) :: x1,y1,z1,x2,y2,z2,r
    integer, intent (in) :: nsk1,nsk2
    call doesnotexist("tube3d")
  end subroutine tube3d

  subroutine txtbgd(i)
    implicit none
    integer, intent (in):: i
    call doesnotexist("txtbgd")
  end subroutine txtbgd

  subroutine txtjus(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("txtjus")
  end subroutine txtjus

  subroutine txture(itmat,nx,ny)
    implicit none
    integer, intent (in) :: nx,ny
    integer, dimension (nx,ny), intent (out) :: itmat
    call doesnotexist("txture")
  end subroutine txture

  subroutine unit(i)
    implicit none
    integer, intent (in) :: i
    call doesnotexist("unit")
  end subroutine unit

  subroutine units(cmod)
    implicit none
    character (len = *), intent (in) :: cmod
    call doesnotexist("units")
  end subroutine units

  subroutine upstr(cstr)
    implicit none
    character (len = *), intent (in out) :: cstr
    call doesnotexist("upstr")
  end subroutine upstr

  subroutine usrpie(iseg,xdat,xper,nrad,noff,ang,nvbox,idrw,iann)
    implicit none
    integer, intent (in out) :: iseg,nrad,noff,nvbox,idrw,iann
    double precision,    intent (in out) :: xdat,xper,ang
    call doesnotexist("usrpie")
  end subroutine usrpie

  subroutine utfint(cstr,iray,n,nl)
    implicit none
    character (len=*), intent (in) :: cstr
    integer, intent (in) :: n
    integer, dimension (n), intent (out) :: iray
    integer, intent (out) :: nl
    call doesnotexist("utfint")
  end subroutine utfint

  subroutine vang3d(a)
    implicit none
    double precision, intent (in) :: a
    call doesnotexist("vang3d")
  end subroutine vang3d

  subroutine vclp3d(x1,x2)
    implicit none
    double precision, intent (in) :: x1,x2
    call doesnotexist("vclp3d")
  end subroutine vclp3d

  subroutine vecclr(iclr)
    implicit none
    integer, intent (in) :: iclr
    call doesnotexist("vecclr")
  end subroutine vecclr

  subroutine vecf3d(xv,yv,zv,xp,yp,zp,n,ivec)
    implicit none
    integer, intent (in) :: n,ivec
    double precision, dimension (n), intent (in) :: xv,yv,zv,xp,yp,zp
    call doesnotexist("vecf3d")
  end subroutine vecf3d

  subroutine vecfld(xv,yv,xp,yp,n,ivec)
    implicit none
    integer, intent (in) :: n,ivec
    double precision, dimension (n), intent (in) :: xv,yv,xp,yp
    call doesnotexist("vecfld")
  end subroutine vecfld

  subroutine vecmat(xmat,ymat,nx,ny,xp,yp,ivec)
    implicit none
    integer, intent (in) :: nx,ny,ivec
    double precision, dimension (nx,ny), intent (in) :: xmat,ymat
    double precision, dimension (nx), intent (in) :: xp
    double precision, dimension (ny), intent (in) :: yp
    call doesnotexist("vecmat")
  end subroutine vecmat

  subroutine vecmat3d(xv,yv,zv,nx,ny,nz,xp,yp,zp,ivec)
    implicit none
    integer, intent (in) :: nx,ny,nz,ivec
    double precision, dimension (nx,ny,nz), intent (in) :: xv,yv,zv
    double precision, dimension (nx), intent (in) :: xp
    double precision, dimension (ny), intent (in) :: yp
    double precision, dimension (nz), intent (in) :: zp
    call doesnotexist("vecmat3d")
  end subroutine vecmat3d

  subroutine vecopt(x,copt)
    implicit none
    double precision, intent (in) :: x
    character (len = *), intent (in) :: copt
    call doesnotexist("vecopt")
  end subroutine vecopt

  subroutine vector(ix1,iy1,ix2,iy2,ivec)
    implicit none
    integer, intent (in) :: ix1,iy1,ix2,iy2,ivec
    call doesnotexist("vector")
  end subroutine vector

  subroutine vectr3(x1,y1,z1,x2,y2,z2,ivec)
    implicit none
    double precision, intent (in) :: x1,y1,z1,x2,y2,z2
    integer, intent (in) :: ivec
    call doesnotexist("vectr3")
  end subroutine vectr3

  subroutine vfoc3d(x,y,z,cview)
    implicit none
    double precision, intent (in) :: x,y,z
    character (len = *), intent (in) :: cview
    call doesnotexist("vfoc3d")
  end subroutine vfoc3d

  subroutine view3d(x,y,z,cview)
    implicit none
    double precision, intent (in) :: x,y,z
    character (len = *), intent (in) :: cview
    call doesnotexist("view3d")
  end subroutine view3d

  subroutine vkxbar(i)
    implicit none
    integer, intent (in) :: i
    call doesnotexist("vkxbar")
  end subroutine vkxbar

  subroutine vkybar(i)
    implicit none
    integer, intent (in) :: i
    call doesnotexist("vkybar")
  end subroutine vkybar

  subroutine vkytit(i)
    implicit none
    integer, intent (in) :: i
    call doesnotexist("vkytit")
  end subroutine vkytit

  subroutine vltfil(cfl, copt)
    implicit none
    character (len = *), intent (in) :: cfl, copt
    call doesnotexist("vltfil")
  end subroutine vltfil

  subroutine vscl3d(x)
    implicit none
    double precision, intent (in) :: x
    call doesnotexist("vscl3d")
  end subroutine vscl3d

  subroutine vtx3d(xray,yray,zray,n,copt)
    implicit none
    integer, intent (in) :: n
    double precision, dimension (n), intent (in) :: xray,yray,zray
    character (len = *), intent (in) :: copt
    call doesnotexist("vtx3d")
  end subroutine vtx3d

  subroutine vtxc3d(xray,yray,zray,ic,n,copt)
    implicit none
    integer, intent (in) :: n
    double precision, dimension (n), intent (in) :: xray,yray,zray
    integer, dimension (n), intent (in) :: ic
    character (len = *), intent (in) :: copt
    call doesnotexist("vtxc3d")
  end subroutine vtxc3d

  subroutine vtxn3d(xray,yray,zray,xn,yn,zn,n,copt)
    implicit none
    integer, intent (in) :: n
    double precision, dimension (n), intent (in) :: xray,yray,zray,xn,yn,zn
    character (len = *), intent (in) :: copt
    call doesnotexist("vtxn3d")
  end subroutine vtxn3d

  subroutine vup3d(a)
    implicit none
    double precision, intent (in) :: a
    call doesnotexist("vup3d")
  end subroutine vup3d

  subroutine wgapp(ip,cstr,id)
    implicit none
    character (len = *), intent (in) :: cstr
    integer, intent (in)  :: ip
    integer, intent (out) :: id
    call doesnotexist("wgapp")
  end subroutine wgapp

  subroutine wgappb(ip,iray,nw,nh,id)
    implicit none
    character (len = 1), intent (in), dimension (*) :: iray
    integer, intent (in)  :: ip,nw,nh
    integer, intent (out) :: id
    call doesnotexist("wgappb")
  end subroutine wgappb

  subroutine wgbas(ip,copt,id)
    implicit none
    character (len = *), intent (in) :: copt
    integer, intent (in)  :: ip
    integer, intent (out) :: id
    call doesnotexist("wgbas")
  end subroutine wgbas

  subroutine wgbox(ip,cstr,isel,id)
    implicit none
    character (len = *), intent (in) :: cstr
    integer, intent (in)  :: ip,isel
    integer, intent (out) :: id
    call doesnotexist("wgbox")
  end subroutine wgbox

  subroutine wgbut(ip,cstr,ival,id)
    implicit none
    character (len = *), intent (in) :: cstr
    integer, intent (in)  :: ip,ival
    integer, intent (out) :: id
    call doesnotexist("wgbut")
  end subroutine wgbut

  subroutine wgcmd(ip,clab,cstr,id)
    implicit none
    character (len = *), intent (in) :: clab,cstr
    integer, intent (in)  :: ip
    integer, intent (out) :: id
    call doesnotexist("wgcmd")
  end subroutine wgcmd

  subroutine wgdlis(ip,cstr,isel,id)
    implicit none
    character (len = *), intent (in) :: cstr
    integer, intent (in)  :: ip,isel
    integer, intent (out) :: id
    call doesnotexist("wgdlis")
  end subroutine wgdlis

  subroutine wgdraw(ip,id)
    implicit none
    integer, intent (in)  :: ip
    integer, intent (out) :: id
    call doesnotexist("wgdraw")
  end subroutine wgdraw

  subroutine wgfil(ip,clab,cstr,cmask,id)
    implicit none
    character (len = *), intent (in) :: clab,cstr,cmask
    integer, intent (in)  :: ip
    integer, intent (out) :: id
    call doesnotexist("wgfil")
  end subroutine wgfil

  subroutine wgfin()
    call doesnotexist("wgfin")
  end subroutine wgfin

  subroutine wgicon(ip,clab,nw,nh,cfl,id)
    implicit none
    character (len = *), intent (in) :: clab,cfl
    integer, intent (in)  :: ip,nw,nh
    integer, intent (out) :: id
    call doesnotexist("wgicon")
  end subroutine wgicon

  subroutine wgimg(ip,clab,iray,nw,nh,id)
    implicit none
    character (len = *), intent (in) :: clab
    character (len = 1), intent (in), dimension (*) :: iray
    integer, intent (in)  :: ip,nw,nh
    integer, intent (out) :: id
    call doesnotexist("wgimg")
  end subroutine wgimg

  subroutine wgini(ctype,id)
    implicit none
    character (len = *), intent (in) :: ctype
    integer, intent (out) :: id
    call doesnotexist("wgini")
  end subroutine wgini

  subroutine wglab(ip,cstr,id)
    implicit none
    character (len = *), intent (in) :: cstr
    integer, intent (in)  :: ip
    integer, intent (out) :: id
    call doesnotexist("wglab")
  end subroutine wglab

  subroutine wglis(ip,cstr,isel,id)
    implicit none
    character (len = *), intent (in) :: cstr
    integer, intent (in)  :: ip,isel
    integer, intent (out) :: id
    call doesnotexist("wglis")
  end subroutine wglis

  subroutine wgltxt(ip,clab,cstr,iper,id)
    implicit none
    character (len = *), intent (in) :: clab,cstr
    integer, intent (in)  :: ip,iper
    integer, intent (out) :: id
    call doesnotexist("wgltxt")
  end subroutine wgltxt

  subroutine wgok(ip,id)
    implicit none
    integer, intent (in)  :: ip
    integer, intent (out) :: id
    call doesnotexist("wgok")
  end subroutine wgok

  subroutine wgpbar(ip,x1,x2,xstp,id)
    implicit none
    integer, intent (in)  :: ip
    double precision, intent (in)    :: x1,x2,xstp
    integer, intent (out) :: id
    call doesnotexist("wgpbar")
  end subroutine wgpbar

  subroutine wgpbut(ip,cstr,id)
    implicit none
    character (len = *), intent (in) :: cstr
    integer, intent (in)  :: ip
    integer, intent (out) :: id
    call doesnotexist("wgpbut")
  end subroutine wgpbut

  subroutine wgpicon(ip,clab,nw,nh,cfl,id)
    implicit none
    character (len = *), intent (in) :: clab,cfl
    integer, intent (in)  :: ip,nw,nh
    integer, intent (out) :: id
    call doesnotexist("wgpicon")
  end subroutine wgpicon

  subroutine wgpimg(ip,clab,iray,nw,nh,id)
    implicit none
    character (len = *), intent (in) :: clab
    character (len = 1), intent (in), dimension (*) :: iray
    integer, intent (in)  :: ip,nw,nh
    integer, intent (out) :: id
    call doesnotexist("wgpimg")
  end subroutine wgpimg

  subroutine wgpop(ip,cstr,id)
    implicit none
    character (len = *), intent (in) :: cstr
    integer, intent (in)  :: ip
    integer, intent (out) :: id
    call doesnotexist("wgpop")
  end subroutine wgpop

  subroutine wgpopb(ip,iray,nw,nh,id)
    implicit none
    character (len = 1), intent (in), dimension (*) :: iray
    integer, intent (in)  :: ip,nw,nh
    integer, intent (out) :: id
    call doesnotexist("wgpopb")
  end subroutine wgpopb

  subroutine wgquit(ip,id)
    implicit none
    integer, intent (in)  :: ip
    integer, intent (out) :: id
    call doesnotexist("wgquit")
  end subroutine wgquit

  subroutine wgscl(ip,cstr,x1,x2,xval,ndez,id)
    implicit none
    character (len = *), intent (in) :: cstr
    integer, intent (in)  :: ip,ndez
    double precision, intent (in)    :: x1,x2,xval
    integer, intent (out) :: id
    call doesnotexist("wgscl")
  end subroutine wgscl

  subroutine wgsep(ip,id)
    implicit none
    integer, intent (in)  :: ip
    integer, intent (out) :: id
    call doesnotexist("wgsep")
  end subroutine wgsep

  subroutine wgstxt(ip,nsize,nmax,id)
    implicit none
    integer, intent (in)  :: ip,nsize,nmax
    integer, intent (out) :: id
    call doesnotexist("wgstxt")
  end subroutine wgstxt

  subroutine wgtxt(ip,cstr,id)
    implicit none
    character (len = *), intent (in) :: cstr
    integer, intent (in)  :: ip
    integer, intent (out) :: id
    call doesnotexist("wgtxt")
  end subroutine wgtxt

  subroutine widbar(i)
    implicit none
    integer, intent (in) :: i
    call doesnotexist("widbar")
  end subroutine widbar

  subroutine wimage(cfl)
    implicit none
    character (len = *), intent (in) :: cfl
    call doesnotexist("wimage")
  end subroutine wimage

  subroutine winapp(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("winapp")
  end subroutine winapp

  subroutine wincbk (callbk,copt)
    implicit none
    character (len = *), intent (in) :: copt

    interface
       subroutine callbk(id,nx,ny,nw,nh)
         implicit none
         integer, intent (in) :: id,nx,ny,nw,nh
       end subroutine callbk
    end interface
    call doesnotexist("wincbk")
  end subroutine wincbk

  subroutine windbr(xk,nx,ny,nw,a)
    implicit none
    double precision, intent (in) :: xk,a
    integer, intent (in) :: nx,ny,nw
    call doesnotexist("windbr")
  end subroutine windbr

  subroutine window(nx,ny,nw,nh)
    implicit none
    integer, intent (in) :: nx,ny,nw,nh
    call doesnotexist("window")
  end subroutine window

  subroutine winfnt(cfnt)
    implicit none
    character (len = *), intent (in) :: cfnt
    call doesnotexist("winfnt")
  end subroutine winfnt

  subroutine winico(cstr)
    implicit none
    character (len = *), intent (in) :: cstr
    call doesnotexist("winico")
  end subroutine winico

  subroutine winid(id)
    implicit none
    integer, intent (out) :: id
    call doesnotexist("winid")
  end subroutine winid

  subroutine winjus(cmod)
    implicit none
    character (len = *), intent (in) :: cmod
    call doesnotexist("winjus")
  end subroutine winjus

  subroutine winkey(cmod)
    implicit none
    character (len = *), intent (in) :: cmod
    call doesnotexist("winkey")
  end subroutine winkey

  subroutine winmod(cmod)
    implicit none
    character (len = *), intent (in) :: cmod
    call doesnotexist("winmod")
  end subroutine winmod

  subroutine winopt(iopt,copt)
    implicit none
    character (len = *), intent (in) :: copt
    integer, intent (in) :: iopt
    call doesnotexist("winopt")
  end subroutine winopt

  subroutine winsiz(nw,nh)
    implicit none
    integer, intent (in) :: nw,nh
    call doesnotexist("winsiz")
  end subroutine winsiz

  subroutine wintit(cstr)
    implicit none
    character (len = *), intent (in) :: cstr
    call doesnotexist("wintit")
  end subroutine wintit

  subroutine wintyp(cmod)
    implicit none
    character (len = *), intent (in) :: cmod
    call doesnotexist("wintyp")
  end subroutine wintyp

  subroutine wmfmod(cmod,ckey)
    implicit none
    character (len = *), intent (in) :: cmod,ckey
    call doesnotexist("wmfmod")
  end subroutine wmfmod

  subroutine world()
    call doesnotexist("world")
  end subroutine world

  subroutine wpixel(ix,iy,n)
    implicit none
    integer, intent (in) :: ix,iy,n
    call doesnotexist("wpixel")
  end subroutine wpixel

  subroutine wpixls(iray,ix,iy,nw,nh)
    implicit none
    integer, intent (in) :: ix,iy,nw,nh
    character (len=1), intent (in), dimension (*) :: iray
    call doesnotexist("wpixls")
  end subroutine wpixls

  subroutine wpxrow(iray,ix,iy,n)
    implicit none
    integer, intent (in) :: ix,iy,n
    character (len=1), intent (in), dimension (n) :: iray
    call doesnotexist("wpxrow")
  end subroutine wpxrow

  subroutine writfl(nlu,iray,nbyt,istat)
    implicit none
    integer, intent (in) :: nlu,nbyt
    character (len=1), intent (in), dimension (nbyt) :: iray
    integer, intent (out) :: istat
    call doesnotexist("writfl")
  end subroutine writfl

  subroutine wtiff(cfl)
    implicit none
    character (len=*), intent (in) :: cfl
    call doesnotexist("wtiff")
  end subroutine wtiff

  subroutine x11fnt(cfnt, copt)
    implicit none
    character (len = *), intent (in) :: cfnt, copt
    call doesnotexist("x11fnt")
  end subroutine x11fnt

  subroutine x11mod(cmod)
    implicit none
    character (len = *), intent (in) :: cmod
    call doesnotexist("x11mod")
  end subroutine x11mod

  function x2dpos(x,y)
    implicit none
    double precision, intent (in) :: x,y
    double precision :: x2dpos
    call doesnotexist("x2dpos")
  end function x2dpos

  function x3dabs(x,y,z)
    implicit none
    double precision, intent (in) :: x,y,z
    double precision :: x3dabs
    call doesnotexist("x3dabs")
  end function x3dabs

  function x3dpos(x,y,z)
    implicit none
    double precision, intent (in) :: x,y,z
    double precision :: x3dpos
    call doesnotexist("x3dpos")
  end function x3dpos

  function x3drel(x,y,z)
    implicit none
    double precision, intent (in) :: x,y,z
    double precision :: x3drel
    call doesnotexist("x3drel")
  end function x3drel

  subroutine xaxgit()
    call doesnotexist("xaxgit")
  end subroutine xaxgit

  subroutine xaxis(a,b,or,step,il,cstr,it,ix,iy)
    implicit none
    double precision, intent (in) :: a,b,or,step
    character (len = *), intent (in) :: cstr
    integer, intent (in) :: il,it,ix,iy
    call doesnotexist("xaxis")
  end subroutine xaxis

  subroutine xaxlg(a,b,or,step,il,cstr,it,ix,iy)
    implicit none
    double precision, intent (in) :: a,b,or,step
    character (len = *), intent (in) :: cstr
    integer, intent (in) :: il,it,ix,iy
    call doesnotexist("xaxlg")
  end subroutine xaxlg

  subroutine xaxmap(a,b,or,step,cstr,it,iy)
    implicit none
    double precision, intent (in) :: a,b,or,step
    integer, intent (in) :: it,iy
    character (len = *), intent (in) :: cstr
    call doesnotexist("xaxmap")
  end subroutine xaxmap

  subroutine xcross()
    call doesnotexist("xcross")
  end subroutine xcross

  subroutine xdraw(xx,yy)
    implicit none
    double precision, intent (in) :: xx,yy
    call doesnotexist("xdraw")
  end subroutine xdraw

  function xinvrs(i)
    implicit none
    integer, intent (in) :: i
    double precision :: xinvrs
    call doesnotexist("xinvrs")
  end function xinvrs

  subroutine xmove(x,y)
    implicit none
    double precision, intent (in) :: x,y
    call doesnotexist("xmove")
  end subroutine xmove

  function xposn(x)
    implicit none
    double precision, intent (in) :: x
    double precision :: xposn
    call doesnotexist("xposn")
  end function xposn

  function y2dpos(x,y)
    implicit none
    double precision, intent (in) :: x,y
    double precision :: y2dpos
    call doesnotexist("y2dpos")
  end function y2dpos

  function y3dabs(x,y,z)
    implicit none
    double precision, intent (in) :: x,y,z
    double precision :: y3dabs
    call doesnotexist("y3dabs")
  end function y3dabs

  function y3dpos(x,y,z)
    implicit none
    double precision, intent (in) :: x,y,z
    double precision :: y3dpos
    call doesnotexist("y3dpos")
  end function y3dpos

  function y3drel(x,y,z)
    implicit none
    double precision, intent (in) :: x,y,z
    double precision :: y3drel
    call doesnotexist("y3drel")
  end function y3drel

  subroutine yaxgit()
    call doesnotexist("yaxgit")
  end subroutine yaxgit

  subroutine yaxis(a,b,or,step,il,cstr,it,ix,iy)
    implicit none
    double precision, intent (in) :: a,b,or,step
    character (len = *), intent (in) :: cstr
    integer, intent (in) :: il,it,ix,iy
    call doesnotexist("yaxis")
  end subroutine yaxis

  subroutine yaxlg(a,b,or,step,il,cstr,it,ix,iy)
    implicit none
    double precision, intent (in) :: a,b,or,step
    character (len = *), intent (in) :: cstr
    integer, intent (in) :: il,it,ix,iy
    call doesnotexist("yaxlg")
  end subroutine yaxlg

  subroutine yaxmap(a,b,or,step,cstr,it,ix)
    implicit none
    double precision, intent (in) :: a,b,or,step
    integer, intent (in) :: it,ix
    character (len = *), intent (in) :: cstr
    call doesnotexist("yaxmap")
  end subroutine yaxmap

  subroutine ycross()
    call doesnotexist("ycross")
  end subroutine ycross

  function yinvrs(i)
    implicit none
    integer, intent (in) :: i
    double precision :: yinvrs
    call doesnotexist("yinvrs")
  end function yinvrs

  function yposn(y)
    implicit none
    double precision, intent (in) :: y
    double precision :: yposn
    call doesnotexist("yposn")
  end function yposn

  function z3dpos(x,y,z)
    implicit none
    double precision, intent (in) :: x,y,z
    double precision :: z3dpos
    call doesnotexist("z3dpos")
  end function z3dpos

  subroutine zaxis(a,b,or,step,il,cstr,it,idir,ix,iy)
    implicit none
    double precision, intent (in) :: a,b,or,step
    integer, intent (in) :: il,it,idir,ix,iy
    character (len = *), intent (in) :: cstr
    call doesnotexist("zaxis")
  end subroutine zaxis

  subroutine zaxlg(a,b,or,step,il,cstr,it,idir,ix,iy)
    implicit none
    double precision, intent (in) :: a,b,or,step
    integer, intent (in) :: il,it,idir,ix,iy
    character (len = *), intent (in) :: cstr
    call doesnotexist("zaxlg")
  end subroutine zaxlg

  subroutine zbfers()
    call doesnotexist("zbfers")
  end subroutine zbfers

  subroutine zbffin()
    call doesnotexist("zbffin")
  end subroutine zbffin

  subroutine zbfini(iret)
    implicit none
    integer, intent (out) :: iret
    call doesnotexist("zbfini")
  end subroutine zbfini

  subroutine zbflin(x1,y1,z1,x2,y2,z2)
    implicit none
    double precision, intent (in) :: x1,y1,z1,x2,y2,z2
    call doesnotexist("zbflin")
  end subroutine zbflin

  subroutine zbfmod(copt)
    implicit none
    character (len = *), intent (in) :: copt
    call doesnotexist("zbfmod")
  end subroutine zbfmod

  subroutine zbfres()
    call doesnotexist("zbfres")
  end subroutine zbfres

  subroutine zbftri(x,y,z,ic)
    implicit none
    double precision, dimension (3), intent (in) :: x,y,z
    integer, dimension (3), intent (in) :: ic
    call doesnotexist("zbftri")
  end subroutine zbftri

  subroutine zscale(a,e)
    implicit none
    double precision, intent (in) :: a,e
    call doesnotexist("zscale")
  end subroutine zscale

  subroutine doesnotexist(funcname)
    character (len = *), intent (in) :: funcname
#ifdef DEBUG
    write (*, *) trim(funcname) // " was called!"
#endif
  end subroutine doesnotexist
