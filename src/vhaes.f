c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     *                                                               *
c     *                  copyright (c) 1998 by UCAR                   *
c     *                                                               *
c     *       University Corporation for Atmospheric Research         *
c     *                                                               *
c     *                      all rights reserved                      *
c     *                                                               *
c     *                      SPHEREPACK version 3.2                   *
c     *                                                               *
c     *       A Package of Fortran77 Subroutines and Programs         *
c     *                                                               *
c     *              for Modeling Geophysical Processes               *
c     *                                                               *
c     *                             by                                *
c     *                                                               *
c     *                  John Adams and Paul Swarztrauber             *
c     *                                                               *
c     *                             of                                *
c     *                                                               *
c     *         the National Center for Atmospheric Research          *
c     *                                                               *
c     *                Boulder, Colorado  (80307)  U.S.A.             *
c     *                                                               *
c     *                   which is sponsored by                       *
c     *                                                               *
c     *              the National Science Foundation                  *
c     *                                                               *
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c
c
c ... file vhaes.f
c
c     this file contains code and documentation for subroutines
c     vhaes and vhaesi
c
c ... files which must be loaded with vhaes.f
c
c     sphcom.f, hrfft.f
c
c                                                                              
c     subroutine vhaes(nlat,nlon,ityp,nt,v,w,idvw,jdvw,br,bi,cr,ci,
c    +                 mdab,ndab,wvhaes,lvhaes,work,lwork,ierror)
c
c     subroutine vhaes performs the vector spherical harmonic analysis
c     on the vector field (v,w) and stores the result in the arrays
c     br, bi, cr, and ci. v(i,j) and w(i,j) are the colatitudinal 
c     (measured from the north pole) and east longitudinal components
c     respectively, located at colatitude theta(i) = (i-1)*pi/(nlat-1)
c     and longitude phi(j) = (j-1)*2*pi/nlon. the spectral
c     representation of (v,w) is given at output parameters v,w in 
c     subroutine vhses.  
c
c     input parameters
c
c     nlat   the number of colatitudes on the full sphere including the
c            poles. for example, nlat = 37 for a five degree grid.
c            nlat determines the grid increment in colatitude as
c            pi/(nlat-1).  if nlat is odd the equator is located at
c            grid point i=(nlat+1)/2. if nlat is even the equator is
c            located half way between points i=nlat/2 and i=nlat/2+1.
c            nlat must be at least 3. note: on the half sphere, the
c            number of grid points in the colatitudinal direction is
c            nlat/2 if nlat is even or (nlat+1)/2 if nlat is odd.
c
c     nlon   the number of distinct londitude points.  nlon determines
c            the grid increment in longitude as 2*pi/nlon. for example
c            nlon = 72 for a five degree grid. nlon must be greater
c            than zero. the axisymmetric case corresponds to nlon=1.
c            the efficiency of the computation is improved when nlon
c            is a product of small prime numbers.
c
c     ityp   = 0  no symmetries exist about the equator. the analysis
c                 is performed on the entire sphere.  i.e. on the
c                 arrays v(i,j),w(i,j) for i=1,...,nlat and 
c                 j=1,...,nlon.   
c
c            = 1  no symmetries exist about the equator. the analysis
c                 is performed on the entire sphere.  i.e. on the
c                 arrays v(i,j),w(i,j) for i=1,...,nlat and 
c                 j=1,...,nlon. the curl of (v,w) is zero. that is, 
c                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0. 
c                 the coefficients cr and ci are zero.
c
c            = 2  no symmetries exist about the equator. the analysis
c                 is performed on the entire sphere.  i.e. on the
c                 arrays v(i,j),w(i,j) for i=1,...,nlat and 
c                 j=1,...,nlon. the divergence of (v,w) is zero. i.e., 
c                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. 
c                 the coefficients br and bi are zero.
c
c            = 3  v is symmetric and w is antisymmetric about the 
c                 equator. the analysis is performed on the northern
c                 hemisphere only.  i.e., if nlat is odd the analysis
c                 is performed on the arrays v(i,j),w(i,j) for 
c                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
c                 even the analysis is performed on the the arrays
c                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c
c            = 4  v is symmetric and w is antisymmetric about the 
c                 equator. the analysis is performed on the northern
c                 hemisphere only.  i.e., if nlat is odd the analysis
c                 is performed on the arrays v(i,j),w(i,j) for 
c                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
c                 even the analysis is performed on the the arrays
c                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c                 the curl of (v,w) is zero. that is, 
c                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0. 
c                 the coefficients cr and ci are zero.
c
c            = 5  v is symmetric and w is antisymmetric about the 
c                 equator. the analysis is performed on the northern
c                 hemisphere only.  i.e., if nlat is odd the analysis
c                 is performed on the arrays v(i,j),w(i,j) for 
c                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
c                 even the analysis is performed on the the arrays
c                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c                 the divergence of (v,w) is zero. i.e., 
c                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. 
c                 the coefficients br and bi are zero.
c
c            = 6  v is antisymmetric and w is symmetric about the 
c                 equator. the analysis is performed on the northern
c                 hemisphere only.  i.e., if nlat is odd the analysis
c                 is performed on the arrays v(i,j),w(i,j) for 
c                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
c                 even the analysis is performed on the the arrays
c                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c
c            = 7  v is antisymmetric and w is symmetric about the 
c                 equator. the analysis is performed on the northern
c                 hemisphere only.  i.e., if nlat is odd the analysis
c                 is performed on the arrays v(i,j),w(i,j) for 
c                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
c                 even the analysis is performed on the the arrays
c                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c                 the curl of (v,w) is zero. that is, 
c                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0. 
c                 the coefficients cr and ci are zero.
c
c            = 8  v is antisymmetric and w is symmetric about the 
c                 equator. the analysis is performed on the northern
c                 hemisphere only.  i.e., if nlat is odd the analysis
c                 is performed on the arrays v(i,j),w(i,j) for 
c                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
c                 even the analysis is performed on the the arrays
c                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
c                 the divergence of (v,w) is zero. i.e., 
c                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. 
c                 the coefficients br and bi are zero.
c
c
c     nt     the number of analyses.  in the program that calls vhaes,
c            the arrays v,w,br,bi,cr, and ci can be three dimensional
c            in which case multiple analyses will be performed.
c            the third index is the analysis index which assumes the 
c            values k=1,...,nt.  for a single analysis set nt=1. the
c            discription of the remaining parameters is simplified
c            by assuming that nt=1 or that all the arrays are two
c            dimensional.
c
c     v,w    two or three dimensional arrays (see input parameter nt)
c            that contain the vector function to be analyzed.
c            v is the colatitudnal component and w is the east 
c            longitudinal component. v(i,j),w(i,j) contain the
c            components at colatitude theta(i) = (i-1)*pi/(nlat-1)
c            and longitude phi(j) = (j-1)*2*pi/nlon. the index ranges
c            are defined above at the input parameter ityp.
c
c     idvw   the first dimension of the arrays v,w as it appears in
c            the program that calls vhaes. if ityp .le. 2 then idvw
c            must be at least nlat.  if ityp .gt. 2 and nlat is
c            even then idvw must be at least nlat/2. if ityp .gt. 2
c            and nlat is odd then idvw must be at least (nlat+1)/2.
c
c     jdvw   the second dimension of the arrays v,w as it appears in
c            the program that calls vhaes. jdvw must be at least nlon.
c
c     mdab   the first dimension of the arrays br,bi,cr, and ci as it
c            appears in the program that calls vhaes. mdab must be at
c            least min0(nlat,nlon/2) if nlon is even or at least
c            min0(nlat,(nlon+1)/2) if nlon is odd.
c
c     ndab   the second dimension of the arrays br,bi,cr, and ci as it
c            appears in the program that calls vhaes. ndab must be at
c            least nlat.
c
c     lvhaes an array which must be initialized by subroutine vhaesi.
c            once initialized, wvhaes can be used repeatedly by vhaes
c            as long as nlon and nlat remain unchanged.  wvhaes must
c            not be altered between calls of vhaes.
c
c     lvhaes the dimension of the array wvhaes as it appears in the
c            program that calls vhaes. define
c
c               l1 = min0(nlat,nlon/2) if nlon is even or
c               l1 = min0(nlat,(nlon+1)/2) if nlon is odd
c
c            and
c
c               l2 = nlat/2        if nlat is even or
c               l2 = (nlat+1)/2    if nlat is odd
c
c            then lvhaes must be at least
c
c            l1*l2(nlat+nlat-l1+1)+nlon+15
c
c
c     work   a work array that does not have to be saved.
c
c     lwork  the dimension of the array work as it appears in the
c            program that calls vhaes. define
c
c               l2 = nlat/2        if nlat is even or
c               l2 = (nlat+1)/2    if nlat is odd
c
c            if ityp .le. 2 then lwork must be at least
c
c                       (2*nt+1)*nlat*nlon
c
c            if ityp .gt. 2 then lwork must be at least
c
c                        (2*nt+1)*l2*nlon   
c
c     **************************************************************
c
c     output parameters
c
c     br,bi  two or three dimensional arrays (see input parameter nt)
c     cr,ci  that contain the vector spherical harmonic coefficients
c            in the spectral representation of v(i,j) and w(i,j) given 
c            in the discription of subroutine vhses. br(mp1,np1),
c            bi(mp1,np1),cr(mp1,np1), and ci(mp1,np1) are computed 
c            for mp1=1,...,mmax and np1=mp1,...,nlat except for np1=nlat
c            and odd mp1. mmax=min0(nlat,nlon/2) if nlon is even or 
c            mmax=min0(nlat,(nlon+1)/2) if nlon is odd. 
c      
c     ierror = 0  no errors
c            = 1  error in the specification of nlat
c            = 2  error in the specification of nlon
c            = 3  error in the specification of ityp
c            = 4  error in the specification of nt
c            = 5  error in the specification of idvw
c            = 6  error in the specification of jdvw
c            = 7  error in the specification of mdab
c            = 8  error in the specification of ndab
c            = 9  error in the specification of lvhaes
c            = 10 error in the specification of lwork
c
c ********************************************************
c
c     subroutine vhaesi(nlat,nlon,wvhaes,lvhaes,work,lwork,dwork,
c    +                  ldwork,ierror)
c
c     subroutine vhaesi initializes the array wvhaes which can then be
c     used repeatedly by subroutine vhaes until nlat or nlon is changed.
c
c     input parameters
c
c     nlat   the number of colatitudes on the full sphere including the
c            poles. for example, nlat = 37 for a five degree grid.
c            nlat determines the grid increment in colatitude as
c            pi/(nlat-1).  if nlat is odd the equator is located at
c            grid point i=(nlat+1)/2. if nlat is even the equator is
c            located half way between points i=nlat/2 and i=nlat/2+1.
c            nlat must be at least 3. note: on the half sphere, the
c            number of grid points in the colatitudinal direction is
c            nlat/2 if nlat is even or (nlat+1)/2 if nlat is odd.
c
c     nlon   the number of distinct londitude points.  nlon determines
c            the grid increment in longitude as 2*pi/nlon. for example
c            nlon = 72 for a five degree grid. nlon must be greater
c            than zero. the axisymmetric case corresponds to nlon=1.
c            the efficiency of the computation is improved when nlon
c            is a product of small prime numbers.
c
c     lvhaes the dimension of the array wvhaes as it appears in the
c            program that calls vhaes. define
c
c               l1 = min0(nlat,nlon/2) if nlon is even or
c               l1 = min0(nlat,(nlon+1)/2) if nlon is odd
c
c            and
c
c               l2 = nlat/2        if nlat is even or
c               l2 = (nlat+1)/2    if nlat is odd
c
c            then lvhaes must be at least
c
c               l1*l2*(nlat+nlat-l1+1)+nlon+15
c
c
c     work   a work array that does not have to be saved.
c
c     lwork  the dimension of the array work as it appears in the
c            program that calls vhaes. lwork must be at least
c
c              3*(max0(l1-2,0)*(nlat+nlat-l1-1))/2+5*l2*nlat
c
c     dwork  an unsaved double precision work space
c
c     ldwork the length of the array dwork as it appears in the
c            program that calls vhaesi.  ldwork must be at least
c            2*(nlat+1)
c
c
c     **************************************************************
c
c     output parameters
c
c     wvhaes an array which is initialized for use by subroutine vhaes.
c            once initialized, wvhaes can be used repeatedly by vhaes
c            as long as nlat or nlon remain unchanged.  wvhaes must not
c            be altered between calls of vhaes.
c
c
c     ierror = 0  no errors
c            = 1  error in the specification of nlat
c            = 2  error in the specification of nlon
c            = 3  error in the specification of lvhaes
c            = 4  error in the specification of lwork
c            = 5  error in the specification of ldwork
c
c
      subroutine vhaes(nlat,nlon,ityp,nt,v,w,idvw,jdvw,br,bi,cr,ci,
     1           mdab,ndab,wvhaes,lvhaes,work,lwork,ierror)
      implicit double precision (a-h,o-z)
      dimension v(idvw,jdvw,*),w(idvw,jdvw,*),br(mdab,ndab,*),
     1          bi(mdab,ndab,*),cr(mdab,ndab,*),ci(mdab,ndab,*),
     2          work(*),wvhaes(*)
      ierror = 1
      if(nlat .lt. 3) return
      ierror = 2
      if(nlon .lt. 1) return
      ierror = 3
      if(ityp.lt.0 .or. ityp.gt.8) return
      ierror = 4
      if(nt .lt. 0) return
      ierror = 5
      imid = (nlat+1)/2
      if((ityp.le.2 .and. idvw.lt.nlat) .or.
     1   (ityp.gt.2 .and. idvw.lt.imid)) return
      ierror = 6
      if(jdvw .lt. nlon) return
      ierror = 7
      mmax = min0(nlat,(nlon+1)/2)
      if(mdab .lt. mmax) return
      ierror = 8
      if(ndab .lt. nlat) return
      ierror = 9
      idz = (mmax*(nlat+nlat-mmax+1))/2
      lzimn = idz*imid
      if(lvhaes .lt. lzimn+lzimn+nlon+15) return
      ierror = 10
      idv = nlat
      if(ityp .gt. 2) idv = imid
      lnl = nt*idv*nlon
      if(lwork .lt. lnl+lnl+idv*nlon) return
      ierror = 0
      ist = 0
      if(ityp .le. 2) ist = imid
      iw1 = ist+1
      iw2 = lnl+1
      iw3 = iw2+ist
      iw4 = iw2+lnl
      jw1 = lzimn+1
      jw2 = jw1+lzimn
      call vhaes1(nlat,nlon,ityp,nt,imid,idvw,jdvw,v,w,mdab,ndab,
     1     br,bi,cr,ci,idv,work,work(iw1),work(iw2),work(iw3),
     2     work(iw4),idz,wvhaes,wvhaes(jw1),wvhaes(jw2))
      return
      end

      subroutine vhaes1(nlat,nlon,ityp,nt,imid,idvw,jdvw,v,w,mdab,
     1   ndab,br,bi,cr,ci,idv,ve,vo,we,wo,work,idz,zv,zw,wrfft)
      implicit double precision (a-h,o-z)
      dimension v(idvw,jdvw,*),w(idvw,jdvw,*),br(mdab,ndab,*),
     1          bi(mdab,ndab,*),cr(mdab,ndab,*),ci(mdab,ndab,*),
     2          ve(idv,nlon,*),vo(idv,nlon,*),we(idv,nlon,*),
     3          wo(idv,nlon,*),work(*),wrfft(*),
     4          zv(idz,*),zw(idz,*)
      nlp1 = nlat+1
      tsn = 2.d0/nlon
      fsn = 4.d0/nlon
      mlat = mod(nlat,2)
      mlon = mod(nlon,2)
      mmax = min0(nlat,(nlon+1)/2)
      imm1 = imid
      if(mlat .ne. 0) imm1 = imid-1
      if(ityp .gt. 2) go to 3  
      do k=1,nt 
      do i=1,imm1
      do 5 j=1,nlon
      ve(i,j,k) = tsn*(v(i,j,k)+v(nlp1-i,j,k))
      vo(i,j,k) = tsn*(v(i,j,k)-v(nlp1-i,j,k))
      we(i,j,k) = tsn*(w(i,j,k)+w(nlp1-i,j,k))
      wo(i,j,k) = tsn*(w(i,j,k)-w(nlp1-i,j,k))
    5 continue
      end do
      end do
      go to 2
    3 do k=1,nt
      do i=1,imm1 
      do 8 j=1,nlon
      ve(i,j,k) = fsn*v(i,j,k)
      vo(i,j,k) = fsn*v(i,j,k)
      we(i,j,k) = fsn*w(i,j,k)
      wo(i,j,k) = fsn*w(i,j,k)
    8 continue
      end do
      end do
    2 if(mlat .eq. 0) go to 7
      do k=1,nt 
      do 6 j=1,nlon
      ve(imid,j,k) = tsn*v(imid,j,k)
      we(imid,j,k) = tsn*w(imid,j,k)
    6 continue
      end do
    7 do 9 k=1,nt
      call hrfftf(idv,nlon,ve(1,1,k),idv,wrfft,work)
      call hrfftf(idv,nlon,we(1,1,k),idv,wrfft,work)
    9 continue 
      ndo1 = nlat
      ndo2 = nlat
      if(mlat .ne. 0) ndo1 = nlat-1
      if(mlat .eq. 0) ndo2 = nlat-1
      if(ityp.eq.2 .or. ityp.eq.5 .or. ityp.eq.8) go to 11 
      do k=1,nt
      do mp1=1,mmax
      do 10 np1=mp1,nlat
      br(mp1,np1,k)=0.d0
      bi(mp1,np1,k)=0.d0
   10 continue
      end do
      end do
   11 if(ityp.eq.1 .or. ityp.eq.4 .or. ityp.eq.7) go to 13 
      do k=1,nt
      do mp1=1,mmax
      do 12 np1=mp1,nlat
      cr(mp1,np1,k)=0.d0
      ci(mp1,np1,k)=0.d0
   12 continue
      end do
      end do
   13 itypp = ityp+1
      go to (1,100,200,300,400,500,600,700,800),itypp
c
c     case ityp=0 ,  no symmetries
c
c     case m=0
c
    1 do k=1,nt
      do i=1,imid
      do 15 np1=2,ndo2,2
      br(1,np1,k) = br(1,np1,k)+zv(np1,i)*ve(i,1,k)
      cr(1,np1,k) = cr(1,np1,k)-zv(np1,i)*we(i,1,k)
   15 continue
      end do
      end do
      do k=1,nt
      do i=1,imm1
      do 16 np1=3,ndo1,2
      br(1,np1,k) = br(1,np1,k)+zv(np1,i)*vo(i,1,k)
      cr(1,np1,k) = cr(1,np1,k)-zv(np1,i)*wo(i,1,k)
   16 continue
      end do
      end do
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) return
      do 20 mp1=2,mmax
      m = mp1-1
      mb = m*(nlat-1)-(m*(m-1))/2
      mp2 = mp1+1
      if(mp1 .gt. ndo1) go to 17
      do k=1,nt
      do i=1,imm1
      do 23 np1=mp1,ndo1,2
      br(mp1,np1,k) = br(mp1,np1,k)+zv(np1+mb,i)*vo(i,2*mp1-2,k)
     1                             +zw(np1+mb,i)*we(i,2*mp1-1,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)+zv(np1+mb,i)*vo(i,2*mp1-1,k)
     1                             -zw(np1+mb,i)*we(i,2*mp1-2,k)
      cr(mp1,np1,k) = cr(mp1,np1,k)-zv(np1+mb,i)*wo(i,2*mp1-2,k)
     1                             +zw(np1+mb,i)*ve(i,2*mp1-1,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-zv(np1+mb,i)*wo(i,2*mp1-1,k)
     1                             -zw(np1+mb,i)*ve(i,2*mp1-2,k)
   23 continue
      end do
      end do
      if(mlat .eq. 0) go to 17
      do k=1,nt
      do 24 np1=mp1,ndo1,2
      br(mp1,np1,k) = br(mp1,np1,k)+zw(np1+mb,imid)*we(imid,2*mp1-1,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)-zw(np1+mb,imid)*we(imid,2*mp1-2,k)
      cr(mp1,np1,k) = cr(mp1,np1,k)+zw(np1+mb,imid)*ve(imid,2*mp1-1,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-zw(np1+mb,imid)*ve(imid,2*mp1-2,k)
   24 continue
      end do
   17 if(mp2 .gt. ndo2) go to 20
      do k=1,nt
      do i=1,imm1
      do 21 np1=mp2,ndo2,2
      br(mp1,np1,k) = br(mp1,np1,k)+zv(np1+mb,i)*ve(i,2*mp1-2,k)
     1                             +zw(np1+mb,i)*wo(i,2*mp1-1,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)+zv(np1+mb,i)*ve(i,2*mp1-1,k)
     1                             -zw(np1+mb,i)*wo(i,2*mp1-2,k)
      cr(mp1,np1,k) = cr(mp1,np1,k)-zv(np1+mb,i)*we(i,2*mp1-2,k)
     1                             +zw(np1+mb,i)*vo(i,2*mp1-1,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-zv(np1+mb,i)*we(i,2*mp1-1,k)
     1                             -zw(np1+mb,i)*vo(i,2*mp1-2,k)
   21 continue
      end do
      end do
      if(mlat .eq. 0) go to 20
      do k=1,nt
      do 22 np1=mp2,ndo2,2
      br(mp1,np1,k) = br(mp1,np1,k)+zv(np1+mb,imid)*ve(imid,2*mp1-2,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)+zv(np1+mb,imid)*ve(imid,2*mp1-1,k)
      cr(mp1,np1,k) = cr(mp1,np1,k)-zv(np1+mb,imid)*we(imid,2*mp1-2,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-zv(np1+mb,imid)*we(imid,2*mp1-1,k)
   22 continue
      end do
   20 continue
      return
c
c     case ityp=1 ,  no symmetries but cr and ci equal zero
c
c     case m=0
c
  100 do k=1,nt
      do i=1,imid
      do 115 np1=2,ndo2,2
      br(1,np1,k) = br(1,np1,k)+zv(np1,i)*ve(i,1,k)
  115 continue
      end do
      end do
      do k=1,nt
      do i=1,imm1
      do 116 np1=3,ndo1,2
      br(1,np1,k) = br(1,np1,k)+zv(np1,i)*vo(i,1,k)
  116 continue
      end do
      end do
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) return
      do 120 mp1=2,mmax
      m = mp1-1
      mb = m*(nlat-1)-(m*(m-1))/2
      mp2 = mp1+1
      if(mp1 .gt. ndo1) go to 117
      do k=1,nt
      do i=1,imm1
      do 123 np1=mp1,ndo1,2
      br(mp1,np1,k) = br(mp1,np1,k)+zv(np1+mb,i)*vo(i,2*mp1-2,k)
     1                             +zw(np1+mb,i)*we(i,2*mp1-1,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)+zv(np1+mb,i)*vo(i,2*mp1-1,k)
     1                             -zw(np1+mb,i)*we(i,2*mp1-2,k)
  123 continue
      end do
      end do
      if(mlat .eq. 0) go to 117
      do k=1,nt
      do 124 np1=mp1,ndo1,2
      br(mp1,np1,k) = br(mp1,np1,k)+zw(np1+mb,imid)*we(imid,2*mp1-1,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)-zw(np1+mb,imid)*we(imid,2*mp1-2,k)
  124 continue
      end do
  117 if(mp2 .gt. ndo2) go to 120
      do k=1,nt
      do i=1,imm1
      do 121 np1=mp2,ndo2,2
      br(mp1,np1,k) = br(mp1,np1,k)+zv(np1+mb,i)*ve(i,2*mp1-2,k)
     1                             +zw(np1+mb,i)*wo(i,2*mp1-1,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)+zv(np1+mb,i)*ve(i,2*mp1-1,k)
     1                             -zw(np1+mb,i)*wo(i,2*mp1-2,k)
  121 continue
      end do
      end do
      if(mlat .eq. 0) go to 120
      do k=1,nt
      do 122 np1=mp2,ndo2,2
      br(mp1,np1,k) = br(mp1,np1,k)+zv(np1+mb,imid)*ve(imid,2*mp1-2,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)+zv(np1+mb,imid)*ve(imid,2*mp1-1,k)
  122 continue
      end do
  120 continue
      return
c
c     case ityp=2 ,  no symmetries but br and bi equal zero   
c
c     case m=0
c
  200 do k=1,nt
      do i=1,imid
      do 215 np1=2,ndo2,2
      cr(1,np1,k) = cr(1,np1,k)-zv(np1,i)*we(i,1,k)
  215 continue
      end do
      end do
      do k=1,nt
      do i=1,imm1
      do 216 np1=3,ndo1,2
      cr(1,np1,k) = cr(1,np1,k)-zv(np1,i)*wo(i,1,k)
  216 continue
      end do
      end do
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) return
      do 220 mp1=2,mmax
      m = mp1-1
      mb = m*(nlat-1)-(m*(m-1))/2
      mp2 = mp1+1
      if(mp1 .gt. ndo1) go to 217
      do k=1,nt
      do i=1,imm1
      do 223 np1=mp1,ndo1,2
      cr(mp1,np1,k) = cr(mp1,np1,k)-zv(np1+mb,i)*wo(i,2*mp1-2,k)
     1                             +zw(np1+mb,i)*ve(i,2*mp1-1,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-zv(np1+mb,i)*wo(i,2*mp1-1,k)
     1                             -zw(np1+mb,i)*ve(i,2*mp1-2,k)
  223 continue
      end do
      end do
      if(mlat .eq. 0) go to 217
      do k=1,nt
      do 224 np1=mp1,ndo1,2
      cr(mp1,np1,k) = cr(mp1,np1,k)+zw(np1+mb,imid)*ve(imid,2*mp1-1,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-zw(np1+mb,imid)*ve(imid,2*mp1-2,k)
  224 continue
      end do
  217 if(mp2 .gt. ndo2) go to 220
      do k=1,nt
      do i=1,imm1
      do 221 np1=mp2,ndo2,2
      cr(mp1,np1,k) = cr(mp1,np1,k)-zv(np1+mb,i)*we(i,2*mp1-2,k)
     1                             +zw(np1+mb,i)*vo(i,2*mp1-1,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-zv(np1+mb,i)*we(i,2*mp1-1,k)
     1                             -zw(np1+mb,i)*vo(i,2*mp1-2,k)
  221 continue
      end do
      end do
      if(mlat .eq. 0) go to 220
      do k=1,nt
      do 222 np1=mp2,ndo2,2
      cr(mp1,np1,k) = cr(mp1,np1,k)-zv(np1+mb,imid)*we(imid,2*mp1-2,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-zv(np1+mb,imid)*we(imid,2*mp1-1,k)
  222 continue
      end do
  220 continue
      return
c
c     case ityp=3 ,  v even , w odd
c
c     case m=0
c
  300 do k=1,nt
      do i=1,imid
      do 315 np1=2,ndo2,2
      br(1,np1,k) = br(1,np1,k)+zv(np1,i)*ve(i,1,k)
  315 continue
      end do
      end do
      do k=1,nt
      do i=1,imm1
      do 316 np1=3,ndo1,2
      cr(1,np1,k) = cr(1,np1,k)-zv(np1,i)*wo(i,1,k)
  316 continue
      end do
      end do
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) return
      do 320 mp1=2,mmax
      m = mp1-1
      mb = m*(nlat-1)-(m*(m-1))/2
      mp2 = mp1+1
      if(mp1 .gt. ndo1) go to 317
      do k=1,nt
      do i=1,imm1
      do 323 np1=mp1,ndo1,2
      cr(mp1,np1,k) = cr(mp1,np1,k)-zv(np1+mb,i)*wo(i,2*mp1-2,k)
     1                             +zw(np1+mb,i)*ve(i,2*mp1-1,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-zv(np1+mb,i)*wo(i,2*mp1-1,k)
     1                             -zw(np1+mb,i)*ve(i,2*mp1-2,k)
  323 continue
      end do
      end do
      if(mlat .eq. 0) go to 317
      do k=1,nt
      do 324 np1=mp1,ndo1,2
      cr(mp1,np1,k) = cr(mp1,np1,k)+zw(np1+mb,imid)*ve(imid,2*mp1-1,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-zw(np1+mb,imid)*ve(imid,2*mp1-2,k)
  324 continue
      end do
  317 if(mp2 .gt. ndo2) go to 320
      do k=1,nt
      do i=1,imm1
      do 321 np1=mp2,ndo2,2
      br(mp1,np1,k) = br(mp1,np1,k)+zv(np1+mb,i)*ve(i,2*mp1-2,k)
     1                             +zw(np1+mb,i)*wo(i,2*mp1-1,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)+zv(np1+mb,i)*ve(i,2*mp1-1,k)
     1                             -zw(np1+mb,i)*wo(i,2*mp1-2,k)
  321 continue
      end do
      end do
      if(mlat .eq. 0) go to 320
      do k=1,nt
      do 322 np1=mp2,ndo2,2
      br(mp1,np1,k) = br(mp1,np1,k)+zv(np1+mb,imid)*ve(imid,2*mp1-2,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)+zv(np1+mb,imid)*ve(imid,2*mp1-1,k)
  322 continue
      end do
  320 continue
      return
c
c     case ityp=4 ,  v even, w odd, and cr and ci equal 0. 
c
c     case m=0
c
  400 do k=1,nt
      do i=1,imid
      do 415 np1=2,ndo2,2
      br(1,np1,k) = br(1,np1,k)+zv(np1,i)*ve(i,1,k)
  415 continue
      end do
      end do
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) return
      do 420 mp1=2,mmax
      m = mp1-1
      mb = m*(nlat-1)-(m*(m-1))/2
      mp2 = mp1+1
      if(mp2 .gt. ndo2) go to 420
      do k=1,nt
      do i=1,imm1
      do 421 np1=mp2,ndo2,2
      br(mp1,np1,k) = br(mp1,np1,k)+zv(np1+mb,i)*ve(i,2*mp1-2,k)
     1                             +zw(np1+mb,i)*wo(i,2*mp1-1,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)+zv(np1+mb,i)*ve(i,2*mp1-1,k)
     1                             -zw(np1+mb,i)*wo(i,2*mp1-2,k)
  421 continue
      end do
      end do
      if(mlat .eq. 0) go to 420
      do k=1,nt
      do 422 np1=mp2,ndo2,2
      br(mp1,np1,k) = br(mp1,np1,k)+zv(np1+mb,imid)*ve(imid,2*mp1-2,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)+zv(np1+mb,imid)*ve(imid,2*mp1-1,k)
  422 continue
      end do
  420 continue
      return
c
c     case ityp=5   v even, w odd, and br and bi equal zero
c
c     case m=0
c
  500 do k=1,nt
      do i=1,imm1
      do 516 np1=3,ndo1,2
      cr(1,np1,k) = cr(1,np1,k)-zv(np1,i)*wo(i,1,k)
  516 continue
      end do
      end do
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) return
      do 520 mp1=2,mmax
      m = mp1-1
      mb = m*(nlat-1)-(m*(m-1))/2
      mp2 = mp1+1
      if(mp1 .gt. ndo1) go to 520
      do k=1,nt
      do i=1,imm1
      do 523 np1=mp1,ndo1,2
      cr(mp1,np1,k) = cr(mp1,np1,k)-zv(np1+mb,i)*wo(i,2*mp1-2,k)
     1                             +zw(np1+mb,i)*ve(i,2*mp1-1,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-zv(np1+mb,i)*wo(i,2*mp1-1,k)
     1                             -zw(np1+mb,i)*ve(i,2*mp1-2,k)
  523 continue
      end do
      end do
      if(mlat .eq. 0) go to 520
      do k=1,nt
      do 524 np1=mp1,ndo1,2
      cr(mp1,np1,k) = cr(mp1,np1,k)+zw(np1+mb,imid)*ve(imid,2*mp1-1,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-zw(np1+mb,imid)*ve(imid,2*mp1-2,k)
  524 continue
      end do
  520 continue
      return
c
c     case ityp=6 ,  v odd , w even
c
c     case m=0
c
  600 do k=1,nt
      do i=1,imid
      do 615 np1=2,ndo2,2
      cr(1,np1,k) = cr(1,np1,k)-zv(np1,i)*we(i,1,k)
  615 continue
      end do
      end do
      do k=1,nt
      do i=1,imm1
      do 616 np1=3,ndo1,2
      br(1,np1,k) = br(1,np1,k)+zv(np1,i)*vo(i,1,k)
  616 continue
      end do
      end do
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) return
      do 620 mp1=2,mmax
      m = mp1-1
      mb = m*(nlat-1)-(m*(m-1))/2
      mp2 = mp1+1
      if(mp1 .gt. ndo1) go to 617
      do k=1,nt
      do i=1,imm1
      do 623 np1=mp1,ndo1,2
      br(mp1,np1,k) = br(mp1,np1,k)+zv(np1+mb,i)*vo(i,2*mp1-2,k)
     1                             +zw(np1+mb,i)*we(i,2*mp1-1,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)+zv(np1+mb,i)*vo(i,2*mp1-1,k)
     1                             -zw(np1+mb,i)*we(i,2*mp1-2,k)
  623 continue
      end do
      end do
      if(mlat .eq. 0) go to 617
      do k=1,nt
      do 624 np1=mp1,ndo1,2
      br(mp1,np1,k) = br(mp1,np1,k)+zw(np1+mb,imid)*we(imid,2*mp1-1,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)-zw(np1+mb,imid)*we(imid,2*mp1-2,k)
  624 continue
      end do
  617 if(mp2 .gt. ndo2) go to 620
      do k=1,nt
      do i=1,imm1
      do 621 np1=mp2,ndo2,2
      cr(mp1,np1,k) = cr(mp1,np1,k)-zv(np1+mb,i)*we(i,2*mp1-2,k)
     1                             +zw(np1+mb,i)*vo(i,2*mp1-1,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-zv(np1+mb,i)*we(i,2*mp1-1,k)
     1                             -zw(np1+mb,i)*vo(i,2*mp1-2,k)
  621 continue
      end do
      end do
      if(mlat .eq. 0) go to 620
      do k=1,nt
      do 622 np1=mp2,ndo2,2
      cr(mp1,np1,k) = cr(mp1,np1,k)-zv(np1+mb,imid)*we(imid,2*mp1-2,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-zv(np1+mb,imid)*we(imid,2*mp1-1,k)
  622 continue
      end do
  620 continue
      return
c
c     case ityp=7   v odd, w even, and cr and ci equal zero
c
c     case m=0
c
  700 do k=1,nt
      do i=1,imm1
      do 716 np1=3,ndo1,2
      br(1,np1,k) = br(1,np1,k)+zv(np1,i)*vo(i,1,k)
  716 continue
      end do
      end do
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) return
      do 720 mp1=2,mmax
      m = mp1-1
      mb = m*(nlat-1)-(m*(m-1))/2
      mp2 = mp1+1
      if(mp1 .gt. ndo1) go to 720
      do k=1,nt
      do i=1,imm1
      do 723 np1=mp1,ndo1,2
      br(mp1,np1,k) = br(mp1,np1,k)+zv(np1+mb,i)*vo(i,2*mp1-2,k)
     1                             +zw(np1+mb,i)*we(i,2*mp1-1,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)+zv(np1+mb,i)*vo(i,2*mp1-1,k)
     1                             -zw(np1+mb,i)*we(i,2*mp1-2,k)
  723 continue
      end do
      end do
      if(mlat .eq. 0) go to 720
      do k=1,nt
      do 724 np1=mp1,ndo1,2
      br(mp1,np1,k) = br(mp1,np1,k)+zw(np1+mb,imid)*we(imid,2*mp1-1,k)
      bi(mp1,np1,k) = bi(mp1,np1,k)-zw(np1+mb,imid)*we(imid,2*mp1-2,k)
  724 continue
      end do
  720 continue
      return
c
c     case ityp=8   v odd, w even, and both br and bi equal zero
c
c     case m=0
c
  800 do k=1,nt
      do i=1,imid
      do 815 np1=2,ndo2,2
      cr(1,np1,k) = cr(1,np1,k)-zv(np1,i)*we(i,1,k)
  815 continue
      end do
      end do
c
c     case m = 1 through nlat-1
c
      if(mmax .lt. 2) return
      do 820 mp1=2,mmax
      m = mp1-1
      mb = m*(nlat-1)-(m*(m-1))/2
      mp2 = mp1+1
      if(mp2 .gt. ndo2) go to 820
      do k=1,nt
      do i=1,imm1
      do 821 np1=mp2,ndo2,2
      cr(mp1,np1,k) = cr(mp1,np1,k)-zv(np1+mb,i)*we(i,2*mp1-2,k)
     1                             +zw(np1+mb,i)*vo(i,2*mp1-1,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-zv(np1+mb,i)*we(i,2*mp1-1,k)
     1                             -zw(np1+mb,i)*vo(i,2*mp1-2,k)
  821 continue
      end do
      end do
      if(mlat .eq. 0) go to 820
      do k=1,nt
      do 822 np1=mp2,ndo2,2
      cr(mp1,np1,k) = cr(mp1,np1,k)-zv(np1+mb,imid)*we(imid,2*mp1-2,k)
      ci(mp1,np1,k) = ci(mp1,np1,k)-zv(np1+mb,imid)*we(imid,2*mp1-1,k)
  822 continue
      end do
  820 continue
      return
      end
c
c     dwork must be of length at least 2*(nlat+1)
c
      subroutine vhaesi(nlat,nlon,wvhaes,lvhaes,work,lwork,dwork,
     +                  ldwork,ierror)
      implicit double precision (a-h,o-z)
      dimension wvhaes(lvhaes),work(lwork)
      double precision dwork(ldwork)
      ierror = 1
      if(nlat .lt. 3) return
      ierror = 2
      if(nlon .lt. 1) return
      ierror = 3
      mmax = min0(nlat,(nlon+1)/2)
      imid = (nlat+1)/2
      lzimn = (imid*mmax*(nlat+nlat-mmax+1))/2
      if(lvhaes .lt. lzimn+lzimn+nlon+15) return
      ierror = 4
      labc = 3*(max0(mmax-2,0)*(nlat+nlat-mmax-1))/2
      if(lwork .lt. 5*nlat*imid+labc) return
      ierror = 5
      if (ldwork .lt. 2*(nlat+1)) return
      ierror = 0
      iw1 = 3*nlat*imid+1
      idz = (mmax*(nlat+nlat-mmax+1))/2
      CALL VEA1(NLAT,NLON,IMID,WVHAES,WVHAES(LZIMN+1),IDZ,
     +          WORK,WORK(IW1),DWORK)
      call hrffti(nlon,wvhaes(2*lzimn+1))
      return
      end
      subroutine vea1(nlat,nlon,imid,zv,zw,idz,zin,wzvin,dwork)
      implicit double precision (a-h,o-z)
      dimension zv(idz,*),zw(idz,*),zin(imid,nlat,3),wzvin(*)
      double precision dwork(*)
      mmax = min0(nlat,(nlon+1)/2)
      call zvinit (nlat,nlon,wzvin,dwork)
      do mp1=1,mmax
      m = mp1-1
      call zvin (0,nlat,nlon,m,zin,i3,wzvin)
      do np1=mp1,nlat
      mn = m*(nlat-1)-(m*(m-1))/2+np1
      do 33 i=1,imid
      zv(mn,i) = zin(i,np1,i3)
   33 continue
      end do
      end do
      call zwinit (nlat,nlon,wzvin,dwork)
      do mp1=1,mmax
      m = mp1-1
      call zwin (0,nlat,nlon,m,zin,i3,wzvin)
      do np1=mp1,nlat
      mn = m*(nlat-1)-(m*(m-1))/2+np1
      do 34 i=1,imid
      zw(mn,i) = zin(i,np1,i3)
   34 continue
      end do
      end do
      return
      end
