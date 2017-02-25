      subroutine bcup(in,out)
c
c  phi(k,l,m,in) is updated to phi(k,l,m,out) for the converging
c  magnetic field term by Barton's monotonic transport method.  
c
c  The points abs(mu) = 1 can be solved for exactly and are in 
c  this subroutine.
c
      include 'fkrplk_par.h'


      common /ce/ de(nemax), beta2(0:nemax), e(-1:nemax+1)
      common /cmu/ dmu(-nmumax+1:nmumax), onmu2(-nmumax:nmumax),
     +     onmu22(-nmumax+1:nmumax),rb3g2(0:nemax), bconv(0:ntaumax)
      integer in,out
      real tphi(0:nemax,-nmumax:nmumax,0:ntaumax),
     &     tphi1(0:nemax,-nmumax:nmumax,0:ntaumax)
c
       tphi=phi(in,:,:,:)
c
      do k=0,ne
         do m=0,ntau
            a = dy*beta(k)
            aa = a*bconv(m)
c
            if(aa .eq. 0.)then
               do l=-nmu,nmu
                  tphi1(k,l,m) = tphi(k,l,m)
               enddo
            else
               psf = exp(-aa)
               tphi1(k,nmu,m) = psf*tphi(k,nmu,m)
               tphi1(k,-nmu,m) = tphi(k,-nmu,m)/psf
               if(aa .gt. 0.)then
                  fl = tphi(k,nmu,m)*onmu22(nmu)
                  do l=nmu-1,-nmu+1,-1
                     lm1 = l-1
                     lp1 = l+1
                     flp = fl
                     d3 = tphi(k,l,m)
                     d2 = 0.5*(d3+tphi(k,lm1,m))
                     delt = 0.5*dmu(l)/dmu(lp1)
                     d1 = d3*(1.+delt)-tphi(k,lp1,m)*delt
c
c  this replaces the if statement below
c
                     si=sign(1.,d3-tphi(k,lm1,m))
                     dk=-amax1(d3,amin1(d1,d2))*(si-1.)/2.
     &                  +amin1(d3,amax1(d1,d2))*(si+1.)/2.


c                     if(d3 .le. tphi(k,lm1,m)) then
c                        dk = amax1(d3,amin1(d1,d2))
c                     else
c                        dk = amin1(d3,amax1(d1,d2))
c                     endif


                     fl = dk*onmu22(l)
                     tphi1(k,l,m)=
     +                    d3+aa*(flp-fl)/(amu(lp1)-amu(lm1))
                  enddo
               else
                  fl = tphi(k,-nmu,m)*onmu22(-nmu+1)
                  do l=-nmu+1,nmu-1
                     lm1 = l-1
                     lp1 = l+1
                     flp = fl
                     d3 = tphi(k,l,m)
                     d2 = 0.5*(tphi(k,lp1,m)+d3)
                     delt = 0.5*dmu(lp1)/dmu(l)
                     d1 = d3*(1.+delt)-tphi(k,lm1,m)*delt
c
c
c  this replaces the if statement below
c
                     si=sign(1.,d3-tphi(k,lm1,m))
                     dk=-amax1(d3,amin1(d1,d2))*(si-1.)/2.
     &                  +amin1(d3,amax1(d1,d2))*(si+1.)/2.


c                     if(d3 .le. tphi(k,lm1,m)) then
c                        dk = amin1(d3,amax1(d1,d2))
c                     else
c                        dk = amax1(d3,amin1(d1,d2))
c                     endif
c
                     fl = dk*onmu22(lp1)
                     tphi1(k,l,m)=
     +                    d3+aa*(fl-flp)/(amu(lp1)-amu(lm1))
                  enddo
               endif
            endif


        enddo
      enddo


      phi(out,:,:,:)=tphi1


      return
      end