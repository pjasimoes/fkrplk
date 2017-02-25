      subroutine eup(in,out)
c
c  phi(k,l,m,in) is updated to phi(k,l,m,out) for the energy term by
c  Barton's method of monotonic transport.
c
      include 'fkrplk_par.h'
      common /ce/ de(nemax), beta2(0:nemax), e(-1:nemax+1)
      integer in,out
      real tphi(0:nemax,-nmumax:nmumax,0:ntaumax),
     &     tphi1(0:nemax,-nmumax:nmumax,0:ntaumax)
c
      tphi(0:ne,-nmu:nmu,0:ntau)=phi(in,0:ne,-nmu:nmu,0:ntau)
c
      do m=0,ntau
         do l=-nmu,nmu
            a = 2.*dy*rlmbda(m)
c     
c     k=ne:  For this point the best value of d2 or d3
c            is used as the transported value; d1 is not 
c            computable.
c
            k = ne
            km1 = k-1
            d3 = tphi(k,l,m)
            d2 = 0.5*(tphi(k,l,m)+tphi(km1,l,m))
c
c  this replaces the if statement below
c
            si0=sign(1.,tphi(k,l,m)-tphi(km1,l,m))
            dk=amax1(d3,d2)*(1.-si0)/2.+amin1(d3,d2)*(si0+1.)/2.
c
c           if(tphi(k,l,m) .le. tphi(km1,l,m)) then
c              dk = amax1(d3,d2)
c           else
c              dk = amin1(d3,d2)
c           endif
c
            f1 = dk/beta2(k)
            tphi1(k,l,m) = d3 - a*f1/(e(ne+1)-e(km1))
            do k=ne-1,1,-1
               kp1 = k+1
               km1 = k-1
               f1p = f1
               d3 = tphi(k,l,m)
               d2 = 0.5*(d3 + tphi(km1,l,m))
               delt = 0.5*(e(k)-e(km1))/de(kp1)
               d1 = d3*(1.+delt)-tphi(kp1,l,m)*delt
c
c  this replaces the if statement below
c
               si=sign(1.,d3-tphi(km1,l,m))
               dk=-amax1(d3,amin1(d1,d2))*(si-1.)/2.
     &           +amin1(d3,amax1(d1,d2))*(si+1.)/2.
c
c               if(d3 .le. tphi(km1,l,m)) then
c                  dk=amax1(d3,amin1(d1,d2))
c               else
c                  dk=amin1(d3,amax1(d1,d2))
c               endif
c
               f1 = dk/beta2(k)
               tphi1(k,l,m) = d3 + a*(f1p-f1)/(e(kp1)-e(km1))
            enddo
c
c     k=0:  For this point the best value between the upwind 
c           and d1 value is used.  Note that phi(0) is always
c           greater than phi(-1).
c
            k = 0
            kp1 = k+1
            km1 = k-1
            f1p = f1
            d3 = tphi(k,l,m)
            delt = 0.5*(e(k)-e(km1))/de(kp1)
            d1 = d3*(1.+delt)-tphi(kp1,l,m)*delt
            dk = amin1(d3,d1)
            f1 = dk/beta2(0)
            tphi1(k,l,m) = d3 + a*(f1p-f1)/(e(kp1)-e(km1))
        enddo
      enddo
c
      phi(out,0:ne,-nmu:nmu,0:ntau)=tphi1(0:ne,-nmu:nmu,0:ntau)
c
      return 
      end