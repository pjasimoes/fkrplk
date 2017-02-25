      subroutine muup(in,out)


c     phi(k,l,m,in) is updated to phi(k,l,m,out) for the pitch angle 
c     diffusion term.  The Crank-Nicholson method is employed and 
c     subroutine cntridag is called to solve the tridiagonal matrix.  


      include 'fkrplk_par.h'
      common /cmu/ dmu(-nmumax+1:nmumax), onmu2(-nmumax:nmumax),
     +     onmu22(-nmumax+1:nmumax), rb3g2(0:nemax), bconv(0:ntaumax)


      real a(0:nemax,-nmumax:nmumax,0:ntaumax),
     &   b(0:nemax,-nmumax:nmumax,0:ntaumax),
     &   c(0:nemax,-nmumax:nmumax,0:ntaumax),
     &   r(0:nemax,-nmumax:nmumax,0:ntaumax),
     &   tphi(0:nemax,-nmumax:nmumax,0:ntaumax),
     &   tphi1(0:nemax,-nmumax:nmumax,0:ntaumax),
     &   gam(0:nemax,-nmumax:nmumax,0:ntaumax)


      integer in,out
c
      tphi(0:ne,-nmu:nmu,0:ntau)=phi(in,0:ne,-nmu:nmu,0:ntau)
c
      do k=0,ne
         do m=0,ntau
            do l=-nmu+1,nmu-1
               b1=beta(k)*.25
               b4=rlmbda(m)*rb3g2(k)
               lp1=l+1
               lm1=l-1
               dydsum=dy/(dmu(lp1)+dmu(l))
               b5=b4*onmu2(l)
               aa=b5*dydsum
               bb1=(-b4*amu(l))*dydsum
c     a(l), b(l), c(l) are the elements of the tridiagonal matrix. 
c     r(l) is the vector product of the matrix and the vector 
c     of updated phi(k,l,m,out).  
               a(k,l,m)=bb1-aa/dmu(l)
               b(k,l,m)=1.+b5*dy/(dmu(l)*dmu(lp1))
               c(k,l,m)=-bb1-aa/dmu(lp1)
               r(k,l,m)=aa*((tphi(k,lp1,m)-tphi(k,l,m))/dmu(lp1)-
     +              (tphi(k,l,m)-tphi(k,lm1,m))/dmu(l))+bb1*
     +              (tphi(k,lp1,m)-tphi(k,lm1,m))+tphi(k,l,m)
            enddo
          enddo
        enddo
c     at amu=-1, the only term remaining is the first derivative 
c     from the pitch angle diffusion term which we approximate by 
c     a forward difference


      do k=0,ne
         do m=0,ntau
            b4=rlmbda(m)*rb3g2(k)
            c(k,-nmu,m)=-b4*dy/dmu(-nmu+1)
            b(k,-nmu,m)=1.-c(k,-nmu,m)
            r(k,-nmu,m)=-c(k,-nmu,m)*(tphi(k,-nmu+1,m)-tphi(k,-nmu,m))
     +           +tphi(k,-nmu,m)
c     at amu=1, use the backward difference
            a(k,nmu,m)=-b4*dy/dmu(nmu)
            b(k,nmu,m)=1.-a(k,nmu,m)
            r(k,nmu,m)=a(k,nmu,m)*(tphi(k,nmu,m)-tphi(k,nmu-1,m))
     +           +tphi(k,nmu,m)
c
c           call cntridag(a,b,c,r,nmu,tphi1)
c inline cntridag
c     subroutine cntridag(a,b,c,r,n,tphi1)
c
            n=nmu
c-          if (b(k,-n,m).eq.0) stop
            bet=b(k,-n,m)
            tphi1(k,-n,m)=r(k,-n,m)/bet
            do j=-n+1,n
              gam(k,j,m)=c(k,j-1,m)/bet
              bet=b(k,j,m)-a(k,j,m)*gam(k,j,m)
c--           if (bet.eq.0) pause
              tphi1(k,j,m)=(r(k,j,m)-a(k,j,m)*tphi1(k,j-1,m))/bet
            enddo


            do j=n-1, -n,-1
              tphi1(k,j,m)=tphi1(k,j,m)-gam(k,j+1,m)*tphi1(k,j+1,m)
           enddo
        enddo
      enddo
c
      phi(out,0:ne,-nmu:nmu,0:ntau)=tphi1(0:ne,-nmu:nmu,0:ntau)
c
      return
      end