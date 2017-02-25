      subroutine tauup(in,out,symloop)


c
c  phi(k,l,m,in) is updated to phi(k,l,m,out) for the tau term using
c  monotonic transport.  Boundary conditions for positive
c  velocity at m=0 are either symetric or no flux comes in from m < 0;
c  for negative velocity the boundary condition at m=ntau is that
c  no flux comes in.
c  Note:  Nonuniform step sizes in tau are not allowed. 
c
      include 'fkrplk_par.h'
      common /ctau/ dtau(ntaumax), rlmbda2(ntaumax)
      real tphi(-nmumax:nmumax,0:ntaumax),
     +     tphi1(-nmumax:nmumax,0:ntaumax)
      logical symloop
      integer in,out
      eps=1e-20
      a0 = dy/dtau(1)
      do 1 k=0,ne
         a1 = a0*beta(k)
         do 7 m=0,ntau
            do 6 l=-nmu,nmu
               tphi(l,m) = phi(in,k,l,m)
               tphi1(l,m) = tphi(l,m)
 6          continue
 7       continue
c
c  NEGATIVE VELOCITIES
c
         do 2 l=-nmu,-1
            a = a1*amu(l)
            a1pa = a*(1.+a)
c
c  m=ntau:  For this point we the boundary condition that 
c           phi = 0. for m > ntau
c
            phic = tphi(l,ntau)
            phin = tphi(l,ntau-1)
            dphi = phic - phin
            delt = -dphi*phic
c
c replace if statement below
c
             si=sign(1.,delt)
             delt=-delt/(phin+eps)*(si+1.)/2.
             tphi1(l,ntau) = phic*(1.+a) - a1pa*delt*(si+1.)/2.
c
c           if(delt .gt. 0.)then
c              delt = -delt/phin
c              tphi1(l,ntau) = phic*(1.+a) - a1pa*delt
c           else
c              delt = 0.
c              tphi1(l,ntau) = phic*(1.+a) 
c           endif


            do 3 m=ntau-1,1,-1
               phip = phic
               phic = phin
               phin = tphi(l,m-1)
               dphip = dphi
               dphi = phic - phin
               deltp = delt
               delt = dphi*dphip
c
c   replace if statement below
c
               si=sign(1.,delt)
               delt=delt/(phip - phin+eps)*(si+1.)/2.
c
c              if(delt .gt. 0.)then
c                 delt=delt/(phip - phin)
c              else
c                 delt = 0.
c              endif
               tphi1(l,m) = phic - a*dphip + a1pa*(deltp - delt)
 3          continue
c
c  m=0:  no problem for negative velocities, but the correction
c        to the upwind method depends on the symmetry of the loop.
c
            phip = phic
            phic = phin
            dphip = dphi
            deltp = delt
c
c ass. symloop=.true.
c
c           if(symloop)then
               dphi = phic - tphi(-l,1)
               delt = dphi*dphip
               si=sign(1.,delt)
               delt=delt/(phip - tphi(-l,1)+eps)*(si+1.)/2.
c              if(delt .gt. 0.)then
c                 delt = delt/(phip - tphi(-l,1))
c              else
c                 delt = 0.
c              endif
c            else
c               delt = tphi(l,0)*dphip
c               si=sign(1.,delt)
c               delt= delt/(tphi(l,1)+eps)*(si+1.)/2.
cc              if(delt .gt. 0.)then
cc                 delt = delt/tphi(l,1)
cc              else
cc                 delt = 0.
cc              endif
c            endif
            tphi1(l,0) = phic - a*dphip + a1pa*(deltp - delt)
 2       continue


c
c  POSITIVE VELOCITIES
c
         do 8 l=1,nmu
            a = a1*amu(l)
            a1ma = a*(1.-a)
c
c  m=0:  The treatment of this point depends on the loop symmety.
c        If symmetric, we use phi(k,l,-m)=phi(k,-l,m) for m > 0.
c        If not symmetric, then the we use phi = 0. for m < 0.
c
            phic = tphi(l,0)
            phin = tphi(l,1)
c
c assume symloop=.true.
c
c            if(symloop)then
               dphi = phic - tphi(-l,1)
               dphip = phin - phic
               delt = dphi*dphip
               si=sign(1.,delt)
               delt= delt/(phin - tphi(-l,1)+eps)*(si+1.)/2.
c              if(delt .gt. 0.)then
c                 delt = delt/(phin - tphi(-l,1))
c              else
c                 delt = 0.
c              endif
               delto = dphi*(tphi(-l,1)-tphi(-l,2))
               si=sign(1.,delto)
               delto=delto/(phic - tphi(-l,2)+eps)*(si+1.)/2. 
c              if(delto .gt. 0.)then
c                 delto = delto/(phic - tphi(-l,2))
c              else
c                 delto = 0.
c              endif
               tphi1(l,0) = phic - a*dphi - a1ma*(delt - delto)
c            else
c               dphip = phin - phic
c               delt = phic*dphip
c               si=sign(1.,delt)
c               delt= delt/(phin+eps)*(si+1.)/2.
cc              if(delt .gt. 0.)then
cc                 delt = delt/phin
cc              else
cc                 delt = 0.
cc              endif
c               tphi1(l,0) = phic*(1.-a) - a1ma*delt
c            endif


            do 9 m=1,ntau-1
               phip = phic
               phic = phin
               phin = tphi(l,m+1)
               delto = delt
               dphi = dphip
               dphip = phin - phic
               delt = dphi*dphip
               si=sign(1.,delt)
               delt=delt/(phin - phip+eps)*(si+1.)/2.
c              if(delt .gt. 0.)then
c                 delt=delt/(phin - phip)
c              else
c                 delt = 0.
c              endif
               tphi1(l,m) = phic - a*dphi - a1ma*(delt - delto)
 9          continue
            phip = phic
            phic = phin
            delto = delt
            dphi = dphip
            delt = -phic*dphi
            si=sign(1.,delt)
            delt= -delt/(phip+eps)*(si+1.)/2.
c           if(delt .gt. 0.)then
c              delt = -delt/phip
c           else
c              delt = 0.
c           endif
            tphi1(l,ntau)=phic- a*dphi - a1ma*(delt-delto)
 8       continue
         do 17 m=0,ntau
            do 16 l=-nmu,nmu
               phi(out,k,l,m) = tphi1(l,m)
 16         continue
 17      continue
 1    continue
      return
      end 