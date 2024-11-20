
*****************************************************************************
*          Molecular Dynamics code to simulate at NVE collectivity 
*          a system of N atoms of Ar inside a cubic box, interacting
*          through a truncated Lennard-Jones force field at 2.5 sigma. 
*          Leap-frog Verlet integration algorithm.
*
*****************************************************************************
      PROGRAM leapfroglj
      implicit double precision(a-h,o-z)
      double precision mass
      REAL*8 :: nid
    
      INTEGER, DIMENSION(25) :: time0 
      REAL*8, DIMENSION(500,25) :: x0, y0, z0, vx0, vy0,vz0 
      
      REAL*8, DIMENSION (101) :: Integral
      REAL*8, DIMENSION (250) :: vacfcut
      REAL*8 :: simpson, h
      INTEGER :: kk, omega

c     1. Defining dimensions
 
      dimension r(3,1000),raux(3,1000),vinf(3,1000),accel(3,1000),
     &          vinfaux(3,1000),accelaux(3,1000)
      dimension g(100)
      dimension ntime(1000),vacf(1000),r2t(1000)
      dimension tspectral(1000)

c     2. Reading data and computing related quantities

      open(1,file='leap-lj.data',status='old')
      read(1,*) nconf
      read(1,*) natoms
      read(1,*) mass
      read(1,*) sigma,epsil
      read(1,*) deltat 
      close(1)
      
      pi = 3.14159265359
      nf = 3*natoms-3 !number of degrees of freedom, since we impose the center of masses not to leave the box, hence we have -3 degrees of freedom
      rc = 2.5d0 !sigma = range of the potential in reduced units
      rc1 = rc*sigma !! Not in reduced units, do not use it
      

c     3. Reading initial configuration (positions, velocities) in A and A/ps

      open(2,file='leap-conf.data',status='old')
      do is = 1,natoms
         read(2,*) (r(l,is),l=1,3)
         read(2,*) (vinf(l,is),l=1,3)
      end do         
      read(2,*) boxlength
      close(2)
c
c     Opening files to write results
c
      open(3,file='energy-leap.dat',status='unknown')
      open(4,file='temp-leap.dat',status='unknown')
      open(5,file='gr.dat',status='unknown')
      open(7,file='pressure.dat',status='unknown')
      open(8,file='energy_with_tail.dat',status='unknown')
      open(9,file='Vacf_and_MSD.dat',status='unknown')
      open(10,FILE="C(w).txt")

c     4. Change to reduced units

      call reduced(natoms,r,vinf,boxlength,deltat,epsil,sigma,
     &mass,uvel)

c     5. Initialization of some variable and vectors.

      ntt0 = 0
      nt0 = 0
      ndelt = 0
      
      rho = DFLOAT(natoms)/(boxlength**3) !! Converts to real, already in reduced units
      delg = boxlength/(2*DFLOAT(100)) ! bin size. Over 2 because we take half of the box's distance
      do i = 1,100 ! 100 = nhis = total number of bins
        g(i) = 0.d0
      enddo
     
      !! initialization for the calculation of VACF and MSD
      dtime = deltat !(time interval between two samples)
      do i = 1, 1000 !(1000 = tmax = total number of time steps of the VACF)
        ntime(i) = 0 !(ntime = number of functions to average for time ‘i’)
        vacf(i) = 0
        r2t(i) = 0
      enddo
      do i = 1, 25 
        time0(i) = 0
      enddo
      do i = 1, natoms
        do j = 1, 25
          x0(i,j) = 0 
          y0(i,j) = 0
          z0(i,j) = 0
          vx0(i,j) = 0 
          vy0(i,j) = 0
          vz0(i,j) = 0
        enddo 
      enddo
      
      raux = r 
      vinfaux =vinf ! raux and vinfaux are needed to avoid applying boundary conditions, r and vinf have BC
      
c     6. Calculate tail correction for energy and pressure and pressures-energy   
      
      !! Calculate formula in reduced units!!!!!!!!!!!   
      deltappot = 16.d0/3.d0*pi*(rho**2)*
     &            (2.d0/3.d0*(1/rc)**9-(1/rc)**3)  
   
      deltaE = (8.d0*pi/3.d0)*natoms*rho*
     &         (1.d0/3.d0*(1/rc)**9-(1/rc)**3) !! Tail correction for the energy
      
    
c     7. Start the loop to generate new configurations    

      do i = 1, nconf
      ppot = 0.d0
      pkin = 0.d0
         call forces(natoms,r,boxlength,accel,rc,epot,g,delg,
     &              ppot)
         call forces(natoms,raux,boxlength,accelaux,rc,epot,g,delg,
     &              ppot) ! Second call of forces function because for the vacf and msd calculation the forces to be applied shouldn't have boundary conditions
         call velpos(natoms,vinf,accel,deltat,r,nf,ecin,temp,
     &              boxlength,pkin,i,ntime,vacf,r2t,x0,y0,z0,
     &              vx0,vy0,vz0,time0,nt0,ntt0,ndelt,raux,
     &              vinfaux,accelaux)

         etot = ecin + epot + deltaE !! We add tail correction for the energy
         ppot = ppot + deltappot
         ptot = ppot + pkin 
        
         write(3,*) i*deltat, etot 
         write(4,*) i*deltat, temp
         write(7,*) i*deltat, ppot, pkin, ptot, deltappot
         write(8,*) i*deltat, epot, ecin, etot, deltaE
      enddo
      
c     8. Make the average of certain magnitudes 
      
      do i=1,100  !100 = nhis = total number of bins
          vb = ((i+1)**3-i**3)*delg**3 !volume
          nid = (4.0/3.0)*(4.d0*DATAN(1.d0))*vb*rho  !number of ideal gas particles in vb
          g(i) =g(i) /(dfloat(nconf*natoms)*nid)  !normalize g(r) and divide over nconf since we are averaging over all the times that we have "added up" g(r)
      enddo

      do ii = 1, 99 !! to save all the components of g(r)
         write(5,*) g(ii), delg*(ii-1) 
      enddo
     
      do i = 1,1000
       time = dtime*(i-1) !! (new time)
       tspectral(i) = time
       vacf(i) = vacf(i)/(DFLOAT(natoms*ntime(i))) !! (average of vacf)
       r2t(i) = r2t(i)/(DFLOAT(natoms*ntime(i))) !! (average of msd)
      enddo
      
cc    9. Save data an write it   
      do i = 1,1000
        write(9,*) dtime*(i-1), vacf(i), r2t(i), ntime(i)
      enddo
      
      do i = 1,250
        vacfcut(i) = vacf(i)
      enddo
      
      !! Calculate S(w) using simpson's integration
      h=deltat
      DO omega = 0,100
      	  !!!!!!! omega = (100/1000)*(kk-1)
          Integral(omega+1) = simpson(vacfcut,tspectral,h,omega)
          WRITE(10,*) Integral(omega+1)
      ENDDO
 
      close(3)
      close(4)
      close(5)
      close(7)
      close(8)
      close(9)
      close(10)

      open(11,file='newconf.data',status='unknown')
      do is = 1,natoms
         write(11,*) (r(l,is)*sigma,l=1,3)
         write(11,*) (vinf(l,is)*uvel,l=1,3)
      end do         
      write(11,*) boxlength*sigma
      close(11)

      stop
      end

*********************************************************
*********************************************************
c              subroutine reduced
*********************************************************
*********************************************************

      subroutine reduced(natoms,r,vinf,boxlength,deltat,
     &epsil,sigma,mass,uvel)
      implicit double precision(a-h,o-z)
      double precision mass
      dimension r(3,1000),vinf(3,1000)

      rgas = 8.314472673d0 !J/(mol*K) 
      utime = sigma*dsqrt(mass/epsil)*dsqrt(10.d0/rgas)
      uvel = sigma/utime !unit of velocity, expressed in A/ps

      boxlength = boxlength/sigma
      deltat = deltat/utime
      do is = 1,natoms
         do l = 1,3
            r(l,is) = r(l,is)/sigma
            vinf(l,is) = vinf(l,is)/uvel 
         end do
      end do

      return
      end

*********************************************************
*********************************************************
c              subroutine forces
*********************************************************
*********************************************************

      subroutine forces(natoms,r,boxlength,accel,rc,epot,g,
     &delg,ppot)
      implicit double precision(a-h,o-z)
      dimension r(3,1000),accel(3,1000)
      dimension g(100)

      do is = 1,natoms
         do l = 1,3
            accel(l,is) = 0.d0 !sets accelerations to 0
         end do 
      end do 
      epot = 0.d0 

c      atom-atom interactions

      do is = 1,natoms-1
         do js = is+1,natoms
            call lj(is,js,r,boxlength,accel,rc,pot,g,delg, ppot)
            epot = epot + pot
         end do
      end do

      return
      end

*********************************************************
*********************************************************
c              subroutine Lennard-Jones
*********************************************************
*********************************************************

      subroutine lj(is,js,r,boxlength,accel,rc,pot,g,delg, ppot)
      implicit double precision(a-h,o-z)
      dimension r(3,1000),accel(3,1000)
      dimension rij(3),g(100)

      rr2 = 0.d0
      pot = 0.d0
      do l = 1,3
         rijl = r(l,js) - r(l,is)
         rij(l) = rijl - boxlength*dnint(rijl/boxlength)
         rr2 = rr2 + rij(l)*rij(l)
      end do

      rr = dsqrt(rr2)
      if (rr.lt.boxlength/2) then ! only within half the box length
        ig = int(rr/delg)
        g(ig) = g(ig) +2 ! contribution for particle i and j. We add two because we are only looking at each pair once, namely particle 1 and 4. But 4 also has 1 at a ditance r. Since we'll divide by num of particles later on, we add the contribution of all the pairs, at once
      endif

      if (rr.lt.rc) then
         ynvrr2 = 1.d0/rr2
         ynvrr6 = ynvrr2*ynvrr2*ynvrr2
         ynvrr12 = ynvrr6*ynvrr6
         forcedist = 24.d0*(2.d0*ynvrr12-ynvrr6)*ynvrr2
         pot = 4.d0*(ynvrr12-ynvrr6)  
         do l = 1,3 
            accel(l,is) = accel(l,is) - forcedist*rij(l)
            accel(l,js) = accel(l,js) + forcedist*rij(l)
            ppot = ppot + 1.0/(3.d0*boxlength**3)*(rij(l)*rij(l)
     &      *forcedist)
         end do

      end if

      return
      end
      

*********************************************************
*********************************************************
c              subroutine velpos
*********************************************************
*********************************************************

c       Calculating velocity at instants t+delta/2 and t 
c       and getting temporal position at instant t + deltat.

      subroutine velpos(natoms,vinf,accel,deltat,r,nf,ecin,temp,
     &                  boxlength,pkin,iitcounter,ntime,vacf,r2t,
     &                  x0,y0,z0,vx0,vy0,vz0,time0,nt0,ntt0,ndelt,
     &                  raux,vinfaux,accelaux)

      implicit double precision(a-h,o-z)
      REAL*8, DIMENSION(natoms, 25) :: x0,y0,z0,vx0,vy0,vz0 
      INTEGER, DIMENSION(25) :: time0
      dimension vinf(3,1000),accel(3,1000),r(3,1000),raux(3,1000),
     &          vinfaux(3,1000),accelaux(3,1000)
      dimension ntime(1000),vacf(1000),r2t(1000)
      it0 = 50
      nt0max = 25
      ntmax = 1000
      

      ecin = 0.d0
      do is = 1,natoms
         v2 = 0.d0
         do l = 1,3 ! We update the positions with and without BC
            vsup = vinf(l,is) + accel(l,is)*deltat
            vsupaux = vinfaux(l,is) + accelaux(l,is)*deltat
            r(l,is) = r(l,is) + vsup*deltat
            raux(l,is) = raux(l,is) + vsupaux*deltat
            v = (vsup+vinf(l,is))/2.d0
            v2 = v2 + v*v
            vinf(l,is) = vsup
            vinfaux(l,is) = vsupaux
         end do
         ecin = ecin + 0.5d0*v2
         pkin = pkin + 1/(3.d0*boxlength**3)*v2
      end do
      
      temp = 2.d0*ecin/dfloat(nf)
      
      !! calculate a sample of data
        if (mod(iitcounter-1,it0).eq.0) then !tcounter-1 because otherwise you wouldn't count the first "window"
          nt0 = nt0 + 1 ! (update number of initial time)
          ntt0 = mod((nt0)-1,(nt0max))+ 1 
          time0((ntt0)) = iitcounter !! (store time of t=0)
         do i = 1, natoms
           x0(i,(ntt0)) = raux(1,i) !!(save position at t=0, see #3)
           vx0(i,(ntt0)) = vinfaux(1,i) !!(save velocity at t=0, see #3)
           y0(i,(ntt0)) = raux(2,i) !!(save position at t=0, see #3)
           vy0(i,(ntt0)) = vinfaux(2,i) !!(save velocity at t=0, see #3)
           z0(i,(ntt0)) = raux(3,i) !!(save position at t=0, see #3)
           vz0(i,(ntt0)) = vinfaux(3,i) !!(save velocity at t=0, see #3)
         enddo
        endif
        do iit=1,min((nt0),(nt0max)) !(update vacf and r2 for t=0)
         ndelt = iitcounter-time0(iit)+1 !(shift time: actual time - t=0)
         if ((ndelt).lt.ntmax) then
          ntime((ndelt)) = ntime((ndelt)) + 1
          do i = 1,natoms
           vacf((ndelt)) = vacf((ndelt)) + vinfaux(1,i)*vx0(i,iit) 
     &        +vinfaux(2,i)*vy0(i,iit) + vinfaux(3,i)*vz0(i,iit)!(update vacf)
           r2t((ndelt)) = r2t((ndelt)) + (raux(1,i)-x0(i,iit))**2+ 
     &        (raux(2,i)-y0(i,iit))**2 + (raux(3,i)-z0(i,iit))**2 !(update msd)
          enddo
         endif
        enddo
ccccc
ccccc       Applying  periodic boundary conditions. IMPORTANT: Only to r not raux.
ccccc
       do is=1,natoms
          do l=1,3
             if (r(l,is).lt.0) r(l,is) = r(l,is) + boxlength
             if (r(l,is).gt.boxlength) r(l,is) = r(l,is) - boxlength
          enddo
       enddo
ccccc
        
      return
      end  
 
*********************************************************
*********************************************************
c              simpson
*********************************************************
*********************************************************
      
      REAL*8 FUNCTION simpson(corr,t,h,omega)
      REAL*8, DIMENSION (4) :: y
      REAL*8 :: term, h
      REAL*8, DIMENSION (250) :: corr
      REAL*8, DIMENSION (1000) :: t
      INTEGER :: omega
      
          DO ii = 1, FLOOR((250-1)/3.) ! 250-1/3 since it's the number of summing terms
              y(1) = corr(3*(ii-1)+1)*cos(omega*t(3*(ii-1)+1))
              y(2) = corr(3*(ii-1)+2)*cos(omega*t(3*(ii-1)+2))
              y(3) = corr(3*(ii-1)+3)*cos(omega*t(3*(ii-1)+3))
              y(4) = corr(3*(ii-1)+4)*cos(omega*t(3*(ii-1)+4))
              term=term+h/8.*(3*y(1)+9*y(2)+9*y(3)+3*y(4))
          ENDDO
      simpson=term
      RETURN
      END
