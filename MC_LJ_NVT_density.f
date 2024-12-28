*****************************************************************************
*          Molecular Dynamics code to simulate at NVE collectivity 
*          a system of N atoms of Ar inside a cubic box, interacting
*          through a truncated Lennard-Jones force field at 2.5 sigma. 
*          Leap-frog Verlet integration algorithm.
*
*****************************************************************************
      PROGRAM MC_LJ_NVT
      implicit double precision(a-h,o-z)
      double precision mass
      REAL*8 :: nid
    

c     1. Defining dimensions
 
      dimension r(3,1000),raux(3,1000),vinf(3,1000),accel(3,1000),
     &          vinfaux(3,1000),accelaux(3,1000)
      dimension g(100)

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

      open(2,file='leap-conf.data',status='old') ! We read the last configuration of the last rho (if we want to calculate rho=0.7, we use last conf of rho=0.8)
      do is = 1,natoms
         read(2,*) (r(l,is),l=1,3)
         read(2,*) (vinf(l,is),l=1,3) !OBS: Now we read the velocities for completenness, even though it's not necessary.
      end do         
      read(2,*) boxlength
      close(2)
      
c     Opening files to write results
c
      open(3,file='thermodynamics_rho_08.dat',status='unknown')
      open(4,file='gr_rho_08.dat',status='unknown')
      open(5,file='last_config_rho_08.dat',status='unknown')


c     4. Change to reduced units
      
      
      call reduced(natoms,r,vinf,boxlength,deltat,epsil,sigma,
     &mass,uvel)

c   


cc    MC STARTS HERE
      ncycl = 1000000
      nsamp = 1000
      ndumped = 0
      delta= 0.25
      beta=1/(2.d0) !2 is the temperature of the isotherm. Since we are working on the (N,V,T) collectivity, T remains constant, hence the MC method will make the system evolve towards thr equilibibriuum temperaature.
      rhoo=0.83
      rhof=0.8
      
      boxlength=boxlength*((rhoo/rhof)**(1./3.))
      ! rc=rc*((0.83/rhof)**(1./3.)) ! rc must not be changed, it's a parameter of the LJ potential
      rho = DFLOAT(natoms)/(boxlength**3) !! Converts to real, already in reduced units
      delg = boxlength/(2*DFLOAT(100)) ! bin size. Over 2 because we take half of the box's distance
      do i = 1,100 ! 100 = nhis = total number of bins
        g(i) = 0.d0
      enddo
      
      do is = 1,natoms      
        do l = 1,3
          r(l,is)=r(l,is)*((rhoo/rhof)**(1./3.))
        end do
      end do
      
c     Calculate tail correction for energy and pressure and pressures-energy   
      
      !! Calculate formula in reduced units !!   
      deltappot = 16.d0/3.d0*pi*(rho**2)*
     &            (2.d0/3.d0*(1/rc)**9-(1/rc)**3)  
   
      deltaE = (8.d0*pi/3.d0)*natoms*rho*
     &         (1.d0/3.d0*(1/rc)**9-(1/rc)**3) !! Tail correction for the energy
      !! NOTE: When calcuating Thermo properties from the computed g(r), this tail has to be changed rc--> rlim (the value of r where g(r) stabilized to one) 
      
      do icycl=1,ncycl !! perform ncycl MC cycles
        call mcmove(r,rc,boxlength,natoms,beta,ndumped,delta) !! displace a particle
        if(mod(icycl,nsamp).eq.0) then !! it makes no sense to calculate the properties every one particle movement because the properties will be almost equal, so every nsamp movements we take the configuration.
          call sample(natoms,r,boxlength,rc,epot,g,
     &                delg,ppot) !! sample averages
          epot = epot + deltaE !! We add tail correction for the energy
          ptot = ppot + rho/beta + deltappot
        
          write(3,*) icycl/nsamp, ptot, epot
        endif
      enddo
      naccepted = dfloat(ncycl - ndumped)/dfloat(ncycl)*100.d0 ! We compute the percentage of accepted trial moves. Important to readjust delta
      write(3,*) naccepted, 0, rho
      
      ! We only change/ compute g when we recalculate termo properties, since it doesnt change much on every iteration.
      ! Hence nconf has to be ncycl/nsamp
      nconf = ncycl/nsamp 
      do i=1,100  !100 = nhis = total number of bins
          vb = ((i+1)**3-i**3)*delg**3 !volume
          nid = (4.0/3.0)*(4.d0*DATAN(1.d0))*vb*rho  !number of ideal gas particles in vb
          g(i) =g(i) /(dfloat(nconf*natoms)*nid)  !normalize g(r) and divide over nconf since we are averaging over all the times that we have "added up" g(r)
      enddo

      do ii = 1, 99 !! to save all the components of g(r)
         write(4,*) g(ii), delg*(ii-0.5) 
      enddo
      
      do is = 1,natoms
        do l = 1,3
          r(l,is)=sigma*r(l,is)
          vinf(l,is)=uvel*vinf(l,is)
        end do   
      end do
      
      do is = 1,natoms
         write(5,*) (r(l,is),l=1,3)
         write(5,*) (vinf(l,is),l=1,3)
      end do   
      boxlength = boxlength*sigma
      write(5,*) boxlength
      
      close(3)
      close(4)
      close(5)
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
c              ener
*********************************************************
*********************************************************

      subroutine ener(r,rc,boxlength,natoms,iio,epot)
      implicit double precision(a-h,o-z)
      dimension r(3,1000)
      dimension rij(3),g(100)
      epot = 0.d0 
      do is = 1,natoms ! since we only need the diference between old and new configurations this means we compute only the enrgy of the affected particle.
        if (is.ne.iio) then
            rr2 = 0.d0
            do l = 1,3
              rijl = r(l,iio) - r(l,is)
              rij(l) = rijl - boxlength*dnint(rijl/boxlength)
              rr2 = rr2 + rij(l)*rij(l)
            end do
            rr = dsqrt(rr2)
            if (rr.lt.rc) then
              ynvrr2 = 1.d0/rr2
              ynvrr6 = ynvrr2*ynvrr2*ynvrr2
              ynvrr12 = ynvrr6*ynvrr6
              forcedist = 24.d0*(2.d0*ynvrr12-ynvrr6)*ynvrr2
              epot = epot + 4.d0*(ynvrr12-ynvrr6)  
            end if
        end if
      end do
      return
      end
      


*********************************************************
*********************************************************
c              mcmove
*********************************************************
*********************************************************
      
      
      SUBROUTINE mcmove(r,rc,boxlength,natoms,beta,ndumped,delta) !! Attempts to displace a particle
      implicit double precision(a-h,o-z)
      dimension ro(3)
      dimension r(3,1000)
      iio=int(rand()*natoms)+1 !! select a particle at random
      
      call ener(r,rc,boxlength,natoms,iio,epoto) !! energy old configuration
      do l=1,3 
        ro(l)=r(l,iio)
        r(l,iio) = ro(l)+ (rand()-0.5)*delta !! give particle random displacement
      !!boundary conditions
         if (r(l,iio).lt.0) r(l,iio) = r(l,iio) + boxlength
         if (r(l,iio).gt.boxlength) r(l,iio) = r(l,iio) - boxlength
      enddo
      call ener(r,rc,boxlength,natoms,iio,epotn) !! energy new configuration
      if(rand().gt.exp(-beta*(epotn-epoto))) then !! inverse acceptance law
        ndumped = ndumped +1 !!COntador de cuantas configuraciones recghazamos. Lo necesitamos para saber si la delta elegida es adecuado.
        do l=1,3 
          r(l,iio) = ro(l) !! give particle random displacement !! reverse the replacement r(o) by rn
        enddo
      endif
      
      return
      end
      
*********************************************************
*********************************************************
c              subroutine sample
*********************************************************
*********************************************************

      subroutine sample(natoms,r,boxlength,rc,epot,g,
     &                  delg,ppot)
      implicit double precision(a-h,o-z)
      dimension r(3,1000)
      dimension g(100)
      epot = 0.d0 
      ppot = 0.d0
c      atom-atom interactions
      do is = 1,natoms-1
         do js = is+1,natoms
            call ljMC(is,js,r,boxlength,rc,pot,g,delg,ppot)
            epot = epot + pot
         end do
      end do
      return
      end

*********************************************************
*********************************************************
c        subroutine Lennard-Jones-Montecarlo
*********************************************************
*********************************************************

      subroutine ljMC(is,js,r,boxlength,rc,pot,g,delg,ppot)
      implicit double precision(a-h,o-z)
      dimension r(3,1000)
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
            ppot = ppot + 1.0/(3.d0*boxlength**3)*(rij(l)*rij(l)
     &      *forcedist)
         end do

      end if

      return
      end



