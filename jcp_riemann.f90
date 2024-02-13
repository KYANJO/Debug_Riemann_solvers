program main
    implicit none
    
    ! Declarations
    integer :: ixy, maxm, meqn, mwaves, maux, mbc, mx, mcapa,i
    integer coordinate_system, imp,M,n,j,mw
    double precision dx, maxcfl, cfl, dt_est, dt,t
    double precision currenttime, maxspeed
    double precision t0, Tfinal, ax, bx, ay, by
    double precision :: drytol, earth_radius, deg2rad, g
    double precision, allocatable :: ql(:,:), qr(:,:), auxl(:,:), auxr(:,:),aux1(:,:), aux2(:,:),aux3(:,:)
    double precision, allocatable :: fwave(:,:,:), s(:,:), amdq(:,:), apdq(:,:)
    double precision, allocatable :: asdq(:,:), bmasdq(:,:), bpasdq(:,:)
    double precision, allocatable :: xe(:), x(:)
    double precision, allocatable :: h(:), hu(:), hv(:)
    double precision, allocatable :: qoldl(:,:), qoldr(:,:)

    logical :: cfl_violated

    cfl_violated =.false.

    ! Initialize variables - Replace with actual values
    ixy = 1
    maxm = 16
    meqn = 3
    mwaves = 3
    maux = 1
    mbc = 2
   !  ax = 953236
   !  bx = 959554
   !  ay = 1832407.25
   !  by = 1848572.75
   ax = -4
   bx = 4
   ay = -2
   by = 2
    cfl = 0.75d0
    t0 = 0.0d0
    Tfinal = 10d0
    
    ! Compute dx and dt_est
    dx = (bx - ax) / maxm
    
    ! Allocate and initialize arrays
    allocate(ql(meqn, 1-mbc:maxm+mbc))
    allocate(qr(meqn, 1-mbc:maxm+mbc))
    allocate(qoldl(meqn, 1-mbc:maxm+mbc))
    allocate(qoldr(meqn, 1-mbc:maxm+mbc))
    allocate(auxl(maux, 1-mbc:maxm+mbc))
    allocate(auxr(maux, 1-mbc:maxm+mbc))
    allocate(fwave(meqn, mwaves, 1-mbc:maxm+mbc))
    allocate(s(mwaves, 1-mbc:maxm+mbc))
    allocate(amdq(meqn, 1-mbc:maxm+mbc))
    allocate(apdq(meqn, 1-mbc:maxm+mbc))

    allocate(aux1(maux, 1-mbc:maxm+mbc))
    allocate(aux2(maux, 1-mbc:maxm+mbc))
    allocate(aux3(maux, 1-mbc:maxm+mbc))
    allocate(asdq(meqn, 1-mbc:maxm+mbc))
    allocate(bmasdq(meqn, 1-mbc:maxm+mbc))
    allocate(bpasdq(meqn, 1-mbc:maxm+mbc))
    allocate(h(1-mbc:maxm+mbc), hu(1-mbc:maxm+mbc), hv(1-mbc:maxm+mbc))
  
    
    ! Initialize arrays with dummy values for demonstration purposes
    mx = maxm
    mcapa = 0
    drytol = 0.001d0
    earth_radius = 6371000.0d0  ! Approximate radius of earth in meters
    deg2rad = 4.0d0*atan(1.0d0) / 180.0d0
    g = 9.81d0  ! Acceleration due to gravity
    coordinate_system = 1
    imp = 1

   ! ql(1,:) = 0.4473496224909771D2
   ! ql(2,:) = -0.1018126050977528D-13
   ! ql(3,:) = -0.6489769399271391D-9

   ! qr(1,:) = 0.2519771920562447D2
   ! qr(2,:) = -0.2016857953898484D-13
   ! qr(3,:) = -0.3325890167957528D-9

   ! auxl(1,:) = 0.9999999999996810d2
   ! ! auxl(2,:) = 0.0d0

   ! auxr(1,:) = 0.9999999999996810d4
   ! auxr(2,:) = 0.0d0

   ql(1,:) = 2.0d0
   ql(2,:) = 0.0d0
   ql(3,:) = 0.0d0

   qr(1,:) = 1.0d0
   qr(2,:) = 0.0d0
   qr(3,:) = 0.0d0

   auxl(1,:) = 0.0d0;
   ! auxl(2,:) = 0.0d0

   auxr(1,:) = 0.0d0
   ! auxr(2,:) = 0.0d0

   ! transverse solver
    aux1(1,:) = 0.0d0
    aux2(1,:) = 0.0d0
    aux3(1,:) = 0.0d0
    
    aux1(2,:) = 0.0d0
    aux2(2,:) = 0.0d0
    aux3(2,:) = 0.0d0
    
    
    
    qoldl = ql
    qoldr = qr
    maxcfl = 0.0d0
    currenttime = 0.0d0
    dt = 0.01d0

   !  do n = 1,M
   do while (currenttime < Tfinal)

        h(1) = qoldl(1,1)
        hu(1) = qoldl(2,1)
        hv(1) = qoldl(3,1)
        h(2) = qoldr(1,2)
        hu(2) = qoldr(2,2)
        hv(2) = qoldr(3,2)

        cfl_violated =.false.  ! Reset cfl_violated flag each time through loop

        do j = 2,maxm
            ! Normal solver
            call fc2d_geoclaw_rpn2(ixy, maxm, meqn, mwaves, maux, mbc, mx, &
                           qoldl, qoldr, auxl, auxr, fwave, s, amdq, apdq, g, drytol, &
                           earth_radius, deg2rad, mcapa)

             i = j
            do mw=1,mwaves
               maxspeed = abs(s(mw,i))
               maxcfl = max(maxcfl,maxspeed*dt/dx)
            enddo

            t = currenttime + dt

            if (maxcfl > cfl) then
               dt = dt*0.5d0
               cfl_violated =.true.
               return
            end if

            ! if (imp == 1 .and. ixy == 1) asdq = amdq
            ! asdq = apdq
            ! call fc2d_geoclaw_rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx,&
            !     qoldl,qoldr,aux1,aux2,aux3,asdq,bmasdq,bpasdq,drytol,g,coordinate_system,&
            !     earth_radius,deg2rad)
      
            ! write out 
            ! do i=2-mbc,mx+mbc
           
           write(*,*) 's = ', s(1:3,i)
           write(*,*) 'fwave = ', fwave(1,1,i), fwave(2,1,i), fwave(3,1,i)
           write(*,*) 'fwave = ', fwave(1,2,i), fwave(2,2,i), fwave(3,2,i)
           write(*,*) 'fwave = ', fwave(1,3,i), fwave(2,3,i), fwave(3,3,i)
           write(*,*) 'amdq = ', amdq(1,i), amdq(2,i), amdq(3,i)
           write(*,*) 'apdq = ', apdq(1,i), apdq(2,i), apdq(3,i)
           write(*,*) ' '

            ! write(*,*) 'asdq = ', asdq(1,i), asdq(2,i), asdq(3,i)
            ! write(*,*) 'bmasdq = ', bmasdq(1,i), bmasdq(2,i), bmasdq(3,i)
            ! write(*,*) 'bpasdq = ', bpasdq(1,i), bpasdq(2,i), bpasdq(3,i)
            ! write(*,*) ' '
         !   stop
           h(j) = h(j) - (dt/dx)*(amdq(1,j) + apdq(1,j))
           hu(j) = hu(j) - (dt/dx)*(amdq(2,j) + apdq(2,j))
           hv(j) = hv(j) - (dt/dx)*(amdq(3,j) + apdq(3,j))

            ! h(j) = h(j) - (dt/dx)*(amdq(1,j) + apdq(1,j)) - (dt/dx)*(bmasdq(1,j) + bpasdq(1,j))
            ! hu(j) = hu(j) - (dt/dx)*(amdq(2,j) + apdq(2,j)) - (dt/dx)*(bmasdq(2,j) + bpasdq(2,j))
            ! hv(j) = hv(j) - (dt/dx)*(amdq(3,j) + apdq(3,j)) - (dt/dx)*(bmasdq(3,j) + bpasdq(3,j))
           
          ! update
          qoldl(1,:) = h(j-1)
          qoldl(2,:) = hu(j-1)
          qoldl(3,:) = hv(j-1)
          
          qoldr(1,:) = h(j)
          qoldr(2,:) = hu(j)
          qoldr(3,:) = hv(j)
          
          write(*,*) "dt = ", dt, "maxcfl = ", maxcfl, "t = ", t
          write(*,*) "qoldl = ", qoldl(1,j), qoldl(2,j), qoldl(3,j)
          write(*,*) "qoldr = ", qoldr(1,j), qoldr(2,j), qoldr(3,j)
        !   write(*,*) "dt = ", dt, " dx = ", dx, "dt/dx = ", dt/dx, " amdq = ", amdq(2,j), " apdq = ", apdq(2,j) 
        !   stop
        enddo
   if (cfl_violated) then
      if (t<Tfinal) then
          continue
      else
         write(*,*) "CFL violation at time ", t
         stop
      end if
   end if   

   ! Compute new time
   currenttime = currenttime + dt
   dt_est = maxcfl * dx/ (maxspeed + 1.d-10)
   dt = min(dt_est, Tfinal - currenttime)

enddo

    
!     ! Normal solver
!     call fc2d_geoclaw_rpn2(ixy, maxm, meqn, mwaves, maux, mbc, mx, &
!                           ql, qr, auxl, auxr, fwave, s, amdq, apdq, g, drytol, &
!                           earth_radius, deg2rad, mcapa)
    
!     ! write out 
!     ! do i=2-mbc,mx+mbc
!   i = 1
!       write(*,*) 's = ', s(1:3,i)
!       write(*,*) 'fwave = ', fwave(1,1,i), fwave(2,1,i), fwave(3,1,i)
!       write(*,*) 'fwave = ', fwave(1,2,i), fwave(2,2,i), fwave(3,2,i)
!       write(*,*) 'fwave = ', fwave(1,3,i), fwave(2,3,i), fwave(3,3,i)
!       write(*,*) 'amdq = ', amdq(1,i), amdq(2,i), amdq(3,i)
!       write(*,*) 'apdq = ', apdq(1,i), apdq(2,i), apdq(3,i)
!       write(*,*) ' '
    ! enddo

    
    
    ! Clean up
    deallocate(ql, qr, auxl, auxr, aux1, aux2, aux3, fwave, s, amdq, apdq, asdq, &
                bpasdq, bmasdq)
    ! deallocate(t, xe, x, h, hu, hv, qoldl, qoldr)
    deallocate(qoldl,qoldr)
end program main


! ======================================================================
       subroutine fc2d_geoclaw_rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx, &
                     ql,qr,auxl,auxr,fwave,s,amdq,apdq,g,drytol, &
                     earth_radius,deg2rad,mcapa) 
! ======================================================================

      ! use geoclaw_module, only: g => grav, drytol => dry_tolerance, rho
      ! use geoclaw_module, only: earth_radius, deg2rad
      ! use amr_module, only: mcapa

      implicit none

      !input
      integer maxm,meqn,maux,mwaves,mbc,mx,ixy

      double precision  fwave(meqn, mwaves, 1-mbc:maxm+mbc)
      double precision  s(mwaves, 1-mbc:maxm+mbc)
      double precision  ql(meqn, 1-mbc:maxm+mbc)
      double precision  qr(meqn, 1-mbc:maxm+mbc)
      double precision  apdq(meqn,1-mbc:maxm+mbc)
      double precision  amdq(meqn,1-mbc:maxm+mbc)
      double precision  auxl(maux,1-mbc:maxm+mbc)
      double precision  auxr(maux,1-mbc:maxm+mbc)

      !local only
      integer m,i,mw,maxiter,mu,nv
      double precision wall(3)
      double precision fw(3,3)
      double precision sw(3)

      double precision hR,hL,huR,huL,uR,uL,hvR,hvL,vR,vL,phiR,phiL,pL,pR
      double precision bR,bL,sL,sR,sRoe1,sRoe2,sE1,sE2,uhat,chat
      double precision s1m,s2m
      double precision hstar,hstartest,hstarHLL,sLtest,sRtest
      double precision tw,dxdc

      logical rare1,rare2
      
      double precision g, drytol, earth_radius, deg2rad
      integer mcapa

      logical debug
      double precision :: dtcom, dxcom, dycom, tcom
      integer :: icom, jcom
      common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom   

      if (ixy == 1) then
         debug = .true.
      else
         debug  = .false.
      endif

      ! In case there is no pressure forcing
      pL = 0.d0
      pR = 0.d0

      !loop through Riemann problems at each grid cell
      do i=2-mbc,mx+mbc

!-----------------------Initializing-----------------------------------
         !inform of a bad riemann problem from the start
         if((qr(1,i-1).lt.0.d0).or.(ql(1,i) .lt. 0.d0)) then
            write(*,*) 'Negative input: hl,hr,i=',qr(1,i-1),ql(1,i),i
         endif

         !Initialize Riemann problem for grid interface
         do mw=1,mwaves
              s(mw,i)=0.d0
                 fwave(1,mw,i)=0.d0
                 fwave(2,mw,i)=0.d0
                 fwave(3,mw,i)=0.d0
         enddo

!        !set normal direction
         if (ixy.eq.1) then
            mu=2
            nv=3
         else
            mu=3
            nv=2
         endif

         !zero (small) negative values if they exist
         if (qr(1,i-1).lt.0.d0) then
               qr(1,i-1)=0.d0
               qr(2,i-1)=0.d0
               qr(3,i-1)=0.d0
         endif

         if (ql(1,i).lt.0.d0) then
               ql(1,i)=0.d0
               ql(2,i)=0.d0
               ql(3,i)=0.d0
         endif

!        if (debug) then
!           write(6,*) 'i = ', i, ' j = ' , jcom
!           write(*,*) 'qr = ', qr(1,i-1), ' ql = ', ql(1,i)
!           write(*,*) 'qr = ', qr(2,i-1), ' ql = ', ql(2,i)
!           write(*,*) 'qr = ', qr(3,i-1), ' ql = ', ql(3,i)
!           write(6,*) ' '
!        endif 

         !skip problem if in a completely dry area
         if (qr(1,i-1) <= drytol .and. ql(1,i) <= drytol) then
            go to 30
         endif

         !Riemann problem variables
         hL = qr(1,i-1) 
         hR = ql(1,i) 
         huL = qr(mu,i-1) 
         huR = ql(mu,i) 
         bL = auxr(1,i-1)
         bR = auxl(1,i)
         
         hvL=qr(nv,i-1) 
         hvR=ql(nv,i)

         ! hL = ql(1,i-1) 
         ! hR = qr(1,i) 
         ! huL = ql(mu,i-1) 
         ! huR = qr(mu,i) 
         ! bL = auxl(1,i-1)
         ! bR = auxr(1,i)
         
         ! hvL=ql(nv,i-1) 
         ! hvR=qr(nv,i)
        

         !check for wet/dry boundary
         if (hR.gt.drytol) then
            uR=huR/hR
            vR=hvR/hR
            phiR = 0.5d0*g*hR**2 + huR**2/hR
         else
            hR = 0.d0
            huR = 0.d0
            hvR = 0.d0
            uR = 0.d0
            vR = 0.d0
            phiR = 0.d0
         endif

         if (hL.gt.drytol) then
            uL=huL/hL
            vL=hvL/hL
            phiL = 0.5d0*g*hL**2 + huL**2/hL
         else
            hL=0.d0
            huL=0.d0
            hvL=0.d0
            uL=0.d0
            vL=0.d0
            phiL = 0.d0
         endif

!        if (debug) then
!              write(6,*) 'i = ', i, ' j = ' , jcom
!              write(*,*) 'hL = ', hL, ' hR = ', hR
!              write(*,*) 'huL = ', huL, ' huR = ', huR
!              write(*,*) 'hvL = ', hvL, ' hvR = ', hvR
!              write(*,*) 'uL = ', uL, ' uR = ', uR
!              write(*,*) 'vL = ', vL, ' vR q= ', vR
!              write(*,*) 'phiL = ', phiL, ' phiR = ', phiR
!              write(*,*) 'bL = ', bL, ' bR = ', bR
!              write(6,*) ' '
!           endif

         wall(1) = 1.d0
         wall(2) = 1.d0
         wall(3) = 1.d0
         if (hR.le.drytol) then
            call riemanntype(hL,hL,uL,-uL,hstar,s1m,s2m, &
                      rare1,rare2,1,drytol,g)
!
!       if (debug) then
!           write(6,*) 'i = ', i, ' j = ' , jcom
!           write(*,*) 'hL = ', hL, ' uL = ', uL
!             write(*,*) 'hstar = ', hstar
!             write(*,*) 's1m = ', s1m, ' s2m = ', s2m
!             write(*,*) 'rare1 = ', rare1, ' rare2 = ', rare2
!             write(6,*) ' '
!          endif
!
            hstartest=max(hL,hstar)
            if (hstartest+bL.lt.bR) then !right state should become ghost values that mirror left for wall problem
!                bR=hstartest+bL
               wall(2)=0.d0
               wall(3)=0.d0
               hR=hL
               huR=-huL
               bR=bL
               phiR=phiL
               uR=-uL
               vR=vL
            elseif (hL+bL.lt.bR) then
               bR=hL+bL
            endif
         elseif (hL.le.drytol) then ! right surface is lower than left topo
            call riemanntype(hR,hR,-uR,uR,hstar,s1m,s2m, &
                                  rare1,rare2,1,drytol,g)

!            if (debug) then
!             write(6,*) 'i = ', i, ' j = ' , jcom
!             write(*,*) 'hL = ', hL, ' uL = ', uL
!               write(*,*) 'hstar = ', hstar
!               write(*,*) 's1m = ', s1m, ' s2m = ', s2m
!               write(*,*) 'rare1 = ', rare1, ' rare2 = ', rare2
!               write(6,*) ' '
!            endif
!
            hstartest=max(hR,hstar)
            if (hstartest+bR.lt.bL) then  !left state should become ghost values that mirror right
!               bL=hstartest+bR
               wall(1)=0.d0
               wall(2)=0.d0
               hL=hR
               huL=-huR
               bL=bR
               phiL=phiR
               uL=-uR
               vL=vR
            elseif (hR+bR.lt.bL) then
               bL=hR+bR
            endif
         endif

!        if (debug) then
!          write(6,*) 'i = ', i, ' j = ' , jcom
!          write(*,*) 'hL = ', hL, ' hR = ', hR
!          write(*,*) 'huL = ', huL, ' huR = ', huR
!          write(*,*) 'hvL = ', hvL, ' hvR = ', hvR
!          write(*,*) 'uL = ', uL, ' uR = ', uR
!          write(*,*) 'vL = ', vL, ' vR q= ', vR
!          write(*,*) 'phiL = ', phiL, ' phiR = ', phiR
!          write(*,*) 'bL = ', bL, ' bR = ', bR
!          write(6,*) ' '
!       endif

         !determine wave speeds
         sL=uL-sqrt(g*hL) ! 1 wave speed of left state
         sR=uR+sqrt(g*hR) ! 2 wave speed of right state

         uhat=(sqrt(g*hL)*uL + sqrt(g*hR)*uR)/(sqrt(g*hR)+sqrt(g*hL)) ! Roe average
         chat=sqrt(g*0.5d0*(hR+hL)) ! Roe average
         sRoe1=uhat-chat ! Roe wave speed 1 wave
         sRoe2=uhat+chat ! Roe wave speed 2 wave

         sE1 = min(sL,sRoe1) ! Eindfeldt speed 1 wave
         sE2 = max(sR,sRoe2) ! Eindfeldt speed 2 wave


         !--------------------end initializing...finally----------
         !solve Riemann problem.
!       if (debug) then
!          write(6,*) 'i = ', i, ' j = ' , jcom
!          write(*,*) 'sL = ', sL, ' sR = ', sR
!          write(*,*) 'sRoe1 = ', sRoe1, ' sRoe2 = ', sRoe2
!          write(*,*) 'sE1 = ', sE1, ' sE2 = ', sE2
!          write(6,*) ' '
!       endif
!         write(6,*) "i = ",icom, "j = ", jcom
         maxiter = 1

         call riemann_aug_JCP(maxiter,3,3,hL,hR,huL, &  
               huR,hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,pL,pR,sE1,sE2, &
               drytol,g,sw,fw,i,jcom,ixy)

!          call riemann_ssqfwave(maxiter,meqn,mwaves,hL,hR,huL,huR,
!      &     hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,pL,pR,sE1,sE2,drytol,g,
!      &     rho,sw,fw)

!          call riemann_fwave(meqn,mwaves,hL,hR,huL,huR,hvL,hvR,
!      &      bL,bR,uL,uR,vL,vR,phiL,phiR,pL,pR,sE1,sE2,drytol,g,rho,sw,
!      &      fw)

!         write (*,*) 'hL = ', hL, ' hR = ', hR
!        write (*,*) 'huL = ', huL, ' huR = ', huR
!        write (*,*) 'hvL = ', hvL, ' hvR = ', hvR
!        write (*,*) 'uL = ', uL, ' uR = ', uR
!        write (*,*) 'vL = ', vL, ' vR q= ', vR
!        write (*,*) 'phiL = ', phiL, ' phiR = ', phiR
!        write (*,*) 'bL = ', bL, ' bR = ', bR
!        stop
!
!        !eliminate ghost fluxes for wall
         do mw=1,3
            sw(mw)=sw(mw)*wall(mw)

               fw(1,mw)=fw(1,mw)*wall(mw) 
               fw(2,mw)=fw(2,mw)*wall(mw)
               fw(3,mw)=fw(3,mw)*wall(mw)
         enddo

         do mw=1,mwaves
            s(mw,i)=sw(mw)
            fwave(1,mw,i)=fw(1,mw)
            fwave(mu,mw,i)=fw(2,mw)
            fwave(nv,mw,i)=fw(3,mw)
!            write(51,515) sw(mw),fw(1,mw),fw(2,mw),fw(3,mw)
!515         format("++sw",4e25.16)
         enddo

!         write out speeds and fwaves
        !  if (debug) then
        !     write(*,*) 'i = ', i, ' j = ' , jcom
        !   write(*,*) 's = ', sw(1), sw(2), sw(3)
        !   write(*,*) 'fwave = ', fw(1,1), fw(2,1), fw(3,1)
        !   write(*,*) 'fwave = ', fw(1,2), fw(2,2), fw(3,2)
        !   write(*,*) 'fwave = ', fw(1,3), fw(2,3), fw(3,3)
        !   write(*,*) ' '
        !  endif
         
!         stop


 30      continue
      enddo


!==========Capacity for mapping from latitude longitude to physical space====
        if (mcapa.gt.0) then
         do i=2-mbc,mx+mbc
          if (ixy.eq.1) then
             dxdc=(earth_radius*deg2rad)
          else
             dxdc=earth_radius*cos(auxr(3,i))*deg2rad
          endif
         !  write(*,*) 'dxdc = ', dxdc, 'earth_radius = ', earth_radius, 'deg2rad = ', deg2rad
         !  write(*,*) ' '
         !  stop

          do mw=1,mwaves
!             if (s(mw,i) .gt. 316.d0) then
!               # shouldn't happen unless h > 10 km!
!                write(6,*) 'speed > 316: i,mw,s(mw,i): ',i,mw,s(mw,i)
!                endif
               s(mw,i)=dxdc*s(mw,i)
               fwave(1,mw,i)=dxdc*fwave(1,mw,i)
               fwave(2,mw,i)=dxdc*fwave(2,mw,i)
               fwave(3,mw,i)=dxdc*fwave(3,mw,i)
          enddo
         enddo
        endif

!===============================================================================


!============= compute fluctuations=============================================
         amdq(1:3,:) = 0.d0
         apdq(1:3,:) = 0.d0
         do i=2-mbc,mx+mbc
            do  mw=1,mwaves
               if (s(mw,i) < 0.d0) then
                     amdq(1:3,i) = amdq(1:3,i) + fwave(1:3,mw,i)
               else if (s(mw,i) > 0.d0) then
                  apdq(1:3,i)  = apdq(1:3,i) + fwave(1:3,mw,i)
               else
                 amdq(1:3,i) = amdq(1:3,i) + 0.5d0 * fwave(1:3,mw,i)
                 apdq(1:3,i) = apdq(1:3,i) + 0.5d0 * fwave(1:3,mw,i)
               endif
            enddo
!        write out fluctuations
!         write(*,*) 'amdq = ', amdq(1:3,i)
!         write(*,*) 'apdq = ', apdq(1:3,i)
         
         enddo
!--       do i=2-mbc,mx+mbc
!--            do m=1,meqn
!--                write(51,151) m,i,amdq(m,i),apdq(m,i)
!--                write(51,152) fwave(m,1,i),fwave(m,2,i),fwave(m,3,i)
!--151             format("++3 ampdq ",2i4,2e25.15)
!--152             format("++3 fwave ",8x,3e25.15)
!--            enddo
!--        enddo

      return
      end subroutine
      
! =====================================================
      subroutine fc2d_geoclaw_rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx,&
      ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq,tol,g,coordinate_system,&
      earth_radius,deg2rad)
! =====================================================
!
!     Riemann solver in the transverse direction using 
!     Jacobian matrix from left cell (if imp==1) or right cell (if imp==2).
!
!     Note this has been modified from the version used in v5.7.x and
!     earlier, where Roe averages based on ql and qr were used, which
!     is not correct.  In addition:
!      - a bug in the second component of the eigenvectors was fixed.
!      - when s(2) is close to zero this component of flux difference
!        is split equally between bmasdq and bpasdq to improve symmetry.
!   
!     Further modified to clean up and avoid a lot of work in dry cells.

!-----------------------last modified October 2020 ----------------------

    !   use geoclaw_module, only: g => grav, tol => dry_tolerance
    !   use geoclaw_module, only: coordinate_system,earth_radius,deg2rad

      implicit none
      
      double precision g, tol, earth_radius,deg2rad
      integer coordinate_system
      integer, intent(in) :: ixy,maxm,meqn,maux,mwaves,mbc,mx,imp

      double precision, intent(in) ::  ql(meqn,1-mbc:maxm+mbc)
      double precision, intent(in) ::  qr(meqn,1-mbc:maxm+mbc)
      double precision, intent(in) ::  asdq(meqn,1-mbc:maxm+mbc)
      double precision, intent(in) ::  aux1(maux,1-mbc:maxm+mbc)
      double precision, intent(in) ::  aux2(maux,1-mbc:maxm+mbc)
      double precision, intent(in) ::  aux3(maux,1-mbc:maxm+mbc)

      double precision, intent(out) ::  bmasdq(meqn,1-mbc:maxm+mbc)
      double precision, intent(out) ::  bpasdq(meqn,1-mbc:maxm+mbc)

      ! local:
      double precision ::  s(mwaves), r(meqn,mwaves), beta(mwaves)
      double precision ::  h,u,v
      double precision ::  delf1,delf2,delf3
      double precision ::  dxdcm,dxdcp,topo1,topo3,eta

      integer :: i,mw,mu,mv
      

      if (ixy == 1) then
         ! normal solve was in x-direction
         mu = 2
         mv = 3
      else
         ! normal solve was in y-direction
         mu = 3
         mv = 2
      endif

      ! initialize all components of result to 0:
      bmasdq(:,:) = 0.d0
      bpasdq(:,:) = 0.d0


      do i=2-mbc,mx+mbc

         if (imp==1) then
            h = ql(1,i-1)
         else
            h = qr(1,i)
         endif

         if (h <= tol) then
             ! fluctuation going into a dry cell, don't know how to split,
             ! so leave bmadsq(:,i)=bpasdq(:,i)=0 and go on to next i:
             cycle  
         endif

         ! compute velocities in relevant cell, and other quantities:

         if (imp==1) then
              ! fluctuation being split is left-going
              u = ql(mu,i-1)/h
              v = ql(mv,i-1)/h
              eta = h + aux2(1,i-1)
              topo1 = aux1(1,i-1)
              topo3 = aux3(1,i-1)
         else
              ! fluctuation being split is right-going
              u = qr(mu,i)/h
              v = qr(mv,i)/h
              eta = h + aux2(1,i)
              topo1 = aux1(1,i)
              topo3 = aux3(1,i)
         endif


         ! check if cell that transverse waves go into are both too high:
         ! Note: prior to v5.8.0 this checked against max rather than min
         if (eta < min(topo1,topo3)) cycle  ! go to next i

         ! if we get here, we want to do the splitting (no dry cells),
         ! so compute the necessary quantities:

         if (coordinate_system == 2) then
            ! on the sphere:
            if (ixy == 2) then
               dxdcp=(earth_radius*deg2rad)
               dxdcm = dxdcp
            else
               if (imp == 1) then
                  dxdcp = earth_radius*cos(aux3(3,i-1))*deg2rad
                  dxdcm = earth_radius*cos(aux1(3,i-1))*deg2rad
               else
                  dxdcp = earth_radius*cos(aux3(3,i))*deg2rad
                  dxdcm = earth_radius*cos(aux1(3,i))*deg2rad
               endif
            endif
         else
            ! coordinate_system == 1 means Cartesian:
            dxdcp = 1.d0
            dxdcm = 1.d0
         endif

        ! Determine some speeds necessary for the Jacobian
            
         ! In v5.7.x and prior versions,
         ! we used left right states to define Roe averages,
         ! which is consistent with those used in rpn2.
         ! But now we are computing upgoing, downgoing waves either in
         ! cell on left (if imp==1) or on right (if imp==2) so we
         ! should possibly use q values in cells above/below,
         ! but these aren't available (only aux values).
         ! At any rate, there is no clear justification for using cells
         ! on the other side of the normal-solve interface.

         ! v5.8.0: modified to use left or right state alone in defining
         ! Jacobian, based on imp:

         s(1) = v - dsqrt(g*h)
         s(2) = v
         s(3) = v + dsqrt(g*h)

        !  write(*,*) 's=',s(1:3)

! c        Determine asdq decomposition (beta)

         delf1 = asdq(1,i)
         delf2 = asdq(mu,i)
         delf3 = asdq(mv, i)

! c         write(*,*) 'delf=',delf1,delf2,delf3

         ! v5.8.0: fixed bug in beta(2): u in place of s(2)=v
         beta(1) = (s(3)*delf1 - delf3) / (s(3) - s(1))
         beta(2) = -u*delf1 + delf2
         beta(3) = (delf3 - s(1)*delf1) / (s(3) - s(1))

! c        Set-up eigenvectors
         r(1,1) = 1.d0
         r(2,1) = u    ! v5.8.0: fixed bug, u not s(2)=v
         r(3,1) = s(1)

         r(1,2) = 0.d0
         r(2,2) = 1.d0
         r(3,2) = 0.d0

         r(1,3) = 1.d0
         r(2,3) = u    ! v5.8.0: fixed bug, u not s(2)=v
         r(3,3) = s(3)

        write(*,*) 's=',s(1:3)
         write (*,*) ' beta=',beta(1:3)
            write (*,*) ' r=',r(1:3,1:3)
            write (*,*) ' '

         ! compute fluctuations

         do  mw=1,3
            if ((s(mw) < 0.d0) .and. (eta >= topo1)) then
                 bmasdq(1,i) =bmasdq(1,i) + dxdcm*s(mw)*beta(mw)*r(1,mw)
                 bmasdq(mu,i)=bmasdq(mu,i)+ dxdcm*s(mw)*beta(mw)*r(2,mw)
                 bmasdq(mv,i)=bmasdq(mv,i)+ dxdcm*s(mw)*beta(mw)*r(3,mw)
            elseif ((s(mw) > 0.d0) .and. (eta >= topo3)) then
                 bpasdq(1,i) =bpasdq(1,i) + dxdcp*s(mw)*beta(mw)*r(1,mw)
                 bpasdq(mu,i)=bpasdq(mu,i)+ dxdcp*s(mw)*beta(mw)*r(2,mw)
                 bpasdq(mv,i)=bpasdq(mv,i)+ dxdcp*s(mw)*beta(mw)*r(3,mw)
            endif
         enddo  ! loop on mw

! c         write (*,*) ' bmasdq=',bmasdq(1:3,i)
! c         write(*,*) ' bpasdq=',bpasdq(1:3,i)      
         
      enddo  ! loop on i

      return
      end
      
      
!----------------------------------------------------------------------
      subroutine riemann_aug_JCP(maxiter,meqn,mwaves,hL,hR,huL,huR,&
         hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,pL,pR,sE1,sE2,drytol,g, &
            sw,fw,i,j,ixy)

      ! solve shallow water equations given single left and right states
      ! This solver is described in J. Comput. Phys. (6): 3089-3113, March 2008
      ! Augmented Riemann Solvers for the Shallow Equations with Steady States and Inundation

      ! To use the original solver call with maxiter=1.

      ! This solver allows iteration when maxiter > 1. The iteration seems to help with
      ! instabilities that arise (with any solver) as flow becomes transcritical over variable topo
      ! due to loss of hyperbolicity.

      implicit none

      !input
      integer meqn,mwaves,maxiter
      double precision fw(meqn,mwaves)
      double precision sw(mwaves)
      double precision hL,hR,huL,huR,bL,bR,uL,uR,phiL,phiR,sE1,sE2
      double precision hvL,hvR,vL,vR,pL,pR
      double precision drytol,g,rho


      !local
      integer m,mw,k,iter
      double precision A(3,3)
      double precision r(3,3)
      double precision lambda(3)
      double precision del(3)
      double precision beta(3)

      double precision delh,delhu,delphi,delb,delnorm
      double precision rare1st,rare2st,sdelta,raremin,raremax
      double precision criticaltol,convergencetol,raretol
      double precision criticaltol_2, hustar_interface
      double precision s1s2bar,s1s2tilde,hbar,hLstar,hRstar,hustar
      double precision huRstar,huLstar,uRstar,uLstar,hstarHLL
      double precision deldelh,deldelphi,delP
      double precision s1m,s2m,hm
      double precision det1,det2,det3,determinant

      logical rare1,rare2,rarecorrector,rarecorrectortest,sonic

      logical debug
      integer ixy, i, j

      if (ixy == 1) then
         debug = .true.
      else
         debug  = .false.
      endif

      !determine del vectors
      delh = hR-hL
      delhu = huR-huL
      delphi = phiR-phiL
      delb = bR-bL
      delP = pR - pL
      delnorm = delh**2 + delphi**2

      call riemanntype(hL,hR,uL,uR,hm,s1m,s2m,rare1,rare2,&
      1,drytol,g)

!     if ( (i == 7) .and. (j == 15)) then
        !  if (debug) then
            write(6,*) 'i = ', i, ' j = ' , j
            write(6,*) 'hL = ', hL, ' hR = ', hR
            write(6,*) 'uL = ', uL, ' uR = ', uR
            write(6,*) 'hm = ', hm
            write(6,*) 's1m = ', s1m, ' s2m = ', s2m
            write(6,*) 'rare1 = ', rare1, ' rare2 = ', rare2
            ! write(6,*) ' '
        !  endif
!     endif

      lambda(1)= min(sE1,s2m) !Modified Einfeldt speed
      lambda(3)= max(sE2,s1m) !Modified Eindfeldt speed
      sE1=lambda(1)
      sE2=lambda(3)
      lambda(2) = 0.d0  ! ### Fix to avoid uninitialized value in loop on mw -- Correct?? ###

      
      hstarHLL = max((huL-huR+sE2*hR-sE1*hL)/(sE2-sE1),0.d0) ! middle state in an HLL solve

!    !determine the middle entropy corrector wave------------------------
      rarecorrectortest=.false.
      rarecorrector=.false.
      if (rarecorrectortest) then
         sdelta=lambda(3)-lambda(1)
         raremin = 0.5d0
         raremax = 0.9d0
         if (rare1.and.sE1*s1m.lt.0.d0) raremin=0.2d0
         if (rare2.and.sE2*s2m.lt.0.d0) raremin=0.2d0
         if (rare1.or.rare2) then
            !see which rarefaction is larger
            rare1st=3.d0*(sqrt(g*hL)-sqrt(g*hm))
            rare2st=3.d0*(sqrt(g*hR)-sqrt(g*hm))
            if (max(rare1st,rare2st).gt.raremin*sdelta.and. &
            max(rare1st,rare2st).lt.raremax*sdelta) then
                  rarecorrector=.true.
               if (rare1st.gt.rare2st) then
                  lambda(2)=s1m
               elseif (rare2st.gt.rare1st) then
                  lambda(2)=s2m
               else
                  lambda(2)=0.5d0*(s1m+s2m)
               endif
            endif
         endif
         if (hstarHLL.lt.min(hL,hR)/5.d0) rarecorrector=.false.
      endif

!    ## Is this correct 2-wave when rarecorrector == .true. ??
      do mw=1,mwaves
         r(1,mw)=1.d0
         r(2,mw)=lambda(mw)
         r(3,mw)=(lambda(mw))**2
      enddo
      if (.not.rarecorrector) then
         lambda(2) = 0.5d0*(lambda(1)+lambda(3))
!        lambda(2) = max(min(0.5d0*(s1m+s2m),sE2),sE1)
         r(1,2)=0.d0
         r(2,2)=0.d0
         r(3,2)=1.d0
      endif
!    !---------------------------------------------------


!    !determine the steady state wave -------------------
      !criticaltol = 1.d-6
      ! MODIFIED:
      criticaltol = max(drytol*g, 1d-6)
      criticaltol_2 = sqrt(criticaltol)
      deldelh = -delb
!     deldelphi = -0.5d0 * (hR + hL) * (g * delb + delp / rho)
      deldelphi = -0.5d0 * (hR + hL) * (g * delb)

!    !determine a few quanitites needed for steady state wave if iterated
      hLstar=hL
      hRstar=hR
      uLstar=uL
      uRstar=uR
      huLstar=uLstar*hLstar
      huRstar=uRstar*hRstar

      !iterate to better determine the steady state wave
      convergencetol=1.d-6
      do iter=1,maxiter
         !determine steady state wave (this will be subtracted from the delta vectors)
         if (min(hLstar,hRstar).lt.drytol.and.rarecorrector) then
            rarecorrector=.false.
            hLstar=hL
            hRstar=hR
            uLstar=uL
            uRstar=uR
            huLstar=uLstar*hLstar
            huRstar=uRstar*hRstar
            lambda(2) = 0.5d0*(lambda(1)+lambda(3))
!          lambda(2) = max(min(0.5d0*(s1m+s2m),sE2),sE1)
            r(1,2)=0.d0
            r(2,2)=0.d0
            r(3,2)=1.d0
         endif

         hbar =  max(0.5d0*(hLstar+hRstar),0.d0)
         s1s2bar = 0.25d0*(uLstar+uRstar)**2 - g*hbar
         s1s2tilde= max(0.d0,uLstar*uRstar) - g*hbar

!       !find if sonicproblem
         ! MODIFIED from 5.3.1 version
         sonic= .false.
         if (abs(s1s2bar) <= criticaltol) then
            sonic= .true.
         else if (s1s2bar*s1s2tilde <= criticaltol**2) then
            sonic= .true.
         else if (s1s2bar*sE1*sE2 <= criticaltol**2) then
            sonic= .true.
         else if (min(abs(sE1),abs(sE2)) < criticaltol_2) then
            sonic= .true.
         else if (sE1 <  criticaltol_2 .and. s1m > -criticaltol_2) then
            sonic= .true.
         else if (sE2 > -criticaltol_2 .and. s2m <  criticaltol_2) then
            sonic= .true.
         else if ((uL+dsqrt(g*hL))*(uR+dsqrt(g*hR)) < 0.d0) then
            sonic= .true.
         else if ((uL- dsqrt(g*hL))*(uR- dsqrt(g*hR)) < 0.d0) then
            sonic= .true.
         end if

!       !find jump in h, deldelh
         if (sonic) then
            deldelh =  -delb
         else
            deldelh = delb*g*hbar/s1s2bar
         endif
!       !find bounds in case of critical state resonance, or negative states
         if (sE1.lt.-criticaltol.and.sE2.gt.criticaltol) then
            deldelh = min(deldelh,hstarHLL*(sE2-sE1)/sE2)
            deldelh = max(deldelh,hstarHLL*(sE2-sE1)/sE1)
         elseif (sE1.ge.criticaltol) then
            deldelh = min(deldelh,hstarHLL*(sE2-sE1)/sE1)
            deldelh = max(deldelh,-hL)
         elseif (sE2.le.-criticaltol) then
            deldelh = min(deldelh,hR)
            deldelh = max(deldelh,hstarHLL*(sE2-sE1)/sE2)
         endif

!       ! adjust deldelh for well-balancing of atmospheri!pressure difference 
!        deldelh = deldelh - delP/(rho*g)

!       !find jump in phi, deldelphi
         if (sonic) then
            deldelphi = -g*hbar*delb
         else
            deldelphi = -delb*g*hbar*s1s2tilde/s1s2bar
         endif
!       !find bounds in case of critical state resonance, or negative states
         deldelphi=min(deldelphi,g*max(-hLstar*delb,-hRstar*delb))
         deldelphi=max(deldelphi,g*min(-hLstar*delb,-hRstar*delb))
!        deldelphi = deldelphi - hbar * delp / rho

         del(1)=delh-deldelh
         del(2)=delhu
         del(3)=delphi-deldelphi

!       !Determine determinant of eigenvector matrix========
         det1=r(1,1)*(r(2,2)*r(3,3)-r(2,3)*r(3,2))
         det2=r(1,2)*(r(2,1)*r(3,3)-r(2,3)*r(3,1))
         det3=r(1,3)*(r(2,1)*r(3,2)-r(2,2)*r(3,1))
         determinant=det1-det2+det3

!       !solve for beta(k) using Cramers Rule=================
         do k=1,3
            do mw=1,3
                  A(1,mw)=r(1,mw)
                  A(2,mw)=r(2,mw)
                  A(3,mw)=r(3,mw)
            enddo
            A(1,k)=del(1)
            A(2,k)=del(2)
            A(3,k)=del(3)
            det1=A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))
            det2=A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
            det3=A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
            beta(k)=(det1-det2+det3)/determinant
         enddo

         !exit if things aren't changing
         if (abs(del(1)**2+del(3)**2-delnorm).lt.convergencetol) exit
         delnorm = del(1)**2+del(3)**2
         !find new states qLstar and qRstar on either side of interface
         hLstar=hL
         hRstar=hR
         uLstar=uL
         uRstar=uR
         huLstar=uLstar*hLstar
         huRstar=uRstar*hRstar
         do mw=1,mwaves
            if (lambda(mw).lt.0.d0) then
               hLstar= hLstar + beta(mw)*r(1,mw)
               huLstar= huLstar + beta(mw)*r(2,mw)
            endif
         enddo
         do mw=mwaves,1,-1
            if (lambda(mw).gt.0.d0) then
               hRstar= hRstar - beta(mw)*r(1,mw)
               huRstar= huRstar - beta(mw)*r(2,mw)
            endif
         enddo

         if (hLstar.gt.drytol) then
            uLstar=huLstar/hLstar
         else
            hLstar=max(hLstar,0.d0)
            uLstar=0.d0
         endif
         if (hRstar.gt.drytol) then
            uRstar=huRstar/hRstar
         else
            hRstar=max(hRstar,0.d0)
            uRstar=0.d0
         endif

      enddo ! end iteration on Riemann problem

      do mw=1,mwaves
         sw(mw)=lambda(mw)
         fw(1,mw)=beta(mw)*r(2,mw)
         fw(2,mw)=beta(mw)*r(3,mw)
         fw(3,mw)=beta(mw)*r(2,mw)
      enddo
      !find transverse components (ie huv jumps).
      ! MODIFIED from 5.3.1 version
      fw(3,1)=fw(3,1)*vL
      fw(3,3)=fw(3,3)*vR
      fw(3,2)= 0.d0
 
      hustar_interface = huL + fw(1,1)   ! = huR - fw(1,3)
      if (hustar_interface <= 0.0d0) then
          fw(3,1) = fw(3,1) + (hR*uR*vR - hL*uL*vL - fw(3,1)- fw(3,3))
        else
          fw(3,3) = fw(3,3) + (hR*uR*vR - hL*uL*vL - fw(3,1)- fw(3,3))
        end if


      return

      end !subroutine riemann_aug_JCP-------------------------------------------------

!============================================================================
      subroutine riemanntype(hL,hR,uL,uR,hm,s1m,s2m,rare1,rare2, &
                   maxiter,drytol,g)

      !determine the Riemann structure (wave-type in each family)


      implicit none

      !input
      double precision hL,hR,uL,uR,drytol,g
      integer maxiter

      !output
      double precision s1m,s2m
      logical rare1,rare2

      !local
      double precision hm,u1m,u2m,um,delu
      double precision h_max,h_min,h0,F_max,F_min,dfdh,F0,slope,gL,gR
      integer iter



!    !Test for Riemann structure

      h_min=min(hR,hL)
      h_max=max(hR,hL)
      delu=uR-uL

      if (h_min.le.drytol) then
         hm=0.d0
         um=0.d0
         s1m=uR+uL-2.d0*sqrt(g*hR)+2.d0*sqrt(g*hL)
         s2m=uR+uL-2.d0*sqrt(g*hR)+2.d0*sqrt(g*hL)
         if (hL.le.0.d0) then
            rare2=.true.
            rare1=.false.
         else
            rare1=.true.
            rare2=.false.
         endif

      else
         F_min= delu+2.d0*(sqrt(g*h_min)-sqrt(g*h_max))
         F_max= delu + &    
                 (h_max-h_min)*(sqrt(.5d0*g*(h_max+h_min)/(h_max*h_min)))

         if (F_min.gt.0.d0) then !2-rarefactions

            hm=(1.d0/(16.d0*g))* &      
                        max(0.d0,-delu+2.d0*(sqrt(g*hL)+sqrt(g*hR)))**2
            um=sign(1.d0,hm)*(uL+2.d0*(sqrt(g*hL)-sqrt(g*hm)))

            s1m=uL+2.d0*sqrt(g*hL)-3.d0*sqrt(g*hm)
            s2m=uR-2.d0*sqrt(g*hR)+3.d0*sqrt(g*hm)

            rare1=.true.
            rare2=.true.

         elseif (F_max.le.0.d0) then !2 shocks

!          !root finding using a Newton iteration on sqrt(h)===
            h0=h_max
            do iter=1,maxiter
               gL=sqrt(.5d0*g*(1/h0 + 1/hL))
               gR=sqrt(.5d0*g*(1/h0 + 1/hR))
               F0=delu+(h0-hL)*gL + (h0-hR)*gR
               dfdh=gL-g*(h0-hL)/(4.d0*(h0**2)*gL)+ &
                        gR-g*(h0-hR)/(4.d0*(h0**2)*gR)
               slope=2.d0*sqrt(h0)*dfdh
               h0=(sqrt(h0)-F0/slope)**2
            enddo
               hm=h0
               u1m=uL-(hm-hL)*sqrt((.5d0*g)*(1/hm + 1/hL))
               u2m=uR+(hm-hR)*sqrt((.5d0*g)*(1/hm + 1/hR))
               um=.5d0*(u1m+u2m)

               s1m=u1m-sqrt(g*hm)
               s2m=u2m+sqrt(g*hm)
               rare1=.false.
               rare2=.false.

         else !one shock one rarefaction
            h0=h_min

            do iter=1,maxiter
               F0=delu + 2.d0*(sqrt(g*h0)-sqrt(g*h_max)) &
               + (h0-h_min)*sqrt(.5d0*g*(1/h0+1/h_min))
               slope=(F_max-F0)/(h_max-h_min)
               h0=h0-F0/slope
            enddo

            hm=h0
            if (hL.gt.hR) then
               um=uL+2.d0*sqrt(g*hL)-2.d0*sqrt(g*hm)
               s1m=uL+2.d0*sqrt(g*hL)-3.d0*sqrt(g*hm)
               s2m=uL+2.d0*sqrt(g*hL)-sqrt(g*hm)
               rare1=.true.
               rare2=.false.
            else
               s2m=uR-2.d0*sqrt(g*hR)+3.d0*sqrt(g*hm)
               s1m=uR-2.d0*sqrt(g*hR)+sqrt(g*hm)
               um=uR-2.d0*sqrt(g*hR)+2.d0*sqrt(g*hm)
               rare2=.true.
               rare1=.false.
            endif
         endif
      endif

      return

      end ! subroutine riemanntype----------------------------------------------------------------
