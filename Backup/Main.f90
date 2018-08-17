    !  Console1.f90
    !
    !  FUNCTIONS:
    !  Console1 - Entry point of console application.
    !

    !****************************************************************************
    !
    !  PROGRAM: Console1
    !
    !  PURPOSE:  Entry point for the console application.
    !
    !****************************************************************************
    !#include "fintrf.h"
    program Console1

    use Header
    use routines

    implicit none

#ifdef mpi
    include 'mpif.h'
#endif        

!#include "../globalMacros.txt"
#include "globalMacros.txt"

    ! Variables
    type (structparamstype) :: params
    type (gridsType) :: grids
    type (modelObjectsType), target :: modelObjects

    real (kind=rk) :: ypath(Tperiods, numSims) !income
    real (kind=rk) :: cpath(Tperiods, numSims)  !consumption
    integer :: lpath(Tperiods, numSims) !labour supply
    real (kind=rk) :: vpath(Tperiods, numSims)  !value
    integer :: apath(Tperiods + 1,numSims) !this is the path at the start of each period, so we include the 'start' of death
    real (kind=rk) :: AIME(Tperiods + 1,numSims)

    real (kind=rk) :: start, finish, moments(2,24),weights(2,24), y(dimEstimation+1), p(dimEstimation+1,dimEstimation)
    integer :: action, ios, requiredl, typeSim
    INTEGER(kind=4) :: iter

#ifdef mpi
    integer :: provided
#endif     

    CHARACTER(len=255) :: cwd
    CALL getcwd(cwd)
    WRITE(*,*) TRIM(cwd)

    !call cpu_time(start)

#ifdef mpi    
    call MPI_Init_thread(MPI_THREAD_MULTIPLE,provided,ierror)!mpi_init
    if (ierror.ne.0) stop 'mpi problem171'
    call mpi_comm_rank(mpi_comm_world, rank, ierror)
    call mpi_comm_size(mpi_comm_world, procSize, ierror)
    if (ierror.ne.0) stop 'mpi problem172'
    if (rank.eq.0) write (*, *) 'Using MPI in solution. Using', procSize, 'cores'
    if (rank.eq.0) write(*,*)
    call mpi_barrier(mpi_comm_world, ierror)
    if (ierror.ne.0) stop 'mpi problem173'
#else
    rank = 0
    procSize = 1
#endif

    params%r= 0.02

    params%mu = 0
    params%sigma = 0.0922
    params%rho = 0.96
    params%hrsWrk = 0.3159
    params%p = 0.06
    !params%delta(3) = 9.083298
    !params%delta(2)=0.029333
    !params%delta(1)=-0.00033

    !Earning type 1 - uneducated no db
    params%delta(1,3) = 9.083298
    params%delta(1,2)=0.029333
    params%delta(1,1)=-0.00023

    params%spouseInc(1) = 5121!6235.8998884204320

    !Earning type 2 - uneducated db
#if TYPES_SIZE > 1
    !if (numPointsType .GE. 2) then
    params%delta(2,3) = 8.5256042
    params%delta(2,2)= 0.061582208
    params%delta(2,1)= -0.0005901
    !end if
    params%spouseInc(2) = 5282
#endif

#if TYPES_SIZE > 2
    !Earning type 3 - educated no db
    !if (numPointsType .GE. 3) then
    params%delta(3,3) = 7.906901
    params%delta(3,2)= 0.094386
    params%delta(3,1)= -0.00091
    !end if
    params%spouseInc(3) = 6840
#endif

#if TYPES_SIZE > 3    
    !Earning type 4 - educated db
    !if (numPointsType .GE. 4) then
    params%delta(4,3) = 7.392151
    params%delta(4,2)= 0.121426
    params%delta(4,1)= -0.00117
    !end if
    params%spouseInc(4) = 7684
#endif


    params%tol = 1e-10
    params%minCons = 1e-5

    params%nu = 0.290307522834662 !0.8 !0.5   !0.290307522834662 ! 0.287177078339264!
    params%beta = 0.985451393874388 ! 0.999 !0.999 !0.985451393874388 !0.985724939013559!
    params%gamma = 2.34092202506161!1.0 !2.34092202506161! 2.31567939133319!
    params%db(1) = 0.596815295617949 !0.796815!0.596815!0.596815295617949 ! 0.590379716068728!
    params%db(2) =  -4.270610627570502E-006!4.270610627570502E-005!-4.270610627570502E-007!-4.270610627570502E-006 !-4.224559772943942E-006 !
    params%thetab =  2.894599354187538E-002 !2.926152647601828E-001!2.926152647601828E-002 !2.894599354187538E-002 !
    params%k=650000

    params%lambda= 1.0!100! 0.001! 0.001!100.0!1.0 !1000000 !0.001!0.0000001!

    call setupMisc(params,grids)
    if (fullLifeCycle) then
        params%startA =0
    else
        params%startA = maxval(grids%initialAssets)
    end if

    action =1
    if (action .EQ. 0) then
        do typeSim = 1, numPointsType
            call getassetgrid( params, grids%maxInc(typeSim,:), grids%Agrid(typeSim,:,:))
        end do
 
        open (unit = 1001,file='..\\..\\moments\\moments.txt', action='read', IOSTAT = ios)
        open (unit = 1002,file='..\\..\\moments\\assetmom.txt', action='read', IOSTAT = ios)
        open (unit = 1003,file='..\\..\\moments\\Weight.txt', action='read', IOSTAT = ios)
        open (unit = 1004,file='..\\..\\moments\\weightAsset.txt', action='read', IOSTAT = ios) 
        read (1001, *) moments(1,:)
        read (1002, *) moments(2,:)
        read (1003, *) weights(1,:)
        read (1004, *) weights(2,:)
        close (unit=1001)
        close (unit=1002)
        close (unit=1003)
        close (unit=1004)      
        p(1,1)=   0.287177126203365 !0.287177126203365 ! 0.290307522834662 ! 0.287177078339264!
        p(1,2) =  0.985724876594556     ! 0.985724876594556 ! 0.985451393874388 !0.985724939013559!
        p(1,3) =   2.31567977728987  !2.31567977728987 ! 2.34092202506161! 2.31567939133319!
        p(1,4) =   0.590379814467926  !0.590379814467926 !0.596815295617949 ! 0.590379716068728!
        p(1,5) = -4.224560477055699E-006 !-4.224560477055699E-006! -4.270610627570502E-006 !-4.224559772943942E-006 !
        p(1,6) =  2.894599836633410E-002 !2.894599836633410E-002! 2.926152647601828E-002 !2.894599354187538E-002 !        
        write (*,*), gmm_criteria(p(1,:))
       
    else if (action .EQ. 1) then
        params%nu =   0.287177126203365 ! 0.290307522834662 ! 0.287177078339264!
        params%beta = 0.985724876594556 ! 0.985451393874388 !0.985724939013559!
        params%gamma =  2.31567977728987 ! 2.34092202506161! 2.31567939133319!
        params%db(1) =  0.590379814467926 !0.596815295617949 ! 0.590379716068728!
        params%db(2) = -4.224560477055699E-006! -4.270610627570502E-006 !-4.224559772943942E-006 !
        params%thetab =  2.894599836633410E-002! 2.926152647601828E-002 !2.894599354187538E-002 !
        params%k=650000   
    
        do typeSim = 1, numPointsType
            call getassetgrid( params, grids%maxInc(typeSim,:), grids%Agrid(typeSim,:,:))
        end do
        
        call solveValueFunction( params, grids, modelObjects, .TRUE. )
        !simulate
        if (rank == 0) then
            !call cpu_time(starsim)
            call simWithUncer(params, grids, modelObjects, grids%Simy, cpath, apath, vpath, lpath, ypath, AIME )
            !call cpu_time(finish)
            !print '("Time = ",f11.3," seconds.")',finish-starsim
            call writetofile(grids,ypath, cpath, apath, vpath, lpath, grids%Simy,AIME)
            !call writetofileByType(grids,ypath, cpath, apath, vpath, lpath, grids%Simy,AIME)
        end if
    else
        if (Tretire+startAge .NE. 60) then
            write (*,*) "Only estimate with SPA = 60"
        end if
#ifdef win
        open (unit = 1001,file='..\\..\\moments\\moments.txt', action='read', IOSTAT = ios)
        open (unit = 1002,file='..\\..\\moments\\assetmom.txt', action='read', IOSTAT = ios)
        open (unit = 1003,file='..\\..\\moments\\Weight.txt', action='read', IOSTAT = ios)
        open (unit = 1004,file='..\\..\\moments\\weightAsset.txt', action='read', IOSTAT = ios)
#else
        open (unit = 1001,file='../moments/moments.txt', action='read', IOSTAT = ios)
        open (unit = 1002,file='../moments/assetmom.txt', action='read', IOSTAT = ios)
        open (unit = 1003,file='../moments/Weight.txt', action='read', IOSTAT = ios)
        open (unit = 1004,file='../moments/weightAsset.txt', action='read', IOSTAT = ios)
#endif
        read (1001, *) moments(1,:)
        read (1002, *) moments(2,:)
        read (1003, *) weights(1,:)
        read (1004, *) weights(2,:)
        close (unit=1001)
        close (unit=1002)
        close (unit=1003)
        close (unit=1004)

        if (rank==0) then
            print '("Setting up initial guess for hilling climbing algorithm")'
        end if

        call initialGuess(rank,params,grids,moments,weights,p,y)

        if (rank==0) open (unit=211, file='..\\out\guess2.txt', status='unknown', action='write')
        call amoeba(p,y, 0.001_rk,gmm_criteria,iter,.TRUE.,0.0_rk) !0.0001_rk !0.001_rk!0.07_rk!0.0001054 !0.000000001_rk !
        if (rank==0) close (unit=211)

        if (rank==0) then
            print '("P = ",f6.3)',P(1,:)
            print '("Y = ",f16.3)',Y
#ifdef win
            inquire (iolength=requiredl)  P(1,:)
            open (unit=201, file='..\\out\params.txt', status='unknown',recl=requiredl, action='write')
            write (201, * ) P(1,:)
#else
            inquire (iolength=requiredl)  P
            open (unit=201, file='./out/params.txt', status='unknown',recl=requiredl, action='write')
            write (201, * ) P(1,1)
            write (201, * ) P(1,2)
            write (201, * ) P(1,3)
            write (201, * ) P(1,4)
            write (201, * ) P(1,5)
#endif
            close (unit=201)

            print '("Generating files")'
        end if

        call solveValueFunction( params, grids, modelObjects, .FALSE. )

        if (rank==0) then
            !call simWithUncer(params, grids, modelObjects, grids%Simy, cpath, apath, vpath, lpath, ypath ,AIME)
            !call writetofile(grids,ypath, cpath, apath, vpath, lpath, grids%Simy,AIME)
            !call writetofileByType(grids,ypath, cpath, apath, vpath, lpath, grids%Simy,AIME)
        end if
#ifdef mpi 
        call mpi_barrier(mpi_comm_world, ierror)
        if (ierror.ne.0) stop 'mpi problem180'
#endif        
    end if

    if (rank==0) then
        call cpu_time(finish)
        print '("Time = ",f11.3," seconds.")',finish-start
    end if
#ifdef mpi
    call mpi_barrier(mpi_comm_world, ierror)
    if (ierror.ne.0) stop 'mpi problem180'
    if (rank.eq.0) then
        call mpi_finalize(ierror)
    end if
    if (ierror.ne.0) stop 'mpi problem190'
#endif

    contains
    function gmm_criteria(control)
    implicit none
    !inputs
    real (kind=rk), intent(in) ::control(:)
    !output
    real (kind=rk) :: gmm_criteria
    if (maxval(control(1:2)) > 1)  then
        gmm_criteria = huge(gmm_criteria)
        return
    end if
    if (minval(control(1:4)) < 0 .or.  control(5)> 0 .OR. control(6)< 0 )  then !.OR. control(6)< 0
        gmm_criteria = huge(gmm_criteria)
        return
    end if
    params%nu = control(1)
    params%beta = control(2)
    params%gamma = control(3)
    params%db(1)= control(4)
    params%db(2)= control(5)
    params%thetab = control(6)
    if (rank==0) then
        write (211,*) "Calcualting GMM for",  control(1),  control(2),  control(3),  control(4),  control(5),  control(6)
    end if
    gmm_criteria = gmm(params,grids,moments,weights) !*-1.0

    end function

    end program Console1

