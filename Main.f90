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
    !type (modelObjectsType), target :: modelObjects

    real (kind=rk) :: ypath(Tperiods, numSims) !income
    real (kind=rk) :: cpath(Tperiods, numSims)  !consumption
    integer :: lpath(Tperiods, numSims) !labour supply
    real (kind=rk) :: vpath(Tperiods, numSims)  !value
    integer :: apath(Tperiods + 1,numSims) !this is the path at the start of each period, so we include the 'start' of death
    integer :: combinedData(3,24,numSims)
    real (kind=rk) :: AIME(Tperiods + 1,numSims)


    real (kind=rk) :: start, finish, moments(2,24),weights(2,24)
    real (kind=rk), allocatable :: y(:), p(:,:)
    integer :: action, ios, requiredl, typeSim
    integer ::SPA
    logical :: delCache
    INTEGER(kind=4) :: iter
    real (kind=rk) :: gmmScorce
    real (kind=rk) :: indxrk
    real (kind=rk) :: optLambda

#ifdef mpi
    integer :: provided
#endif     

    !CHARACTER(len=255) :: cwd
    !call chdir( '..' )
    !CALL getcwd(cwd)
    !WRITE(*,*) TRIM(cwd)

    call cpu_time(start)

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

    call setModel
    allocate(y(dimEstimation+1), p(dimEstimation+1,dimEstimation))

    params%r= 0.02

    params%percentCons = 0.0 !0.1
    params%mu = 0
    params%sigma = 0.0922
    params%rho = 0.96
    params%hrsWrk = 0.3159
    params%p = 0.06

    !Earning type 1 - uneducated no db
    params%delta(1,3) = 9.083298
    params%delta(1,2)=0.029333
    params%delta(1,1)=-0.00023

    params%spouseInc(1) = 5121!6235.8998884204320

    !Earning type 2 - uneducated db
#if TYPES_SIZE > 1
    params%delta(2,3) = 8.5256042
    params%delta(2,2)= 0.061582208
    params%delta(2,1)= -0.0005901
    params%spouseInc(2) = 5282
#endif

#if TYPES_SIZE > 2
    params%delta(3,3) = 7.906901
    params%delta(3,2)= 0.094386
    params%delta(3,1)= -0.00091
    params%spouseInc(3) = 6840
#endif

#if TYPES_SIZE > 3    
    params%delta(4,3) = 7.392151
    params%delta(4,2)= 0.121426
    params%delta(4,1)= -0.00117
    params%spouseInc(4) = 7684
#endif

    params%tol = 1e-10
    params%minCons = 1e-5

    params%nu =   0.287177126203365 ! 0.290307522834662 ! 0.287177078339264!
    params%beta = 0.985724876594556 ! 0.985451393874388 !0.985724939013559!
    params%gamma =  2.31567977728987 ! 2.34092202506161! 2.31567939133319!
    params%db(1) =  0.590379814467926 !0.596815295617949 ! 0.590379716068728!
    params%db(2) = -4.224560477055699E-006! -4.270610627570502E-006 !-4.224559772943942E-006 !
    params%thetab =  2.894599836633410E-002! 2.926152647601828E-002 !2.894599354187538E-002 !
    params%k=650000

    
    !params%nu = 0.233985260541000       
    !params%beta = 0.999690531106298      
    !params%gamma = 1.55555894964362     
    !params%db(1) = 0.878735699727095      
    !params%db(2) = -2.837850432880502E-006 
    !params%thetab = 1.944448764319950E-002

         
        params%nu = 0.236302115158612!  0.236281327423160  !0.287177126203365 !0.866530785705649 ! 0.290307522834662 ! 0.287177078339264!
        params%beta =  0.999961762412824   !0.999986702072559 !0.985724876594556 !0.990655741245757 ! 0.985451393874388 !0.985724939013559!
        params%gamma =  1.57096164606611   !1.57082344698532 ! 2.31567977728987 !1.78141736803550  ! 2.34092202506161! 2.31567939133319!
        params%db(1) = 0.887436687382750     !0.887358618682115! 0.596815295617949 !0.590379814467926 ! 0.590379716068728!
        params%db(2) = -2.865950010026133E-006 !-2.865697889512581E-006!-4.224560477055699E-006!-1.274722676211178E-005!  -4.270610627570502E-006 !-4.224559772943942E-006 !
        params%thetab = 1.963702135613103E-002 !1.963529386755248E-002 !2.894599836633410E-002!8.734191576465597E-002!  2.926152647601828E-002 !2.894599354187538E-002 !
        params%lambda =  1.000087978748379E-003 ! 7.804665289405568E-002
        params%k=650000      
        
 
       
    
     
  
 
  
       
 
  
    
    
    !params%lambda= 0.001!0.0000000001! 1.0!0.01!!!

    call setupMisc(params,grids)
    if (fullLifeCycle) then
        params%startA =0
    else
        params%startA = maxval(grids%initialAssets)
    end if

#ifdef debugMPI
    if (rank==0) pause
#endif    
    action =0
    if (action .EQ. 0) then
        do typeSim = 1, numPointsType
            call getassetgrid( params, grids%maxInc(typeSim,:), grids%Agrid(typeSim,:,:))
        end do

        open (unit = 1001,file=trim(pathMoments) // 'moments.txt', action='read', IOSTAT = ios)
        open (unit = 1002,file=trim(pathMoments) // 'assetmom.txt', action='read', IOSTAT = ios)
        open (unit = 1003,file=trim(pathMoments) // 'Weight.txt', action='read', IOSTAT = ios)
        open (unit = 1004,file=trim(pathMoments) // 'weightAsset.txt', action='read', IOSTAT = ios)
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

        !p(1,1) =  2.241321348790637E-002 !0.287177126203365 !0.866530785705649 ! 0.290307522834662 ! 0.287177078339264!
        !p(1,2)=  0.975935121023040 ! 0.985724876594556 !0.990655741245757 ! 0.985451393874388 !0.985724939013559!
        !p(1,3) =  0.180731055791926 ! 2.31567977728987 !1.78141736803550  ! 2.34092202506161! 2.31567939133319!
        !p(1,4) =  4.607716845543490E-002 ! 0.596815295617949 !0.590379814467926 ! 0.590379716068728!
        !p(1,5) =  -3.297128051827123E-007 !-4.224560477055699E-006!-1.274722676211178E-005!  -4.270610627570502E-006 !-4.224559772943942E-006 !
        !p(1,6) =  2.259138287169178E-003 !2.894599836633410E-002!8.734191576465597E-002!  2.926152647601828E-002 !2.894599354187538E-002 !
        !p(1,7) =  0.001 ! 7.804665289405568E-002

        p(1,1) = 0.236200045532760 !0.240191353041616
        p(1,2) = 0.999995559030882  !0.999477774429760
        p(1,3) = 1.57028307631512 ! 1.93680557456545
        p(1,4) = 0.887053363134520  !0.493786285559129
        p(1,5) = -2.864712075930465E-006 ! -3.533369493611974E-006
        p(1,6) = 1.962853923390662E-002 !2.421007064408893E-002
        optLambda =  9.996559952862711E-004 !0.836387480497043 
                   
         open (unit = 666, form="unformatted", file=trim(path) // 'guessP.txt', status='unknown', ACCESS="STREAM", action='read', IOSTAT = ios)
         read (666) p
         close (unit=666)
         
         open (unit = 667, form="unformatted", file=trim(path) // 'guessY.txt', status='unknown', ACCESS="STREAM", action='read', IOSTAT = ios)
         read (667) Y
         close (unit=667)         
 
        optLambda = p(1,7)
        !Lambda = 0.436387480497043 is best!!!!
        do indxrk = -0.3, 0.3, 0.1
            p(1,7) = optLambda + indxrk*optLambda
            if (rank==0) write (*,*), 'For lambda = ', p(1,7), indxrk
            gmmScorce =  gmm_criteria(p(1,:))
            if (rank==0) write (*,*), 'GMM Criteria is ', gmmScorce
        end do
    else if (action .EQ. 1) then
        params%nu =   0.236281327423160  !0.287177126203365 !0.866530785705649 ! 0.290307522834662 ! 0.287177078339264!
        params%beta = 0.999986702072559 !0.985724876594556 !0.990655741245757 ! 0.985451393874388 !0.985724939013559!
        params%gamma = 1.57082344698532 ! 2.31567977728987 !1.78141736803550  ! 2.34092202506161! 2.31567939133319!
        params%db(1) = 0.887358618682115! 0.596815295617949 !0.590379814467926 ! 0.590379716068728!
        params%db(2) = -2.865697889512581E-006!-4.224560477055699E-006!-1.274722676211178E-005!  -4.270610627570502E-006 !-4.224559772943942E-006 !
        params%thetab = 1.963529386755248E-002 !2.894599836633410E-002!8.734191576465597E-002!  2.926152647601828E-002 !2.894599354187538E-002 !
        !params%lambda =  7.804665289405568E-002
        params%k=650000
            
       

        do typeSim = 1, numPointsType
            call getassetgrid( params, grids%maxInc(typeSim,:), grids%Agrid(typeSim,:,:))
        end do

        call solveValueFunction( params, grids, .TRUE., .TRUE. )
        !simulate
#ifdef mpi
        call mpi_barrier(mpi_comm_world, ierror)
        if (ierror.ne.0) stop 'mpi problem180'
#endif         
        if (rank == 0) then
            delCache = .FALSE.
            if (modelChoice>1) then
                if (counterFact) then
                    call simWithUncer(params, grids, grids%Simy, cpath, apath, vpath, lpath, ypath, AIME, 1, .TRUE. ) !5
                    call writetofile(grids,ypath, cpath, apath, vpath, lpath, grids%Simy,AIME)
                    call writetofileByType(grids,ypath, cpath, apath, vpath, lpath, grids%Simy,AIME)
                else
                    delCache = .FALSE.
                    do SPA =1,3
                        write (*,*) "Sim True SPA", 59+SPA
                        if (SPA==3) delCache = .TRUE.
                        call simWithUncer(params, grids, grids%Simy, cpath, apath, vpath, lpath, ypath, AIME, SPA, .FALSE. )
                        combinedData(SPA,:,:) = apath(52-startAge+1:75-startAge+1,:)*(2*lpath(52-startAge+1:75-startAge+1,:)-1)
                        call writetofileAge(grids,ypath, cpath, apath, vpath, lpath, grids%Simy,AIME, 59+SPA)
                    end do
                    inquire (iolength=requiredl) combinedData
                    !print *, path
                    !print *, trim(path) // 'CombinedData'
                    open (unit=201, form="unformatted", file=trim(path) // 'CombinedData', ACCESS="STREAM", action='write', IOSTAT = ios)
                    write (201)  combinedData
                    close( unit=201)
                end if
            else
                call simWithUncer(params, grids, grids%Simy, cpath, apath, vpath, lpath, ypath, AIME, TrueSPA, .TRUE. )
                call writetofile(grids,ypath, cpath, apath, vpath, lpath, grids%Simy,AIME)
                call writetofileByType(grids,ypath, cpath, apath, vpath, lpath, grids%Simy,AIME)
            end if
        end if
    else
        if (TrueSpa .NE. 1) then
            write (*,*) "Only estimate with SPA = 60"
        end if
        open (unit = 1001,file=trim(pathMoments) // 'moments.txt', action='read', IOSTAT = ios)
        open (unit = 1002,file=trim(pathMoments) // 'assetmom.txt', action='read', IOSTAT = ios)
        open (unit = 1003,file=trim(pathMoments) // 'Weight.txt', action='read', IOSTAT = ios)
        open (unit = 1004,file=trim(pathMoments) // 'weightAsset.txt', action='read', IOSTAT = ios)

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

        !call initialGuess(rank,params,grids,moments,weights,p,y)
         open (unit = 666, form="unformatted", file=trim(path) // 'guessP.txt', status='unknown', ACCESS="STREAM", action='read', IOSTAT = ios)
         read (666) p
         close (unit=666)
         
         open (unit = 667, form="unformatted", file=trim(path) // 'guessY.txt', status='unknown', ACCESS="STREAM", action='read', IOSTAT = ios)
         read (667) Y
         close (unit=667)
#ifdef mpi
        call mpi_barrier(mpi_comm_world, ierror)
        if (ierror.ne.0) stop 'mpi problem180'
#endif          

        !open (unit=201, form="unformatted", file=trim(path) // 'CombinedData', ACCESS="STREAM", action='write', IOSTAT = ios)
        !open (unit=211, file='..\\out\guess2.txt', status='unknown', action='write')
        call amoeba(p,y, 0.00001_rk,gmm_criteria,iter,.TRUE.,0.0_rk) !0.0001_rk!0.001_rk!0.0001_rk !0.001_rk!0.07_rk!0.0001054 !0.000000001_rk !
        !if (rank==0) close (unit=211)

        if (rank==0) then
            print '("P = ",f6.3)',P(1,:)
            print '("Y = ",f16.3)',Y
            
            inquire (iolength=requiredl)  P(1,:)
            open (unit = 201,file=trim(path) // 'params.txt',  status='unknown',recl=requiredl, action='write', IOSTAT = ios)
            !open (unit=201, file='..\\out\params.txt', status='unknown',recl=requiredl, action='write')

            write (201, * ) P(1,:)

            close (unit=201)

            print '("Generating files")'
        end if

        call solveValueFunction( params, grids, .FALSE., .FALSE. )

        if (rank==0) then
            call simWithUncer(params, grids, grids%Simy, cpath, apath, vpath, lpath, ypath ,AIME, TrueSPA, .TRUE.)
            call writetofile(grids,ypath, cpath, apath, vpath, lpath, grids%Simy,AIME)
            call writetofileByType(grids,ypath, cpath, apath, vpath, lpath, grids%Simy,AIME)
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
    if (dimEstimation == 7) params%lambda = control(7)
    if (rank==0) then
        !write (211,*) "Calcualting GMM for",  control(1),  control(2),  control(3),  control(4),  control(5),  control(6)
    end if
    gmm_criteria = gmm(params,grids,moments,weights) !*-1.0

    end function

    end program Console1