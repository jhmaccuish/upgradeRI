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
    use bobyqa_module

    implicit none

#ifdef mpiBuild
    include 'mpif.h'
#endif        

#include "globalMacros.txt"

    ! Variables
    type (structparamstype) :: params
    type (gridsType) :: grids

    real (kind=rk) :: ypath(Tperiods, numSims) !income
    real (kind=rk) :: cpath(Tperiods, numSims)  !consumption
    integer :: lpath(Tperiods, numSims) !labour supply
    real (kind=rk) :: vpath(Tperiods, numSims)  !value
    integer :: apath(Tperiods + 1,numSims) !this is the path at the start of each period, so we include the 'start' of death
    integer :: combinedData(3,24,numSims)
    real (kind=rk) :: AIME(Tperiods + 1,numSims)

    !For regression
    real (kind=rk):: yreg(24, 3*numsims)
    real (kind=rk) :: xreg(24, 3*numsims, 24+2)
    logical :: mask(24,3*numsims), ols
    real (kind=rk) :: beta(24+2, 1)
    real (kind=rk), allocatable :: yreg2(:, :)
    real (kind=rk), allocatable :: xreg2(:, :, :)
    logical, allocatable :: mask2(:,:)

    real (kind=rk) :: start, finish, moments(2,24),weights(2,24)
    real (kind=rk), allocatable :: y(:), p(:,:)
    integer :: action, ios, requiredl, typeSim
    integer ::SPA,i, lamb
    logical :: delCache
    INTEGER(kind=4) :: iter
    real (kind=rk) :: gmmScorce
    real (kind=rk) :: indxrk
    real (kind=rk) :: optLambda
    logical:: recal = .true.
    logical :: file_exists
    real (kind=rk), allocatable :: x(:), xl(:), xu(:)
    real (kind=rk) :: rhobeg, rhoend
    real (kind=rk) :: importStructure(6,24*3*numsims), amax(Tperiods), medianA
    real (kind=rk), allocatable :: sterrs(:)
    character(len=1024) :: outFile
    integer :: j, dimAboveMed, counter

#ifdef mpiBuild
    integer :: provided
#endif     

    !CHARACTER(len=255) :: cwd
    !call chdir( '..' )
    !CALL getcwd(cwd)
    !WRITE(*,*) TRIM(cwd)
    !call cpu_time(timeHack)
    call cpu_time(start)

#ifdef mpiBuild    
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
    !if (NM) then
    !    allocate(y(dimEstimation+1), p(dimEstimation+1,dimEstimation))
    !else
    allocate(x(dimEstimation), xl(dimEstimation), xu(dimEstimation),sterrs(dimEstimation))
    rhobeg = 0.4
    !end if


    !    !Earning type 1 - uneducated couple
    !    params%delta(1,3) = 9.083298
    !    params%delta(1,2)=0.029333
    !    params%delta(1,1)=-0.00023
    !
    !    params%spouseInc(1) = 5121!6235.8998884204320
    !
    !    !Earning type 2 - uneducated single
    !#if TYPES_SIZE > 1
    !    params%delta(2,3) = 8.5256042
    !    params%delta(2,2)= 0.061582208
    !    params%delta(2,1)= -0.0005901
    !    params%spouseInc(2) = 5282
    !#endif
    !
    !    !Earning type 3 - educated single
    !#if TYPES_SIZE > 2
    !    params%delta(3,3) = 7.906901
    !    params%delta(3,2)= 0.094386
    !    params%delta(3,1)= -0.00091
    !    params%spouseInc(3) = 6840
    !#endif
    !
    !#if TYPES_SIZE > 3
    !    params%delta(4,3) = 7.392151
    !    params%delta(4,2)= 0.121426
    !    params%delta(4,1)= -0.00117
    !    params%spouseInc(4) = 7684
    !#endif


    !params%thetab = 700000000 !2.899476882946300E-002
    !historic excahnge rate 1.657347 taken from https://www.ofx.com/en-gb/forex-news/historical-exchange-rates/yearly-average-rates/
    !params%k= 215000*1.657347*(101.9/162.9)!from De Nardi, French, Jones
    params%k = 215*0.647491/0.8131991


    params%nu = 0.7 !0.465961054093721 !0.489447013499279! 0.455688293518331 ! 0.518965309049116 !  0.509238327116098  !0.503632895230470 !0.525701440155937 ! 0.330639470815664
    params%beta = 1.0 !0.995215700427888 ! 0.973148936962498  !0.986727981493142! 0.996592226255899  !  0.986357800876955  ! 1.00000000000000! 0.999965608083641  ! 0.986701722470645
    params%gamma = 7 !  4.00794815142518  !3.17936997498581! 4.25299731750991 ! 5.49316011629655 !6.32167143715174! 4.31646852872327 !6.59553977038513 !  1.25523751463278
    !params%thetab = 100.0*700000000 ! 28072812873.2562 ! 28319331859.8928 !942211880.074711 ! 28778043030.9356 !37857189401.5407 ! 43461745124.7349 ! 8380790861.40518


    allocate(y(4))
    open (unit = 304,file=trim(path) // 'paramsBin',  ACCESS="STREAM", action='read')
    read (304) y
    close(unit=304)

    params%nu = y(1)
    params%beta = y(2)
    params%gamma = y(3)
    params%thetab = y(4)
    !params%thetab = 100.00

    params%lambda = 0.00001
    !params%lambda =  0.000001 !1.0 !1.0E-8 !1.0 !100.0 !9.996559952862711E-004
    params%lambda = 0.0000001
    !params%lambda = 0.000001 !0.0001
    if (dimEstimation == 4) then
        params%lambda = 1.0
    else
        params%lambda =  1.0E-6 !1.0E-7 !100 !1.0E-5 !9.0E-8!
    end if
    !params%thetab = 100 !39991832534.4503! 700000000 !  30000.00 !0.0 !params%thetab - 0.25*params%thetab

    params%BNDnu(1) =  0.1
    params%BNDnu(2) =  0.7
    params%BNDbeta(1) =  0.85
    params%BNDbeta(2) =  1.05
    params%BNDgamma(1) = 1.1
    params%BNDgamma(2) = 9
    params%BNDthetab(1) = 1.0 !0.01*700000000
    params%BNDthetab(2) = 30000.00 !100.0*700000000

    x(1) = (params%nu - params%BNDnu(1))/(params%BNDnu(2)-params%BNDnu(1))
    x(2) = (params%beta - params%BNDbeta(1))/(params%BNDbeta(2)-params%BNDbeta(1))
    x(3) = (params%gamma - params%BNDgamma(1))/(params%BNDgamma(2)-params%BNDgamma(1))
    x(4) = (params%thetab -  params%BNDthetab(1))/(params%BNDthetab(2)-params%BNDthetab(1))

    call setupMisc(params,grids)
    params%startA = maxval(grids%initialAssets)

#ifdef debugMPI
    if (rank==0) pause
#endif

#ifdef _WIN64
    if (rank == 0) then
        write (*,*) "Press 1 to check GMM values, 2 for single run, 3 to estimate, 4 repeat regressions, 5 Calculate SE"
        read (*,*) action
    end if
#ifdef mpiBuild
    call MPI_Bcast( action, 1, MPI_INTEGER ,0, mpi_comm_world, ierror)
    if (ierror.ne.0) stop 'mpi problem171'
#endif       
#else
    action =3
#endif  
    open (unit = 1001,file=trim(pathMoments) // 'moments.txt', action='read', IOSTAT = ios)
    open (unit = 1002,file=trim(pathMoments) // 'assetmom.txt', action='read', IOSTAT = ios)
    open (unit = 1003,file=trim(pathMoments) // 'Weight.txt', action='read', IOSTAT = ios)
    open (unit = 1004,file=trim(pathMoments) // 'weightAsset.txt', action='read', IOSTAT = ios)
    read (1001, *) moments(1,:)
    read (1002, *) moments(2,:)
    read (1003, *) weights(1,:)
    read (1004, *) weights(2,:)
    weights = 1/weights
    close (unit=1001)
    close (unit=1002)
    close (unit=1003)
    close (unit=1004)
    !Set asset grid
    do typeSim = 1, numPointsType
        call getassetgrid( params, grids%maxInc(typeSim,:), grids%Agrid(typeSim,:,:))
    end do

    if (action .EQ. 1) then
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
        weights = 1/weights
        close (unit=1001)
        close (unit=1002)
        close (unit=1003)
        close (unit=1004)
        !Set asset grid
        do typeSim = 1, numPointsType
            call getassetgrid( params, grids%maxInc(typeSim,:), grids%Agrid(typeSim,:,:))
        end do
        !gmmScorce =  gmm_criteria(x)
        gmmScorce = gmm(params,grids,moments,weights, .false.)
        !
        if (rank==0) write (*,*), 'GMM Criteria is ', gmmScorce
        gmmScorce = gmm(params,grids,moments,weights, .false.,.false. )
        if (rank==0) write (*,*), 'GMM Criteria is ', gmmScorce
        stop
        gmmScorce = gmm(params,grids,moments,weights, .false.,.false. )
        if (rank==0) write (*,*), 'GMM Criteria is ', gmmScorce
        gmmScorce = gmm(params,grids,moments,weights, .false.,.false. )
        if (rank==0) write (*,*), 'GMM Criteria is ', gmmScorce
        gmmScorce = gmm(params,grids,moments,weights, .false.,.false. )
        if (rank==0) write (*,*), 'GMM Criteria is ', gmmScorce

        !call func(4, x, gmmScorce)
        !if (rank==0) write (*,*), 'GMM Criteria is ', gmmScorce
        !p(1,1) = 0.236200045532760 !0.240191353041616
        !p(1,2) = 0.999995559030882  !0.999477774429760
        !p(1,3) = 1.57028307631512 ! 1.93680557456545
        !p(1,4) = 0.887053363134520  !0.493786285559129
        !p(1,5) = -2.864712075930465E-006 ! -3.533369493611974E-006
        !p(1,6) = 1.962853923390662E-002 !2.421007064408893E-002
        !optLambda =  9.996559952862711E-004 !0.836387480497043
        !
        !open (unit = 666, form="unformatted", file=trim(path) // 'guessP.txt', status='unknown', ACCESS="STREAM", action='read', IOSTAT = ios)
        !read (666) p
        !close (unit=666)
        !
        !open (unit = 667, form="unformatted", file=trim(path) // 'guessY.txt', status='unknown', ACCESS="STREAM", action='read', IOSTAT = ios)
        !read (667) Y
        !close (unit=667)
        !
        !optLambda = p(1,7)
        !!Lambda = 0.436387480497043 is best!!!!
        !do indxrk = -0.3, 0.3, 0.1
        !    p(1,7) = optLambda + indxrk*optLambda
        !    if (rank==0) write (*,*), 'For lambda = ', p(1,7), indxrk
        !    gmmScorce =  gmm_criteria(p(1,:))
        !    if (rank==0) write (*,*), 'GMM Criteria is ', gmmScorce
        !end do
    else if (action .EQ. 2) then

        do typeSim = 1, numPointsType
            call getassetgrid( params, grids%maxInc(typeSim,:), grids%Agrid(typeSim,:,:))
        end do
        do lamb=6,6 !3 !1 !2
            params%lambda=10.0**-lamb
            if (rank==0) write (*,*) "Lambda is ", params%lambda
            if (rank==0) write (*,*) params%nu, params%beta, params%gamma, params%thetab
            if (recal) call solveValueFunction( params, grids, .TRUE., .TRUE. )
            !simulate
#ifdef mpiBuild
            call mpi_barrier(mpi_comm_world, ierror)
            if (ierror.ne.0) stop 'mpi problem180'
#endif         
            if (rank == 0) then
                delCache = .FALSE.
                if (modelChoice>1) then
                    if (counterFact) then
                        call simWithUncer(params, grids, grids%Simy, cpath, apath, vpath, lpath, ypath, AIME, 1, .FALSE. ) !5
                        call writetofile(grids, ypath, cpath, apath, vpath, lpath, grids%Simy,AIME)
                        call writetofileByType(grids,params,ypath, cpath, apath, vpath, lpath, grids%Simy,AIME)
                    else
                        delCache = .FALSE.
                        do SPA =1,3
                            write (*,*) "Sim True SPA", 59+SPA
                            if (SPA==3) delCache = .TRUE.
                            call simWithUncer(params, grids, grids%Simy, cpath, apath, vpath, lpath, ypath, AIME, SPA, .FALSE. )
                            combinedData(SPA,:,:) = apath(52-startAge+1:75-startAge+1,:)*(2*lpath(52-startAge+1:75-startAge+1,:)-1)
                            call writetofileAge(grids,params,ypath, cpath, apath, vpath, lpath, grids%Simy,AIME, 59+SPA)
                        end do
                        inquire (iolength=requiredl) combinedData
                        !print *, path
                        !print *, trim(path) // 'CombinedData'
                        open (unit=201, form="unformatted", file=trim(path) // 'CombinedData', ACCESS="STREAM", action='write', IOSTAT = ios)
                        write (201)  combinedData
                        close( unit=201)
                    end if
                else
                    call simWithUncer(params, grids, grids%Simy, cpath, apath, vpath, lpath, ypath, AIME, TrueSPA, .FALSE. )
                    !do i = 1, Tperiods
                    !    amax(i) =  maxval( grids%Agrid(1,1, apath(i,:) ))
                    !end do
                    !write (*,*) maxval(amax)
                    call writetofile(grids,ypath, cpath, apath, vpath, lpath, grids%Simy,AIME)
                    call writetofileByType(grids,params,ypath, cpath, apath, vpath, lpath, grids%Simy,AIME)
                end if
            end if
        end do
    elseif (action .EQ. 3) then
        if (TrueSpa .NE. 1) then
            write (*,*) "Only estimate with SPA = 60"
        end if
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
        weights = 1/weights
        !        !if (NM)   then
        !            if (rank==0) then
        !                print '("Setting up initial guess for hilling climbing algorithm")'
        !            end if
        !
        !#ifdef _WIN64
        !            if (rank == 0) then
        !                write (*,*) "Press 1 to load, 2 to do fresh"
        !                read (*,*) action
        !            end if
        !#ifdef mpiBuild
        !            call MPI_Bcast( action, 1, MPI_INTEGER ,0, mpi_comm_world, ierror)
        !            if (ierror.ne.0) stop 'mpi problem171'
        !#endif
        !
        !            if (action == 1) then
        !                open (unit = 666, form="unformatted", file=trim(path) // 'guessP', status='unknown', ACCESS="STREAM", action='read', IOSTAT = ios)
        !                read (666) p
        !                close (unit=666)
        !
        !                open (unit = 667, form="unformatted", file=trim(path) // 'guessY', status='unknown', ACCESS="STREAM", action='read', IOSTAT = ios)
        !                read (667) Y
        !                close (unit=667)
        !            else
        !                call initialGuess(rank,params,grids,moments,weights,p,y)
        !            end if
        !#else
        !            call initialGuess(rank,params,grids,moments,weights,p,y)
        !#endif
        !
        !
        !#ifdef mpiBuild
        !            call mpi_barrier(mpi_comm_world, ierror)
        !            if (ierror.ne.0) stop 'mpi problem180'
        !#endif
        !
        !            !open (unit=201, form="unformatted", file=trim(path) // 'CombinedData', ACCESS="STREAM", action='write', IOSTAT = ios)
        !            !open (unit=211, file='..\\out\guess2.txt', status='unknown', action='write')
        !            call amoeba(p,y, 0.00001_rk,gmm_criteria,iter,.TRUE.,0.0_rk) !0.0001_rk!0.001_rk!0.0001_rk !0.001_rk!0.07_rk!0.0001054 !0.000000001_rk !
        !            !if (rank==0) close (unit=211)
        !
        !            if (rank==0) then
        !                print '("P = ",f6.3)',P(1,:)
        !                print '("Y = ",f16.3)',Y
        !
        !                inquire (iolength=requiredl)  P(1,:)
        !                open (unit = 201,file=trim(path) // 'params.txt',  status='unknown',recl=requiredl, action='write', IOSTAT = ios)
        !                !open (unit=201, file='..\\out\params.txt', status='unknown',recl=requiredl, action='write')
        !
        !                write (201, * ) P(1,:)
        !
        !                close (unit=201)
        !
        !                print '("Generating files")'
        !            end if
        !        else
        xl = 0.0
        xu = 1.0
        rhoend =  0.00001_rk !0.0000001_rk !
        INQUIRE(file=trim(path_bobyqa) // 'location', EXIST=file_exists)
        !if (rank==0) write(*,*) file_exists
        call bobyqa (dimEstimation, 2*dimEstimation+1, x, xl, xu, rhobeg, rhoend, 3, 2000, func, .TRUE., file_exists)
        !       end if
        if (rank==0) write (*,*) x
        if (rank==0) write (*,*) params%nu, params%beta, params%gamma, params%thetab
        if (rank==0) then
            !allocate(y(dimEstimation))
            y(1) = params%nu
            y(2) = params%beta
            y(3) = params%gamma
            y(4) = params%thetab
            inquire (iolength=requiredl)  y
            open (unit = 302,file=trim(path) // 'params.txt', recl=requiredl, action='write', IOSTAT = ios)
            write (*,*) ios
            write (302, * ) y
            close(unit=302)

            inquire (iolength=requiredl)  x
            open (unit = 303,file=trim(path) // 'x.txt', recl=requiredl, action='write', IOSTAT = ios)
            write (303, * ) x
            close(unit=303)

            open (unit = 304,file=trim(path) // 'paramsBin', status='replace', ACCESS="STREAM", action='write')
            write (*,*) ios
            write (304) y
            close(unit=304)

            open (unit = 305,file=trim(path) // 'xBin',   status='replace', ACCESS="STREAM", action='write')
            write (305) x
            close(unit=305)

        end if
        call solveValueFunction( params, grids, .FALSE., .FALSE. )

        if (rank==0) then
            call simWithUncer(params, grids, grids%Simy, cpath, apath, vpath, lpath, ypath ,AIME, TrueSPA, .TRUE.)
            call writetofile(grids,ypath, cpath, apath, vpath, lpath, grids%Simy,AIME)
            !call writetofileByType(grids,params,ypath, cpath, apath, vpath, lpath, grids%Simy,AIME)
        end if
#ifdef mpiBuild 
        call mpi_barrier(mpi_comm_world, ierror)
        if (ierror.ne.0) stop 'mpi problem180'
#endif   
    elseif (action .EQ. 4) then
        do lamb=1,1
            params%lambda=10.0**-lamb
            write (*,*) "Lambda is ", params%lambda
            do typeSim = 1, numPointsType
                call getassetgrid( params, grids%maxInc(typeSim,:), grids%Agrid(typeSim,:,:))
            end do

            xreg = 0.0
            do SPA =1,3
                if (rank == 0) write (*,*) ' SPA ', spa
                call setSPA(SPA)
                if (modelChoice==1 .OR. SPA==1) then
                    if (recal) call solveValueFunction( params, grids, .TRUE., .TRUE. )
                end if
                if (rank == 0) then
                    call simWithUncer(params, grids, grids%Simy, cpath, apath, vpath, lpath, ypath, AIME, TrueSPA, .false. )
                    combinedData(SPA,:,:) = apath(52-startAge+1:75-startAge+1,:)*(2*lpath(52-startAge+1:75-startAge+1,:)-1)
                    call writetofileAge(grids,params,ypath, cpath, apath, vpath, lpath, grids%Simy,AIME, 59+SPA)
                    do i = 1, 24
                        yreg(i,1+(SPA-1)*numsims:SPA*numsims) = lpath(i,:)
                        xreg(i,1+(SPA-1)*numsims:SPA*numsims,1) = abs(i<Tretire+SPA-1)
                        xreg(i,1+(SPA-1)*numsims:SPA*numsims,2) = grids%Agrid(1,1, apath(i,:) )
                        xreg(i,1+(SPA-1)*numsims:SPA*numsims,3) = 1.0
                        if (i>1) xreg(i,1+(SPA-1)*numsims:SPA*numsims,i+2) = 1.0
                    end do
                end if
            end do
            !real (kind=rk) :: ypath(Tperiods, numSims)
            !        real (kind=rk), intent(in) :: y(t, n)
            !real (kind=rk), intent(in) :: x(t, n, k)
            !logical, intent(in) :: mask(t,n), ols
            !real (kind=rk), intent(out) :: beta(k, 1)
            if (rank == 0) then
                mask = .true.
                call doReg(yreg, xreg, 3*numsims, 24, 24+2, mask, .false., beta)
                write (*,*) 'Treatment effect', beta(1,1)
                inquire (iolength=requiredl) combinedData
                !print *, path
                !print *, trim(path) // 'CombinedData'
                !open (unit=201, form="unformatted", file=trim(path) // 'CombinedData', ACCESS="STREAM", action='write', IOSTAT = ios)
                !write (201)  combinedData
                !close( unit=201)

                importStructure = 0.0
                medianA = median(xreg(9,:,2))
                dimAboveMed = count(xreg(9,:,2)>medianA)
                write (*,*) dimAboveMed
                counter = 0
                allocate(yreg2(24, dimAboveMed), xreg2(24, dimAboveMed, 24+2), mask2(24,dimAboveMed))
                do j = 1,3*numsims
                    if (xreg(9,j,2) > medianA) then
                        counter = counter + 1
                        !write(*,*) counter
                        yreg2(:, counter) = yreg(:, j)
                        xreg2(:, counter, :) =  xreg(:, j, :)
                        mask2(:,counter) = mask(:,j)
                    end if
                    do i=1,24

                        importStructure(1,(j-1)*24+i) = j
                        importStructure(2,(j-1)*24+i) = 51 + i
                        importStructure(3,(j-1)*24+i) = yreg(i,j)
                        if (j <= numsims ) importStructure(4,(j-1)*24+i) = 60
                        if (j <= 2*numsims .AND. j> numsims ) importStructure(4,(j-1)*24+i) = 61
                        if ( j> 2*numsims ) importStructure(4,(j-1)*24+i) = 62
                        !if (i >=11) importStructure(4,(j-1)*24+i) = 1.0 !xreg(i,j,2)
                        !!xreg(i,1+(SPA-1)*numsims:SPA*numsims,3) = 1.0
                        !if (i >=10 .AND. j <= 2*numsims ) importStructure(4,(j-1)*24+i) = 1.0
                        !if (i >=9 .AND. j <= numsims ) importStructure(4,(j-1)*24+i) = 1.0
                        importStructure(5,(j-1)*24+i) = xreg(i,j,2)
                        importStructure(6,(j-1)*24+i) = xreg(i,j,1)
                        !if (i>1) xreg(i,j,1+i) =1.0
                        !xreg(i, j, 24+2) = 1.0
                    end do
                end do
                write (outfile,'(I1)') lamb
                outfile = trim(path) // 'testcase_output_fort' // trim(outfile) // '.txt'
                open (unit = 1001,file=outfile, action='write', IOSTAT = ios)
                write (1001,  '(6F15.2)',  IOSTAT = ios) importStructure
                close (unit=1001)

                call doReg(yreg2, xreg2, dimAboveMed, 24, 24+2, mask2, .false., beta)
                write (*,*) 'Treatment effect above median', beta(1,1)
                deallocate(yreg2, xreg2, mask2)
            end if
        end do
    elseif (action .EQ. 5) then
        call getStandardErrors(dimEstimation, params, grids, moments, weights, sterrs)

    end if

    if (rank==0) then
        call cpu_time(finish)
        print '("Time = ",f11.3," seconds.")',finish-start
    end if
#ifdef mpiBuild
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
    if (minval(control(1:4)) < 0 )  then !.OR. control(6)< 0
        gmm_criteria = huge(gmm_criteria)
        return
    end if
    params%nu = control(1)
    params%beta = control(2)
    params%gamma = control(3)
    params%thetab = control(4)
    if (dimEstimation == 5) params%lambda = control(5)
    if (rank==0) then
        !write (211,*) "Calcualting GMM for",  control(1),  control(2),  control(3),  control(4),  control(5),  control(6)
    end if
    gmm_criteria = gmm(params,grids,moments,weights) !*-1.0

    end function

    subroutine func (n, x, f)  !! calfun interface
    implicit none
    integer,intent(in)               :: n
    real(rk),dimension(:),intent(in) :: x
    real(rk),intent(out)             :: f
    params%nu = params%BNDnu(1)+x(1)*(params%BNDnu(2)-params%BNDnu(1))
    params%beta = params%BNDbeta(1)+x(2)*(params%BNDbeta(2)-params%BNDbeta(1))
    params%gamma = params%BNDgamma(1)+ x(3)*(params%BNDgamma(2)-params%BNDgamma(1))
    params%thetab = params%BNDthetab(1)+x(4)*(params%BNDthetab(2)-params%BNDthetab(1))
    if (n == 5) params%lambda = x(5)

    f = gmm(params,grids,moments,weights)
    end subroutine func

    end program Console1