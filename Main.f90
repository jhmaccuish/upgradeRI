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
    real (kind=rk) :: agridslice(Tperiods+1,numPointsA), maxIncslice(Tperiods)

    !For regression
    real (kind=rk):: yreg(regPeriods, 3*numsims)
    real (kind=rk) :: xreg(regPeriods, 3*numsims, regPeriods-1+controls), assetstore( 3*numsims)
    logical :: mask(regPeriods,3*numsims) !, ols
    real (kind=rk) :: beta(regPeriods+2, 1)
    real (kind=rk), allocatable :: yreg2(:, :)
    real (kind=rk), allocatable :: xreg2(:, :, :)
    logical, allocatable :: mask2(:,:)

    real (kind=rk) :: start, finish
    real (kind=rk), allocatable :: moments(:,:),weights(:,:)
    INTEGER (KIND=8) :: c1loc, c2loc
    real (kind=rk), allocatable :: y(:) !, p(:,:)
    integer :: action, ios, requiredl, typeSim
    integer ::SPA,i, lamb
    logical :: delCache
    !INTEGER(kind=4) :: iter
    real (kind=rk) :: gmmScorce
    !real (kind=rk) :: indxrk
    !real (kind=rk) :: optLambda
    logical:: recal = .true.
    logical :: file_exists
    real (kind=rk), allocatable :: x(:), xl(:), xu(:)
    real (kind=rk) :: rhobeg, rhoend, medianA
    !real (kind=rk) :: importStructure(6,24*3*numsims),   !,amax(Tperiods)
    real (kind=rk), allocatable :: sterrs(:)
    character(len=1024) :: outFile
    integer :: j, dimAboveMed, counter, total
    character(len=2) :: temp2
    type (modelObjectsType) :: modelObjects(numPointsType,Tperiods)

#ifdef mpiBuild
    integer :: provided
#endif     

    !CHARACTER(len=255) :: cwd
    !call chdir( '..' )
    !CALL getcwd(cwd)
    !WRITE(*,*) TRIM(cwd)
    !call cpu_time(timeHack)
    call cpu_time(start)
    CALL SYSTEM_CLOCK(c1loc)

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
    write (*,*) "Done"
    allocate(x(dimEstimation), xl(dimEstimation), xu(dimEstimation),sterrs(dimEstimation))
    rhobeg = 0.4

    !historic excahnge rate 1.657347 taken from https://www.ofx.com/en-gb/forex-news/historical-exchange-rates/yearly-average-rates/
    !params%k= 215000*1.657347*(101.9/162.9)!from De Nardi, French, Jones
    params%k = 215*0.647491/0.8131991

    allocate(y(dimEstimation))
    open (unit = 304,file=trim(path) // 'paramsBin',  ACCESS="STREAM", action='read')
    read (304) y(1:4)
    close(unit=304)

    params%nu = y(1)
    params%beta =  y(2)
    params%gamma = y(3)
    params%thetab =  y(4)



    if (dimEstimation == 4) then
        params%lambda = 1.0
    else
        params%lambda =  1.0E-6 !1.0E-7 !100 !1.0E-5 !9.0E-8!
    end if

    params%BNDnu(1) =  0.1
    params%BNDnu(2) =  0.7
    params%BNDbeta(1) =  0.85
    params%BNDbeta(2) =  1.05
    params%BNDgamma(1) = 1.5
    params%BNDgamma(2) = 9
    params%BNDthetab(1) = 1.0 !0.01*700000000
    params%BNDthetab(2) = 30000.00 !100.0*700000000
    params%BNDlambda(1) = 1.0D-9 !0.01*700000000
    params%BNDlambda(2) = 0.1 !100.0*700000000

    x(1) = (params%nu - params%BNDnu(1))/(params%BNDnu(2)-params%BNDnu(1))
    x(2) = (params%beta - params%BNDbeta(1))/(params%BNDbeta(2)-params%BNDbeta(1))
    x(3) = (params%gamma - params%BNDgamma(1))/(params%BNDgamma(2)-params%BNDgamma(1))
    x(4) = (params%thetab -  params%BNDthetab(1))/(params%BNDthetab(2)-params%BNDthetab(1))
    if (dimEstimation == 5) then
        x(5) = (params%lambda -  params%BNDlambda(1))/(params%BNDlambda(2)-params%BNDlambda(1))
    end if

    call setupMisc(params,grids)
    params%startA = maxval(grids%initialAssets)

#ifdef debugMPI
    if (rank==0) pause
#endif

#ifdef _WIN64
    !if (rank == 0) then
    !    write (*,*) "Press 1 to check GMM values, 2 for single run, 3 to estimate, 4 repeat regressions, 5 Calculate SE"
    !    read (*,*) action
    !end if
    action =2
#ifdef mpiBuild
    call MPI_Bcast( action, 1, MPI_INTEGER ,0, mpi_comm_world, ierror)
    if (ierror.ne.0) stop 'mpi problem171'
#endif       
#else
    action =3
#endif 

    !Select the moments to match
    !if (action /= 2) then
    !        #ifdef _WIN64
    !        if (rank == 0) then
    !            write (*,*) "Press 1 to check GMM values, 2 for single run, 3 to estimate, 4 repeat regressions, 5 Calculate SE"
    !            read (*,*) action
    !        end if
    !#ifdef mpiBuild
    !        call MPI_Bcast( action, 1, MPI_INTEGER ,0, mpi_comm_world, ierror)
    !        if (ierror.ne.0) stop 'mpi problem171'
    !#endif
    !#else
    !        action =3
    !#endif
    select case (momentsToUse)
    case (1)
        allocate( moments(3,24),weights(3,24))
        open (unit = 1001,file=trim(pathMoments) // 'LSMoments_60.txt', action='read', IOSTAT = ios)
        open (unit = 1002,file=trim(pathMoments) // 'AssetMoments_60.txt', action='read', IOSTAT = ios)
        open (unit = 1003,file=trim(pathMoments) // 'WeightLS_60.txt', action='read', IOSTAT = ios)
        open (unit = 1004,file=trim(pathMoments) // 'weightAsset_60.txt', action='read', IOSTAT = ios)
        open (unit = 1005,file=trim(pathMoments) // 'Belief_60.txt', action='read', IOSTAT = ios)
        open (unit = 1006,file=trim(pathMoments) // 'WeightBelief_60.txt', action='read', IOSTAT = ios)
        read (1001, *) moments(1,:)
        read (1002, *) moments(2,:)
        read (1003, *) weights(1,:)
        read (1004, *) weights(2,:)
        if (dimEstimation >= 5) then
            read (1005, *) moments(3,1)
            read (1006, *) weights(3,1)
        end if
        weights = 1/max(weights,0.0000001)
        close (unit=1001)
        close (unit=1002)
        close (unit=1003)
        close (unit=1004)
        close (unit=1005)
        close (unit=1006)
    case (2)
        allocate( moments(3,1),weights(3,1))
        open (unit = 1001,file=trim(pathMoments) // 'RegA.txt', action='read', IOSTAT = ios)
        open (unit = 1002,file=trim(pathMoments) // 'RegAboveMed.txt', action='read', IOSTAT = ios)
        open (unit = 1003,file=trim(pathMoments) // 'WeightReg.txt', action='read', IOSTAT = ios)
        open (unit = 1004,file=trim(pathMoments) // 'WeightRegAboveMed.txt', action='read', IOSTAT = ios)
        open (unit = 1005,file=trim(pathMoments) // 'Belief_60.txt', action='read', IOSTAT = ios)
        open (unit = 1006,file=trim(pathMoments) // 'WeightBelief_60.txt', action='read', IOSTAT = ios)
        read (1001, *) moments(1,1)
        read (1002, *) moments(2,1)
        read (1003, *) weights(1,1)
        read (1004, *) weights(2,1)
        if (dimEstimation >= 5) then
            read (1005, *) moments(3,1)
            read (1006, *) weights(3,1)
        end if
        weights = 1/weights
        close (unit=1001)
        close (unit=1002)
        close (unit=1003)
        close (unit=1004)
        close (unit=1005)
        close (unit=1006)
    case(3)
        allocate( moments(3*3,24),weights(3*3,24))
        do spa=1,3
            write (temp2,'(I2)') 59+spa

            open (unit = 1001,file=trim(trim(pathMoments) // 'LSMoments_' // trim(adjustl(temp2)) ) // '.txt', action='read', IOSTAT = ios)
            open (unit = 1002,file=trim(trim(pathMoments) // 'AssetMoments_' // trim(adjustl(temp2)) ) // '.txt', action='read', IOSTAT = ios)
            open (unit = 1003,file=trim(trim(pathMoments) // 'WeightLS_' // trim(adjustl(temp2)) ) // '.txt', action='read', IOSTAT = ios)
            open (unit = 1004,file=trim(trim(pathMoments) // 'weightAsset_' // trim(adjustl(temp2)) ) // '.txt', action='read', IOSTAT = ios)
            open (unit = 1005,file=trim(trim(pathMoments) // 'Belief_' // trim(adjustl(temp2)) ) // '.txt', action='read', IOSTAT = ios)
            open (unit = 1006,file=trim(trim(pathMoments) // 'WeightBelief_' // trim(adjustl(temp2)) ) // '.txt', action='read', IOSTAT = ios)
            read (1001, *) moments((spa-1)*3+1,:)
            read (1002, *) moments((spa-1)*3+2,:)
            read (1003, *) weights((spa-1)*3+1,:)
            read (1004, *) weights((spa-1)*3+2,:)
            if (dimEstimation >= 5) then
                read (1005, *) moments((spa-1)*3+3,1)
                read (1006, *) weights((spa-1)*3+3,1)
            end if
            weights = 1/weights
            close (unit=1001)
            close (unit=1002)
            close (unit=1003)
            close (unit=1004)
            close (unit=1005)
            close (unit=1006)

        end do
    end select
    !end if
    !Set asset grid
    !do typeSim = 1, numPointsType
    !    call getassetgrid( params, grids%maxInc(typeSim,:), grids%Agrid(typeSim,:,:))
    !end do
    do typeSim = 1, numPointsType
        call getassetgrid( params,maxIncslice ,agridslice)
        grids%Agrid(typeSim,:,:) = agridslice
        grids%maxInc(typeSim,:) = maxIncslice
    end do

    if (action .EQ. 1) then
        !Set asset grid
        do typeSim = 1, numPointsType
            call getassetgrid( params, grids%maxInc(typeSim,:), grids%Agrid(typeSim,:,:))
        end do
        gmmScorce = gmm(params,grids,moments,weights, .false.)
        !
        if (rank==0) write (*,*), 'GMM Criteria is ', gmmScorce
        gmmScorce = gmm(params,grids,moments,weights, .false.,.false. )
        if (rank==0) write (*,*), 'GMM Criteria is ', gmmScorce
    else if (action .EQ. 2) then
        do typeSim = 1, numPointsType
            call getassetgrid( params,maxIncslice ,agridslice)
            grids%Agrid(typeSim,:,:) = agridslice
            grids%maxInc(typeSim,:) = maxIncslice
        end do
        !upper limit 10?
        lamb=  8 !6 !7
        !lamb = 1
        ! 9 !8 !0!  6 ! 5 !3 ! 0! 2! 4 !8 !3 !1 !2
        params%lambda=10.0**-lamb
        !params%gamma = 5
        !x(1) = 6.522674D-01
        !x(2) = 5.626104D-01
        !x(3) =5.224251D-01
        !x(4) = 6.383511D-01
        !x(5) = 4.051395D-01
        !params%nu = params%BNDnu(1)+x(1)*(params%BNDnu(2)-params%BNDnu(1))
        !params%beta = params%BNDbeta(1)+x(2)*(params%BNDbeta(2)-params%BNDbeta(1))
        !params%gamma = params%BNDgamma(1)+ x(3)*(params%BNDgamma(2)-params%BNDgamma(1))
        !params%thetab = params%BNDthetab(1)+x(4)*(params%BNDthetab(2)-params%BNDthetab(1))
        if (dimEstimation == 5) then
        params%lambda = params%BNDlambda(1)+x(5)*(params%BNDlambda(2)-params%BNDlambda(1))
        ! = locx5/abs(((0.6**(1-params%nu)*5000**params%nu)**(1-params%gamma))/(1-params%gamma))
        !write (*,*) params%lambda
        end if

        write (*,*)     params%nu, params%beta,  params%gamma,  params%thetab
         if (dimEstimation == 5)  write (*,*) params%lambda
        call solveOuterLoop(params, grids )
        !        do typeSim = 1, numPointsType
        !            call getassetgrid( params, grids%maxInc(typeSim,:), grids%Agrid(typeSim,:,:))
        !        end do
        !        do lamb=6,6 !3 !1 !2
        !            params%lambda=10.0**-lamb
        !            if (rank==0) write (*,*) "Lambda is ", params%lambda
        !            if (rank==0) write (*,*) params%nu, params%beta, params%gamma, params%thetab
        !            if (recal) call solveValueFunction( params, grids, .TRUE., .TRUE. )
        !            !simulate
        !#ifdef mpiBuild
        !            call mpi_barrier(mpi_comm_world, ierror)
        !            if (ierror.ne.0) stop 'mpi problem180'
        !#endif
        !            if (rank == 0) then
        !                delCache = .FALSE.
        !                if (modelChoice>1) then
        !                    if (counterFact) then
        !                        call simWithUncer(params, grids, grids%Simy, cpath, apath, lpath, ypath, AIME, 1, .FALSE. ) !5 !vpath,
        !                        call writetofile(grids, ypath, cpath, apath, vpath, lpath, grids%Simy,AIME)
        !                        call writetofileByType(grids,params,ypath, cpath, apath, vpath, lpath, grids%Simy,AIME)
        !                    else
        !                        delCache = .FALSE.
        !                        do SPA =1,3
        !                            write (*,*) "Sim True SPA", 59+SPA
        !                            if (SPA==3) delCache = .TRUE.
        !                            call simWithUncer(params, grids, grids%Simy, cpath, apath, lpath, ypath, AIME, SPA, .FALSE. ) !vpath,
        !                            combinedData(SPA,:,:) = apath(52-startAge+1:75-startAge+1,:)*(2*lpath(52-startAge+1:75-startAge+1,:)-1)
        !                            call writetofileAge(grids,params,ypath, cpath, apath, vpath, lpath, grids%Simy,AIME, 59+SPA)
        !                        end do
        !                        inquire (iolength=requiredl) combinedData
        !                        !print *, path
        !                        !print *, trim(path) // 'CombinedData'
        !                        open (unit=201, form="unformatted", file=trim(path) // 'CombinedData', ACCESS="STREAM", action='write', IOSTAT = ios)
        !                        write (201)  combinedData
        !                        close( unit=201)
        !                    end if
        !                else
        !                    call simWithUncer(params, grids, grids%Simy, cpath, apath,  lpath, ypath, AIME, TrueSPA, .FALSE. ) !vpath,
        !                    !do i = 1, Tperiods
        !                    !    amax(i) =  maxval( grids%Agrid(1,1, apath(i,:) ))
        !                    !end do
        !                    !write (*,*) maxval(amax)
        !                    call writetofile(grids,ypath, cpath, apath, vpath, lpath, grids%Simy,AIME)
        !                    call writetofileByType(grids,params,ypath, cpath, apath, vpath, lpath, grids%Simy,AIME)
        !                end if
        !            end if
        !        end do
    elseif (action .EQ. 3) then
        if (TrueSpa .NE. 1) then
            write (*,*) "Only estimate with SPA = 60"
        end if
        do typeSim = 1, numPointsType
            call getassetgrid( params, grids%maxInc(typeSim,:), grids%Agrid(typeSim,:,:))
        end do

        xl = 0.0
        xu = 1.0
        rhoend =  0.000001_rk!0.00001_rk
        INQUIRE(file=trim(path_bobyqa) // 'location', EXIST=file_exists)
        !if (rank==0) write(*,*) file_exists
        call bobyqa(dimEstimation, 2*dimEstimation+1, x, xl, xu, rhobeg, rhoend, 3, 2000, func, .TRUE., file_exists)
        !       end if
        if (rank==0) write (*,*) x
        if (rank==0) write (*,*) params%nu, params%beta, params%gamma, params%thetab
        if (rank==0) then
            y(1) = params%nu
            y(2) = params%beta
            y(3) = params%gamma
            y(4) = params%thetab
            if (dimEstimation==5) then
                y(5) = params%lambda
            end if

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
        call solveValueFunction( params, grids, modelObjects, .FALSE., .FALSE. )

        if (rank==0) then
            call simWithUncer(params, grids, grids%Simy, cpath, apath, lpath, ypath ,AIME, TrueSPA, modelObjects, .TRUE.) !vpath,
            call writetofile(grids,ypath, cpath, apath, vpath, lpath, grids%Simy,AIME)
            !call writetofileByType(grids,params,ypath, cpath, apath, vpath, lpath, grids%Simy,AIME)
        end if
#ifdef mpiBuild 
        call mpi_barrier(mpi_comm_world, ierror)
        if (ierror.ne.0) stop 'mpi problem180'
#endif   
    elseif (action .EQ. 4) then
        write (*,*)     params%nu, params%beta,  params%gamma,  params%thetab
        write (*,*) params%lambda
        do typeSim = 1, numPointsType
            call getassetgrid( params, grids%maxInc(typeSim,:), grids%Agrid(typeSim,:,:))
        end do

        xreg = 0.0
        do SPA =1,3
            if (rank == 0) write (*,*) ' SPA ', spa
            call setSPA(SPA)
            if (modelChoice==1 .OR. SPA==1) then
                if (recal) call solveValueFunction( params, grids, modelObjects, .TRUE., .TRUE. )
            end if
            if (rank == 0) then
                call simWithUncer(params, grids, grids%Simy, cpath, apath, lpath, ypath, AIME, SPA, modelObjects, .false. ) ! vpath,
                !combinedData(SPA,:,:) = apath(52-startAge+1:75-startAge+1,:)*(2*lpath(52-startAge+1:75-startAge+1,:)-1)
                call writetofileAge(grids,params,ypath, cpath, apath, vpath, lpath, grids%Simy,AIME, 59+SPA)
                do i = firstReg, lastReg
                    !y variable is labour supply
                    yreg(i-firstReg+1,1+(SPA-1)*numsims:SPA*numsims) = lpath(i,:)
                    !Indicator below spa
                    xreg(i-firstReg+1,1+(SPA-1)*numsims:SPA*numsims,1) = abs(i<Tretire+SPA-1)
                    !Wealth
                    !xreg(i-firstReg+1,1+(SPA-1)*numsims:SPA*numsims,2) = grids%Agrid(1,1, apath(i,:) )
                    if (i==9) then
                        assetstore(1+(SPA-1)*numsims:SPA*numsims) = grids%Agrid(1,1, apath(i,:) )
                    end if
                    !Constant
                    xreg(i-firstReg+1,1+(SPA-1)*numsims:SPA*numsims,controls) = 1.0
                    !Age dummy with first period exclude category
                    if (i>1) xreg(i,1+(SPA-1)*numsims:SPA*numsims,i-firstReg+controls) = 1.0
                end do
            end if
            if (modelChoice==1 ) then
                do i=1,Tperiods
                    do j=1,4
                        if (allocated(modelObjects(j,i)%EV) ) deallocate(modelObjects(j,i)%EV )
                        if (allocated(modelObjects(j,i)%policy) ) deallocate(modelObjects(j,i)%policy)
                        if (allocated(modelObjects(j,i)%V) ) deallocate(modelObjects(j,i)%V)
                    end do
                end do
            end if
        end do

        if (rank == 0) then
            mask = .true.
            call doReg(yreg, xreg, 3*numsims, regPeriods, regPeriods-1+controls, mask, .false., beta)
            write (*,*) 'Treatment effect', beta(1,1)
            !inquire (iolength=requiredl) combinedData
            !print *, path
            !print *, trim(path) // 'CombinedData'
            !open (unit=201, form="unformatted", file=trim(path) // 'CombinedData', ACCESS="STREAM", action='write', IOSTAT = ios)
            !write (201)  combinedData
            !close( unit=201)

            !importStructure = 0.0
            !medianA = median(xreg(9,:,2))
            medianA = 34869.00 ! 24846.00
            write (*,*) "median assets", medianA
            dimAboveMed = count(assetstore>medianA)
            Total = count(xreg(9,:,2)>0.0)
            write (*,*) dimAboveMed, 'of', total, 'above med'
            if (dimAboveMed > 10) then
                counter = 0
                allocate(yreg2(regPeriods, dimAboveMed), xreg2(regPeriods, dimAboveMed, regPeriods-1+controls), mask2(regPeriods,dimAboveMed))
                do j = 1,3*numsims
                    if (assetstore(j) > medianA) then
                        counter = counter + 1
                        !write(*,*) counter
                        yreg2(:, counter) = yreg(:, j)
                        xreg2(:, counter, :) =  xreg(:, j, :)
                        mask2(:,counter) = mask(:,j)
                    end if
                    do i= firstReg, lastReg

                        !importStructure(1,(j-1)*24+i) = j
                        !importStructure(2,(j-1)*24+i) = 51 + i
                        !importStructure(3,(j-1)*24+i) = yreg(i,j)
                        !if (j <= numsims ) importStructure(4,(j-1)*24+i) = 60
                        !if (j <= 2*numsims .AND. j> numsims ) importStructure(4,(j-1)*24+i) = 61
                        !if ( j> 2*numsims ) importStructure(4,(j-1)*24+i) = 62
                        !!if (i >=11) importStructure(4,(j-1)*24+i) = 1.0 !xreg(i,j,2)
                        !!!xreg(i,1+(SPA-1)*numsims:SPA*numsims,3) = 1.0
                        !!if (i >=10 .AND. j <= 2*numsims ) importStructure(4,(j-1)*24+i) = 1.0
                        !!if (i >=9 .AND. j <= numsims ) importStructure(4,(j-1)*24+i) = 1.0
                        !importStructure(5,(j-1)*24+i) = xreg(i,j,2)
                        !importStructure(6,(j-1)*24+i) = xreg(i,j,1)
                        !!if (i>1) xreg(i,j,1+i) =1.0
                        !!xreg(i, j, 24+2) = 1.0
                    end do
                end do
                write (*,*) "Mean wealth of those above median: ", sum(xreg2(9, :, 2))/dimAboveMed
                !write (outfile,'(I1)') lamb
                !outfile = trim(path) // 'testcase_output_fort' // trim(outfile) // '.txt'
                !open (unit = 1001,file=outfile, action='write', IOSTAT = ios)
                !write (1001,  '(6F15.2)',  IOSTAT = ios) importStructure
                !close (unit=1001)

                call doReg(yreg2, xreg2, dimAboveMed, regPeriods, regPeriods-1+controls, mask2, .false., beta)
                write (*,*) 'Treatment effect above median', beta(1,1)
                deallocate(yreg2, xreg2, mask2)
            end if
        end if

    elseif (action .EQ. 5) then
        call getStandardErrors(dimEstimation, params, grids, moments, weights, sterrs)

    end if

    if (rank==0) then
        call cpu_time(finish)
        CALL SYSTEM_CLOCK(c2loc)

        print '("Time = ",f11.3," seconds.")', (c2loc-c1loc)/rate
        print '("Tot CPU time = ",f11.3," seconds.")', finish-start
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
    real(rk) :: locx5

    params%nu = params%BNDnu(1)+x(1)*(params%BNDnu(2)-params%BNDnu(1))
    params%beta = params%BNDbeta(1)+x(2)*(params%BNDbeta(2)-params%BNDbeta(1))
    params%gamma = params%BNDgamma(1)+ x(3)*(params%BNDgamma(2)-params%BNDgamma(1))
    params%thetab = params%BNDthetab(1)+x(4)*(params%BNDthetab(2)-params%BNDthetab(1))
    if (n == 5) then
        params%lambda = params%BNDlambda(1)+x(5)*(params%BNDlambda(2)-params%BNDlambda(1))
        ! = locx5/abs(((0.6**(1-params%nu)*5000**params%nu)**(1-params%gamma))/(1-params%gamma))
        !write (*,*) params%lambda
    end if

    f = gmm(params,grids,moments,weights)
    end subroutine func

    end program Console1
