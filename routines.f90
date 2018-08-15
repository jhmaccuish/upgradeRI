    module routines
    use Header
    use routines_generic

    implicit none
!#include "../globalMacros.txt"
#include "globalMacros.txt"

    interface unpackArray
    module procedure unpackArrayInt
    module procedure unpackArrayReal
    module procedure unpackArraySmall
    end interface unpackArray

    contains

    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    !!Get Income grid
    subroutine getIncomeGrid(params, grids)
    implicit none

    !Changing
    type (structparamstype), intent(inout) :: params
    type (gridsType), intent(inout) :: grids

    !outputs
    !real (kind=rk) :: YgridOut(:,:), incTransitionMrxOut(:,:), maxInc(:), AIMEgrid(:,:), benefit(:)

    !local
    real (kind=rk) :: sig_inc, ly(numPointsProd), upper(numPointsType,Tperiods+1), a !Q(numPointsY,numPointsY)
    real (kind=rk) :: Ygrid(numPointsType,Tperiods,numPointsProd), incTransitionMrx(numPointsProd,numPointsProd)
    integer :: t, i, workAge, typeSim

    sig_inc = params%sigma/((1-params%rho**2)**0.5)

    call tauchen(numPointsProd,params%mu,params%rho,params%sigma,normBnd,ly,incTransitionMrx)

    params%pension = 107.45*52
    upper(:,1) = 0
    do typeSim = 1, numPointsType
        do t=1 , Tperiods
            workAge = startAge - 20 + t
            Ygrid(typeSim,t,:)= exp(ly+params%delta(typeSim,1)*workAge**2+params%delta(typeSim,2)*workAge+params%delta(typeSim,3)-grids%fc(t))
            !if (typeSim==numPointsType) then
            grids%maxInc(typeSim,t) = exp((normBnd * sig_inc)+params%delta(typeSim,1)*workAge**2+params%delta(typeSim,2)*workAge+params%delta(typeSim,3)-grids%fc(t))
            !This should be discounted to average not sum
            upper(typeSim,t+1) = upper(typeSim,t) + grids%maxInc(typeSim,t)
            if (t <=spouseretire) then
                a = (upper(typeSim,t+1)/t)/(numAIME-1)
                grids%AIMEgrid(typeSim,t,:) = a*(/(i,i=0,numAIME-1)/)
            else
                grids%AIMEgrid(typeSim,t,:) = grids%AIMEgrid(typeSim,t-1,:)
            end if
            !if (t <=5) then
            !    grids%benefit(t) = 57.90*52
            !else
            grids%benefit(t) = 0!73.10*52
            !end if
            !end if
            !grids%YgridExtrap(typeSim,t,1) = 0.0
            grids%YgridExtrap(typeSim,t,:) = Ygrid(typeSim,t,:)
        end do
        call addUnemploymentshock(Ygrid(typeSim,:,:), incTransitionMrx,grids,typeSim)
        grids%AIMEgrid(typeSim,Tperiods+1,:) = grids%AIMEgrid(typeSim,Tperiods,:)
    end do


    end subroutine

    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    !!Inpliments tauchen to discretise AR1
    subroutine tauchen(N,mu,rho,sigma,m,Z,Zprob)

    implicit none
    !Function TAUCHEN
    !
    !Purpose:    Finds a Markov chain whose sample paths
    !            approximate those of the AR(1) process
    !                z(t+1) = (1-rho)*mu + rho * z(t) + eps(t+1)
    !            where eps are normal with stddev sigma
    !
    !Format:     {Z, Zprob} = Tauchen(N,mu,rho,sigma,m)
    !
    !Input:      N       scalar, number of nodes for Z
    !            mu      scalar, unconditional mean of process
    !            rho     scalar
    !            sigma   scalar, std. dev. of epsilons
    !            m       max +- std. devs.
    !
    !Output:     Z       N*1 vector, nodes for Z
    !            Zprob   N*N matrix, transition probabilities
    !
    !    Martin Floden
    !    Fall 1996
    !
    !    This procedure is an implementation of George Tauchen's algorithm
    !    described in Ec. Letters 20 (1986) 177-181.
    !

    !inputs
    integer, intent(in) :: N, m
    real (kind=rk), intent(in):: mu, rho, sigma

    !outputs
    real (kind=rk), intent(out) :: Z(N), Zprob(N,N)

    !local
    real (kind=rk) :: a, zstep
    integer :: i, j, k

    a     = (1-rho)*mu;

    Z(N)  = m * sqrt(sigma**2 / (1 - rho**2))
    Z(1)  = -Z(N)
    zstep = (Z(N) - Z(1)) / (N - 1)

    do i = 2, (N-1)
        Z(i) = Z(1) + zstep * (i - 1)
    end do

    Z = Z + a / (1-rho);

    do j = 1, N
        do k = 1, N
            if (k == 1) then
                Zprob(j,k) = cdf_normal((Z(1) - a - rho * Z(j) + zstep / 2) / sigma)
            elseif (k == N) then
                Zprob(j,k) = 1 - cdf_normal((Z(N) - a - rho * Z(j) - zstep / 2) / sigma)
            else
                Zprob(j,k) = cdf_normal((Z(k) - a - rho * Z(j) + zstep / 2) / sigma) - &
                    cdf_normal((Z(k) - a - rho * Z(j) - zstep / 2) / sigma);
            end if
        end do
    end do

    end subroutine

    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    !!Normal CDF
    function cdf_normal(x)
    REAL(KIND=rk), INTENT(in)  :: X
    REAL(KIND=rk) :: cdf_normal

    cdf_normal =  0.5d0*( 1+erf(x/SQRT(2.d0)) )

    end function

    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    !!Solve value function by backward induction
    subroutine solveValueFunction( params, grids, modelObjects, show )
    implicit none

#ifdef mpi
    include 'mpif.h'
#endif     
    !inputs
    type (structparamstype), intent(in) :: params
    type (gridsType), intent(in) :: grids
    logical :: show

    !outputs
    type (modelObjectsType), intent(out) :: modelObjects
    !local
    integer :: ixt, ixAIME, ixA, ixY, ixL, ixA1
    real (kind=rk) :: negV, A,A1, Y, lbA1, ubA1, AIME, EV1(numPointsA,numAIME), EV1SPA(numPointsA,numAIME,numPointsSPA), Va2 ! Agrid1(numPointsA)
    real (kind=rk) :: AIME1grid(numAIME), policyA1temp, negVtemp, realisedV(numPointsY), va, check, mat(numPointsSPA,numPointsL,numPointsA), checkA(numPointsA),  check2, checkA2(numPointsA)
    real (kind=rk) :: temp(numPointsA, numAIME, numPointsY),  tempA(numPointsA, numAIME, numPointsY,numPointsA)
    integer:: indexBigN(2), indexSmalln(2), singleIndex, typeSim, AssetChoice(numPointsType, numPointsA, numAIME, numPointsY), lchoice(numPointsType, numPointsA, numAIME, numPointsY)
    integer :: numTasks, tasksPerCore, leftOverTasks, thisCoreStart, thisCoreEnd, requiredl, spa, minSPA, ixl2, ixA12
    character(len=1024) :: outFile
#ifdef mpi
    integer:: thisCoreStartStore, thisCoreEndStore, count

    real (kind=rk), allocatable :: VecV(:), Vecpolicy(:), VecEV(:)
    real (kind=rk) :: tempV(mpiDim*numPointsY*numPointsSPA), temppolicy(mpiDim*numPointsSPA*numPointsY*numPointsL*numPointsA), tempEV(mpiDim*numPointsY*numPointsSPA)
    real (kind=rk), allocatable :: LocV(:), Locpolicy(:), LocEV(:)
    integer :: locVecSize, otherDimP, testRank, thisCoreStartTest, thisCoreEndTest, start, finish, age, otherDimV

    !Test
    real (kind=rk) :: testC

    allocate(VecV(mpiDim*numPointsY*numPointsSPA), Vecpolicy(mpiDim*numPointsSPA*numPointsY*numPointsL*numPointsA), VecEV(mpiDim*numPointsY*numPointsSPA))
#endif   

    !To large to be static object
    allocate(modelObjects%EV(Tperiods+1, numPointsType, numPointsA, numAIME, numPointsSPA,numPointsY))
    allocate(modelObjects%policy(Tperiods, numPointsType, numPointsA, numAIME, numPointsSPA, numPointsY,numPointsL,numPointsA))
    allocate(modelObjects%V(Tperiods+1, numPointsType,numPointsA, numAIME, numPointsSPA,numPointsY))

    !----------------------------------------------------------------------------------!
    ! Initialize MPI
    !----------------------------------------------------------------------------------!

    indexBigN(1) = numPointsA
    indexBigN(2) = numAIME
    numTasks = mpiDim

    tasksPerCore = int(numTasks/procSize)
    leftOverTasks = numTasks - tasksPerCore*procSize
    if (leftOverTasks.gt.0) tasksPerCore = tasksPerCore + 1
    thisCoreStart = rank*tasksPerCore + 1
    thisCoreEnd = min((rank+1)*tasksPerCore, numTasks)
#ifdef mpi
    call mpi_barrier(mpi_comm_world, ierror)
    if (ierror.ne.0) stop 'mpi problem180'
#endif 

    !Set the terminal value function and expected value function to 0
    modelObjects%EV(Tperiods + 1,:,:,:,:,:)  = 0          ! continuation value at T-1
    modelObjects%V(Tperiods + 1,:,:,:,:,:) = 0

    !Initialise everything
    modelObjects%EV(:,:, :,:,:,:)  = 0          ! continuation value at T-1
    modelObjects%V(:,:,:,:,:,:) = 0
    modelObjects%policy(:,:,:,:,:,:,:,:) = 0.0

#ifdef debugMPI
    if (rank==0) pause
#endif

    do typeSim = 1, numPointsType
        if (show .AND. rank==0) WRITE(*,*)  'Solving for type ', typeSim, ' of ',  numPointsType
        do ixt=Tperiods,1, -1                               ! Loop from time T-1 to 1
            AIME1grid = grids%AIMEgrid(typeSim,ixt + 1, :)
            do ixAIME = 1, numAIME
                do ixA = 1, numPointsA                   ! points on asset grid

                    indexSmalln(1) = ixa
                    indexSmalln(2) = ixAIME
                    !indexSmalln(3) = spa WHAT TO DO WITH MPI????
                    singleIndex = sumindex(indexBigN, indexSmalln, .FALSE.)

                    if ((singleIndex .ge. thisCoreStart) .and. (singleIndex .le. thisCoreEnd)) then
                        if (ixt < stopwrok) then
                            !Although doesn't recieve income still need to loop around
                            !hypothetical income because participation doesn't effect
                            !earning potential
                            ! STEP 1. solve problem at grid points in assets, income + labour choices
                            ! ---------------------------------------------------------
                            do ixY = 1, numPointsY               ! points on income grid
                                if (ixt >= TendRI) then
                                    EV1  = modelObjects%EV(ixt + 1,typeSim,:,:,1,ixY)  ! relevant section of EV matrix (in assets tomorrow)
                                    !include arbirtary SPA and copy results as you are definitely above it by this point and so don't care
                                    call solvePeriodRE(params, grids, ixy, ixA, ixAIME ,ixt,typeSim, Tretire-1+TrueSPA, EV1, ixA1,ixl,Va)

                                    modelObjects%policy(ixt,typeSim,ixA,ixAIME, :,ixY,ixl,ixA1) = 1.0
                                    if (ixt==1) then
                                        AssetChoice(typeSim,ixA,ixAIME,ixY) = ixA1
                                        lChoice(typeSim,ixA,ixAIME,ixY) = ixl
                                    end if 
                                    !!modelObjects%policyC(ixt,typeSim, ixA, ixAIME,:, ixY) = modelObjects%policyC(ixt,typeSim, ixA, ixAIME,1, ixY)
                                    !modelObjects%policyL(ixt,typeSim,ixA,ixAIME,:,ixY,ixl) = 1.0
                                    modelObjects%V(ixt,typeSim, ixA, ixAIME,:,ixY) = Va
                                elseif (ixt == TendRI-1) then
                                    EV1  = modelObjects%EV(ixt + 1,typeSim,:,:,1,ixY)  ! relevant section of EV matrix (in assets tomorrow)
                                    !If you SPA is max possible in the period before you know the SPA but haven't recieved it yet so it is ratonal
                                    !but you do care that you haven't yet received the SPA
                                    call solvePeriodRE(params, grids, ixy, ixA, ixAIME ,ixt,typeSim, TendRI, EV1, ixA1,ixl,Va)
                                    modelObjects%policy(ixt,typeSim,ixA,ixAIME, numPointsSPA,ixY,ixl,ixA1) = 1.0
                                    !modelObjects%policyL(ixt,typeSim,ixA,ixAIME,numPointsSPA,ixY,ixl) = 1.0
                                    modelObjects%V(ixt,typeSim, ixA, ixAIME,numPointsSPA,ixY) = Va
                                    !All other possible SPA are equivalent as you have recieved the SPA and so know what it is and don't care
                                    call solvePeriodRE(params, grids, ixy, ixA, ixAIME ,ixt,typeSim, TendRI-1, EV1, ixA1,ixl,Va)
                                    modelObjects%policy(ixt,typeSim,ixA,ixAIME, (/(spa,spa=1,numPointsSPA-1)/),ixY,ixl,ixA1) = 1.0
                                    !modelObjects%policyL(ixt,typeSim,ixA,ixAIME,(/(spa,spa=1,numPointsSPA-1)/),ixY,ixl) = 1.0
                                    modelObjects%V(ixt,typeSim, ixA, ixAIME,(/(spa,spa=1,numPointsSPA-1)/),ixY) = Va

                                else
                                    EV1SPA  = modelObjects%EV(ixt + 1,typeSim,:,:,:,ixY)  ! relevant section of EV matrix (in assets tomorrow)
                                    !Set SPA to be TendRI as just need something greater than current age
                                    minSPA = maxval((/ixt-Tretire+2,1/))
                                    mat=0.0
                                    call solvePeriodRI(params, grids, ixy, ixA, ixAIME ,ixt,typeSim, TendRI, EV1SPA, &
                                        mat(minSPA:numPointsSPA,:,:) ,modelObjects%V(ixt,typeSim, ixA, ixAIME,minSPA:numPointsSPA,ixY) )
                                    modelObjects%policy(ixt,typeSim,ixA,ixAIME, :,ixY,:,:) = mat

                                    call solvePeriodRE(params, grids, ixy, ixA, ixAIME ,ixt,typeSim, TendRI, EV1SPA(:,:,numPointsSPA-1), ixA12,ixl2,Va2)
                                    if (ixy==5 .AND. ixa ==6 .AND. ixAIME==1) then
                                        call solvePeriodRE(params, grids, ixy, ixA, ixAIME ,ixt,typeSim, TendRI, EV1SPA(:,:,numPointsSPA), ixA1,ixl,Va)
                                        check = sum(mat(numPointsSPA,2,:))!sum(modelObjects%policy(ixt,typeSim,ixA,ixAIME, 11,ixY,:,:))
                                        write (*,*) check, ixl
                                    end if
                                    !call solvePeriodRE(params, grids, ixy, ixA, ixAIME ,ixt,typeSim, TendRI, EV1SPA(:,:,numPointsSPA), ixA1,ixl,Va)
                                    !check = sum(mat(numPointsSPA,2,:))!sum(modelObjects%policy(ixt,typeSim,ixA,ixAIME, 11,ixY,:,:))
                                    !checkA = sum(mat(numPointsSPA,:,:),1)
                                    !!check2 = sum(mat(numPointsSPA-1,2,:))!sum(modelObjects%policy(ixt,typeSim,ixA,ixAIME, 11,ixY,:,:))
                                    !!checkA2 = sum(mat(numPointsSPA-1,:,:),1)
                                    if (ixt>=Tretire) then
                                        !solve once with SPA = current period
                                        !I know SPA and recieve it so continaution value is the same for all possible spa < ixt
                                        EV1  = modelObjects%EV(ixt + 1,typeSim,:,:,1,ixY)! relevant section of EV matrix (in assets tomorrow)
                                        call solvePeriodRE(params, grids, ixy, ixA, ixAIME ,ixt,typeSim, ixt, EV1, ixA1,ixl,Va)
                                        modelObjects%policy(ixt,typeSim,ixA,ixAIME, (/(spa,spa=1,ixt-Tretire+1)/),ixY,ixl,ixA1) = 1.0
                                        !modelObjects%policyL(ixt,typeSim,ixA,ixAIME,(/(spa,spa=Tretire,ixt-Tretire+1)/),ixY,ixl) = 1.0
                                        modelObjects%V(ixt,typeSim, ixA, ixAIME,(/(spa,spa=1,ixt-Tretire+1)/),ixY) = Va
                                    end if
                                end if

                            end do
                        else
                            !Can't work at this stage and you are past any possible SPA so pass in arbitrary value for ixy (= 1) and SPA9=60) and
                            !hardcode ixl=0(not work)
                            call maxutility(params, grids,ixt,typeSim,ixa,ixAIME,0,  1,Tretire,  modelObjects%EV(ixt + 1,typeSim,:,:,1,1),ixa1,va)

                            !Copy to all values of SPA and ixy as don't affect decsion at these ages

                            modelObjects%policy(ixt,typeSim,ixA,ixAIME,:,:,1,ixA1)=1.0
                            modelObjects%V(ixt,typeSim, ixA,ixAIME,:,:)       = va
                        end if

                        ! STEP 2. integrate out income today conditional on income
                        ! yesterday to get EV and EdU
                        ! --------------------------------------------------------
                        if (ixt >= TendRI) then
                            realisedV(:) = modelObjects%V(ixt,typeSim, ixA,ixAIME,1,:)
                            do ixY = 1,numPointsY,1
                                modelObjects%EV(ixt,typeSim, ixA,ixAIME,1,ixY)  = dot_product(grids%incTransitionMrx(ixY,:),realisedV)
                                modelObjects%EV(ixt,typeSim, ixA,ixAIME,:,ixY) = modelObjects%EV(ixt,typeSim, ixA,ixAIME,1,ixY)
                            end do !ixY
                        else
                            do spa = 1,numPointsSPA
                                realisedV(:) = modelObjects%V(ixt,typeSim, ixA,ixAIME,spa,:)
                                do ixY = 1,numPointsY,1
                                    modelObjects%EV(ixt,typeSim, ixA,ixAIME,spa,ixY)  = dot_product(grids%incTransitionMrx(ixY,:),realisedV)
                                end do !ixY
                            end do
                        end if
                    end if
                end do
            end do

#ifdef mpi 
            call mpi_barrier(mpi_comm_world, ierror)
            if (ierror.ne.0) stop 'mpi problem180'

            locVecSize = thisCoreEnd - thisCoreStart + 1
            otherDimP = numPointsSPA*numPointsY*numPointsL*numPointsA
            otherDimV = numPointsSPA*numPointsY
            allocate(Locpolicy(locVecSize*otherDimP),LocEV(locVecSize*otherDimV),locV(locVecSize*otherDimV))

            call unpackArrays(modelObjects%policy(ixt,typeSim, :, :,:, :,:,:), modelObjects%V(ixt,typeSim, :, :, :,:), modelObjects%EV(ixt,typeSim, :, :,:, :), &
                locpolicy, locV, locEV, mpiDim,otherDimP,otherDimV,thisCoreStart,thisCoreEnd)

            call mpi_allgather(LocV, locVecSize*otherDimV, mpi_double_precision, VecV, locVecSize*otherDimV, mpi_double_precision, mpi_comm_world, ierror)

            call mpi_allgather(Locpolicy, locVecSize*otherDimP, mpi_double_precision, Vecpolicy, locVecSize*otherDimP, mpi_double_precision, mpi_comm_world, ierror)

            call mpi_allgather(LocEV(1), locVecSize*otherDimV, mpi_double_precision, VecEV(1), locVecSize*otherDimV, mpi_double_precision, mpi_comm_world, ierror)

            deallocate(LocV,Locpolicy,LocEV)
            tempV = VecV
            temppolicy = Vecpolicy
            tempEV = VecEV
            do testRank = 0, (procSize-1)
                thisCoreStartTest = testrank*tasksPerCore + 1
                thisCoreEndTest = min((testrank+1)*tasksPerCore, numTasks)
                locVecSize = thisCoreEndTest - thisCoreStartTest + 1
                allocate(LocV(locVecSize*otherDimV),Locpolicy(locVecSize*otherDimP),LocEV(locVecSize*otherDimV))

                start = (thisCoreStartTest-1)*otherDimP+1
                finish = thisCoreEndTest*otherDimP

                locpolicy = temppolicy(start:finish)

                !Distribute the contigous section between the contigous section in each column corresponding to the non-splitting dimesnions of the array
                do count=1,otherDimP
                    start = (count-1)*mpiDim+thisCoreStartTest
                    finish = (count-1)*mpiDim+thisCoreEndTest
                    Vecpolicy(start:finish) = locpolicy((count-1)*locVecSize+1:(count-1)*locVecSize+locVecSize)
                end do

                start = (thisCoreStartTest-1)*otherDimV+1
                finish = thisCoreEndTest*otherDimV

                locV= tempV(start:finish)
                locEV = tempEV(start:finish)

                !Distribute the contigous section between the contigous section in each column corresponding to the non-splitting dimesnions of the array
                do count=1,otherDimV
                    start = (count-1)*mpiDim+thisCoreStartTest
                    finish = (count-1)*mpiDim+thisCoreEndTest
                    VecV(start:finish) = locV((count-1)*locVecSize+1:(count-1)*locVecSize+locVecSize)
                    VecEV(start:finish) = locEV((count-1)*locVecSize+1:(count-1)*locVecSize+locVecSize)
                end do

                deallocate(LocV,Locpolicy,LocEV)
            end do!
            modelObjects%V(ixt,typeSim, :, :, :,:) = reshape(VecV, (/numPointsA, numAIME, numPointsSPA,numPointsY/))
            modelObjects%policy(ixt,typeSim, :, :, :,:,:,:) = reshape(Vecpolicy, (/numPointsA, numAIME,numPointsSPA, numPointsY,numPointsL,numPointsA/))
            modelObjects%EV(ixt,typeSim, :, :, :,:) = reshape(VecEV, (/numPointsA, numAIME, numPointsSPA,numPointsY/))

            call mpi_barrier(mpi_comm_world, ierror)
            if (ierror.ne.0) stop 'mpi problem180'
#endif  

            if (show .AND. rank==0) WRITE(*,*)  'Passed period ', ixt, ' of ',  Tperiods
            !if (show .AND. rank==0) WRITE(*,*)  'Average Value at ', ixt, ' is ',   params%lambda*sum(modelObjects%V(ixt,typeSim,:,:, TrueSPA,:))/(size(modelObjects%V(ixt,typeSim,:,:, TrueSPA,:)))
            !if (show .AND. rank==0 .AND. ixt<= spouseretire) WRITE(*,*)  'Sum Policy L', sum(policyL(ixt,typeSim, :, :, :))
            if (rank==0 .AND. intermediateToFile .AND. ixt <10) then
                !!L
                write (outFile, *), trim("..\\temp\polL"),ixt,".txt"
                !check = sum(mat(numPointsSPA,2,:))!sum(modelObjects%policy(ixt,typeSim,ixA,ixAIME, 11,ixY,:,:))
                outfile=ADJUSTL(outfile)
                temp =  sum(modelObjects%policy(ixt,typeSim,:,:, TrueSPA,:,2,:),4 )
                inquire (iolength=requiredl) temp
                open (unit=201, file=outfile, status='unknown',recl=requiredl, action='write')
                write (201, *) temp
                close( unit=201)

                write (outFile, *), trim("..\\temp\polA"),ixt,".txt"
                !check = sum(mat(numPointsSPA,2,:))!sum(modelObjects%policy(ixt,typeSim,ixA,ixAIME, 11,ixY,:,:))
                outfile=ADJUSTL(outfile)
                tempA =  sum(modelObjects%policy(ixt,typeSim,:,:, TrueSPA,:,:,:),4 )
                inquire (iolength=requiredl) tempA
                open (unit=201, file=outfile, status='unknown',recl=requiredl, action='write')
                write (201, *) tempA
                close( unit=201)

            end if
        end do !ixt

    end do !type

    end subroutine

    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    !!Objective function
    function objectivefunc(params, grids, A1, A0, Y,L,ixP,ixType, AIME, EV1)
    implicit none
    !inputs
    type (structparamstype), intent(in) :: params
    type (gridsType), intent(in) :: grids
    real (kind=rk), intent(in) :: A1, A0, Y, AIME, EV1(:)
    integer, intent(in) :: ixP, L, ixType
    !real (kind=rk), intent(in), optional :: lambda
    !ouptut
    real (kind=rk) :: objectivefunc
    !local
    real (kind=rk) :: cons, VA1, VB1 !, locLambda

    !Get tomorrow's consumption (cons), the value of left over assets (VA1) and
    !total value (u(c) + b * VA1
    cons = A0  + Y - (A1)/(1+params%r)

    !!THIS IS WRONG SHOULD BE NEXT PERIODS AIME NOT THIS PERIODS
    call  linearinterp1(grids%AIMEgrid(ixType,ixP + 1, :),EV1, numAIME, AIME, VA1, 1, 1 )
    VB1 = params%thetab*((A1+params%K)**(1-params%gamma))/(1-params%gamma)
    !if (present(lambda)) then
    !    locLambda = lambda
    !else
    !    locLambda = 1.0
    !end if
    objectivefunc = utility(params,cons,L)/params%lambda + params%beta * ((1- grids%mortal(ixP))* VA1+grids%mortal(ixP)*VB1)

    !! ------------------------------------------------------------------------
    !The optimisation routine that we will use searches for the minimum of the
    !function. We want the maximum. So we multiply out function here by -1 so
    !that the optimiser will fill the minimum of the negative of our function,
    !i.e. the maximum of our functino

    !objectivefunc = - objectivefunc;

    end function

    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    !!Utility Function
    function utility(params,cons,L)
    implicit none
    !This function takes consumption as an argument and returns utility. The
    !utility functin is CRRA except that we add a very small number (eps) to
    !consumption so that the computer can deal wi
    !inputs
    type (structparamstype), intent(in) :: params
    integer, intent(in) :: L
    real (kind=rk), intent(in) :: cons

    !outpus
    real (kind=rk) :: utility, les

    if (cons<=0) then
        print *, 'Error in utility! Consumption is LE 0'
        stop
    end if
    !10/112 comuniting time -10/112
    les=(L)*(1-params%hrsWrk )+(1-L);
    if (params%gamma == 1) then
        utility = log(cons**params%nu*les**(1-params%nu));
    else
        utility= ((cons**params%nu*les**(1-params%nu))**(1-params%gamma)  )/(1-params%gamma);
    end if

    end function

    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    !!Defined Benefite Pension function
    function dbPension(params,AIME)
    implicit none
    !calibrated to give 25,000 max from 60,000 AIME upwards
    !db2 = level/(point^2-2point), db1 = 6.9447e-06*points
    !point

    !inputs
    type (structparamstype), intent(in) :: params
    real (kind=rk), intent(in) :: AIME
    !outputs
    real (kind=rk) :: dbPension
    !local
    real (kind=rk):: bound
    bound = -params%db(1) /(2.0*params%db(2))
    if (AIME < bound) then
        dbPension = params%db(2)*AIME**2 + params%db(1)*AIME;
    else
        dbPension = params%db(2)*bound**2 + params%db(1)*bound !25000
    endif
    end function

    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    !!Simulation Subroutine
    subroutine simWithUncer(params, grids, modelObjects,  yex, c, a, v, l, y, AIME )
    use Header
    implicit none

    !inputs
    type (structparamstype), intent(in) :: params
    type (gridsType), intent(in) :: grids
    type (modelObjectsType), intent(inout) :: modelObjects
    !real (kind=rk), intent(in) :: policyA1(Tperiods,  numPointsType, numPointsA, numAIME, numPointsY)
    !integer, intent(in) :: policyL(Tperiods,  numPointsType, numPointsA, numAIME, numPointsY)
    !real (kind=rk), intent(in) :: EV(Tperiods+1,  numPointsType, numPointsA, numAIME, numPointsY)
    real (kind=rk), intent(in) :: yex(Tperiods, numSims)

    !outputs
    real (kind=rk), intent(out) :: y(Tperiods, numSims) !income
    real (kind=rk), intent(out) :: c(Tperiods, numSims)  !consumption
    integer, intent(out) :: l(Tperiods, numSims) !labour supply
    real (kind=rk), intent(out) :: v(Tperiods, numSims)  !value
    integer, intent(out) :: a(Tperiods + 1,numSims) !this is the path at the start of each period, so we include the 'start' of death
    real (kind=rk), intent(out) :: AIME(Tperiods + 1,numSims)

    !local
    real (kind=rk) :: startingA(numSims), startAIME, Aval,AIMEval, check, check1, check2, checkMax
    integer :: s, t, idxA(1), idxAIME(1), idxY(1), workAge
    integer :: seedIn, Lcube(8), typeSim
    real (kind=rk) ::  uniformRand(numSims), ltemp, lbA1, EV1(numPointsA ,numAIME)
    INTEGER :: n, i, uniformInt(numSims), j
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
    integer :: unemployed !1 employed, 2 unemploeyed
    type (modelObjectsType):: LocmodelObjects(2)
    !real (kind=rk):: LocpolicyA1(Tperiods,  numPointsType, numPointsA, numAIME, numPointsSPA,numPointsProd,2)
    !integer :: LocpolicyL(Tperiods,  numPointsType, numPointsA, numAIME, numPointsSPA,numPointsProd,2)
    real (kind=rk) :: LocEV(Tperiods+1,  numPointsType, numPointsA, numAIME, numPointsSPA,numPointsProd,2);
    real (kind=rk) :: gridY(numPointsProd), Yval, logitShock(numSims,Tperiods), ccp(numPointsL*numPointsA), policy(numPointsA, numAIME,numPointsL*numPointsA), ccpTemp(numPointsL*numPointsA)
    real (kind=rk) :: mat(numPointsL,numPointsA), cumsum
    integer :: idxa1
    
    !deallocate(modelObjects%v)
    LocmodelObjects(1)%policy =  modelObjects%policy(:,:,:,:,:,(/(i,i=1,numPointsY-1,2)/),:,:)
    LocmodelObjects(2)%policy =  modelObjects%policy(:,:,:,:,:,(/(i,i=2,numPointsY,2)/),:,:)
    !deallocate(modelObjects%policy)
    LocmodelObjects(1)%EV =  modelObjects%EV(:,:,:,:,:,(/(i,i=1,numPointsY-1,2)/))
    LocmodelObjects(2)%EV =  modelObjects%EV(:,:,:,:,:,(/(i,i=2,numPointsY,2)/))
    !deallocate(modelObjects%ev)
    
    !LocmodelObjects(1)%policy =>  modelObjects%policy(:,:,:,:,:,(/(i,i=1,numPointsY-1,2)/),:,:)
    !LocmodelObjects(2)%policy =>  modelObjects%policy(:,:,:,:,:,(/(i,i=2,numPointsY,2)/),:,:)
    !LocmodelObjects(1)%EV =>  modelObjects%EV(:,:,:,:,:,(/(i,i=1,numPointsY-1,2)/))
    !LocmodelObjects(2)%EV =>  modelObjects%EV(:,:,:,:,:,(/(i,i=2,numPointsY,2)/))
    
    !call random_seed(size = n)
    !allocate(seed(n))
    !seed(n) = 16101988
    !call random_seed(put = seed)
    !!get uniform random number
    call random_number(uniformrand)
    do i = 1, numsims
        startinga(i) = grids%initialassets(i) !uniformrand(i)*numsims)+1
        if (startinga(i) < 0.0) startinga(i) = 0.0
    end do
    startaime = 15500

    call RANDOM_NUMBER(logitShock)
    checkMax = 0.0
    do s = 1,numsims,1
        if (numpointstype == 1) then
            typesim = 1
        else
            if (real(s,rk)/numsims <0.16) then
                typesim = 1
            else if (real(s,rk)/numsims <0.48) then
                typesim = 2
            else if (real(s,rk)/numsims <0.59) then
                typesim =3
            else
                typesim =4
            end if
        end if
        !spa = 60 + mod(s,6)

        a(1, s) = minloc(abs(startinga(s)-grids%Agrid(typesim,1,:)),1)
        aime(1,s)=startaime
        do t = 1,tperiods,1                              ! loop through time periods for a particular individual
            if (yex(t,s) <= 0.0001 ) then
                unemployed = 2
#if PROD_SIZE == 5
                !if (numpointsprod==5) then
                gridy = (/1,2,3,4,5/)
#elif  PROD_SIZE == 10                   
                !else if (numpointsprod==10) then
                gridy = (/1,2,3,4,5,6,7,8,9,10/)
                !end if
#endif
                yval = idxy(1)
            else
                unemployed = 1
                !everyone start employed so this is ok
                idxy=minloc(abs(yex(t,s)-grids%ygridextrap(typesim,t,:)))
                gridy =grids%ygrid(typesim,t,(/(i,i=1,numpointsy-1,2)/))
                yval = yex(t,s)
                y(t, s) =  yex(t,s)
            end if

            idxaime=minloc(abs(aime(t, s)-grids%aimegrid(typesim,t,:)))
            !check = sum(modelObjects%policy(t,typeSim,a(1, s),idxaime(1), 1,idxy(1),:,:))
            !=locpolicyl(t,typesim,a(1, s),idxaime(1),spa,idxy(1),unemployed)
            do i=1,numPointsL
                do j=1,numPointsA
                    policy(:,:,j+(i-1)*numPointsA) = LocmodelObjects(unemployed)%policy(t,typesim,:,:,TrueSPA,idxy(1),i,j)
                end do
            end do
            ccp =  reshape(LocmodelObjects(unemployed)%policy(t,typesim,a(1, s),idxaime(1),TrueSPA,idxy(1),:,:),(/numPointsL*numPointsA/))
            check1 = abs(sum(ccp))
            y(t, s) = yex(t,s)
            !if (maxval(lcube)==minval(lcube)) then
            call linearinterp2_vec(grids%aimegrid(typesim,t,:),gridy, numaime, numpointsprod,(numPointsL*numPointsA), aime(t, s), yval, ccpTemp, policy)  !,real(locpolicyl(t, typesim, :, :,spa,:,unemployed),rk)
            ccpTemp = ccpTemp/sum(ccpTemp)
            check2 = abs(sum(ccpTemp))
            check = abs(sum(ccp - ccpTemp))
            if (check > checkMax) checkMax = check
            if (abs(sum(ccp - ccpTemp)) > 0.00001) then
                write (*,*) "Larege Error"
                !lba1 = grids%agrid(typesim,t + 1, 1);          ! lower bound: assets tomorrow
                !ev1 = modelobjects%ev(t + 1,typesim,:,:,spa,idxy(1))
                !!we can only end up here for ages below 80 as for >80 all entries in poll are 0
                !!this also takes account of uneployment shock as the wage offer they have is 0 so won't work
                !call solveperiod(params, grids, y(t, s), a(t, s), aime(t, s) ,t, typesim, lba1, ev1, &
                !    grids%benefit(t),spa,a(t+1, s),c(t, s),l(t,s),v(t  , s))
            else
                cumsum = 0.0
                do i = 1,size(ccp)
                    cumsum=cumsum+ccp(i)
                    if (cumsum > logitShock(s,t)) exit
                end do
                l(t,s) = 1- mod(i,2)
                idxa1= (i-1)/2+1
                a(t+1, s) = idxa1
                mat = LocmodelObjects(unemployed)%policy(t,typesim,a(1, s),idxaime(1),TrueSPA,idxy(1),:,:)
                !LocmodelObjects(unemployed)%policy(t,typesim,a(1, s),idxaime(1),TrueSPA,idxy(1),:,:)
                !call linearinterp3_withextrap(grids%agrid(typesim,t,:), grids%aimegrid(typesim,t,:), gridy, &
                !    numpointsa, numaime, numpointsprod, a(t, s), aime(t, s), yval, a(t+1, s),  locpolicya1(t,typesim, :, :,spa,:,unemployed))
                !call linearinterp3_withextrap(grids%agrid(typesim,t,:), grids%aimegrid(typesim,t,:), gridy,&
                !    numpointsa, numaime, numpointsprod+1,  a(t, s), aime(t, s), yval, v(t  , s),  locev(t,typesim, :, :,spa,:,unemployed))
            end if

            if (l(t,s) .eq. 0) then
                y(t, s)=grids%benefit(t)
            end if

            !    workage = startage - 20 + t
            !    if  (t < tretire) then ! spouseretire
            !        if (l(t,s) .eq. 0) then
            !            aime(t+1, s) =   aime(t, s) * (workage-1)/workage
            !        else
            !            aime(t+1, s) =  y(t,s)/workage + aime(t, s) * (workage-1)/workage
            !        end if
            !    else
            !        aime(t+1, s) =   aime(t, s)
            !    end if
            !
            !    call gross2net(params,y(t, s),t,l(t,s),typesim, aime(t, s),spa)
            !    ! tomorrow's optimal assets
            !    !! this is wrong because above is wrong!
            !    c(t, s) = a(t, s)  + y(t, s) - (a(t+1, s)/(1+params%r))
            !
        end do   !t
    end do! s
    if (checkMax > 0.00001) write (*,*), "Larege Error = ", checkMax
    !write (*,*) "Starting A ", sum(startinga)/size(startinga)
    !write (*,*) "Sim Start A ", sum(grids%Agrid(1,1,a(1,:)))/real(numSims,rk)
    end subroutine
    ! ---------------------------------------------------------------------------------------------------------!
    !!truncate obsrvtion
    !---------------------------------------------------------------------------------------------------------!
    function truncate(y, negtrunc, postrunc)
    implicit none
    !input
    real (kind=rk), intent(in) :: y, negtrunc, postrunc

    !output
    real (kind=rk) :: truncate

    ! Truncate input if its too big or too small
    if (y < negtrunc) then ! If y is less than the value negtrunc
        truncate = negtrunc;
    elseif (y > postrunc) then
        truncate = postrunc; ! If y is greater than the value postrunc
    else
        truncate = y;
    end if
    end function
    ! ---------------------------------------------------------------------------------------------------------!
    !!Convert Groos to net income
    !---------------------------------------------------------------------------------------------------------!
    subroutine gross2net(params,grids,Y,ixt,ixl,ixType, AIME, spa)
    implicit none

    !inputs
    type (structparamstype), intent(in) :: params
    type (gridsType), intent(in) :: grids

    integer, intent(in) :: ixt,ixl, ixType, spa
    real (kind=rk), intent(in) :: AIME

    !changing
    real (kind =rk), intent(inout) :: Y

    !!You can't work over certain age so whatever input income most be benefits until add on pension
    !if (ixt >= stopwrok) then
    !    Y    = grids%benefit(ixt)
    !end if

    !This previously assumed spouseretire > spa
    !Add own and husbands pension

    if (ixt >= spa) then
        Y = Y + params%pension
        if  (ixt >= spouseretire) then
            Y = Y + params%pension
            if (ixL==0 .AND. mod(ixType,2)==0) then !
                Y =   Y + dbPension(params,AIME)
            end if
        else
            Y = Y +params%spouseInc(ixType)
        end if
    else
        if  (ixt < spouseretire) then
            Y = Y +params%spouseInc(ixType)
        else
            Y = Y + params%pension
        end if
    end if

    end subroutine

    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    !!Returns GMM Cretieria
    function gmm(params,grids,target,weights)
    implicit none

#ifdef mpi
    include 'mpif.h'
#endif  

    !inputs
    type (structparamstype), intent(in) :: params
    type (gridsType), intent(inout) :: grids
    real (kind=rk), intent(in) :: target(:,:), weights(:,:)

    !output
    real (kind=rk) :: gmm

    !local
    type (modelObjectsType) :: modelObjects
    !real (kind=rk) :: V(numPointsType,Tperiods+1, numPointsA, numPointsY, numAIME)
    !real (kind=rk) :: policyA1(numPointsType,Tperiods, numPointsA, numPointsY, numAIME)
    !real (kind=rk) :: policyC(numPointsType,Tperiods, numPointsA, numPointsY, numAIME)
    !integer :: policyL(numPointsType,Tperiods, numPointsA, numPointsY, numAIME)
    !real (kind=rk) :: EV(numPointsType,Tperiods+1, numPointsA, numPointsY, numAIME)
    real (kind=rk) :: ypath(Tperiods, numSims) !income
    real (kind=rk) :: cpath(Tperiods, numSims)  !consumption
    integer :: lpath(Tperiods, numSims) !labour supply
    real (kind=rk) :: vpath(Tperiods, numSims)  !value
    integer :: apath(Tperiods + 1,numSims) !this is the path at the start of each period, so we include the 'start' of death
    real (kind=rk) ::  meanL(Tperiods), meanA(Tperiods)
    real (kind=rk) :: AIME(Tperiods + 1,numSims)

    integer :: n, typeSim
    integer :: i, k, start, finish
    integer, allocatable:: rangeSims(:)
    !Set asset grid
    do typeSim = 1, numPointsType
        call getassetgrid( params, grids%maxInc(typeSim,:), grids%Agrid(typeSim,:,:))
    end do
    !solve
    call solveValueFunction( params, grids, modelObjects, .FALSE. )
    if (rank==0) then
        !simulate
        call simWithUncer(params, grids, modelObjects, grids%Simy, cpath, apath, vpath, lpath, ypath, AIME )
        do n=1,Tperiods
            meanL(n)=sum(real(lpath(n,:),rk))/real(numSims,rk) !size(lpath(n,:))
            !meanA(n)=sum(real(Apath(n,:),rk))/real(numSims,rk) !size(lpath(n,:))
            meanA(n)=0.0
            do k=1,numPointsType
                if (k == 1)  then
                    start = 1
                    if (numPointsType>1) then
                        finish = 0.16*numSims
                    else
                        finish = numSims
                    end if
                    if (n>1) then
                        deallocate(rangeSims)
                    end if
                    allocate(rangeSims(finish))
                    rangeSims = (/(i,i=1,finish)/)
                else if (k==2) then
                    deallocate(rangeSims)
                    start = finish +1
                    finish = 0.48*numSims
                    allocate(rangeSims(finish-start+1))
                    rangeSims = (/(i,i=start,finish)/)
                    !real(s,rk)/numSims <0.48
                else if (k==3) then
                    deallocate(rangeSims)
                    start = finish +1
                    finish = 0.59*numSims
                    allocate(rangeSims(finish-start+1))
                    rangeSims = (/(i,i=start,finish)/)
                    !real(s,rk)/numSims <0.59
                else
                    deallocate(rangeSims)
                    start = finish +1
                    finish = numSims
                    allocate(rangeSims(finish-start+1))
                    rangeSims = (/(i,i=start,finish)/)
                end if
                meanA(n)=meanA(n)+sum(real(grids%Agrid(k,n,apath(n,rangeSims)),rk))
            end do
            meanA(n)= meanA(n)/real(numSims,rk)
        end do

        if (fullLifeCycle) then
            gmm = dot_product(abs(meanL(32:32+23)-target(1,:)),abs(meanL(32:32+23)-target(1,:)))! + &
                !dot_product(abs(meanA(32:32+23)-target(2,:)),abs(meanA(32:32+23)-target(2,:)))
            else
        gmm = dot_product(abs(meanL(33 - (StartAge-20):33 - (StartAge-20)+23)-target(1,:)),weights(1,:)*abs(meanL(33 - (StartAge-20):33 - (StartAge-20)+23)-target(1,:))) + &
            dot_product(abs(meanA(33 - (StartAge-20):33 - (StartAge-20)+23)-target(2,:)),weights(2,:)*abs(meanA(33 - (StartAge-20):33 - (StartAge-20)+23)-target(2,:)))
        end if
    end if

#ifdef mpi 
    call MPI_Bcast( gmm,   1,   mpi_double_precision, 0, mpi_comm_world, ierror)
    call mpi_barrier(mpi_comm_world, ierror)
    if (ierror.ne.0) stop 'mpi problem180'
#endif

    end function
    ! ---------------------------------------------------------------------------------------------------------!
    ! !Write to file
    !---------------------------------------------------------------------------------------------------------!
    subroutine writetofile(grids, ypath, cpath, apath, vpath, lpath, yemp, AIME)
    implicit none

    !inputs
    type (gridsType), intent(inout) :: grids
    !type (structparamstype), intent(in) :: params
    real (kind=rk), intent(in) :: ypath(Tperiods, numSims) !income
    real (kind=rk), intent(in) :: cpath(Tperiods, numSims)  !consumption
    integer, intent(in) :: lpath(Tperiods, numSims) !labour supply
    real (kind=rk), intent(in) :: vpath(Tperiods, numSims)  !value
    integer, intent(in) :: apath(Tperiods + 1,numSims) !this is the path at the start of each period, so we include the 'start' of death
    real (kind=rk), intent(in)  :: yemp(Tperiods, numSims)
    real (kind=rk), intent(in)  :: AIME(Tperiods + 1,numSims)

    !local
    integer :: n, requiredl , i
    real(kind=rk) :: meanL(Tperiods), meanV(Tperiods), meanA(Tperiods), meanC(Tperiods), meanY(Tperiods), medianA
    real(kind=rk) :: meanYemp(Tperiods), meanAIME(Tperiods), meanPoor(Tperiods), meanRich(Tperiods), numLC(Tperiods)
    integer :: lrich(Tperiods,numSims/2), lpoor(Tperiods,numSims/2)
    integer :: rich, poor, k, start, finish
    integer, allocatable:: rangeSims(:)

#ifdef win   


    do n=1,Tperiods
        meanL(n)=sum(real(lpath(n,:),rk))/real(numSims,rk) !size(lpath(n,:))
        !write (201, * ) meanL(n)
        meanV(n)=sum(real(vpath(n,:),rk))/real(numSims,rk)
        !write (202, * ) meanV(n)
        meanA(n)=0.0
        do k=1,numPointsType
            if (k == 1)  then
                start = 1
                if (numPointsType>1) then
                    finish = 0.16*numSims
                else
                    finish = numSims
                end if
                if (n>1) then
                    deallocate(rangeSims)
                end if
                allocate(rangeSims(finish))
                rangeSims = (/(i,i=1,finish)/)
            else if (k==2) then
                deallocate(rangeSims)
                start = finish +1
                finish = 0.48*numSims
                allocate(rangeSims(finish-start+1))
                rangeSims = (/(i,i=start,finish)/)
                !real(s,rk)/numSims <0.48
            else if (k==3) then
                deallocate(rangeSims)
                start = finish +1
                finish = 0.59*numSims
                allocate(rangeSims(finish-start+1))
                rangeSims = (/(i,i=start,finish)/)
                !real(s,rk)/numSims <0.59
            else
                deallocate(rangeSims)
                start = finish +1
                finish = numSims
                allocate(rangeSims(finish-start+1))
                rangeSims = (/(i,i=start,finish)/)
            end if
            meanA(n)=meanA(n)+sum(real(grids%Agrid(k,n,apath(n,rangeSims)),rk))

        end do
        meanA(n)= meanA(n)/real(numSims,rk)
        !write (203, * ) meanA(n)
        meanC(n)=sum(real(cpath(n,:),rk))/real(numSims,rk)
        !write (204, * ) meanC(n)
        meanY(n)=sum(real(ypath(n,:),rk))/real(numSims,rk)
        !write (205, * ) meanY(n)
        meanYemp(n)=sum(real(yemp(n,:),rk))/real(numSims,rk)
        !write (205, * ) meanY(n)
        meanAIME(n)=sum(real(AIME(n,:),rk))/real(numSims,rk)
        numLC(n) = count(apath(n,:).eq.0)
    end do
    !L
    inquire (iolength=requiredl)  meanL
    open (unit=201, file='..\\out\lpath.txt', status='unknown',recl=requiredl, action='write')
    write (201, '(6E15.3)' ) meanL
    close( unit=201)
    !V
    inquire (iolength=requiredl)  meanV
    open (unit=202, file='..\\out\Vpath.txt', status='unknown',recl=requiredl, action='write')
    write (202, '(6E15.3)' ) meanV
    close( unit=202)
    !A
    inquire (iolength=requiredl)  meanA
    open (unit=203, file='..\\out\Apath.txt', status='unknown',recl=requiredl, action='write')
    write (203, '(6E15.3)' ) meanA
    close( unit=203)
    !C
    inquire (iolength=requiredl)  meanC
    open (unit=204, file='..\\out\Cpath.txt', status='unknown',recl=requiredl, action='write')
    write (204, '(6E15.3)' ) meanC
    close( unit=204)
    !Y
    inquire (iolength=requiredl)  meanY
    open (unit=205, file='..\\out\Ypath.txt', status='unknown',recl=requiredl, action='write')
    write (205, '(6E15.3)' ) meanY
    close( unit=205)

    inquire (iolength=requiredl)  meanYemp
    open (unit=206, file='..\\out\YempPath.txt', status='unknown',recl=requiredl, action='write')
    write (206, '(6E15.3)' ) meanYemp
    close( unit=206)

    inquire (iolength=requiredl)  meanAIME
    open (unit=207, file='..\\out\AIMEPath.txt', status='unknown',recl=requiredl, action='write')
    write (207, '(6E15.3)' ) meanAIME
    close( unit=207)

    inquire (iolength=requiredl)  numLC
    open (unit=207, file='..\\out\numLC.txt', status='unknown',recl=requiredl, action='write')
    write (207, '(6E15.3)' ) numLC
    close( unit=207)

#else
    !medianA = median(apath(Tretire,:))
    rich = 0
    poor = 0
    do n=1,numSims
        if (apath(Tretire,n) <= medianA) then
            poor = poor + 1
            lpoor(:,poor) = lpath(:,n)
        else
            rich = rich + 1
            lrich(:,rich) = lpath(:,n)
        end if
    end do

    inquire (iolength=requiredl)  meanL
    open (unit=201, file='./out/lpath', status='unknown',recl=requiredl, action='write')
    write (201, * ) 'Header'

    inquire (iolength=requiredl)  meanV
    open (unit=202, file='./out/Vpath', status='unknown',recl=requiredl, action='write')
    write (202, * ) 'Header'

    inquire (iolength=requiredl)  meanA
    open (unit=203, file='./out/Apath', status='unknown',recl=requiredl, action='write')
    write (203, * ) 'Header'

    inquire (iolength=requiredl)  meanC
    open (unit=204, file='./out/Cpath', status='unknown',recl=requiredl, action='write')
    write (204, * ) 'Header'

    inquire (iolength=requiredl)  meany
    open (unit=205, file='./out/Ypath', status='unknown',recl=requiredl, action='write')
    write (205, * ) 'Header'

    inquire (iolength=requiredl)  meanYemp
    open (unit=206, file='./out/YempPath', status='unknown',recl=requiredl, action='write')
    write (206, * ) 'Header'

    inquire (iolength=requiredl)  meanAIME
    open (unit=207, file='./out/AIMEPath', status='unknown',recl=requiredl, action='write')
    write (207, * ) 'Header'

    inquire (iolength=requiredl)  meanPoor
    open (unit=208, file='./out/PoorPath', status='unknown',recl=requiredl, action='write')
    write (208, * ) 'Header'

    inquire (iolength=requiredl)  meanRich
    open (unit=209, file='./out/RichPath', status='unknown',recl=requiredl, action='write')
    write (209, * ) 'Header'

    inquire (iolength=requiredl)  meanRich
    open (unit=210, file='./out/ldata', status='unknown',recl=requiredl, action='write')


    inquire (iolength=requiredl)  meanRich
    open (unit=211, file='./out/adata', status='unknown',recl=requiredl, action='write')

    inquire (iolength=requiredl)  meanRich
    open (unit=212, file='./out/ydata', status='unknown',recl=requiredl, action='write')

    do n=1,Tperiods
        meanL(n)=sum(real(lpath(n,:),rk))/real(numSims,rk)
        write (201, * ) meanL(n)
        meanV(n)=sum(real(vpath(n,:),rk))/real(numSims,rk)
        write (202, * ) meanV(n)
        meanA(n)=sum(real(apath(n,:),rk))/real(numSims,rk)
        write (203, * ) meanA(n)
        meanC(n)=sum(real(cpath(n,:),rk))/real(numSims,rk)
        write (204, * ) meanC(n)
        meanY(n)=sum(real(ypath(n,:),rk))/real(numSims,rk)
        write (205, * ) meanY(n)
        meanYemp(n)=sum(real(yemp(n,:),rk))/real(numSims,rk)
        write (206, * ) meanYemp(n)
        meanAIME(n)=sum(real(AIME(n,:),rk))/real(numSims,rk)
        write (207, * ) meanAIME(n)
        meanRich(n)=sum(real(lrich(n,:),rk))/real(numSims/2,rk)
        write (209, * ) meanRich(n)
        meanPoor(n)=sum(real(lpoor(n,:),rk))/real(numSims/2,rk)
        write (208, * ) meanPoor(n)
    end do
    do i=1,numsims
        do n=32,55
            write (210, * ) lpath(n,i)
            write (211, * ) apath(n,i)
            write (212, * ) lpath(n,i)*ypath(n,i)
        end do
    end do

    close( unit=201)
    close( unit=202)
    close( unit=203)
    close( unit=204)
    close( unit=205)
    close( unit=206)
    close( unit=207)
    close( unit=208)
    close( unit=209)
    close( unit=210)
    close( unit=211)
    close( unit=212)
#endif

    end subroutine
    ! ---------------------------------------------------------------------------------------------------------!
    ! !Write to file
    !---------------------------------------------------------------------------------------------------------!
    subroutine writetofileByType(grids,ypath, cpath, apath, vpath, lpath, yemp, AIME)
    implicit none

    !inputs
    type (gridsType), intent(inout) :: grids
    !type (structparamstype), intent(in) :: params
    real (kind=rk), intent(in) :: ypath(Tperiods, numSims) !income
    real (kind=rk), intent(in) :: cpath(Tperiods, numSims)  !consumption
    integer, intent(in) :: lpath(Tperiods, numSims) !labour supply
    real (kind=rk), intent(in) :: vpath(Tperiods, numSims)  !value
    integer, intent(in) :: apath(Tperiods + 1,numSims) !this is the path at the start of each period, so we include the 'start' of death
    real (kind=rk), intent(in)  :: yemp(Tperiods, numSims)
    real (kind=rk), intent(in)  :: AIME(Tperiods + 1,numSims)

    !local
    integer :: n, requiredl , i, finish, start
    integer(1) :: k
    real(kind=rk) :: meanL(Tperiods), meanV(Tperiods), meanA(Tperiods), meanC(Tperiods), meanY(Tperiods), medianA
    real(kind=rk) :: meanYemp(Tperiods), meanAIME(Tperiods), meanPoor(Tperiods), meanRich(Tperiods), numLC(Tperiods)
    integer :: lrich(Tperiods,numSims/2), lpoor(Tperiods,numSims/2)
    integer :: rich, poor
    character(len=1024) :: outFile
    integer, allocatable :: rangeSims(:)

#ifdef win
    do k=1,numPointsType
        if (k == 1)  then
            start = 1
            finish = 0.16*numSims
            allocate(rangeSims(finish))
            rangeSims = (/(i,i=1,finish)/)
        else if (k==2) then
            deallocate(rangeSims)
            start = finish +1
            finish = 0.48*numSims
            allocate(rangeSims(finish-start+1))
            rangeSims = (/(i,i=start,finish)/)
            !real(s,rk)/numSims <0.48
        else if (k==3) then
            deallocate(rangeSims)
            start = finish +1
            finish = 0.59*numSims
            allocate(rangeSims(finish-start+1))
            rangeSims = (/(i,i=start,finish)/)
            !real(s,rk)/numSims <0.59
        else
            deallocate(rangeSims)
            start = finish +1
            finish = numSims
            allocate(rangeSims(finish-start+1))
            rangeSims = (/(i,i=start,finish)/)
        end if


        do n=1,Tperiods
            meanL(n)=sum(real(lpath(n,rangeSims),rk))/real(finish-start+1,rk) !size(lpath(n,:))
            !write (201, * ) meanL(n)
            meanV(n)=sum(real(vpath(n,rangeSims),rk))/real(finish-start+1,rk)
            !write (202, * ) meanV(n)
            meanA(n)=sum(real(grids%Agrid(k,n,apath(n,rangeSims)),rk))/real(numSims,rk)!um(real(apath(n,rangeSims),rk))/real(finish-start+1,rk)
            !write (203, * ) meanA(n)
            meanC(n)=sum(real(cpath(n,rangeSims),rk))/real(finish-start+1,rk)
            !write (204, * ) meanC(n)
            meanY(n)=sum(real(ypath(n,rangeSims),rk))/real(finish-start+1,rk)
            !write (205, * ) meanY(n)
            meanYemp(n)=sum(real(yemp(n,rangeSims),rk))/real(finish-start+1,rk)
            !write (205, * ) meanY(n)
            meanAIME(n)=sum(real(AIME(n,rangeSims),rk))/real(finish-start+1,rk)
            numLC(n) = count(apath(n,:).eq.0)
        end do

        !L
        write (outFile, *), trim("..\\out\lpath"),k,".txt"
        outfile=ADJUSTL(outfile)
        inquire (iolength=requiredl)  meanL
        open (unit=201, file=outfile, status='unknown',recl=requiredl, action='write')
        write (201, '(6E15.3)' ) meanL
        close( unit=201)

        !V
        write (outFile, *), trim("..\\out\Vpath"),k,".txt"
        outfile=ADJUSTL(outfile)
        inquire (iolength=requiredl)  meanV
        open (unit=202, file=outFile, status='unknown',recl=requiredl, action='write')
        write (202, '(6E15.3)' ) meanV
        close( unit=202)

        !A
        write (outFile, *), trim("..\\out\Apath"),k,".txt"
        outfile=ADJUSTL(outfile)
        inquire (iolength=requiredl)  meanA
        open (unit=203, file=outFile, status='unknown',recl=requiredl, action='write')
        write (203, '(6E15.3)' ) meanA
        close( unit=203)

        !C
        write (outFile, *), trim("..\\out\Cpath"),k,".txt"
        outfile=ADJUSTL(outfile)
        inquire (iolength=requiredl)  meanC
        open (unit=204, file=outFile, status='unknown',recl=requiredl, action='write')
        write (204, '(6E15.3)' ) meanC
        close( unit=204)

        !Y
        write (outFile, *), trim("..\\out\Ypath"),k,".txt"
        outfile=ADJUSTL(outfile)
        inquire (iolength=requiredl)  meanY
        open (unit=205, file=outFile, status='unknown',recl=requiredl, action='write')
        write (205, '(6E15.3)' ) meanY
        close( unit=205)

        write (outFile, *), trim("..\\out\YempPath"),k,".txt"
        outfile=ADJUSTL(outfile)
        inquire (iolength=requiredl)  meanYemp
        open (unit=206, file=outFile, status='unknown',recl=requiredl, action='write')
        write (206, '(6E15.3)' ) meanYemp
        close( unit=206)

        write (outFile, *), trim("..\\out\AIMEPath"),k,".txt"
        outfile=ADJUSTL(outfile)
        inquire (iolength=requiredl)  meanAIME
        open (unit=207, file=outFile, status='unknown',recl=requiredl, action='write')
        write (207, '(6E15.3)' ) meanAIME
        close( unit=207)

        write (outFile, *), trim("..\\out\numLC"),k,".txt"
        outfile=ADJUSTL(outfile)
        inquire (iolength=requiredl)  numLC
        open (unit=207, file=outFile, status='unknown',recl=requiredl, action='write')
        write (207, '(6E15.3)' ) numLC
        close( unit=207)
    end do
#else
    !medianA = median(apath(Tretire,:))
    rich = 0
    poor = 0
    do n=1,numSims
        if (apath(Tretire,n) <= medianA) then
            poor = poor + 1
            lpoor(:,poor) = lpath(:,n)
        else
            rich = rich + 1
            lrich(:,rich) = lpath(:,n)
        end if
    end do

    inquire (iolength=requiredl)  meanL
    open (unit=201, file='./out/lpath', status='unknown',recl=requiredl, action='write')
    write (201, * ) 'Header'

    inquire (iolength=requiredl)  meanV
    open (unit=202, file='./out/Vpath', status='unknown',recl=requiredl, action='write')
    write (202, * ) 'Header'

    inquire (iolength=requiredl)  meanA
    open (unit=203, file='./out/Apath', status='unknown',recl=requiredl, action='write')
    write (203, * ) 'Header'

    inquire (iolength=requiredl)  meanC
    open (unit=204, file='./out/Cpath', status='unknown',recl=requiredl, action='write')
    write (204, * ) 'Header'

    inquire (iolength=requiredl)  meany
    open (unit=205, file='./out/Ypath', status='unknown',recl=requiredl, action='write')
    write (205, * ) 'Header'

    inquire (iolength=requiredl)  meanYemp
    open (unit=206, file='./out/YempPath', status='unknown',recl=requiredl, action='write')
    write (206, * ) 'Header'

    inquire (iolength=requiredl)  meanAIME
    open (unit=207, file='./out/AIMEPath', status='unknown',recl=requiredl, action='write')
    write (207, * ) 'Header'

    inquire (iolength=requiredl)  meanPoor
    open (unit=208, file='./out/PoorPath', status='unknown',recl=requiredl, action='write')
    write (208, * ) 'Header'

    inquire (iolength=requiredl)  meanRich
    open (unit=209, file='./out/RichPath', status='unknown',recl=requiredl, action='write')
    write (209, * ) 'Header'

    inquire (iolength=requiredl)  meanRich
    open (unit=210, file='./out/ldata', status='unknown',recl=requiredl, action='write')


    inquire (iolength=requiredl)  meanRich
    open (unit=211, file='./out/adata', status='unknown',recl=requiredl, action='write')

    inquire (iolength=requiredl)  meanRich
    open (unit=212, file='./out/ydata', status='unknown',recl=requiredl, action='write')

    do n=1,Tperiods
        meanL(n)=sum(real(lpath(n,:),rk))/real(numSims,rk)
        write (201, * ) meanL(n)
        meanV(n)=sum(real(vpath(n,:),rk))/real(numSims,rk)
        write (202, * ) meanV(n)
        meanA(n)=sum(real(apath(n,:),rk))/real(numSims,rk)
        write (203, * ) meanA(n)
        meanC(n)=sum(real(cpath(n,:),rk))/real(numSims,rk)
        write (204, * ) meanC(n)
        meanY(n)=sum(real(ypath(n,:),rk))/real(numSims,rk)
        write (205, * ) meanY(n)
        meanYemp(n)=sum(real(yemp(n,:),rk))/real(numSims,rk)
        write (206, * ) meanYemp(n)
        meanAIME(n)=sum(real(AIME(n,:),rk))/real(numSims,rk)
        write (207, * ) meanAIME(n)
        meanRich(n)=sum(real(lrich(n,:),rk))/real(numSims/2,rk)
        write (209, * ) meanRich(n)
        meanPoor(n)=sum(real(lpoor(n,:),rk))/real(numSims/2,rk)
        write (208, * ) meanPoor(n)
    end do
    do i=1,numsims
        do n=32,55
            write (210, * ) lpath(n,i)
            write (211, * ) apath(n,i)
            write (212, * ) lpath(n,i)*ypath(n,i)
        end do
    end do

    close( unit=201)
    close( unit=202)
    close( unit=203)
    close( unit=204)
    close( unit=205)
    close( unit=206)
    close( unit=207)
    close( unit=208)
    close( unit=209)
    close( unit=210)
    close( unit=211)
    close( unit=212)
#endif

    end subroutine

    ! ---------------------------------------------------------------------------------------------------------!
    ! get asset grid
    !---------------------------------------------------------------------------------------------------------!
    subroutine getassetgrid( params, maxInc, Agrid)
    implicit none

    !inputs
    type (structparamstype), intent(in) :: params
    real (kind=rk), intent(in) ::  maxInc(Tperiods)

    !outputs
    real (kind=rk), intent(out) :: Agrid(Tperiods+1, numPointsA)

    !local
    real (kind=rk) :: maxA(Tperiods+1), loggrid(numPointsA), span, test
    integer :: ixt, i

    !Set maximum assets
    maxA(1) = params%startA
    do ixt = 2, Tperiods+1
        maxA(ixt) = (maxA(ixt - 1) + maxInc(ixt-1) ) * (1+params%r)
    end do

    !Create asset grid
    do ixt = 1, Tperiods+1
        !span = maxA(ixt)/ (numPointsA-1)
        !Agrid(ixt, :) = span*(/(i,i=0,numPointsA-1)/)

        span = (log(1.0+maxA(ixt))-log(1.0))/ (numPointsA-1)
        loggrid = log(1.0)+span*(/(i,i=0,numPointsA-1)/)
        Agrid(ixt, :) = (/((exp(loggrid(i))-1.0),i=1,numPointsA)/)

        !span =  (log(1.0+log(1.0+log(1+maxA(ixt)))) - log(1.0+log(1.0+log(1.0))) )/ (numPointsA-1)
        !loggrid = log(1.0+log(1.0+log(1.0))) + span*(/(i,i=0,numPointsA-1)/)
        !Agrid(ixt, :) = (/(exp(exp(exp(loggrid(i))-1.0)-1.0)-1.0,i=1,numPointsA)/) !exp(exp(exp(loggrid)-1)-1)-1
    end do
    test = sum(Agrid(1, :))/size(Agrid(ixt, :))

    end subroutine
    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    !!solve period
    !! Now solve for discretised problem
    ! ---------------------------------------------------------------------------------------------------------!
    subroutine solvePeriodRE(params, grids, ixy, ixA, ixAIME ,ixt,ixType,spa,  EV1, policyA1,policyL,V)
    implicit none

    !input
    real(kind=rk), intent(in) ::  EV1(:,:)
    type (structparamstype), intent(in) :: params
    type (gridsType), intent(in) :: grids
    integer, intent(in) :: ixt,ixType,ixa,ixAIME,ixy, spa

    !output
    real(kind=rk), intent(out) ::  V
    integer, intent(out) :: policyL, policyA1

    !local
    integer :: ixl, policyA1temp
    real(kind=rk) :: negV, negVtemp


    negV = -huge(negv) !-log(0.0) !inf
    !workAge = startAge - 20 + ixt
    do ixL = 0,(numPointsL-1),1           ! points of income choice
        !! Value of income and information for optimisation
        !Y    = ixL*Yin+ (1-ixL)*benefit
        !
        !!AIME only depends on earned income so add spousal
        !!and pensions after calculating it
        !call gross2net(params,grids,Y,ixt,ixl,ixType,AIME,spa)


        call maxutility(params, grids,ixt,ixType,ixa,ixAIME,ixl,ixy,spa,EV1,policyA1temp,negVtemp)
        if (negVtemp > negV) then
            negV = negVtemp
            policyA1=policyA1temp
            policyL=ixL+1
        end if
    end do

    V  = negV


    end subroutine
    ! ---------------------------------------------------------------------------------------------------------!
    !!solve period
    !! when RI is imporant
    ! ---------------------------------------------------------------------------------------------------------!
    subroutine solvePeriodRI(params, grids, ixy, ixA, ixAIME ,ixt,ixType,spa,  EV1, ccpOut, v)
    implicit none

    !input
    real(kind=rk), intent(in) ::  EV1(:,:,:)
    type (structparamstype), intent(in) :: params
    type (gridsType), intent(in) :: grids
    integer, intent(in) :: ixt,ixType,ixa,ixAIME,ixy, spa

    !!output
    real(kind=rk), intent(out) :: ccpOut(:,:,:), v(:)
    !local
    integer, parameter :: hp = selected_real_kind(16)
    integer :: ixl,  ixA1, A1, i, j, iter, dim, dimD ,K, labourChoices
    real(kind=rk) :: Y,  ubA1, A,AIME, cons, va1, vb1, const(numPointsL,numPointsA,grids%supportSPA(ixt)), check, u(numPointsL,numPointsA),EVloc(numAIME) !EVloc(grids%supportSPA(ixt))
    !real(kind=rk) ::
    real(kind=rk), allocatable :: values(:), parameters(:), ccp(:), ccpMat(:,:,:), eye(:,:),  GrossUtil(:,:), test(:)
    integer, save :: saveDim(2), lastixT = 0, unemp
    real(kind=rk), allocatable, save :: locinitialGuessRI(:,:)

    !If suffered unemployment shock then can't choose to work
    labourChoices=mod(ixy,2)*numPointsL+(1-mod(ixy,2))*1

    do ixL = 0,(labourChoices-1),1           ! points of income choice
        ! Value of income and information for optimisation
        A    = grids%Agrid(ixType,ixt, ixA)            ! assets today
        Y    = grids%Ygrid(ixType,ixt, ixY)
        Y    = ixL*Y+ (1-ixL)*grids%benefit(ixt)
        AIME = grids%AIMEgrid(ixType,ixt,ixAIME)

        !Also doesn't matter if income is arbitrary as in that case they are too old for AIME to increase
        call nextAIME(grids,ixt,ixType,ixl,y,AIME)

        call gross2net(params,grids,Y,ixt,ixl,ixType,AIME,spa)


        ubA1 = (A + Y - params%minCons)*(1+params%r);    ! upper bound: assets tomorrow
        do ixA1 = 1, numPointsA
            if (grids%Agrid(ixType,ixt+1, ixA1) > ubA1) then
                const(ixl+1,ixa1:numPointsA,:) = -huge(const(ixl+1,ixa1,1))
                exit
            end if
            A1 = grids%Agrid(ixType,ixt+1, ixA1)

            cons = A  + Y - (A1)/(1+params%r)
            u(ixl+1,ixa1) = utility(params,cons,ixL)

            do i = 1, grids%supportSPA(ixt)-1
                !(/(j,j=numPointsSPA-grids%supportSPA(ixt),numPointsSPA)/)
                !
                EVloc =  (1-params%p)*EV1(ixA1,:,numPointsSPA-grids%supportSPA(ixt)+i)+params%p*EV1(ixA1,:,numPointsSPA-grids%supportSPA(ixt)+i+1)
                const(ixl+1,ixa1,i) = objectivefunc(params, grids,grids%Agrid(ixType,ixt+1, ixA1), A, Y,ixL,ixt,ixType, AIME,EVloc)
                !call  linearinterp1(grids%AIMEgrid(ixType,ixt + 1, :),EV1(ixA1,:,numPointsSPA-grids%supportSPA(ixt)+i), numAIME, AIME, EVloc(i), 1, 1 )
            end do
            !do i = 1, grids%supportSPA(ixt)-1
            !    va1 = (1-params%p)*EVloc(i)+params%p*EVloc(i+1)!dot_product(grids%posteriorSPA(ixt,(/(j,j=numPointsSPA-grids%supportSPA(ixt),numPointsSPA)/)),EVloc)
            !    VB1 = params%thetab*((A1+params%K)**(1-params%gamma))/(1-params%gamma)
            !    const(ixl+1,ixa1,i) = u(ixl+1,ixa1) + params%beta * ((1- grids%mortal(ixt))*VA1+grids%mortal(ixt)*VB1)
            !end do
            EVloc =  EV1(ixA1,:,numPointsSPA)
            const(ixl+1,ixa1,i) = objectivefunc(params, grids,grids%Agrid(ixType,ixt+1, ixA1), A, Y,ixL,ixt,ixType, AIME,EVloc)

            !va1 =EVloc( grids%supportSPA(ixt))
            !VB1 = params%thetab*((A1+params%K)**(1-params%gamma))/(1-params%gamma)
            !const(ixl+1,ixa1, grids%supportSPA(ixt)) = u(ixl+1,ixa1) + params%beta * ((1- grids%mortal(ixt))*VA1+grids%mortal(ixt)*VB1)

        end do
    end do

    dimD = labourChoices*(ixa1-1)
    dim = grids%supportSPA(ixt)*labourChoices*(ixa1-1)
    allocate(values(dim+1),parameters(dim),ccp(dim),ccpMat(grids%supportSPA(ixt),labourChoices,ixa1-1), GrossUtil(dimD,grids%supportSPA(ixt)),test(dimD))
    !Default ccp to zero for regions that are ignored
    ccp = 0.0
    !Want to save grid to reuse over employment states
    unemp = mod(ixy-1,2)+1
    if (ixt/=lastixT) then
        parameters = grids%initialGuessRI(1,1:dim)
        if (allocated(locinitialGuessRI)) deallocate(locinitialGuessRI)
        allocate(locinitialGuessRI(2,grids%supportSPA(ixt)*labourChoices*numPointsA))
    else if (ixy<=2) then
        parameters = grids%initialGuessRI(unemp,1:dim)
    else
        parameters=0.0
        if (saveDim(unemp) > grids%supportSPA(ixt)*(ixa1-1) ) saveDim(unemp) = grids%supportSPA(ixt)*(ixa1-1)
        do k=1,labourChoices
            parameters((k-1)*saveDim(unemp)+1:k*saveDim(unemp))= &
                locinitialGuessRI(unemp,(k-1)*grids%supportSPA(ixt)*numPointsA+1:(k-1)*grids%supportSPA(ixt)*numPointsA+saveDim(unemp))
        end do
        !Make sure we don't pass zeros across as gets stuck at them if give them as initial conds
        parameters = parameters + 0.1/dim
    end if
    !store last calculated period to know when we change period
    lastixT = ixt

    !parameters = grids%initialGuessRI(1,1:dim)
    !if (mod(ixy,2)==0) then
    !    parameters(grids%supportSPA(ixt)*(ixa1-1)+1:dim)=0.0
    !end if

    !Normalise guess to be prob dist
    do i=1,grids%supportspa(ixt)
        check = sum(parameters((/(i+j*grids%supportspa(ixt),j=0,dimd-1)/)))
        parameters((/(i+j*grids%supportspa(ixt),j=0,dimd-1)/))=parameters((/(i+j*grids%supportspa(ixt),j=0,dimd-1)/))*check**-1
    end do

    check = 1.0
    !ccp=parameters
    iter=0
    !Convergence criteria should depend of cost of attention otherwise utility can get very small and it just assigns equiprobable dist
    do while (check >= 0.0001*params%lambda )!0.001!0.005!0.0001!0.00005 !0.00001 !/params%lambda
        iter=iter+1
        ccp = func(parameters)
        check = sum(abs(ccp-parameters))
        parameters=ccp
    end do

    !Check is good enough solution
    !check = func(ccp)
    if (check/dim>0.02) then
        write (*,*) "CCP don't solve system of equations", check/dim, iter, ixa,ixy
        stop
    end if
    do i=1,grids%supportSPA(ixt)
        check = sum(ccp((/(i+j*grids%supportSPA(ixt),j=0,Dimd-1)/)))
        if (check > 1.01 .OR. check<0.99) then
            write (*,*) "CCP doesn't sum to 1"
            stop
        end if
        !ccp((/(i+j*grids%supportSPA(ixt),j=0,Dimd-1)/))=ccp((/(i+j*grids%supportSPA(ixt),j=0,Dimd-1)/))*check**-1
    end do

    !!Cache for use on next round
    if (ixy < numPointsY) then
        saveDim(unemp) = grids%supportSPA(ixt)*(ixa1-1)
        locinitialGuessRI(unemp,:)=0.0
        do k=1,labourChoices
            locinitialGuessRI(unemp,(k-1)*grids%supportSPA(ixt)*numPointsA+1:(k-1)*grids%supportSPA(ixt)*numPointsA+saveDim(unemp))= ccp((k-1)*saveDim(unemp)+1:k*saveDim(unemp))
        end do
    end if

    ccpOut=0.0
    !ccpMat=reshape(ccp, (/grids%supportSPA(ixt),labourChoices,ixa1-1/))
    do i=1,labourChoices
        do j =1,ixa1-1
            do k=1,grids%supportSPA(ixt)
                ccpMat(k,i,j) =  ccp(k+(j-1)*grids%supportSPA(ixt)+(i-1)*saveDim(unemp))
            end do
        end do
    end do
    ccpOut(1:grids%supportSPA(ixt),1:labourChoices,1:ixa1-1)=ccpMat
    do i=1,grids%supportSPA(ixt)
        check = sum(ccpOut(i,:,:))
        if (check > 1.01 .OR. check<0.99) then
            write (*,*) "CCP doesn't sum to 1"
            stop
        end if
        !ccp((/(i+j*grids%supportSPA(ixt),j=0,Dimd-1)/))=ccp((/(i+j*grids%supportSPA(ixt),j=0,Dimd-1)/))*check**-1
    end do

    !Calculate continuation value
    do i=1,grids%supportSPA(ixt)
        test= exp(GrossUtil(:,i))
        v(i) = log(sum(test))
    end do

    contains
    function func(p) result(res)
    implicit none
    !input
    real(kind=rk), intent(in) :: p(:)
    !output
    real(kind=rk), allocatable :: res(:)
    !local
    real(kind=hp) :: systemEq(dim),test(dim), loc!,grids%supportSPA(ixt)
    integer:: d, y, locl, loca1, j

    allocate(res(dim))

    do d=1,dimD
        do y=1,grids%supportSPA(ixt)
            locl=(d-1)/(ixa1-1)+1
            locA1=mod(d-1,ixa1-1)+1
            GrossUtil(d,y) =  const(locl,locA1,y) + log(dot_product(grids%posteriorSPA(ixt,numPointsSPA-grids%supportSPA(ixt)+1:numPointsSPA),p((d-1)*grids%supportSPA(ixt)+1:d*grids%supportSPA(ixt))))
        end do
    end do

    do d=1,dimD
        do y=1,grids%supportSPA(ixt)
            res((d-1)*grids%supportSPA(ixt)+y) = exp(GrossUtil(d,y)) / sum(exp(GrossUtil(:,y)))
        end do
    end do

    end function
    end subroutine

    !!-----------------------------------------------------------------------------------------------------------!
    ! Unpack all arrays
    !!-----------------------------------------------------------------------------------------------------------!
    subroutine unpackArrays(p, V,eV,vecP,vecV, vecEV,dim1,dimP,dimV,thisCoreStart,thisCoreEnd)
    implicit none

    !inputs
    integer, intent(in) :: dim1,dimP, dimV, thisCoreStart,thisCoreEnd
    real(kind=rk), intent(in) :: p( :, :, :,:,:,:)
    real(kind=rk), intent(in) :: V( :, :, :,:), eV( :, :, :,:)

    !outputs
    real(kind=rk), intent(out) :: vecP(:)
    real(kind=rk), intent(out) :: vecV(:), vecEV(:)

    call unpackArray(p,vecP,dim1,dimP,thisCoreStart,thisCoreEnd)
    !call unpackArray(pA,vecPA,dim1,dim2,thisCoreStart,thisCoreEnd)
    !call unpackArray(pC,vecPC,dim1,dim2,thisCoreStart,thisCoreEnd)
    call unpackArray(v,vecV,dim1,dimV,thisCoreStart,thisCoreEnd)
    call unpackArray(eV,vecEV,dim1,dimV,thisCoreStart,thisCoreEnd)

    end subroutine
    !!-----------------------------------------------------------------------------------------------------------!
    ! Unpack integer array
    !!-----------------------------------------------------------------------------------------------------------!
    subroutine unpackArrayInt(Array,vec,dim1,dim2,thisCoreStart,thisCoreEnd)
    implicit none
    !inputs
    integer, intent(in) :: Array( :, :, :,:), dim1,dim2,thisCoreStart,thisCoreEnd

    !outputs
    integer, intent(out) :: vec(:)

    !local
    integer, allocatable :: mat(:,:)

    allocate( mat(dim1,dim2))

    mat = reshape(Array,(/dim1,dim2/))
    vec = reshape(mat(thisCoreStart:thisCoreEnd,:),(/(thisCoreEnd-thisCoreStart+1)*dim2/))

    end subroutine
    !!-----------------------------------------------------------------------------------------------------------!
    ! Unpack integer array
    !!-----------------------------------------------------------------------------------------------------------!
    subroutine unpackArrayReal(Array,vec,dim1,dim2,thisCoreStart,thisCoreEnd)
    implicit none

    !inputs
    real(kind=rk), intent(in) :: Array( :, :, :,:,:,:)
    integer, intent(in) :: dim1,dim2,thisCoreStart,thisCoreEnd

    !outputs
    real(kind=rk), intent(out) :: vec(:)

    !local
    real(kind=rk), allocatable :: mat(:,:)

    allocate( mat(dim1,dim2))

    mat = reshape(Array,(/dim1,dim2/))
    vec = reshape(mat(thisCoreStart:thisCoreEnd,:),(/(thisCoreEnd-thisCoreStart+1)*dim2/))

    end subroutine
    !!-----------------------------------------------------------------------------------------------------------!
    ! Unpack integer small array
    !!-----------------------------------------------------------------------------------------------------------!
    subroutine unpackArraySmall(Array,vec,dim1,dim2,thisCoreStart,thisCoreEnd)
    implicit none

    !inputs
    real(kind=rk), intent(in) :: Array( :, :, :,:)
    integer, intent(in) :: dim1,dim2,thisCoreStart,thisCoreEnd

    !outputs
    real(kind=rk), intent(out) :: vec(:)

    !local
    real(kind=rk), allocatable :: mat(:,:)

    allocate( mat(dim1,dim2))

    mat = reshape(Array,(/dim1,dim2/))
    vec = reshape(mat(thisCoreStart:thisCoreEnd,:),(/(thisCoreEnd-thisCoreStart+1)*dim2/))

    end subroutine

    !!-----------------------------------------------------------------------------------------------------------!
    !Intialise gues for Nedler-Mead Algorithm
    !!-----------------------------------------------------------------------------------------------------------!
    subroutine initialGuess(rank,params,grids,moments,weights,p,y)
    implicit none

    !inputs
    integer, intent(in) :: rank
    real(kind=rk), intent(in) :: moments(2,24), weights(:,:)

    !changing
    type (structparamstype), intent(inout) :: params
    type (gridsType), intent(inout) :: grids

    !outpus
    real(kind=rk), intent(out) :: y(dimEstimation+1), p(dimEstimation+1,dimEstimation)

    !local
    integer :: i, seedIn, n
    real (kind=rk) ::  uniformRand(dimEstimation+1)
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed

    if (fullLifeCycle) then
        if (rank==0) then
            print '("Guess 1")'
        end if
        params%nu =      0.38022456150504280
        params%beta =  0.97235322545400193
        params%gamma =  2.092041817380061
        params%db(1) = 0.91387622628345655
        params%db(2) = -4.7393420983952571E-005

        p(1,1) = params%nu
        p(1,2) = params%beta
        p(1,3) = params%gamma
        p(1,4) = params%db(1)
        p(1,5) = params%db(2)
        y(1) = gmm_criteria(p(1,:))

        if (rank==0) then
            print '("Guess 2")'
        end if
        p(2,1) = 0.4637
        p(2,2) = 0.970
        P(2,3) = 1
        p(2,4) = params%db(1)*0.9
        p(2,5) = params%db(2)*1.2
        y(2) = gmm_criteria(p(2,:))

        if (rank==0) then
            print '("Guess 3")'
        end if
        p(3,1) = 0.322
        p(3,2) = 0.9843
        P(3,3) = 2
        p(3,4) = params%db(1)*1.1
        p(3,5) = params%db(2)*0.7
        y(3) = gmm_criteria(p(3,:))

        if (rank==0) then
            print '("Guess 4")'
        end if
        p(4,1) = 0.55
        p(4,2) = 0.96
        P(4,3) = 0.5
        p(4,4) = params%db(1)*1.3
        p(4,5) = params%db(2)*0.95
        y(4) = gmm_criteria(p(4,:))

        if (rank==0) then
            print '("Guess 5")'
        end if
        p(5,1) = 0.15
        p(5,2) = 0.9999
        P(5,3) = 4
        p(5,4) = params%db(1)*0.85
        p(5,5) = params%db(2)*1.15
        y(5) = gmm_criteria(p(5,:))

        if (rank==0) then
            print '("Guess 6")'
        end if
        p(6,1) = 0.27
        p(6,2) = 0.986
        P(6,3) = 0.9
        p(6,4) = params%db(1)*0.99
        p(6,5) = params%db(2)*0.87
        y(6) = gmm_criteria(p(6,:))

    else

        seedIn = 16101988
        !Set seed
        CALL RANDOM_SEED(size = n)
        ALLOCATE(seed(n))
        seed = seedIn * (/ (i - 1, i = 1, n) /)
        CALL RANDOM_SEED(PUT = seed)
        DEALLOCATE(seed)

        !!get uniform random number
        CALL RANDOM_NUMBER(uniformRand)
        uniformRand = -0.5+uniformRand
        do i = 1, dimEstimation+1
            if (rank==0) then
                write (*,*) "Guess ", i
            end if
            p(i,1) = params%nu*(1+uniformRand(i))
            p(i,2) = (0.95+0.1*uniformRand(i))!params%beta*(1+uniformRand(i))
            p(i,3) = params%gamma*(1+uniformRand(i))
            p(i,4) = params%db(1)*(1+uniformRand(i))
            p(i,5) = params%db(2)*(1+uniformRand(i))
            p(i,6) = params%thetab*(1+uniformRand(i))
            y(i) = gmm_criteria(p(i,:))
            !write (*,*) uniformRand(i), y(i)
        end do
    end if

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
    if (minval(control(1:4)) < 0 .or.  control(5)> 0 )  then !.OR. control(6)< 0
        gmm_criteria = huge(gmm_criteria)
        return
    end if

    params%nu = control(1)
    params%beta = control(2)
    params%gamma = control(3)
    params%db(1)= control(4)
    params%db(2)= control(5)
    params%thetab = control(6)
    gmm_criteria = gmm(params,grids,moments,weights) !*-1.0

    end function
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Add Unemployment shocks
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine addUnemploymentshock(Ygrid, incTransitionMrx,grids,typeSim)
    implicit none

    !Input
    real (kind=rk), intent(in) :: Ygrid(:,:), incTransitionMrx(:,:)
    integer, intent(in) :: typeSim

    !Changingin
    type (gridsType), intent(inout) :: grids

    !local
    real (kind=rk) :: reemployRaw(4,10), unemployRaw(4,10)
    integer :: i

    unemployRaw(1,:) = (/0.227, 0.049, 0.027, 0.014, 0.004, 0.007, 0.013, 0.017, 0.007, 0.01/)
    unemployRaw(2,:) = (/0.099, 0.026, 0.01, 0.013, 0.008, 0.012, 0.003, 0.003, 0.002, 0.0/)
    unemployRaw(3,:) = (/0.201, 0.024, 0.054, 0.032, 0.029, 0.038, 0.02, 0.005, 0.004, 0.01/)
    unemployRaw(4,:) = (/0.063, 0.016, 0.01, 0.006, 0.003, 0.002, 0.004, 0.001, 0.003, 0.0/)

#if PROD_SIZE == 5
    !if (numPointsProd==5) then
    do i=1,numPointsType
        grids%unemploy(i,1) =unemployRaw(i,1)+unemployRaw(i,2)!(/0.1762, 0.0415,0.0257,0.0165,0.0090/)
        grids%unemploy(i,2) =unemployRaw(i,3)+unemployRaw(i,4)
        grids%unemploy(i,3) =unemployRaw(i,5)+unemployRaw(i,6)
        grids%unemploy(i,4) =unemployRaw(i,7)+unemployRaw(i,8)
        grids%unemploy(i,5) =unemployRaw(i,9)+unemployRaw(i,10)
    end do
#elif  PROD_SIZE == 10          
    !else if (numPointsProd==10) then
    grids%unemploy = unemployRaw
    !end if
#endif
    reemployRaw(1,:) = (/0.13, 0.043, 0.022, 0.005, 0.002, 0.007, 0.01, 0.002, 0.0, 0.007 /)
    reemployRaw(2,:) = (/0.174, 0.005, 0.013, 0.017, 0.013, 0.017, 0.013, 0.003, 0.003, 0.0 /)
    reemployRaw(3,:) = (/0.118, 0.029, 0.037, 0.011, 0.008, 0.013, 0.003, 0.0, 0.005, 0.0 /)
    reemployRaw(4,:) = (/0.238, 0.032, 0.02, 0.008, 0.008, 0.012, 0.008, 0.004, 0.0, 0.004/)

#if PROD_SIZE == 5
    !if (numPointsProd==5) then
    !grids%reemploy = (/ 0.1923, 0.0333, 0.0200, 0.0108,0.0047/)
    do i=1,numPointsType
        grids%reemploy(i,1) =reemployRaw(i,1)+reemployRaw(i,2)!(/0.1762, 0.0415,0.0257,0.0165,0.0090/)
        grids%reemploy(i,2) =reemployRaw(i,3)+reemployRaw(i,4)
        grids%reemploy(i,3) =reemployRaw(i,5)+reemployRaw(i,6)
        grids%reemploy(i,4) =reemployRaw(i,7)+reemployRaw(i,8)
        grids%reemploy(i,5) =reemployRaw(i,9)+reemployRaw(i,10)
    end do
#elif  PROD_SIZE == 10         
    !else if (numPointsProd==10) then
    grids%reemploy = reemployRaw
    !end if
#endif        

    grids%Ygrid(typeSim,:,:) = 0.0
    grids%incTransitionMrx(:,:) = 0.0
    grids%Ygrid(typeSim,:,(/(i,i=1,numPointsY-1,2)/))= Ygrid(:,:)
    do i = 1, numPointsY
        if (mod(i,2)==1) then
            grids%incTransitionMrx(i,(/(i,i=1,numPointsY-1,2)/)) = (1-grids%unemploy(typeSim,i/2+1))*incTransitionMrx(i/2+1,:)
            grids%incTransitionMrx(i,(/(i,i=2,numPointsY,2)/)) = grids%unemploy(typeSim,i/2+1)*incTransitionMrx(i/2+1,:)
        else
            grids%incTransitionMrx(i,i-1) = grids%reemploy(typeSim,i/2)
            grids%incTransitionMrx(i,i) = 1-grids%reemploy(typeSim,i/2)
        end if
    end do

    end subroutine
    !!-----------------------------------------------------------------------------------------------------------!
    !!-----------------------------------------------------------------------------------------------------------!
    !Initalise sstuff
    !!-----------------------------------------------------------------------------------------------------------!
    !!-----------------------------------------------------------------------------------------------------------!
    subroutine setupMisc(params,grids)
    implicit none

    !Changingin
    type (structparamstype), intent(inout) :: params
    type (gridsType), intent(inout) :: grids

    !local
    REAL(kind=rk) :: temp(85),  temp2(86), temp3(numSims, Tperiods)
    integer :: ios, I
    real (kind=rk) ::  sig_inc, sig_initial, conditioningProb
    real (kind=rk) :: e(Tperiods, numSims)   ! the innovations to log income
    real (kind=rk) :: logy1(1,numSims)        ! draws for initial income
    real (kind=rk) :: ly(Tperiods, numSims)           ! log income
    integer :: s, t,  workAge
    real (kind=rk) ::  uniformRand(Tperiods,numSims)
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
    integer :: seedIn, n, productiv(1), typeSim, requiredl, finish
    logical :: unemployed
    character(len=1024) :: outFile
    integer, allocatable :: rangeSims(:)

#ifdef win
    open (unit = 1001,file='..\\..\\moments\\InitialAssets.txt', action='read', IOSTAT = ios)
#else
    open (unit = 1001,file='../moments/InitialAssets.txt', action='read', IOSTAT = ios)
#endif
    read (1001, *,  IOSTAT = ios) grids%initialAssets
    close (unit=1001)

    temp = (/-0.01374, -0.01069, -0.00767, -0.00467, -0.00169, 0.00126, 0.00418,  0.00708, 0.00996, 0.01281, 0.01563, 0.01843, &
        0.02121, 0.02396, 0.02668,  0.02938,  0.03206, 0.03471, 0.03733,  0.03994, 0.04251, 0.04506, 0.04759, 0.05009, &
        0.05256, 0.05501,  0.05744,  0.05984,  0.06221, 0.06457, 0.06689,  0.06919, 0.07147, 0.07372, 0.07594, 0.07814, &
        0.08032, 0.08247,  0.08460,  0.08670,  0.08877, 0.09082, 0.09285,  0.09485, 0.09683, 0.09878, 0.10070, 0.10261, &
        0.10448, 0.10633, 0.10816,  0.10996,  0.11174, 0.11349, 0.11521,  0.11691, 0.11859, 0.12024, 0.12187, 0.12347, &
        0.12505, (I*0.0,I=1,24)/)
    grids%fc = temp(startAge-20+1:85)

    call getIncomeGrid(params, grids)

    temp2 = (/0.0, 0.0, 0.0, 0.0, 0.0, &
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.011071, 0.011907, 0.012807, &
        0.013676, 0.01475, 0.015818, 0.016846, 0.018079, 0.019343, 0.020659, 0.0218, 0.023505, 0.025202, 0.02696, &
        0.028831, 0.031017, 0.033496, 0.036024, 0.038912, 0.042054, 0.045689, 0.049653, 0.054036, 0.05886, 0.064093, &
        0.069636, 0.07533, 0.081069, 0.086912, 0.093067, 0.099807, 0.107483, 0.116125, 0.125196, 0.134361, 0.143881, &
        0.1542, 0.165675, 0.17842, 0.192363, 0.2117, 0.1672,0.1565, 0.1485,0.1459, 1.0/);
    grids%mortal = temp2(startAge-20+1:86)

    ! Obtain time series of incomes for our simulated individuals
    ! Draw random draws for starting income and for innovations
    sig_inc = params%sigma/ ((1-params%rho**2)**0.5)
    sig_initial = 0.2950 !mean([0.073 0.053 0.110 0.112])^0.5;
    ! normally distributed random draws for the innovation
#ifdef win
    open (unit = 1001,file='..\\Errors.txt', action='read', IOSTAT = ios)
#else
    open (unit = 1001,file='./data/Errors.txt', action='read', IOSTAT = ios)
#endif
    read (1001, *) temp3
    close (unit=1001)
    e = transpose(temp3)
#ifdef win
    open (unit = 1001,file='..\\IntProd.txt', action='read', IOSTAT = ios)
#else
    open (unit = 1001,file='./data/IntProd.txt', action='read', IOSTAT = ios)
#endif

    read (1001, *) logy1
    close (unit=1001)

    seedIn = 01042017
    !Set seed
    !CALL RANDOM_SEED(size = n)
    !ALLOCATE(seed(n))
    !seed = seedIn * (/ (i - 1, i = 1, n) /)
    !CALL RANDOM_SEED(PUT = seed)
    !DEALLOCATE(seed)

    !!get uniform random number
    CALL RANDOM_NUMBER(uniformRand)

    do s = 1, numSims, 1                           ! loop through individuals
        if (numPointsType == 1) then
            typeSim = 1
        else
            if (real(s,rk)/numSims <0.16) then
                typeSim = 1
            else if (real(s,rk)/numSims <0.48) then
                typeSim = 2
            else if (real(s,rk)/numSims <0.59) then
                typeSim =3
            else
                typeSim =4
            end if
        end if
        unemployed = .FALSE.
        ! Get all the incomes, recursively
        ly(1, s) = truncate(logy1(1, s), -normBnd*sig_inc,normBnd*sig_inc )
        grids%Simy(1, s) = exp(ly(1, s)+params%delta(typeSim,1)*(startAge - 20 + 1)**2+params%delta(typeSim,2)*(startAge - 20 + 1)+params%delta(typeSim,3)-grids%fc(1))
        do t = 1,Tperiods-1,1                              ! loop through time periods for a particular individual
            workAge = startAge - 20 + t
            if (unemployed) then
                if (uniformRand(t,s) < grids%reemploy(typeSim,productiv(1))) then
                    unemployed = .FALSE.
                    ly(t+1, s) = (1 -params%rho) * params%mu + params%rho * ly(t, s) + e(t + 1, s)
                    ly(t+1, s) = truncate(ly(t+1, s), -normBnd*sig_inc,normBnd*sig_inc )
                    grids%Simy(t+1, s) = exp( ly(t+1, s) + params%delta(typeSim,1)*(workAge+1)**2+params%delta(typeSim,2)*(workAge+1)+params%delta(typeSim,3)-grids%fc(t) )
                else
                    ly(t+1, s)  = ly(t, s)
                    grids%Simy(t+1, s) = 0
                end if
            else
                productiv = minloc((grids%ygrid(typeSim,t,(/(i,i=1,numPointsY-1,2)/))-grids%Simy(t, s)),(grids%ygrid(typeSim,t,(/(i,i=1,numPointsY-1,2)/))-grids%Simy(t, s))>0)
                if (uniformRand(t,s) < grids%unemploy(typeSim,productiv(1))) then
                    unemployed = .TRUE.
                    ly(t+1, s)  = ly(t, s)
                    grids%Simy(t+1, s) = 0
                else
                    ly(t+1, s) = (1 -params%rho) * params%mu + params%rho * ly(t, s) + e(t + 1, s)
                    ly(t+1, s) = truncate(ly(t+1, s), -normBnd*sig_inc,normBnd*sig_inc )
                    grids%Simy(t+1, s) = exp( ly(t+1, s) + params%delta(typeSim,1)*(workAge+1)**2+params%delta(typeSim,2)*(workAge+1)+params%delta(typeSim,3)-grids%fc(t) )
                end if
            end if
        end do ! t

    end do! s

    if (TendRI>Tretire) then
        grids%supportSPA(1:Tretire-1) = numPointsSPA
        do s = Tretire, TendRI
            grids%supportSPA(s) = grids%supportSPA(s-1) - 1
        end do

        !posteriorSPA(TendRI,numPointsSPA)
        grids%posteriorSPA=0.0
        do t = 1, TendRI
            workAge = startAge - 20 + t
            if (t<Tretire) then
                do s=0,(numPointsSPA -2)
                    grids%posteriorSPA(t,s+1) =  choose(workAge,s)*(params%p**s)*((1-params%p)**(workAge-s))
                end do
                grids%posteriorSPA(t,numPointsSPA)=1-sum(grids%posteriorSPA(t,:))
            else
                grids%posteriorSPA(t,1:(t-Tretire+1))=0.0
                conditioningProb=0.0
                do s=0,(t-Tretire)
                    conditioningProb =  choose(workAge,s)*(params%p**s)*((1-params%p)**(workAge-s))
                end do
                do s=(t-Tretire+1),(numPointsSPA -2)
                    grids%posteriorSPA(t,s+1) =  (choose(workAge,s)*(params%p**s)*((1-params%p)**(workAge-s)))/(1-conditioningProb)
                end do
                grids%posteriorSPA(t,numPointsSPA)=1-sum(grids%posteriorSPA(t,:))
            end if
        end do
    end if
    !seedIn =


    !!Set seed

    !ALLOCATE(seed(n))
    !seed = seedIn * (/ (i - 1, i = 1, n) /)
    !CALL RANDOM_SEED(PUT = seed)
    !DEALLOCATE(seed)

    !!get uniform random number
    CALL RANDOM_NUMBER(grids%initialGuessRI)



    !if (rank==0) then
    !    !write (*,*) rank
    !    finish = 0.16*numSims
    !    allocate(rangeSims(finish))
    !    rangeSims = (/(i,i=1,finish)/)
    !    write (outFile, *), trim("..\\temp\ypath"),1,".txt"
    !    outfile=ADJUSTL(outfile)
    !    inquire (iolength=requiredl)  grids%Simy(:, rangeSims)
    !    open (unit=201, file=outfile, status='unknown',recl=requiredl, action='write')
    !    write (201, '(6E15.3)' ) grids%Simy(:, rangeSims)
    !    close( unit=201)
    !end if

    end subroutine

    !!-----------------------------------------------------------------------------------------------------------!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Law motion AIME
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!-----------------------------------------------------------------------------------------------------------!
    subroutine nextAIME(grids,ixt,typeSim,ixl,yin, AIME)
    implicit none

    !inputs
    !type (structparamstype), intent(in) :: params
    type (gridsType), intent(in) :: grids
    integer, intent(in) :: ixt, typeSIm, ixl
    real(kind=rk), intent(in):: yin

    !changing
    real (kind=rk), intent(inout) :: AIME

    !local
    integer :: workAge

    workAge = startAge - 20 + ixt
    ! Next peridos AIME
    if (ixL==1 .AND.  ixt < spouseretire ) then ! spouseretire
        AIME =  Yin/workAge + AIME * (workAge-1)/workAge
    else if ((ixL==0 .AND.  ixt < spouseretire )) then !spouseretire
        AIME = AIME * (workAge-1)/workAge
    end if

    end subroutine nextAIME
    !!-----------------------------------------------------------------------------------------------------------!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !In period without RI max utility
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!-----------------------------------------------------------------------------------------------------------!
    subroutine maxutility(params, grids,ixt,typeSim,ixa,ixAIME,ixl,ixy,spa,EV,policyA1temp,negV)
    implicit none

    !inputs
    type (structparamstype), intent(in) :: params
    type (gridsType), intent(in) :: grids
    integer, intent(in) :: ixt,typeSim,ixa,ixAIME,ixl,ixy, spa
    real (kind=rk), intent(in) :: EV(:,:)

    !output
    real(kind=rk), intent(out) :: negV
    integer, intent(out) :: policyA1temp

    !local
    real (kind=rk) :: A,Y,AIME,ubA1, negVtemp, EV1(numAIME)
    integer :: ixA1
    negV = -huge(negv)

    ! Value of income and information for optimisation
    A    = grids%Agrid(typeSim,ixt, ixA)            ! assets today
    Y    = grids%Ygrid(typeSim,ixt, ixY)
    Y    = ixL*Y+ (1-ixL)*grids%benefit(ixt)
    AIME = grids%AIMEgrid(typeSim,ixt,ixAIME)
    !Also doesn't matter if income is arbitrary as in that case they are too old for AIME to increase
    call nextAIME(grids,ixt,typeSim,ixl,y,AIME)
    !set spa=60 as it shouldn't matter as this subroutine only called if above SPA
    call gross2net(params,grids,Y,ixt,ixl,typeSim,AIME,spa)

    ubA1 = (A + Y - params%minCons)*(1+params%r);    ! upper bound: assets tomorrow
    do ixA1 = 1, numPointsA
        if (grids%Agrid(typeSim,ixt+1, ixA1) > ubA1) then
            exit
        end if
        EV1  = EV(ixA1,:)  ! relevant section of EV matrix (in assets tomorrow) because income is irrelevant
        negVtemp = objectivefunc(params, grids,grids%Agrid(typeSim,ixt+1, ixA1), A, Y,ixL,ixt,typeSim, AIME,EV1)
        if (negVtemp > negV) then
            negV = negVtemp
            policyA1temp = ixa1
        end if
    end do

    end subroutine
    !!-----------------------------------------------------------------------------------------------------------!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !In period without RI max utility
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!-----------------------------------------------------------------------------------------------------------!
    !function integrateSPA(grids,EV,ixt,knowsSPA)
    !implicit none
    !
    !!inputs
    !type (gridsType), intent(in) :: grids
    !real (kind=rk), intent(in) :: EV(:,:,:)
    !logical, intent(in) :: knowsSPA
    !integer, intent(in) :: ixt
    !
    !!output
    !real (kind=rk) :: integrateSPA(:,:)
    !
    !!local
    !integer:: i, j
    !
    !do i=1,size(EV,1)
    !    do j=1,size(EV,2)
    !        if (knowsSPA) then
    !            integrateSPA(i,j) =
    !        else
    !
    !        end if
    !    end do
    !end do
    !
    !end function
    end module routines
