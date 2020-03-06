    module routines
    use Header
    use routines_generic
    use bobyqa_module

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

    !sig_inc = params%sigma/((1-params%rho**2)**0.5)
    !
    !call tauchen(numPointsProd,params%mu,params%rho,params%sigma,normBnd,ly,incTransitionMrx)

    !params%pension = 107.45*52
    upper(:,1) = 0
    do typeSim = 1, numPointsType
        sig_inc = params%sigma(typeSim)/((1-params%rho(typeSim)**2)**0.5)

        call tauchen(numPointsProd,params%mu,params%rho(typeSim),params%sigma(typeSim),normBnd,ly,incTransitionMrx)
        do t=1 , Tperiods
            workAge = startAge + t - 1
            Ygrid(typeSim,t,:)= weeksYear*exp(ly+params%delta(typeSim,1)*workAge**2+params%delta(typeSim,2)*workAge+params%delta(typeSim,3)-grids%fc(t))
            grids%maxInc(typeSim,t) = weeksYear*exp((normBnd * sig_inc)+params%delta(typeSim,1)*workAge**2+params%delta(typeSim,2)*workAge+params%delta(typeSim,3)-grids%fc(t))
            !!This should be discounted to average not sum
            upper(typeSim,t+1) = upper(typeSim,t) + grids%maxInc(typeSim,t)
            if (t <=params%ERA(typeSim)) then
                a = ((upper(typeSim,t+1))/t)/(numAIME-1) !-3700

                grids%AIMEgrid(typeSim,t,:) = a*(/(i,i=0,numAIME-1)/)
                grids%AIMEgrid(typeSim,t,:) = grids%AIMEgrid(typeSim,t,:)! +3700
            else
                !upper(typeSim,t+1) = ((t-1)*upper(typeSim,t) + grids%maxInc(typeSim,t))/t
                !if (t <=params%ERA(typeSim)) then
                !    !upper(typeSim,t+1) = 0.75* upper(typeSim,t+1)
                !    a = ((0.6*upper(typeSim,t+1)-2000.0))/(numAIME-1) !/t
                !
                !    grids%AIMEgrid(typeSim,t,:) = a*(/(i,i=0,numAIME-1)/)
                !    grids%AIMEgrid(typeSim,t,:) = grids%AIMEgrid(typeSim,t,:) + 2000.0! +3700
                !else
                grids%AIMEgrid(typeSim,t,:) = grids%AIMEgrid(typeSim,t-1,:)
            end if
            if (t < stopwrok) then
                grids%benefit(t) = 73.10*52
            else
                grids%benefit(t) = 0.0
            end if
            grids%YgridExtrap(typeSim,t,:) = Ygrid(typeSim,t,:)
        end do
        call addUnemploymentshock(Ygrid(typeSim,:,:), incTransitionMrx,grids,typeSim)
        grids%AIMEgrid(typeSim,Tperiods+1,:) = grids%AIMEgrid(typeSim,Tperiods,:)
        !write(*,*) -params%pension(typeSim,1) /(2.0*params%pension(typeSim,2))
        !write(*,*) -params%db(typeSim,1) /(2.0*params%db(typeSim,2))
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
    subroutine solveValueFunction( params, grids, show, intermediateToFile )
    implicit none

#ifdef mpiBuild
    include 'mpif.h'
#endif     
    !inputs
    type (structparamstype), intent(in) :: params
    type (gridsType), intent(in) :: grids
    logical, intent(in) :: show
    logical, intent(in) :: intermediateToFile

    !local
    type (modelObjectsType) :: modelObjects
    !local
    integer :: ixt, ixAIME, ixA, ixY, ixL, ixA1
    real (kind=rk) :: EV1(numPointsA,numAIME), EV1SPA(numPointsA,numAIME,numPointsSPA), EV2(numPointsA,numAIME)!negV, Va2, A,A1, Y, lbA1, ubA1, AIME, Agrid1(numPointsA)
    real (kind=rk) :: contEV(numPointsA, numAIME, numPointsSPA,numPointsY)
    real (kind=rk) :: AIME1grid(numAIME), realisedV(numPointsY), mat(numPointsSPA,numPointsL,numPointsA), q(numPointsL,numPointsA),va !, checkA(numPointsA), check,  policyA1temp, negVtemp,
    real (kind=rk) :: uConst(numPointsSPA,numPointsL,numPointsA) !, check2, checkA2(numPointsA)
    real (kind=rk) ::poilcyL1(numPointsA,numPointsY), val1(numPointsA,numPointsY) !, temp(numPointsA, numAIME, numPointsSPA, numPointsY,numPointsL,numPointsA), EVPath(stopwrok-1)
    !real (kind=rk) ::tempL(numPointsA, numAIME, numPointsY), poilcyA1(numPointsA,numPointsY), poilcyA1wrk(numPointsA,numPointsY) !,tempA(numPointsA, numAIME, numPointsY,numPointsA), poilcyA1dnt(numPointsA,numPointsY)
    real (kind=hp) :: mutualInfo(numPointsA,numPointsY), signal(numPointsSPA,numPointsA,numPointsY) ,totalP(numPointsA,numPointsY,numpointsA*numPointsL) !, ent(numPointsA,numPointsY), testSUM(5)
    integer:: indexBigN(2), indexSmalln(2), singleIndex, typeSim !, AssetChoice(numPointsType, numPointsA, numAIME, numPointsY), lchoice(numPointsType, numPointsA, numAIME, numPointsY)
    integer :: numTasks, tasksPerCore, leftOverTasks, thisCoreStart, thisCoreEnd, requiredl, spa, minSPA !, ixl2, ixA12, i,j
    character(len=1024) :: outFile !, outFile2
    !character (50) :: format_numpeepcols_int
#ifdef mpiBuild
    integer::  lcount, i !thisCoreStartStore, thisCoreEndStore,

    real (kind=rk), allocatable :: VecV(:), Vecpolicy(:), VecEV(:)
    real (kind=rk) :: tempV(mpiDim*numPointsY*numPointsSPA), temppolicy(mpiDim*numPointsSPA*numPointsY*numPointsL*numPointsA), tempEV(mpiDim*numPointsY*numPointsSPA)
    real (kind=rk), allocatable :: LocV(:), Locpolicy(:), LocEV(:)
    integer :: locVecSize, otherDimP, testRank, thisCoreStartTest, thisCoreEndTest, start, finish,  otherDimV  !, age
    integer(kind=4):: recvcounts(0:procSize-1), displ(0:procSize-1), mpiSIze

    !Test
    !real (kind=rk) :: testC

    allocate(VecV(mpiDim*numPointsY*numPointsSPA), Vecpolicy(mpiDim*numPointsSPA*numPointsY*numPointsL*numPointsA), VecEV(mpiDim*numPointsY*numPointsSPA))
#endif   

    !To large to be static object
    allocate(modelObjects%EV(numPointsA, numAIME, numPointsSPA,numPointsY))
    allocate(modelObjects%policy(numPointsA, numAIME, numPointsSPA, numPointsY,numPointsL,numPointsA))
    allocate(modelObjects%q(numPointsA, numAIME, numPointsY,numPointsL,numPointsA))
    allocate(modelObjects%u(numPointsA, numAIME, numPointsSPA, numPointsY,numPointsL,numPointsA))
    allocate(modelObjects%V(numPointsA, numAIME, numPointsSPA,numPointsY))

    !----------------------------------------------------------------------------------!
    ! Initialize MPI
    !----------------------------------------------------------------------------------!

    indexBigN(1) = numPointsA
    indexBigN(2) = numAIME
    numTasks = mpiDim

    tasksPerCore = int(numTasks/procSize)
    leftOverTasks = numTasks - tasksPerCore*procSize

    thisCoreStart = rank*tasksPerCore + 1 + max(rank + leftOverTasks - procsize,0)
    thisCoreEnd = (rank+1)*tasksPerCore + max(rank + leftOverTasks +1 - procsize,0)
#ifdef mpiBuild
    call mpi_barrier(mpi_comm_world, ierror)
    if (ierror.ne.0) stop 'mpi problem180'
#endif 

    !!Set the terminal value function and expected value function to 0
    !modelObjects%EV(Tperiods + 1,:,:,:,:,:)  = 0          ! continuation value at T-1
    !modelObjects%V(Tperiods + 1,:,:,:,:,:) = 0

    !Initialise everything
    modelObjects%EV  = 0.0          ! continuation value at T-1
    modelObjects%V = 0.0
    modelObjects%policy= 0.0
    !utlityShifter
    do typeSim = 1, numPointsType
        if (show .AND. rank==0) WRITE(*,*)  'Solving for type ', typeSim, ' of ',  numPointsType
        contEV = 0.0
        do ixt=Tperiods,1, -1                               ! Loop from time T-1 to 1
            AIME1grid = grids%AIMEgrid(typeSim,ixt + 1, :)
            modelObjects%EV  = 0.0          ! continuation value at T-1
            modelObjects%V = 0.0
            modelObjects%policy= 0.0
            poilcyL1 = 0.0
            val1 = 0.0
            totalP = 0.0
            mutualInfo = 0.0
            signal = 0.0
            totalP = 0.0
            do ixAIME = 1, numAIME
                do ixA = 1, numPointsA                   ! points on asset grid

                    indexSmalln(1) = ixa
                    indexSmalln(2) = ixAIME
                    singleIndex = sumindex(indexBigN, indexSmalln, .FALSE.)

                    if ((singleIndex .ge. thisCoreStart) .and. (singleIndex .le. thisCoreEnd)) then
                        if (ixt < stopwrok) then
                            !Although doesn't recieve income still need to loop around
                            !hypothetical income because participation doesn't effect
                            !earning potential
                            ! STEP 1. solve problem at grid points in assets, income + labour choices
                            ! ---------------------------------------------------------
                            do ixY = 1, numPointsY               ! points on income grid
                                if (ixt >= EndPeriodRI) then
                                    EV1  = contEV(:,:,1,ixY)  ! relevant section of EV matrix (in assets tomorrow)
                                    !include arbirtary SPA and copy results as you are definitely above it by this point and so don't care
                                    call solvePeriodRE(params, grids, ixy, ixA, ixAIME ,ixt,typeSim, Tretire-1+TrueSPA, EV1, ixA1,ixl,Va)
                                    modelObjects%policy(ixA,ixAIME, :,ixY,ixl,ixA1) = 1.0
                                    modelObjects%V( ixA, ixAIME,:,ixY) = Va
                                    !if (ixt==1) then
                                    !    AssetChoice(typeSim,ixA,ixAIME,ixY) = ixA1
                                    !    lChoice(typeSim,ixA,ixAIME,ixY) = ixl
                                    !end if
                                else
                                    if (uncerRE) then
                                        minSPA = maxval((/ixt-Tretire+2,1/))
                                        EV2 =  contEV(:,:,numPointsSPA,ixY)
                                        do spa = numPointsSPA,minSPA,-1
                                            EV1  = params%p*EV2+ (1-params%p)*contEV(:,:,spa,ixY)  ! relevant section of EV matrix (in assets tomorrow)
                                            !include arbirtary SPA and copy results as you are definitely above it by this point and so don't care
                                            call solvePeriodRE(params, grids, ixy, ixA, ixAIME ,ixt,typeSim, Tretire-1+spa, EV1, ixA1,ixl,Va)
                                            modelObjects%policy(ixA,ixAIME, spa,ixY,ixl,ixA1) = 1.0
                                            modelObjects%V( ixA, ixAIME,spa,ixY) = Va
                                            EV2 = contEV(:,:,spa,ixY)
                                        end do
                                        !if (spa >= 1 ) then
                                        !    modelObjects%policy(ixA,ixAIME, 1:spa,ixY,ixl,ixA1) = 1.0
                                        !    modelObjects%V( ixA, ixAIME,1:spa,ixY) = Va
                                        !end if
                                        if (ixt>=Tretire) then
                                            !solve once with SPA = current period
                                            !I know SPA and recieve it so continaution value is the same for all possible spa < ixt
                                            EV1  = contEV(:,:,1,ixY)! relevant section of EV matrix (in assets tomorrow)
                                            call solvePeriodRE(params, grids, ixy, ixA, ixAIME ,ixt,typeSim, ixt, EV1, ixA1,ixl,Va)
                                            modelObjects%policy(ixA,ixAIME, (/(spa,spa=1,ixt-Tretire+1)/),ixY,ixl,ixA1) = 1.0
                                            modelObjects%V(ixA, ixAIME,(/(spa,spa=1,ixt-Tretire+1)/),ixY) = Va
                                        end if
                                    else
                                        if (ixt == TendRI-1) then
                                            EV1  = contEV(:,:,1,ixY)  ! relevant section of EV matrix (in assets tomorrow)
                                            !If you SPA is max possible in the period before you know the SPA but haven't recieved it yet so it is ratonal
                                            !but you do care that you haven't yet received the SPA
                                            call solvePeriodRE(params, grids, ixy, ixA, ixAIME ,ixt,typeSim, TendRI, EV1, ixA1,ixl,Va)
                                            modelObjects%policy(ixA,ixAIME, numPointsSPA,ixY,ixl,ixA1) = 1.0
                                            !modelObjects%q(ixA,ixAIME, ixY,ixl,ixA1) = 1.0
                                            modelObjects%V( ixA, ixAIME,numPointsSPA,ixY) = Va
                                            !All other possible SPA are equivalent as you have recieved the SPA and so know what it is and don't care
                                            call solvePeriodRE(params, grids, ixy, ixA, ixAIME ,ixt,typeSim, TendRI-1, EV1, ixA1,ixl,Va)
                                            modelObjects%policy(ixA,ixAIME, (/(spa,spa=1,numPointsSPA-1)/),ixY,ixl,ixA1) = 1.0
                                            !modelObjects%q(ixA,ixAIME,ixY,ixl,ixA1) = 1.0
                                            modelObjects%V(ixA, ixAIME,(/(spa,spa=1,numPointsSPA-1)/),ixY) = Va
                                        else
                                            EV1SPA  = contEV(:,:,:,ixY)  ! relevant section of EV matrix (in assets tomorrow)
                                            !Set SPA to be TendRI as just need something greater than current age
                                            minSPA = maxval((/ixt-Tretire+2,1/))
                                            mat=0.0
                                            call solvePeriodRI(params, grids, ixy, ixA, ixAIME ,ixt,typeSim, TendRI, EV1SPA, &
                                                mat(minSPA:numPointsSPA,:,:) ,modelObjects%V(ixA, ixAIME,minSPA:numPointsSPA,ixY), q, uConst )
                                            modelObjects%policy(ixA,ixAIME, :,ixY,:,:) = mat
                                            modelObjects%q(ixA,ixAIME, ixY,:,:) = q
                                            modelObjects%u(ixA,ixAIME, :,ixY,:,:) = uConst
                                            !if (ixy==5 .AND. ixa ==6 .AND. ixAIME==1) then
                                            !    call solvePeriodRE(params, grids, ixy, ixA, ixAIME ,ixt,typeSim, TendRI, EV1SPA(:,:,numPointsSPA), ixA1,ixl,Va)
                                            !    check = sum(mat(numPointsSPA,2,:))!sum(modelObjects%policy(ixt,typeSim,ixA,ixAIME, 11,ixY,:,:))
                                            !    write (*,*) check, ixl
                                            !end if
                                            !call solvePeriodRE(params, grids, ixy, ixA, ixAIME ,ixt,typeSim, TendRI, EV1SPA(:,:,numPointsSPA), ixA1,ixl,Va)
                                            !checkA = sum(mat(numPointsSPA,:,:),1)
                                            !!check2 = sum(mat(numPointsSPA-1,2,:))!sum(modelObjects%policy(ixt,typeSim,ixA,ixAIME, 11,ixY,:,:))
                                            !!checkA2 = sum(mat(numPointsSPA-1,:,:),1)
                                            if (ixt>=Tretire) then
                                                !solve once with SPA = current period
                                                !I know SPA and recieve it so continaution value is the same for all possible spa < ixt
                                                EV1  = contEV(:,:,1,ixY)! relevant section of EV matrix (in assets tomorrow)
                                                call solvePeriodRE(params, grids, ixy, ixA, ixAIME ,ixt,typeSim, ixt, EV1, ixA1,ixl,Va)
                                                modelObjects%policy(ixA,ixAIME, (/(spa,spa=1,ixt-Tretire+1)/),ixY,ixl,ixA1) = 1.0
                                                ! modelObjects%q(ixA,ixAIME, ixY,ixl,ixA1) = 1.0
                                                modelObjects%V(ixA, ixAIME,(/(spa,spa=1,ixt-Tretire+1)/),ixY) = Va
                                            end if
                                        end if
                                    end if
                                end if
                            end do
                        else
                            !Can't work at this stage and you are past any possible SPA so pass in arbitrary value for ixy (= 1) and SPA9=60) and
                            !hardcode ixl=0(not work)
                            call maxutility(params, grids,ixt,typeSim,ixa,ixAIME,0,  1,Tretire,  contEV(:,:,1,1),ixa1,va)

                            !Copy to all values of SPA and ixy as don't affect decsion at these ages
                            modelObjects%policy(ixA,ixAIME,:,:,1,ixA1)=1.0
                            modelObjects%V(ixA,ixAIME,:,:)       = va
                        end if

                        ! STEP 2. integrate out income today conditional on income
                        ! yesterday to get EV and EdU
                        ! --------------------------------------------------------
                        if (ixt >= EndPeriodRI) then
                            realisedV(:) = modelObjects%V(ixA,ixAIME,1,:)
                            do ixY = 1,numPointsY,1
                                modelObjects%EV(ixA,ixAIME,1,ixY)  = dot_product(grids%incTransitionMrx(ixY,:),realisedV)
                                modelObjects%EV(ixA,ixAIME,:,ixY) = modelObjects%EV(ixA,ixAIME,1,ixY)
                            end do !ixY
                        else
                            do spa = 1,numPointsSPA
                                realisedV(:) = modelObjects%V(ixA,ixAIME,spa,:)
                                do ixY = 1,numPointsY,1
                                    modelObjects%EV(ixA,ixAIME,spa,ixY)  = dot_product(grids%incTransitionMrx(ixY,:),realisedV)
                                    if (modelObjects%EV(ixA,ixAIME,spa,ixY) /= modelObjects%EV(ixA,ixAIME,spa,ixY)) then
                                        continue
                                    end if
                                end do !ixY
                            end do
                        end if
                    end if
                end do
            end do

#ifdef mpiBuild 
            call mpi_barrier(mpi_comm_world, ierror)
            if (ierror.ne.0) stop 'mpi problem180'

            locVecSize = thisCoreEnd - thisCoreStart + 1
            otherDimP = numPointsSPA*numPointsY*numPointsL*numPointsA
            otherDimV = numPointsSPA*numPointsY

            allocate(Locpolicy(locVecSize*otherDimP),LocEV(locVecSize*otherDimV),locV(locVecSize*otherDimV))
            call unpackArrays(modelObjects%policy(:, :,:, :,:,:), modelObjects%V(:, :, :,:), modelObjects%EV(:, :,:, :), &
                locpolicy, locV, locEV, mpiDim,otherDimP,otherDimV,thisCoreStart,thisCoreEnd)

            recvcounts(0:procsize - leftOverTasks-1) = tasksPerCore*otherDimP
            recvcounts(procsize - leftOverTasks:procsize-1) = (tasksPerCore+1)*otherDimP
            displ = (/((i*tasksPerCore  + max(i + leftOverTasks - procsize,0))*otherDimP,i=0,procSize-1)/)
            mpiSIze = locVecSize*otherDimP

            call mpi_gatherV(Locpolicy, mpiSIze, mpi_double_precision, Vecpolicy,  recvcounts, displ, mpi_double_precision, 0, mpi_comm_world, ierror)

            recvcounts(0:procsize - leftOverTasks-1) = tasksPerCore*otherDimV
            call mpi_barrier(mpi_comm_world, ierror)

            recvcounts(procsize - leftOverTasks:procsize-1) = (tasksPerCore+1)*otherDimV
            displ = (/((i*tasksPerCore  + max(i + leftOverTasks - procsize,0))*otherDimV,i=0,procSize-1)/)
            mpiSIze = locVecSize*otherDimV
            call mpi_barrier(mpi_comm_world, ierror)

            call mpi_allgatherV(LocV, mpiSIze, mpi_double_precision, VecV, recvcounts, displ, mpi_double_precision, mpi_comm_world, ierror)
            call mpi_allgatherV(LocEV, mpiSIze, mpi_double_precision, VecEV, recvcounts, displ, mpi_double_precision, mpi_comm_world, ierror)

            deallocate(LocV,Locpolicy,LocEV)
            tempV = VecV
            temppolicy = Vecpolicy
            tempEV = VecEV
            do testRank = 0, (procSize-1)
                thisCoreStartTest = testrank*tasksPerCore + 1  + max(testrank + leftOverTasks - procsize,0)
                thisCoreEndTest = (testrank+1)*tasksPerCore + max(testrank + leftOverTasks +1 - procsize,0)
                locVecSize = thisCoreEndTest - thisCoreStartTest + 1
                allocate(LocV(locVecSize*otherDimV),Locpolicy(locVecSize*otherDimP),LocEV(locVecSize*otherDimV))

                start = (thisCoreStartTest-1)*otherDimP+1
                finish = thisCoreEndTest*otherDimP

                locpolicy = temppolicy(start:finish)

                !Distribute the contigous section between the contigous section in each column corresponding to the non-splitting dimesnions of the array
                do lcount=1,otherDimP
                    start = (lcount-1)*mpiDim+thisCoreStartTest
                    finish = (lcount-1)*mpiDim+thisCoreEndTest
                    Vecpolicy(start:finish) = locpolicy((lcount-1)*locVecSize+1:(lcount-1)*locVecSize+locVecSize)
                end do

                start = (thisCoreStartTest-1)*otherDimV+1
                finish = thisCoreEndTest*otherDimV

                locV= tempV(start:finish)
                locEV = tempEV(start:finish)

                !Distribute the contigous section between the contigous section in each column corresponding to the non-splitting dimesnions of the array
                do lcount=1,otherDimV
                    start = (lcount-1)*mpiDim+thisCoreStartTest
                    finish = (lcount-1)*mpiDim+thisCoreEndTest
                    VecV(start:finish) = locV((lcount-1)*locVecSize+1:(lcount-1)*locVecSize+locVecSize)
                    VecEV(start:finish) = locEV((lcount-1)*locVecSize+1:(lcount-1)*locVecSize+locVecSize)
                end do

                deallocate(LocV,Locpolicy,LocEV)
            end do!
            modelObjects%V( :, :, :,:) = reshape(VecV, (/numPointsA, numAIME, numPointsSPA,numPointsY/))
            modelObjects%policy( :, :, :,:,:,:) = reshape(Vecpolicy, (/numPointsA, numAIME,numPointsSPA, numPointsY,numPointsL,numPointsA/))
            modelObjects%EV( :, :, :,:) = reshape(VecEV, (/numPointsA, numAIME, numPointsSPA,numPointsY/))

            call mpi_barrier(mpi_comm_world, ierror)
            if (ierror.ne.0) stop 'mpi problem180'
#endif  
            contEV =  modelObjects%EV
            if (show .AND. rank==0) WRITE(*,*)  'Passed period ', ixt, ' of ',  Tperiods

            if (rank==0) then
                write (outFile, *), trim(trim(pathDataStore) // "polType"),typeSim,trim("Period"),ixt,".txt"
                outfile=ADJUSTL(outfile)
                inquire (iolength=requiredl)  modelObjects%policy
                open (unit=201,form="unformatted", file=outfile, status='unknown',recl=requiredl, action='write')
                write (201)  modelObjects%policy
                close( unit=201)

                write (outFile, *), trim(trim(pathDataStore) // "EVType"),typeSim,trim("Period"),ixt,".txt"
                outfile=ADJUSTL(outfile)
                inquire (iolength=requiredl)  modelObjects%EV
                open (unit=201,form="unformatted", file=outfile, status='unknown',recl=requiredl, action='write')
                write (201)  modelObjects%EV
                close( unit=201)

                write (outFile, *), trim(trim(pathDataStore) // "qType"),typeSim,trim("Period"),ixt,".txt"
                outfile=ADJUSTL(outfile)
                inquire (iolength=requiredl)  modelObjects%q
                open (unit=201,form="unformatted", file=outfile, status='unknown',recl=requiredl, action='write')
                write (201)  modelObjects%q
                close( unit=201)

                write (outFile, *), trim(trim(pathDataStore) // "uType"),typeSim,trim("Period"),ixt,".txt"
                outfile=ADJUSTL(outfile)
                inquire (iolength=requiredl)  modelObjects%u
                open (unit=201,form="unformatted", file=outfile, status='unknown',recl=requiredl, action='write')
                write (201)  modelObjects%u
                close( unit=201)

                if (ixt < stopwrok .AND. intermediateToFile) then
                    !    write (format_numpeepcols_int,*),'(',numPointsY
                    !    format_numpeepcols_int = trim(format_numpeepcols_int) // 'E20.7E4)'
                    !    select case(modelChoice)
                    !    case(1)
                    !        path = "C:\Users\Uctphen\Dropbox\SourceCode\upgradeProject\VSProj - Copy\tempBaseline\"
                    !    case(2)
                    !        path = "C:\Users\Uctphen\Dropbox\SourceCode\upgradeProject\VSProj - Copy\tempRE\"
                    !        case default
                    !        path = "C:\Users\Uctphen\Dropbox\SourceCode\upgradeProject\VSProj - Copy\temp\"
                    !    end select
                    !
                    !    EVPath(ixt) = sum(modelObjects%EV(:,:, 11,:))
                    !
                    !    do spa = 1, numPointsSPA
                    !        !Policy L
                    !        poilcyL1 = sum(modelObjects%policy(:,4, spa,:,2,:),3)
                    !        write (outFile, *), trim(trim(path) // "policyL"),typeSim,trim("Period"),ixt,trim("SPA"),Tretire+startAge-2+spa,".txt"
                    !        outfile=ADJUSTL(outfile)
                    !        inquire (iolength=requiredl)  transpose(poilcyL1)
                    !        open (unit=201, file=outfile, status='unknown',recl=requiredl, action='write')
                    !        write (201,format_numpeepcols_int)  transpose(poilcyL1)
                    !        close( unit=201)
                    !
                    !        !Policy A
                    !        do i =1, numpointsA !grids%Agrid(typesim,1,:))
                    !            do j=1,numpointsY
                    !                poilcyA1(i,j) = dot_product(sum(modelObjects%policy(i,4, spa,j,:,:),1),grids%Agrid(typesim,1,:))
                    !                if (poilcyL1(i,j) > 0) then
                    !                    poilcyA1wrk(i,j) = dot_product(modelObjects%policy(i,4, spa,j,2,:)/poilcyL1(i,j),grids%Agrid(typesim,1,:))
                    !                else
                    !                    poilcyA1wrk(i,j) = 0.0
                    !                end if
                    !                if (poilcyL1(i,j) <1) then
                    !                    poilcyA1dnt(i,j) = dot_product(modelObjects%policy(i,4, spa,j,1,:)/(1-poilcyL1(i,j)),grids%Agrid(typesim,1,:))
                    !                else
                    !                    poilcyA1dnt(i,j) = 0.0
                    !                end if
                    !                if (abs(poilcyA1(i,j) - ((1-poilcyL1(i,j))*poilcyA1dnt(i,j)+poilcyL1(i,j)*poilcyA1wrk(i,j)))>1.0) then
                    !                    write (*,*) "Prob error"
                    !                    stop
                    !                end if
                    !            end do
                    !        end do
                    !        write (outFile, *), trim(trim(path) // "policyA"),typeSim,trim("Period"),ixt,trim("SPA"),Tretire+startAge-2+spa,".txt"
                    !        outfile=ADJUSTL(outfile)
                    !        inquire (iolength=requiredl)  transpose(poilcyA1)
                    !        open (unit=201, file=outfile, status='unknown',recl=requiredl, action='write')
                    !        write (201,format_numpeepcols_int)  transpose(poilcyA1)
                    !        close( unit=201)
                    !
                    !        write (outFile, *), trim(trim(path) // "policyAwrk"),typeSim,trim("Period"),ixt,trim("SPA"),Tretire+startAge-2+spa,".txt"
                    !        outfile=ADJUSTL(outfile)
                    !        inquire (iolength=requiredl)  transpose(poilcyA1wrk)
                    !        open (unit=201, file=outfile, status='unknown',recl=requiredl, action='write')
                    !        write (201,format_numpeepcols_int)  transpose(poilcyA1wrk)
                    !        close( unit=201)
                    !
                    !        write (outFile, *), trim(trim(path) // "policyAdnt"),typeSim,trim("Period"),ixt,trim("SPA"),Tretire+startAge-2+spa,".txt"
                    !        outfile=ADJUSTL(outfile)
                    !        inquire (iolength=requiredl)  transpose(poilcyA1dnt)
                    !        open (unit=201, file=outfile, status='unknown',recl=requiredl, action='write')
                    !        write (201,format_numpeepcols_int)  transpose(poilcyA1dnt)
                    !        close( unit=201)
                    !
                    !        !Value
                    !        val1 = modelObjects%v(:,4, spa,:)
                    !        write (outFile, *), trim(trim(path) // "Val"),typeSim,trim("Period"),ixt,trim("SPA"),Tretire+startAge-2+spa,".txt"
                    !        outfile=ADJUSTL(outfile)
                    !        inquire (iolength=requiredl)  transpose(val1 )
                    !        open (unit=201, file=outfile, status='unknown',recl=requiredl, action='write')
                    !        write (201,format_numpeepcols_int)  transpose(val1 )
                    !        close( unit=201)
                    !
                    !        !Value
                    !        !val1 = modelObjects%ev(:,4, spa,:)
                    !        !write (outFile, *), trim(trim(path) // "EVal"),typeSim,trim("Period"),ixt,trim("SPA"),Tretire+startAge-2+spa,".txt"
                    !        !outfile=ADJUSTL(outfile)
                    !        !inquire (iolength=requiredl)  transpose(val1 )
                    !        !open (unit=201, file=outfile, status='unknown',recl=requiredl, action='write')
                    !        !write (201,format_numpeepcols_int)  transpose(val1 )
                    !        !close( unit=201)
                    !
                    !
                    !        if ( modelChoice==3 .AND. ixt <TendRI-1) then
                    !            FORALL(I = 1:numpointsA, J = 1:numpointsY) signal(spa,i,j) =  entropy(real(reshape(modelObjects%policy(i,4, spa,j,:,:),(/numpointsA*numPointsL/)),hp))
                    !            FORALL(I = 1:numpointsA, J = 1:numpointsY) totalP(i,j,:) = totalP(i,j,:) + grids%posteriorSPA(ixt,spa)*reshape(modelObjects%policy(i,4, spa,j,:,:),(/numpointsA*numPointsL/))
                    !
                    !            if (ixt==8 ) then
                    !
                    !                write (outFile, *), trim("C:\Users\Uctphen\Dropbox\SourceCode\upgradeProject\VSProj - Copy\Entropy\Pol"),trim("SPA"),Tretire+startAge-2+spa,".txt"
                    !                outfile=ADJUSTL(outfile)
                    !                inquire (iolength=requiredl)  modelObjects%policy(5,4, spa,5,:,:)
                    !                open (unit=201, file=outfile, status='unknown',recl=requiredl, action='write')
                    !                write (201,*)  modelObjects%policy(5,4, spa,5,:,:)
                    !                close( unit=201)
                    !
                    !            end if
                    !        end if
                    !    end do
                end if
                if ( modelChoice==3 .AND. ixt <TendRI-1 .AND. intermediateToFile) then
                    !FORALL(I = 1:numpointsA, J = 1:numpointsY)  mutualInfo(i,j) = entropy( totalP(i,j,:)) - dot_product(grids%posteriorSPA(ixt,:),signal(:,i,j))
                    !if (minval(mutualInfo) < -1.0E-16 ) then
                    !    testSUM(1) = sum(grids%posteriorSPA(ixt,:))
                    !    testSUM(2) = sum(modelObjects%policy(8,4, 11,1,:,:))
                    !    testSUM(3) = sum(modelObjects%policy(8,4, 10,1,:,:))
                    !    testSUM(4) = sum(modelObjects%policy(8,4, 9,1,:,:))
                    !    testSUM(5) = sum(totalP(8,1,:))
                    !    FORALL(I = 1:numpointsA, J = 1:numpointsY)  ent(i,j) = entropy( totalP(i,j,:))
                    !    write (*,*) "Neg. mutaul info!", minval(mutualInfo)
                    !end if
                    !write (outFile, *), trim(trim(path) // "MutualInfo"),typeSim,trim("Period"),ixt,".txt"
                    !outfile=ADJUSTL(outfile)
                    !inquire (iolength=requiredl)  transpose(mutualInfo )
                    !open (unit=201, file=outfile, status='unknown',recl=requiredl, action='write')
                    !write (201,format_numpeepcols_int)  transpose(mutualInfo )
                    !close( unit=201)
                end if
            end if


            if (rank==0 .AND. intermediateToFile .AND. ixt <10) then
                !!L
                !write (outFile, *), trim("..\\temp\polL"),ixt,".txt"
                !!check = sum(mat(numPointsSPA,2,:))!sum(modelObjects%policy(ixt,typeSim,ixA,ixAIME, 11,ixY,:,:))
                !outfile=ADJUSTL(outfile)
                !inquire (iolength=requiredl) tempL
                !open (unit=201, file=outfile, status='unknown',recl=requiredl, action='write')
                !write (201, *) tempL
                !close( unit=201)
                !
                !write (outFile, *), trim("..\\temp\polA"),ixt,".txt"
                !!check = sum(mat(numPointsSPA,2,:))!sum(modelObjects%policy(ixt,typeSim,ixA,ixAIME, 11,ixY,:,:))
                !outfile=ADJUSTL(outfile)
                !inquire (iolength=requiredl) tempA
                !open (unit=201, file=outfile, status='unknown',recl=requiredl, action='write')
                !write (201, *) tempA
                !close( unit=201)

            end if
        end do !ixt

        if (rank==0 .AND. intermediateToFile) then
            !write (outFile, *), trim(trim(path) // "EVpath"),typeSim,".txt"
            !outfile=ADJUSTL(outfile)
            !inquire (iolength=requiredl) Evpath
            !open (unit=201, file=outfile, status='unknown',recl=requiredl, action='write')
            !write (201,'(F20.10)')  EVpath
            !close( unit=201)
        end if
    end do !type
#ifdef mpiBuild
    call mpi_barrier(mpi_comm_world, ierror)
    if (ierror.ne.0) stop 'mpi problem180'
    deallocate(VecV, Vecpolicy, VecEV)
#endif 

    !write (*,*) "Done"

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

    !ouptut
    real (kind=rk) :: objectivefunc
    !local
    real (kind=rk) :: cons, VA1, VB1

    !Get tomorrow's consumption (cons), the value of left over assets (VA1) and
    !total value (u(c) + b * VA1
    cons = A0  + Y - (A1)/(1+params%r)

    !!THIS IS WRONG SHOULD BE NEXT PERIODS AIME NOT THIS PERIODS
    call  linearinterp1(grids%AIMEgrid(ixType,ixP + 1, :),EV1, numAIME, AIME, VA1, 1, 1 )
    VB1 = params%thetab*((A1+params%K)**(1-params%gamma))/(1-params%gamma)

    objectivefunc = utility(params,cons,L)/params%lambda + params%beta * ((1- grids%mortal(ixP))* VA1+grids%mortal(ixP+1)*VB1)


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
    real (kind=rk) :: cevChng

    !outpus
    real (kind=rk) :: utility, les

    if (cons<=0) then
        print *, 'Error in utility! Consumption is LE 0'
        stop
    end if
    cevChng = (1.0 + params%percentCons)
    !10/112 comuniting time -10/112
    les=(L)*(1-params%hrsWrk )+(1-L);
    if (params%gamma == 1) then
        utility = log((cons*cevChng)**params%nu*les**(1-params%nu));
    else
        utility= ((((cons*cevChng)**params%nu)*(les**(1-params%nu)))**(1-params%gamma)  )/(1-params%gamma);
    end if

    end function

    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    !!Defined Benefite Pension function
    function dbPension(params,AIME, ixType)
    implicit none
    !inputs
    type (structparamstype), intent(in) :: params
    real (kind=rk), intent(in) :: AIME
    integer, intent(in) :: ixType
    !outputs
    real (kind=rk) :: dbPension
    !local
    real (kind=rk):: bound
    bound = -params%db(ixType,1) /(2.0*params%db(ixType,2))
    if (AIME < bound) then
        dbPension = params%db(ixType,2)*AIME**2 + params%db(ixType,1)*AIME;
    else
        dbPension = params%db(ixType,2)*bound**2 + params%db(ixType,1)*bound
    endif
    dbPension= weeksyear*dbPension !*0.0 !temPchange
    end function
    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    !!State Pension function
    function statePension(params,AIME, ixType)
    implicit none
    !inputs
    type (structparamstype), intent(in) :: params
    real (kind=rk), intent(in) :: AIME
    integer, intent(in) :: ixType
    !outputs
    real (kind=rk) :: statePension
    !local
    real (kind=rk):: bound
    bound = -params%pension(ixType,1) /(2.0*params%pension(ixType,2))
    if (AIME < bound) then
        statePension = params%pension(ixType,2)*AIME**2 + params%pension(ixType,1)*AIME;
    else
        statePension = params%pension(ixType,2)*bound**2 + params%pension(ixType,1)*bound
    endif
    statePension = weeksyear*statePension
    !statePension = weeksyear !*100.0

    end function
    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    !!Simulation Subroutine v,
    subroutine simWithUncer(params, grids,  yex, c, a,  l, y, AIME, SPAin, delCache, error, posOut, prioOut )
    use Header
    implicit none

    !inputs
    type (structparamstype), intent(in) :: params
    type (gridsType), intent(inout) :: grids
    integer, intent(in) :: SPAin
    logical, intent(in) :: delCache

    real (kind=rk), intent(in) :: yex(Tperiods, numSims)

    !outputs
    real (kind=rk), intent(out) :: y(Tperiods, numSims) !income
    real (kind=rk), intent(out) :: c(Tperiods, numSims)  !consumption
    integer, intent(out) :: l(Tperiods, numSims) !labour supply
    !real (kind=rk), intent(out) :: v(Tperiods, numSims)  !value
    integer, intent(out) :: a(Tperiods + 1,numSims) !this is the path at the start of each period, so we include the 'start' of death
    real (kind=rk), intent(out) :: AIME(Tperiods + 1,numSims)
    real (kind=rk), intent(out), optional :: error

    integer :: pos(7,numSims,1), prio(7,numSims,1)
    integer, intent(out), optional :: posOut(7,numSims,1), prioOut(7,numSims,1)

    !local
    type (modelObjectsType) :: modelObjects
    real (kind=rk) :: startingA(numSims), startAIME(numSims), check, check1, check2, checkMax !, Aval, AIMEval
    integer :: s, t, idxAIME(1), idxY(numSims) !, workAge, idxA(1)
    integer :: typeSim !,seedIn, Lcube(8)
    !real (kind=rk) ::  ltemp, lbA1, EV1(numPointsA ,numAIME)
    INTEGER :: i, j !, n, uniformInt(numSims)
    !INTEGER, DIMENSION(:), ALLOCATABLE :: seed
    integer :: unemployed(numSims) !1 employed, 2 unemploeyed
    !type (locmodelObjectsType):: LocmodelObjects(2)
    type (modelObjectsType):: LocmodelObjects(2)
    !real (kind=rk) :: LocEV(Tperiods+1,  numPointsType, numPointsA, numAIME, numPointsSPA,numPointsProd,2);
    real (kind=rk) :: gridY(numPointsProd), Yval, ccp(numPointsL*numPointsA), policy(numAIME, numpointsprod, numPointsL*numPointsA), ccpTemp(numPointsL*numPointsA),  ccpSPA(numPointsL*numPointsA,numPointsSPA)
    real (kind=rk) :: mat(numPointsL,numPointsA), cumsum, posterior(Tperiods + 1,numSims, numPointsSPA), prior(Tperiods + 1, numSims, numPointsSPA)!, denom, q, unifposterior( numPointsSPA)
    integer :: idxa1, typesimOld, lcount, decision
    character(len=1024) :: outFile
    integer :: SPA, requiredl, spaResponse(434) !, mode(1)
    character (20) :: format_numpeepcols_int
    real (kind=rk) :: sppath(numSims), dbpath(numsims), q, weight = 0.8

    real (kind=rk), allocatable :: locposteriorSPA(:,:,:, :,  :,:)

    integer :: maxVec(5), pointState(TendRI, numSims,5), statesVisited(TendRI,numPointsType,numPointsA, numAIME,  numPointsY) = 0
    integer :: aimeInt(Tperiods,numSims), yint(Tperiods, numSims)

    logical :: maks(TendRI,numPointsType,numPointsA, numAIME, numPointsY)

    SPA = SPAin
    !call random_number(uniformrand)
    do lcount = 1, numsims
        i = mod(lcount-1,obsInitDist) + 1
        startinga(lcount) = grids%initialassets(i,2) !uniformrand(i)*numsims)+1
        !if (startinga(lcount) < 0.0) startinga(lcount) = 0.0
        startaime(lcount) = grids%initialassets(i,1)
    end do
    !startaime = 7351 !15500 !

    !To large to be static object
    allocate(modelObjects%EV(numPointsA, numAIME, numPointsSPA,numPointsY))
    allocate(modelObjects%policy(numPointsA, numAIME, numPointsSPA, numPointsY,numPointsL,numPointsA))
    allocate(modelObjects%u(numPointsA, numAIME, numPointsSPA, numPointsY,numPointsL,numPointsA))
    allocate(modelObjects%q(numPointsA, numAIME, numPointsY,numPointsL,numPointsA))
    allocate(modelObjects%V(numPointsA, numAIME, numPointsSPA,numPointsY))
    allocate(locposteriorSPA(TendRI,numPointsType,numPointsA, numAIME,  numPointsY,numPointsSPA))
    !open (unit = 1001,file=trim(pathMoments) // 'prior.txt', action='read')
    !read (1001, *) spaResponse
    !close (unit = 1001)
    !spaResponse = spaResponse - 59
    !do lcount = 1, numsims
    !    i = mod(lcount-1,434) + 1
    !    prior(1,lcount,:) = grids%posteriorSPA(1,1,1,1,1,:)
    !    if (.NOT. spaResponse(i) == 3) then
    !        prior(1,lcount,spaResponse(i)) = prior(1,lcount,spaResponse(i)) +  0.1 + (prior(1,lcount,3)- prior(1,lcount,spaResponse(i)))
    !        prior(1,lcount,:) = prior(1,lcount,:) / sum(prior(1,lcount,:))
    !    end if
    !end do
    checkMax = 0.0
    typesimOld = 4
    !write (*,*) 'b'
    s=0
    !unifposterior( :) = grids%posteriorSPA(1,:)

    do t = 1,tperiods,1                              ! loop through time periods for a particular individual
        if (counterFact == .TRUE. .AND. t==9) then
            !SPA = SPA + 1
        end if
        !mode = MAXLOC(unifposterior)
        !write (*,*) 'Uniformed mode', mode(1)
        !unifposterior = matmul(params%spaTransMat,unifposterior)
        do s = 1,numsims,1
            !Initialise beliefs
            if (modelChoice > 2 .AND. t==1) prior(1,s,:) = grids%posteriorSPA(1,1,1, 1,1,:) !grids%posteriorSPA(1,:)
            !write (*,*) "sim ", s
            if (numpointstype == 1) then
                typesim = 1
            else
                if (real(s,rk)/numsims < params%fracType(1)) then
                    typesim = 1
                else if (real(s,rk)/numsims < params%fracType(2)) then
                    typesim = 2
                else if (real(s,rk)/numsims < params%fracType(3)) then
                    typesim =3
                else
                    typesim =4
                end if
            end if
            if (t==1) then
                a(1, s) = minloc(abs(startinga(s)-grids%Agrid(typesim,1,:)),1)
                aime(1,s)=startaime(s)
            end if
            !write (*,*) 'c'
            if (typesimOld .NE. typesim) then
                write (outFile, *), trim(trim(pathDataStore) // "polType"),typeSim,trim("Period"),t,".txt"
                outfile=ADJUSTL(outfile)
                !write (*,*) 'd', outFile
                open (unit=201,form="unformatted", file=outfile, status='unknown', action='read')
                read (201) modelObjects%policy
                if (delCache) then
                    close( unit=201,  status='delete' )
                else
                    close( unit=201 )
                end if

                write (outFile, *), trim(trim(pathDataStore) // "qType"),typeSim,trim("Period"),t,".txt"
                outfile=ADJUSTL(outfile)
                !write (*,*) 'd', outFile
                open (unit=201,form="unformatted", file=outfile, status='unknown', action='read')
                read (201) modelObjects%q

                !i=19
                write (outFile, *), trim(trim(pathDataStore) // "uType"),typeSim,trim("Period"),t,".txt"
                outfile=ADJUSTL(outfile)
                !write (*,*) 'd', outFile
                open (unit=201,form="unformatted", file=outfile, status='unknown', action='read')
                read (201) modelObjects%u
                !stop

                if (delCache) then
                    close( unit=201,  status='delete' )
                else
                    close( unit=201 )
                end if

                write (outFile, *), trim(trim(pathDataStore) // "EVType"),typeSim,trim("Period"),t,".txt"
                outfile=ADJUSTL(outfile)
                !write (*,*) 'e', outFile
                open (unit=201,form="unformatted", file=outfile, status='unknown',action='read')
                read (201) modelObjects%EV
                if (delCache) then
                    close( unit=201,  status='delete' )
                else
                    close( unit=201 )
                end if

                !write (*,*) 'f'
                LocmodelObjects(1)%policy =  modelObjects%policy(:,:,:,(/(i,i=1,numPointsY-1,2)/),:,:)
                !write (*,*) 'g'
                LocmodelObjects(2)%policy =  modelObjects%policy(:,:,:,(/(i,i=2,numPointsY,2)/),:,:)

                LocmodelObjects(1)%u =  modelObjects%u(:,:,:,(/(i,i=1,numPointsY-1,2)/),:,:)
                !write (*,*) 'g'
                LocmodelObjects(2)%u =  modelObjects%u(:,:,:,(/(i,i=2,numPointsY,2)/),:,:)

                LocmodelObjects(1)%q =  modelObjects%q(:,:,(/(i,i=1,numPointsY-1,2)/),:,:)
                !write (*,*) 'g'
                LocmodelObjects(2)%q =  modelObjects%q(:,:,(/(i,i=2,numPointsY,2)/),:,:)

                LocmodelObjects(1)%EV =  modelObjects%EV(:,:,:,(/(i,i=1,numPointsY-1,2)/))
                LocmodelObjects(2)%EV =  modelObjects%EV(:,:,:,(/(i,i=2,numPointsY,2)/))

            end if
            !write (*,*) 'h'
            typesimOld = typesim
            !write (*,*) 'f'
            if (yex(t,s) <= 0.0001 ) then
                unemployed(s) = 2
#if PROD_SIZE == 5
                gridy = (/1,2,3,4,5/)
#elif  PROD_SIZE == 10                   
                gridy = (/1,2,3,4,5,6,7,8,9,10/)
#endif
                yval = idxy(s)
            else
                unemployed(s) = 1
                !everyone start employed so this is ok
                idxy(s)=minloc(abs(yex(t,s)-grids%ygridextrap(typesim,t,:)), DIM = 1)
                gridy =grids%ygrid(typesim,t,(/(i,i=1,numpointsy-1,2)/))
                yval = yex(t,s)
                y(t, s) =  yex(t,s)
            end if

            idxaime=minloc(abs(aime(t, s)-grids%aimegrid(typesim,t,:)))
            do i=1,numPointsL
                do j=1,numPointsA
                    policy(:,:,j+(i-1)*numPointsA) = LocmodelObjects(unemployed(s))%policy(a(t, s) ,:,SPA,:,i,j)
                end do
            end do
            ccp =  reshape(LocmodelObjects(unemployed(s))%policy(a(t, s),idxaime(1),SPA,idxy(s),:,:),(/numPointsL*numPointsA/))
            check1 = abs(sum(ccp))
            !y(t, s) = yex(t,s)

            if ( t < Tretire - 1 +SPA) then
                statesVisited(t,typesim,a(t, s),idxaime(1),idxy(s)) =  statesVisited(t,typesim,a(t, s),idxaime(1),idxy(s)) + 1
                pointState(t, s,:) = (/t,typesim,a(t, s),idxaime(1),idxy(s)/)
                aimeInt(t,s) = idxaime(1)
                yint(t, s) = idxy(s)
            end if

            call linearinterp2_vec(grids%aimegrid(typesim,t,:),gridy, numaime, numpointsprod,(numPointsL*numPointsA), aime(t, s), yval, ccpTemp, policy)
            ccpTemp = ccpTemp/sum(ccpTemp)
            check2 = abs(sum(ccpTemp))
            check = abs(sum(ccp - ccpTemp))
            if (check > checkMax) checkMax = check
            if (abs(sum(ccp - ccpTemp)) > 0.00001) then
                write (*,*) "Larege Error", abs(sum(ccp - ccpTemp))
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
                    if (cumsum > grids%logitShock(s,t)) exit
                end do
                decision = i
                l(t,s) = 1- mod(i,2)
                idxa1= (i-1)/2+1
                a(t+1, s) = idxa1
                mat = LocmodelObjects(unemployed(s))%policy(a(t, s),idxaime(1),SPA,idxy(s),:,:)
                if (modelChoice>2 .and. t < Tretire - 1 +SPA ) then
                    !ccp(i)
                    !q = dot_product(LocmodelObjects(unemployed(s))%policy(a(t, s),idxaime(1),numPointsSPA-grids%supportSPA(t)+1:numPointsSPA,idxy(s),l(t, s)+1,a(t+1, s)),grids%posteriorSPA(t,numPointsSPA-grids%supportSPA(t)+1:numPointsSPA))
                    posterior(t,s,:) = 0.0
                    do i = 1, numPointsSPA-grids%supportSPA(t)
                        if (i >=SPA) posterior(t,s,SPA) = 1.0
                    end do
                    !if (LocmodelObjects(unemployed(s))%q(a(t, s),idxaime(1),idxy(s),l(t, s)+1,a(t+1, s))<0.95) then
                    !    continue
                    !end if

                    do i =numPointsSPA-grids%supportSPA(t)+1,numPointsSPA
                        ccpSPA(:,i) =  reshape(LocmodelObjects(unemployed(s))%policy(a(t, s),idxaime(1),i,idxy(s),:,:),(/numPointsL*numPointsA/))
                    end do
                    q = dot_product(ccpSPA(decision,:),prior(t,s,:))

                    do i =numPointsSPA-grids%supportSPA(t)+1,numPointsSPA
                        posterior(t,s,i)=prior(t,s,i)*ccpSPA(decision,i)/ q! LocmodelObjects(unemployed(s))%q(a(t, s),idxaime(1),idxy(s),l(t, s)+1,a(t+1, s)) !q

                        !LocmodelObjects(unemployed(s))%q(a(t, s),idxaime(1),idxy(s),l(t, s)+1,a(t+1, s))
                        !!denom = &   !reshape(LocmodelObjects(unemployed(s))%p(a(t, s),idxaime(1),SPA,idxy(s),:,:),(/numPointsL*numPointsA/))
                        !!dot_product(reshape(LocmodelObjects(unemployed(s))%u(a(t, s),idxaime(1),i,idxy(s),:,:),(/numPointsL*numPointsA/)),reshape(LocmodelObjects(unemployed(s))%q(a(t, s),idxaime(1),idxy(s),:,:),(/numPointsL*numPointsA/)))
                        !posterior(s,i)=prior(s,i)*exp(LocmodelObjects(unemployed(s))%u(a(t, s),idxaime(1),i,idxy(s),l(t, s)+1,a(t+1, s)))!-denom !utility
                        !!LocmodelObjects(unemployed(s))%policy(a(t, s),idxaime(1),SPA,idxy(s),:,:)


                    end do
                    if (sum(posterior(t,s,:))>1.01 .or. minval(posterior(t,s,:)) < 0) then
                        write(*,*) "Wrong posterior"
                    end if
                    prior(t+1,s,:) = matmul(params%spaTransMat,posterior(t,s,:))
                    if (t <8) then
                        pos(t,s,:) = MAXLOC(posterior(t,s,:))
                        prio(t,s,:) = MAXLOC(prior(t,s,:))
                    end if
                end if
            end if

            if (l(t,s) .eq. 0) then
                !if (unemployed(s) == 2) then
                !    y(t, s)=grids%benefit(t)
                !else
                y(t, s) = 0.0
                !end if
            end if

            call nextAIME(grids,params,t,typeSim,l(t,s),y(t,s), aime(t, s))
            aime(t+1, s) =   aime(t, s)
            if (unemployed(s) == 2) then
                y(t, s)=grids%benefit(t)
            end if
            call gross2net(params,grids,y(t, s),t,l(t,s),typesim, aime(t, s),Tretire-1+SPA, dbpath(s), sppath(s) )
            c(t, s) = a(t, s)  + y(t, s) - (a(t+1, s)/(1+params%r))
            !
        end do
        !write(*,*) 'Mean inc ', sum(y(t,:))/numsims, 'in ', t
        !error
    end do !t

    !maks = statesVisited > 0
    !write (*,*) COUNT(maks)
    !write (*,*) real(sum(statesVisited),rk)/real(COUNT(maks),rk)
    !write (*,*) sum(statesVisited)
    !write (*,*) Tretire -2 +SPA, numSims, (Tretire -2 +SPA)*numSims
    !write (*,*) real(COUNT(maks),rk)/ ((Tretire -2 +SPA)*numPointsType*numPointsA*numAIME*numPointsY)!real(size(maks),rk)
    !!write (*,*) "Aime"
    !!do i=1,numAIME
    !!    write (*,*) sum(statesVisited(:,:,:,i,:))
    !!end do
    !
    !maxVec = maxloc(statesVisited)
    !write (*,*) statesVisited(maxVec(1),maxVec(2),maxVec(3),maxVec(4),maxVec(5))
    !t = maxVec(1)
    !check = 0.0
    !do s= 1, numsims
    !    if (real(s,rk)/numsims < params%fracType(1)) then
    !        typesim = 1
    !    else if (real(s,rk)/numsims < params%fracType(2)) then
    !        typesim = 2
    !    else if (real(s,rk)/numsims < params%fracType(3)) then
    !        typesim =3
    !    else
    !        typesim =4
    !    end if
    !    if (typesim == maxVec(2) .AND. a(t, s) == maxvec(3) .and. aimeInt(t, s) == maxvec(4) .and. yint(t, s) == maxvec(5) ) then
    !        write (*,*) "*************************************************"
    !        write (*,*) posterior(t,s,:)
    !        write (*,*) "Error = ", 0.5*sum(abs(posterior(t,s,:)-grids%posteriorSPA(t,typesim,a(t, s), aimeInt(t, s),  yint(t, s),:)))
    !        check = check +  0.5*sum(abs(posterior(t,s,:)-grids%posteriorSPA(t,typesim,a(t, s), aimeInt(t, s),  yint(t, s),:)))
    !        grids%posteriorSPA(t,typesim,a(t, s), aimeInt(t, s),  yint(t, s),:) = posterior(t,s,:)
    !        write (*,*) "*************************************************"
    !    end if
    !end do

    check = 0.0
    statesVisited = 0.0
    do t=1,Tretire - 2 +SPA
        do s= 1, numsims
            if (real(s,rk)/numsims < params%fracType(1)) then
                typesim = 1
            else if (real(s,rk)/numsims < params%fracType(2)) then
                typesim = 2
            else if (real(s,rk)/numsims < params%fracType(3)) then
                typesim =3
            else
                typesim =4
            end if
            if (statesVisited(t,typesim,a(t, s),aimeInt(t, s), yint(t, s)) == 0) then
                locposteriorSPA(t,typesim,a(t, s), aimeInt(t, s),  yint(t, s),:) = 0.0
            end if
            statesVisited(t,typesim,a(t, s),aimeInt(t, s), yint(t, s)) = statesVisited(t,typesim,a(t, s),aimeInt(t, s), yint(t, s)) + 1
            locposteriorSPA(t,typesim,a(t, s), aimeInt(t, s),  yint(t, s),:) = &
                ((statesVisited(t,typesim,a(t, s),aimeInt(t, s), yint(t, s))-1)*locposteriorSPA(t,typesim,a(t, s), aimeInt(t, s),  yint(t, s),:) + posterior(t,s,:) ) / statesVisited(t,typesim,a(t, s),aimeInt(t, s), yint(t, s))
            !if (typesim == maxVec(2) .AND. a(t, s) == maxvec(3) .and. aimeInt(t, s) == maxvec(4) .and. yint(t, s) == maxvec(5) ) then
            !write (*,*) "*************************************************"
            !write (*,*) posterior(t,s,:)
            !write (*,*) "Error = ", 0.5*sum(abs(posterior(t,s,:)-grids%posteriorSPA(t,typesim,a(t, s), aimeInt(t, s),  yint(t, s),:)))
            !check = check +  0.5*sum(abs(posterior(t,s,:)-grids%posteriorSPA(t,typesim,a(t, s), aimeInt(t, s),  yint(t, s),:)))
            !grids%posteriorSPA(t,typesim,a(t, s), aimeInt(t, s),  yint(t, s),:) = posterior(t,s,:)
            !write (*,*) "*************************************************"
            !end if
        end do
    end do
    grids%posteriorSPA = weight*locposteriorSPA + (1-weight)*grids%posteriorSPA
    check = 0.0
    check2 = 0.0
    do t=1,Tretire - 2 +SPA
        do s= 1, numsims
            if (real(s,rk)/numsims < params%fracType(1)) then
                typesim = 1
            else if (real(s,rk)/numsims < params%fracType(2)) then
                typesim = 2
            else if (real(s,rk)/numsims < params%fracType(3)) then
                typesim =3
            else
                typesim =4
            end if
            check = check +  0.5*sum(abs(posterior(t,s,:)-grids%posteriorSPA(t,typesim,a(t, s), aimeInt(t, s),  yint(t, s),:)))
            check2 = max(check2,0.5*sum(abs(posterior(t,s,:)-grids%posteriorSPA(t,typesim,a(t, s), aimeInt(t, s),  yint(t, s),:))))
            if (check2> 0.01) then
                continue
            end if
        end do
    end do
    if (present(error)) error = check / sum(statesVisited)
    !write(*,*) "Average error", check/statesVisited(maxVec(1),maxVec(2),maxVec(3),maxVec(4),maxVec(5))
    !write (*,*) 'no'
    if (checkMax > 0.00001) write (*,*), "Larege Error = ", checkMax
    write (*,*), "Laregest beliefe error = ", check2

    if (present(posOut)) posOut = pos
    if (present(prioOut)) prioOut = prio
    !write(*,*) 'Mean db ', sum(dbpath)/size(dbpath)
    !write(*,*) 'Mean sp ', sum(sppath)/size(sppath)
    !write(*,*) 'Mean inc ', sum(y)/size(y)
    !write(*,*) 'Mena emp int asset', sum(startinga)/real(numsims,rk)
    !write(*,*) 'Mena sim int asset', sum(grids%agrid(1,1,a(1,:)))/real(numsims,rk)
    write (outFile, *), trim(trim(pathDataStore) // "EVType"),typeSim,trim("Period"),t,".txt"
    write (outFile, *), trim(trim(path) // "posdata"), spa,  ".txt"
    outfile=ADJUSTL(outfile)
    write (format_numpeepcols_int,*),'(',numSims
    format_numpeepcols_int = trim(format_numpeepcols_int) // 'I2)'
    inquire (iolength=requiredl)  transpose(pos(:,:,1))
    open (unit=212, file=outFile, status='unknown',recl=requiredl, action='write')
    write (212,format_numpeepcols_int  ) transpose(pos(:,:,1))
    close( unit=212)


    write (outFile, *), trim(trim(path) // "pridata.txt")
    outfile=ADJUSTL(outfile)
    write (format_numpeepcols_int,*),'(',numSims
    format_numpeepcols_int = trim(format_numpeepcols_int) // 'I2)'
    inquire (iolength=requiredl)  transpose(prio(:,:,1))
    open (unit=212, file=outFile, status='unknown',recl=requiredl, action='write')
    write (212,format_numpeepcols_int  ) transpose(prio(:,:,1))
    close( unit=212)

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
    subroutine gross2net(params,grids,Y,ixt,ixl,ixType, AIME, spa, dbout,spout)
    implicit none

    !inputs
    type (structparamstype), intent(in) :: params
    type (gridsType), intent(in) :: grids

    integer, intent(in) :: ixt,ixl, ixType, spa
    real (kind=rk), intent(in) :: AIME

    !changing
    real (kind =rk), intent(inout) :: Y

    !output
    real (kind=rk), optional, intent(out) :: dbout, spout

    real (kind=rk) :: db, sp
    !!You can't work over certain age so whatever input income most be benefits until add on pension
    !if (ixt >= stopwrok) then
    !    Y    = grids%benefit(ixt)
    !end if
    db = 0.0
    sp = 0.0
    if (ixt >= spa) then
        sp  = statePension(params,AIME,ixType)
        Y = Y + sp!params%pension
    end if

    if  (ixt >= params%ERA(ixType)) then
        !if (ixL==0 .AND. mod(ixType,2)==0) then !
        db = dbPension(params,AIME,ixType)
        Y =   Y + db
        !end if
    else

    end if

    if (present(dbout)) dbout = db
    if (present(spout)) spout = sp
    Y = Y +params%spouseInc(ixType,ixt)

    end subroutine

    ! ---------------------------------------------------------------------------------------------------------!
    !!Returns GMM Cretieria
    !----------------------------------------------------------------------------------------------------------!
    function gmm(params,grids,target,weights, overridePF, overrideRecal, unweightedunsquare)
    implicit none

#ifdef mpiBuild
    include 'mpif.h'
#endif  

    !inputs
    type (structparamstype), intent(in) :: params
    type (gridsType), intent(inout) :: grids
    real (kind=rk), intent(in) :: target(:,:), weights(:,:)
    logical, intent(in), optional :: overridePF, overrideRecal
    real (kind=rk), intent(out), optional :: unweightedunsquare(nummoments)

    !output
    real (kind=rk) :: gmm

    !local
    real (kind=rk) :: ypath(Tperiods, numSims) !income
    real (kind=rk) :: cpath(Tperiods, numSims)  !consumption
    integer :: lpath(Tperiods, numSims) !labour supply
    real (kind=rk) :: vpath(Tperiods, numSims)  !value
    integer :: apath(Tperiods + 1,numSims) !this is the path at the start of each period, so we include the 'start' of death
    real (kind=rk) ::  meanL(Tperiods), meanA(Tperiods)
    real (kind=rk) :: AIME(Tperiods + 1,numSims)
    real (kind=rk) :: locL(24), locA(24)
    character(len=1024) :: outFile

    integer :: n, requiredl !, typeSim
    !integer :: i, k, start, finish
    !integer, allocatable:: rangeSims(:)
    logical :: delpf, reclloc
    integer, save :: filnum = 0
    character (20) :: format_numpeepcols_int
    !!Set asset grid
    !do typeSim = 1, numPointsType
    !    call getassetgrid( params, grids%maxInc(typeSim,:), grids%Agrid(typeSim,:,:))
    !end do
    !solve
    delpf = .true.
    reclloc = .true.
    if (present(overridePF)) delPF = overridePF
    if (present(overrideRecal)) reclloc = overrideRecal

    if (reclloc) call solveValueFunction( params, grids, .false., .false. )
    meanL = 0.0
    meanA = 0.0
    locL = 0.0
    locA = 0.0
    apath = 0.0
    lpath = 0
    cpath = 0.0
    AIME = 0.0
    vpath = 0.0

    if (rank==0) then
        !simulate
        call simWithUncer(params, grids, grids%Simy, cpath, apath,  lpath, ypath, AIME, TrueSPA, delpf ) !vpath,
        filnum = filnum+1

        !L
        write (outFile, *), trim(trim(path) // "lpath"),filnum,".txt"
        outfile=ADJUSTL(outfile)

        write (format_numpeepcols_int,*),'(',numSims !- 2*16383
        format_numpeepcols_int = trim(format_numpeepcols_int) // 'I2)'
        inquire (iolength=requiredl)  transpose(lpath(1:10,:))
        open (unit=212, file=outFile, status='unknown',recl=requiredl, action='write')
        write (212,format_numpeepcols_int  ) transpose(lpath(1:10,:)) !32767:numSims
        close( unit=212)
        !
        !inquire (iolength=requiredl)  lpath
        !open (unit=201, file=outfile, status='unknown',recl=requiredl, action='write')
        !write (201, '(6E15.7)' ) lpath
        !close( unit=201)

        do n=1,Tperiods
            !meanA(n)=sum(real(Apath(n,:),rk))/real(numSims,rk) !size(lpath(n,:))
            meanL(n)=sum(real(lpath(n,:),rk))/real(numSims,rk)
            meanA(n)=sum(grids%Agrid(1,1, apath(n,:) ))/real(numSims,rk)
            !    do k=1,numPointsType
            !        if (k == 1)  then
            !            start = 1
            !            if (numPointsType>1) then
            !                finish = params%fracType(1)*numSims
            !            else
            !                finish = numSims
            !            end if
            !            if (n>1) then
            !                deallocate(rangeSims)
            !            end if
            !            allocate(rangeSims(finish))
            !            rangeSims = (/(i,i=1,finish)/)
            !        else if (k==2) then
            !            deallocate(rangeSims)
            !            start = finish +1
            !            finish = params%fracType(2)*numSims
            !            allocate(rangeSims(finish-start+1))
            !            rangeSims = (/(i,i=start,finish)/)
            !        else if (k==3) then
            !            deallocate(rangeSims)
            !            start = finish +1
            !            finish = params%fracType(3)*numSims
            !            allocate(rangeSims(finish-start+1))
            !            rangeSims = (/(i,i=start,finish)/)
            !        else
            !            deallocate(rangeSims)
            !            start = finish +1
            !            finish = numSims
            !            allocate(rangeSims(finish-start+1))
            !            rangeSims = (/(i,i=start,finish)/)
            !        end if
            !        meanA(n)=meanA(n)+sum(real(grids%Agrid(k,n,apath(n,rangeSims)),rk))
            !    end do
            !    meanA(n)= meanA(n)/real(numSims,rk)
        end do

        locL = meanL(33 - (StartAge-20):33 - (StartAge-20)+23)-target(1,:)
        locA = meanA(33 - (StartAge-20):33 - (StartAge-20)+23)-target(2,:)
        if (present(unweightedunsquare)) then
            unweightedunsquare(1:24) = locl
            unweightedunsquare(25:48) = locA
        end if
        if (dimEstimation >= 5) then

        end if
        gmm = dot_product(weights(1,:)*locL,weights(1,:)*locL) + dot_product(weights(2,:)*locA,weights(2,:)*locA)
    end if

#ifdef mpiBuild 
    call MPI_Bcast( gmm,   1,   mpi_double_precision, 0, mpi_comm_world, ierror)
    call mpi_barrier(mpi_comm_world, ierror)
    if (ierror.ne.0) stop 'mpi problem180'
    if (present(unweightedunsquare)) then
        call MPI_Bcast( unweightedunsquare,   48,   mpi_double_precision, 0, mpi_comm_world, ierror)
        call mpi_barrier(mpi_comm_world, ierror)
        if (ierror.ne.0) stop 'mpi problem180'
    end if

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
    integer :: n, requiredl !, i
    real(kind=rk) :: meanL(Tperiods), meanV(Tperiods), meanA(Tperiods), meanC(Tperiods), meanY(Tperiods) !, medianA
    real(kind=rk) :: meanYemp(Tperiods), meanAIME(Tperiods), numLC(Tperiods) !, meanPoor(Tperiods), meanRich(Tperiods)
    real(kind=rk) :: adata(24,numsims), ydata(24,numsims)
    integer :: ldata(24,numsims)
    !integer :: lrich(Tperiods,numSims/2), lpoor(Tperiods,numSims/2)
    !integer :: rich, poor, k, start, finish
    !integer, allocatable:: rangeSims(:)
    character (20) :: format_numpeepcols_int
    character(len=1024) :: outFile

    !write (outFile, *), trim(trim(path) // "policyL"),typeSim,trim("Period"),ixt,trim("SPA60"),".txt"
    !outfile=ADJUSTL(outfile)

    do n=1,Tperiods
        meanL(n)=sum(real(lpath(n,:),rk))/real(numSims,rk) !size(lpath(n,:))
        meanL(n)=sum(real(lpath(n,:),rk))/real(numSims,rk)
        meanA(n)=sum(grids%Agrid(1,1, apath(n,:) ))/real(numSims,rk)
        !write (201, * ) meanL(n)
        meanV(n)=sum(real(vpath(n,:),rk))/real(numSims,rk)
        !write (202, * ) meanV(n)
        !meanA(n)=0.0
        !do k=1,numPointsType
        !    if (k == 1)  then
        !        start = 1
        !        if (numPointsType>1) then
        !            finish = 0.16*numSims
        !        else
        !            finish = numSims
        !        end if
        !        if (n>1) then
        !            deallocate(rangeSims)
        !        end if
        !        allocate(rangeSims(finish))
        !        rangeSims = (/(i,i=1,finish)/)
        !    else if (k==2) then
        !        deallocate(rangeSims)
        !        start = finish +1
        !        finish = 0.48*numSims
        !        allocate(rangeSims(finish-start+1))
        !        rangeSims = (/(i,i=start,finish)/)
        !        !real(s,rk)/numSims <0.48
        !    else if (k==3) then
        !        deallocate(rangeSims)
        !        start = finish +1
        !        finish = 0.59*numSims
        !        allocate(rangeSims(finish-start+1))
        !        rangeSims = (/(i,i=start,finish)/)
        !        !real(s,rk)/numSims <0.59
        !    else
        !        deallocate(rangeSims)
        !        start = finish +1
        !        finish = numSims
        !        allocate(rangeSims(finish-start+1))
        !        rangeSims = (/(i,i=start,finish)/)
        !    end if
        !    meanA(n)=meanA(n)+sum(real(grids%Agrid(k,n,apath(n,rangeSims)),rk))
        !    if (n>=52-startAge+1 .AND. n <=75-startAge+1) then
        !        adata(n,rangeSims) = real(grids%Agrid(k,n,apath(n,rangeSims)),rk)
        !    end if
        !end do
        !
        !meanA(n)= meanA(n)/real(numSims,rk)
        if (n>=52-startAge+1 .AND. n <=75-startAge+1) then
            ldata(n,:) =lpath(n,:)
            adata(n,:) =grids%Agrid(1,1, apath(n,:) )
            ydata(n,:) =lpath(n,:)*yemp(n,:)
        end if


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
    write (outFile, *), trim(trim(path) // "lpath.txt")
    outfile=ADJUSTL(outfile)
    inquire (iolength=requiredl)  meanL
    open (unit=201, file=outfile, status='unknown',recl=requiredl, action='write')
    write (201, '(6E15.7)' ) meanL
    close( unit=201)

    !V
    write (outFile, *), trim(trim(path) // "Vpath.txt")
    outfile=ADJUSTL(outfile)
    inquire (iolength=requiredl)  meanV
    open (unit=202, file=outfile, status='unknown',recl=requiredl, action='write')
    write (202, '(6E15.7)' ) meanV
    close( unit=202)

    !A
    write (outFile, *), trim(trim(path) // "Apath.txt")
    outfile=ADJUSTL(outfile)
    inquire (iolength=requiredl)  meanA
    open (unit=203, file=outfile, status='unknown',recl=requiredl, action='write')
    write (203, '(6E15.7)' ) meanA
    close( unit=203)

    !C
    write (outFile, *), trim(trim(path) // "Cpath.txt")
    outfile=ADJUSTL(outfile)
    inquire (iolength=requiredl)  meanC
    open (unit=204, file=outFile, status='unknown',recl=requiredl, action='write')
    write (204, '(6E15.7)' ) meanC
    close( unit=204)

    !Y
    write (outFile, *), trim(trim(path) // "Ypath.txt")
    outfile=ADJUSTL(outfile)
    inquire (iolength=requiredl)  meanY
    open (unit=205, file=outfile, status='unknown',recl=requiredl, action='write')
    write (205, '(6E15.7)' ) meanY
    close( unit=205)

    write (outFile, *), trim(trim(path) // "YempPath.txt")
    outfile=ADJUSTL(outfile)
    inquire (iolength=requiredl)  meanYemp
    open (unit=206, file=outFile, status='unknown',recl=requiredl, action='write')
    write (206, '(6E15.7)' ) meanYemp
    close( unit=206)

    write (outFile, *), trim(trim(path) // "AIMEPath.txt")
    outfile=ADJUSTL(outfile)
    inquire (iolength=requiredl)  meanAIME
    open (unit=207, file=outFile, status='unknown',recl=requiredl, action='write')
    write (207, '(6E15.7)' ) meanAIME
    close( unit=207)

    write (outFile, *), trim(trim(path) // "numLC.txt")
    outfile=ADJUSTL(outfile)
    inquire (iolength=requiredl)  numLC
    open (unit=207, file=outfile, status='unknown',recl=requiredl, action='write')
    write (207, '(6E15.7)' ) numLC
    close( unit=207)

    !!!!!!
    write (format_numpeepcols_int,*),'(',numSims
    format_numpeepcols_int = trim(format_numpeepcols_int) // 'E15.7)'

    write (outFile, *), trim(trim(path) // "ydata.txt")
    outfile=ADJUSTL(outfile)
    inquire (iolength=requiredl)  transpose(ydata)
    open (unit=211, file=outFile, status='unknown',recl=requiredl, action='write')
    write (211, format_numpeepcols_int ) transpose(ydata)
    close( unit=211)

    write (outFile, *), trim(trim(path) // "adata.txt")
    outfile=ADJUSTL(outfile)
    inquire (iolength=requiredl)  transpose(adata)
    open (unit=211, file=outFile, status='unknown',recl=requiredl, action='write')
    write (211, format_numpeepcols_int ) transpose(adata)
    close( unit=211)

    write (outFile, *), trim(trim(path) // "ldata.txt")
    outfile=ADJUSTL(outfile)
    write (format_numpeepcols_int,*),'(',numSims
    format_numpeepcols_int = trim(format_numpeepcols_int) // 'I2)'
    inquire (iolength=requiredl)  transpose(ldata)
    open (unit=212, file=outFile, status='unknown',recl=requiredl, action='write')
    write (212,format_numpeepcols_int  ) transpose(ldata)
    close( unit=212)


    end subroutine
    ! ---------------------------------------------------------------------------------------------------------!
    ! !Write to file
    !---------------------------------------------------------------------------------------------------------!
    subroutine writetofileByType(grids,params,ypath, cpath, apath, vpath, lpath, yemp, AIME)
    implicit none

    !inputs
    type (gridsType), intent(inout) :: grids
    type (structparamstype), intent(in) :: params
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
    real(kind=rk) :: meanL(Tperiods), meanV(Tperiods), meanA(Tperiods), meanC(Tperiods), meanY(Tperiods) !, medianA
    real(kind=rk) :: meanYemp(Tperiods), meanAIME(Tperiods), numLC(Tperiods) !, meanPoor(Tperiods), meanRich(Tperiods)
    !integer :: lrich(Tperiods,numSims/2), lpoor(Tperiods,numSims/2)
    !integer :: rich, poor
    character(len=1024) :: outFile
    integer, allocatable :: rangeSims(:)

#ifdef win
    do k=1,numPointsType
        if (k == 1)  then
            start = 1
            finish = params%fracType(1)*numSims
            allocate(rangeSims(finish))
            rangeSims = (/(i,i=1,finish)/)
        else if (k==2) then
            deallocate(rangeSims)
            start = finish +1
            finish = params%fracType(2)*numSims
            allocate(rangeSims(finish-start+1))
            rangeSims = (/(i,i=start,finish)/)
            !real(s,rk)/numSims <0.48
        else if (k==3) then
            deallocate(rangeSims)
            start = finish +1
            finish = params%fracType(3)*numSims
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
        write (outFile, *), trim(trim(path) // "lpath"),k,".txt"
        outfile=ADJUSTL(outfile)
        inquire (iolength=requiredl)  meanL
        open (unit=201, file=outfile, status='unknown',recl=requiredl, action='write')
        write (201, '(6E15.3)' ) meanL
        close( unit=201)

        !V
        write (outFile, *), trim(trim(path) // "Vpath"),k,".txt"
        !write (outFile, *), trim("..\\out\Vpath"),k,".txt"
        outfile=ADJUSTL(outfile)
        inquire (iolength=requiredl)  meanV
        open (unit=202, file=outFile, status='unknown',recl=requiredl, action='write')
        write (202, '(6E15.3)' ) meanV
        close( unit=202)

        !A
        write (outFile, *), trim(trim(path) // "Apath"),k,".txt"
        !write (outFile, *), trim("..\\out\Apath"),k,".txt"
        outfile=ADJUSTL(outfile)
        inquire (iolength=requiredl)  meanA
        open (unit=203, file=outFile, status='unknown',recl=requiredl, action='write')
        write (203, '(6E15.3)' ) meanA
        close( unit=203)

        !C
        write (outFile, *), trim(trim(path) // "Cpath"),k,".txt"
        !write (outFile, *), trim("..\\out\Cpath"),k,".txt"
        outfile=ADJUSTL(outfile)
        inquire (iolength=requiredl)  meanC
        open (unit=204, file=outFile, status='unknown',recl=requiredl, action='write')
        write (204, '(6E15.3)' ) meanC
        close( unit=204)

        !Y
        write (outFile, *), trim(trim(path) // "Ypath"),k,".txt"
        !write (outFile, *), trim("..\\out\Ypath"),k,".txt"
        outfile=ADJUSTL(outfile)
        inquire (iolength=requiredl)  meanY
        open (unit=205, file=outFile, status='unknown',recl=requiredl, action='write')
        write (205, '(6E15.3)' ) meanY
        close( unit=205)

        write (outFile, *), trim(trim(path) // "Yemppath"),k,".txt"
        !write (outFile, *), trim("..\\out\YempPath"),k,".txt"
        outfile=ADJUSTL(outfile)
        inquire (iolength=requiredl)  meanYemp
        open (unit=206, file=outFile, status='unknown',recl=requiredl, action='write')
        write (206, '(6E15.3)' ) meanYemp
        close( unit=206)

        write (outFile, *), trim(trim(path) // "AIMEpath"),k,".txt"
        !write (outFile, *), trim("..\\out\AIMEPath"),k,".txt"
        outfile=ADJUSTL(outfile)
        inquire (iolength=requiredl)  meanAIME
        open (unit=207, file=outFile, status='unknown',recl=requiredl, action='write')
        write (207, '(6E15.3)' ) meanAIME
        close( unit=207)

        write (outFile, *), trim(trim(path) // "numLC"),k,".txt"
        !write (outFile, *), trim("..\\out\numLC"),k,".txt"
        outfile=ADJUSTL(outfile)
        inquire (iolength=requiredl)  numLC
        open (unit=207, file=outFile, status='unknown',recl=requiredl, action='write')
        write (207, '(6E15.3)' ) numLC
        close( unit=207)
    end do
#else

#endif

    end subroutine
    ! ---------------------------------------------------------------------------------------------------------!
    ! !Write to file
    !---------------------------------------------------------------------------------------------------------!
    subroutine writetofileAge(grids,params, ypath, cpath, apath, vpath, lpath, yemp, AIME, SPA)
    implicit none

    !inputs
    type (gridsType), intent(inout) :: grids
    type (structparamstype), intent(in) :: params
    real (kind=rk), intent(in) :: ypath(Tperiods, numSims) !income
    real (kind=rk), intent(in) :: cpath(Tperiods, numSims)  !consumption
    integer, intent(in) :: lpath(Tperiods, numSims) !labour supply
    real (kind=rk), intent(in) :: vpath(Tperiods, numSims)  !value
    integer, intent(in) :: apath(Tperiods + 1,numSims) !this is the path at the start of each period, so we include the 'start' of death
    real (kind=rk), intent(in)  :: yemp(Tperiods, numSims)
    real (kind=rk), intent(in)  :: AIME(Tperiods + 1,numSims)
    integer, intent(in) :: SPA

    !local
    integer :: n, requiredl , i
    real (kind=rk) :: meanA(Tperiods)
    !real(kind=rk) :: meanL(Tperiods), meanV(Tperiods), meanC(Tperiods), meanY(Tperiods), medianA
    !real(kind=rk) :: meanYemp(Tperiods), meanAIME(Tperiods), meanPoor(Tperiods), meanRich(Tperiods), numLC(Tperiods)
    real(kind=rk) :: adata(24,numsims), ydata(24,numsims)
    integer :: ldata(24,numsims)
    !integer :: lrich(Tperiods,numSims/2), lpoor(Tperiods,numSims/2)
    integer :: k, start, finish !,rich, poor
    integer, allocatable:: rangeSims(:)
    character (20) :: format_numpeepcols_int
    character(len=1024) :: outFile

#ifdef win 

    do n=1,Tperiods
        do k=1,numPointsType
            if (k == 1)  then
                start = 1
                if (numPointsType>1) then
                    finish = params%fracType(1)*numSims
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
                finish = params%fracType(2)*numSims
                allocate(rangeSims(finish-start+1))
                rangeSims = (/(i,i=start,finish)/)
            else if (k==3) then
                deallocate(rangeSims)
                start = finish +1
                finish = params%fracType(3)*numSims
                allocate(rangeSims(finish-start+1))
                rangeSims = (/(i,i=start,finish)/)
            else
                deallocate(rangeSims)
                start = finish +1
                finish = numSims
                allocate(rangeSims(finish-start+1))
                rangeSims = (/(i,i=start,finish)/)
            end if
            meanA(n)=meanA(n)+sum(real(grids%Agrid(k,n,apath(n,rangeSims)),rk))
            if (n>=52-startAge+1 .AND. n <=75-startAge+1) then
                adata(n,rangeSims) = real(grids%Agrid(k,n,apath(n,rangeSims)),rk)
            end if
        end do

        if (n>=52-startAge+1 .AND. n <=75-startAge+1) then
            ldata(n,:) =lpath(n,:)
            ydata(n,:) = lpath(n,:)*yemp(n,:)
        end if
    end do

    !!!!!!
    write (format_numpeepcols_int,*),'(',numSims
    format_numpeepcols_int = trim(format_numpeepcols_int) // 'E15.7)'

    write (outFile, *), trim(trim(path) // "ydata"), SPA
    outFile = trim(outFile) // ".txt"
    outfile=ADJUSTL(outfile)
    inquire (iolength=requiredl)  transpose(ydata)
    open (unit=211, file=outFile, status='unknown',recl=requiredl, action='write')
    write (211, format_numpeepcols_int ) transpose(ydata)
    close( unit=211)

    write (outFile, *), trim(trim(path) // "adata"), SPA
    outFile = trim(outFile) // ".txt"
    outfile=ADJUSTL(outfile)
    inquire (iolength=requiredl)  transpose(adata)
    open (unit=211, file=outFile, status='unknown',recl=requiredl, action='write')
    write (211, format_numpeepcols_int ) transpose(adata)
    close( unit=211)

    write (outFile, *), trim(trim(path) // "ldata"), SPA
    outFile = trim(outFile) // ".txt"
    outfile=ADJUSTL(outfile)
    write (format_numpeepcols_int,*),'(',numSims
    format_numpeepcols_int = trim(format_numpeepcols_int) // 'I2)'
    inquire (iolength=requiredl)  transpose(ldata)
    open (unit=212, file=outFile, status='unknown',recl=requiredl, action='write')
    write (212,format_numpeepcols_int  ) transpose(ldata)
    close( unit=212)
#else

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
    real (kind=rk) ::  maxmaxA, growth, error, a, b, errorgrid(0:100) !, test, loggrid(numPointsA), maxA(Tperiods+1), span
    integer :: ixt, i, i1(1), i2(1)


    !!Set maximum assets
    !maxA(1) = params%startA
    !do ixt = 2, Tperiods+1
    !    maxA(ixt) = (maxA(ixt - 1) + maxInc(ixt-1) ) * (1+params%r)
    !end do
    maxmaxA = 1500000.00 !1250000.00! maxval(maxA)

    do i = 0,100
        growth = real(i,rk)/100_rk
        !if (i==0) growth = 0.001
        errorgrid(i) = func(growth)
    end do
    i1 = maxloc(errorgrid)
    i2(1) = i1(1) - 2
    !growth = real(i1(1),rk)/100_rk
    !i2 = minloc((/errorgrid(i1-1),errorgrid(i1+1)/))
    !if (i1(1)> i2(1)) then
    !    a = i2(1)/10
    !    b = i1(1)/10
    !else
    a = real(i2(1),rk)/100_rk
    b = real(i1(1),rk)/100_rk
    !end if

    error = golden_generic(a, b, growth, func,0.001_rk,.false.)
    call getnodes(Agrid(1, :), 0.0_rk, maxmaxA, numPointsA, growth)
    do ixt = 2, Tperiods+1
        Agrid(ixt, :) = Agrid(1, :)
    end do

    contains
    function func(growth)
    implicit none
    real (kind=rk), intent(in) :: growth
    real (kind=rk) :: func

    call getnodes(Agrid(1, :), 0.0_rk, maxmaxA, numPointsA, growth)
    func = -abs(500.0 - Agrid(1, 2))
    end function
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
    integer :: ixl, policyA1temp !, maxL
    real(kind=rk) :: negV, negVtemp


    negV = -huge(negv) !-log(0.0) !inf

    !if have unemployment shock can't work#
    !maxL = 1 - abs(mod(ixy))
    do ixL = 0,(numPointsL-1),1           ! points of income choice

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
    subroutine solvePeriodRI(params, grids, ixy, ixA, ixAIME ,ixt,ixType,spa,  EV1, ccpOut, v, q, uConst)
    implicit none

    !input
    real(kind=rk), intent(in) ::  EV1(:,:,:)
    type (structparamstype), intent(in) :: params
    type (gridsType), intent(in) :: grids
    integer, intent(in) :: ixt,ixType,ixa,ixAIME,ixy, spa

    !!output
    real(kind=rk), intent(out) :: ccpOut(:,:,:), v(:), q(:,:), uConst(:,:,:)
    !local
    !integer, parameter :: hp = selected_real_kind(16)
    integer :: ixl,  ixA1, i, j, dim, dimD ,K, labourChoices, offset, maxmaxA, sizeBuffer !A1, testInt
    integer, allocatable :: maxA(:)
    real(kind=rk) :: Y,  ubA1, A,AIME, EVloc(numAIME), checkSave(grids%supportSPA(ixt)), checkOther(grids%supportSPA(ixt)) !cons, va1, vb1, check
    real(kind=rk) :: scale !, toputil
    real(kind=rk), allocatable :: values(:), ccpMat(:,:,:) !, eye(:,:)
    real(kind=rk), allocatable :: const(:,:,:), GrossUtil(:,:)
    real(kind=rk), allocatable  ::  parameters(:), ccp(:), test(:), utilVec(:), buffer(:,:)
    integer, save :: saveDim(2,numPointsL), unemp, locl, locA1, morethan(11) = 0 !,choices = 0!, div, lastixT = 0
    integer, allocatable :: rang(:), bufferInd(:)
    !real(kind=rk), allocatable, save :: locinitialGuessRI(:,:)
    logical :: converged !, unchanged
    !character (len=2000) :: text(2)

    !temp
    !real(kind=rk), allocatable :: tempC1(:,:), tempC2(:,:)
    real(kind=rk) ::  lowestUtil, psum, temp !diff,
    !integer :: iter, D
    !REAL(kind=rk) :: error
    !LOGICAL :: logcheck
    integer, allocatable          :: kwa(:)
    integer                       :: n, me, m, np, mnn2, maxit, maxfun, iprint, &
        maxnm, iout, mode, ifail, lwa, lkwa, lactiv
    double precision, allocatable :: x(:), g(:), df(:), dg(:,:), u(:), &
        xl(:), xu(:), c(:,:),  wa(:), d(:), payoff(:,:), denom(:)
    double precision                 f(1), acc, accqp, stpmin, rho
    logical, allocatable          :: active(:)
    logical                       :: lql
    external                      :: ql
    integer :: iter !meq, ierr, nact,

    ixa1 = 0
    lowestUtil = 0.0
    !If suffered unemployment shock then can't choose to work
    labourChoices=mod(ixy,2)*numPointsL+(1-mod(ixy,2))*1
    allocate(maxA(labourChoices), const(labourChoices,numPointsA,grids%supportSPA(ixt)))
    maxA(1) = numPointsA !default to all asset level possible
    !if (ixy == 4 .AND. ixt ==14 .AND. ixAIME == 2) then
    !    write(*,*) 'stop'
    !end if
    !if ( ixy == 4 .and. ixAIME == 2 .and. ixt==8 .and. ixType == 4 .and. ixa == 4) then
    if ( ixy == 1 .and. ixAIME == 1 .and. ixt==9 .and. ixType == 4 .and. ixa == 3) then
        continue
    end if
    do ixL = 0,(labourChoices-1),1           ! points of income choice
        ! Value of income and information for optimisation
        if (ixa1 > numPointsA) maxA(ixl) = numPointsA
        A    = grids%Agrid(ixType,ixt, ixA)            ! assets today

        Y    = grids%Ygrid(ixType,ixt, ixY)
        Y    = ixL*Y

        AIME = grids%AIMEgrid(ixType,ixt,ixAIME)

        call nextAIME(grids,params,ixt,ixType,ixl,y,AIME)

        Y    = ixL*Y+ abs(mod(ixy,2)==0 )*grids%benefit(ixt)

        call gross2net(params,grids,Y,ixt,ixl,ixType,AIME,spa)

        ubA1 = (A + Y - params%minCons)*(1+params%r);    ! upper bound: assets tomorrow
        do ixA1 = 1, numPointsA
            if (grids%Agrid(ixType,ixt+1, ixA1) > ubA1) then
                const(ixl+1,ixa1:numPointsA,:) = -huge(const(ixl+1,ixa1,1))
                uConst(:,ixl+1,ixa1:numPointsA) = -huge(const(ixl+1,ixa1,1))
                maxA(ixl+1) = ixa1 - 1
                exit
            end if

            do i = 1, grids%supportSPA(ixt)-1
                EVloc =  (1-params%p)*EV1(ixA1,:,numPointsSPA-grids%supportSPA(ixt)+i)+params%p*EV1(ixA1,:,numPointsSPA-grids%supportSPA(ixt)+i+1)
                const(ixl+1,ixa1,i) = objectivefunc(params, grids,grids%Agrid(ixType,ixt+1, ixA1), A, Y,ixL,ixt,ixType, AIME,EVloc)
                !if (const(ixl+1,ixa1,i) <lowestUtil) lowestUtil = const(ixl+1,ixa1,i)
                uConst(numPointsSPA-grids%supportSPA(ixt)+i,ixl+1,ixa1) = const(ixl+1,ixa1,i)
            end do
            EVloc =  EV1(ixA1,:,numPointsSPA)
            const(ixl+1,ixa1,i) = objectivefunc(params, grids,grids%Agrid(ixType,ixt+1, ixA1), A, Y,ixL,ixt,ixType, AIME,EVloc)
            uConst(numPointsSPA-grids%supportSPA(ixt)+i,ixl+1,ixa1) = const(ixl+1,ixa1,i)
            !if (const(ixl+1,ixa1,i) < lowestUtil) lowestUtil = const(ixl+1,ixa1,i)
        end do

    end do
    if (ixa1 > numPointsA) maxA(labourChoices) = numPointsA
    maxmaxA = maxval(maxA)

    !allocate(tempC1(maxA(1),grids%supportSPA(ixt)))
    !tempC1 = const(1,1:maxA(1),:)
    !if (labourChoices==2) then
    !    allocate(tempC2(maxA(2),grids%supportSPA(ixt)))
    !    tempC2 = const(2,1:maxA(2),:)
    !    if (rank==0) write(*,*) max(maxval(tempC1),maxval(tempC2)) - min(minval(tempC1),minval(tempC2))
    !    diff = max(maxval(tempC1),maxval(tempC2)) - min(minval(tempC1),minval(tempC2))
    !else
    !    if (rank==0) write(*,*) maxval(tempC1) - minval(tempC1)
    !    diff = maxval(tempC1) - minval(tempC1)
    !end if

    scale = maxval(const)
    !temp = log(1.0d300)
    scale= scale!-temp
    const = const-scale

    const=exp(const)


    dimD = sum(maxA)!labourChoices*(ixa1-1)
    dim = grids%supportSPA(ixt)*dimD !*labourChoices*(ixa1-1)
    allocate(values(dim+1),parameters(dim),ccp(dim),ccpMat(grids%supportSPA(ixt),labourChoices,maxmaxA), buffer(dimD,grids%supportSPA(ixt)),bufferInd(dimD), GrossUtil(dimD,grids%supportSPA(ixt)),test(dimD),utilVec(dim),rang(dimd))
    unemp = mod(ixy-1,2)+1

    locl =1
    locA1 = 0
    do i=1,dimD
        locA1 = locA1 + 1
        if (locA1 > maxa(locl)) then
            locl = locl+1
            locA1=1
        end if
        do j = 1,grids%supportSPA(ixt)
            GrossUtil(i,j) = const(locl,locA1,j)
        end do
    end do

    !if ( ixy == 1 .and. ixAIME == 2 .and. ixt==16 .and. ixType == 1 .and. ixa == 2) then
    !    continue
    !end if
    call skylineBNL(GrossUtil,dimD,grids%supportSPA(ixt),buffer,bufferInd,sizeBuffer )
    if (minval(buffer(1:sizeBuffer,:))< 0.0) then
        write (*,*) 'break'
    end if
    morethan(11) = morethan(11) + 1
    do i=1,10
        if (sizebuffer > i ) then
            morethan(i) = morethan(i) + 1
        else
            exit
        end if
    end do
    !if (sizebuffer >= 3 ) then

    !end if
    !if (sizebuffer > choices) then
    !    choices = sizebuffer
    !    !write (*,*) 'num chouices', choices
    !    !write (*,*) 'dim reduc ' , real(choices,rk)/real(dimD,rk)
    !end if
    !return

    if (ixt==9) then
        if (ixa == 20) then
            if (ixAIME == 2) then
                if (ixy == 3) then
                    iout = 6
                end if
            end if
        end if
    end if
    if (sizebuffer> 1) then
        if (sizebuffer>grids%supportSPA(ixt)) then
            !write(*,*) "undefined!"
        end if
        !   Set some constants and initial values

        iout   = 6
        acc    = 1.0d-14 !1.0d-9
        if (params%lambda < 1.0D-11 ) acc = 0.01
        accqp  = 0 !1.0d-14
        stpmin = 1.0d-5
        maxit  = 100
        if (params%lambda < 1.0D-11 ) maxit = 1000 
        maxfun = 20
        if (params%lambda < 1.0D-9 ) maxfun = 50
        maxnm  = 10
        rho    = 0
        lql    = .true.
        iprint = 0
        n      = sizebuffer
        np     = 1
        m      = 1
        me     = 1
        mnn2   = m + n + n + 2
        mode   = 0
        ifail  = 0
        lwa    = 3*(n+1)*(n+1)/2 + 33*(n+1) + 9*m + 150
        lkwa   = n + 30
        lactiv = 2*m + 10

        !   Allocate arrays

        allocate ( x(n+1), xl(n+1), xu(n+1), df(n+1), g(m), &
            dg(m,n+1), u(mnn2), c(n+1,n+1), d(n+1), &
            wa(lwa), kwa(lkwa), active(lactiv), payoff(sizebuffer,grids%supportSPA(ixt)), denom(grids%supportSPA(ixt)))


        !   Set starting values

        xl = 0.0
        x  = 0.5
        xu = 1.0

        do i=1, sizebuffer
            do j =1,  grids%supportSPA(ixt)
                payoff(i,j) = buffer(i,j)
                if (payoff(i,j) <= 0.0) then
                    payoff(i,j) = tiny(payoff)
                end if
            end do
        end do

        denom = 0.0
        do j =1,  grids%supportSPA(ixt)
            do i=1, sizebuffer
                denom(j) = denom(j) + x(i)*payoff(i,j)
            end do
        end do
        f = 0
        do j =1,  grids%supportSPA(ixt)
            f(1) = f(1) - grids%posteriorSPA(ixt,ixtype,ixa, ixAIME,  ixy,numPointsSPA-grids%supportSPA(ixt)+j)*log(denom(j))
        end do
        g(1) =  1.0- sum(x(1:n))

        !============================================================
        !        if (ifail.EQ.-1) goto 4
        !2       continue
        !============================================================
        !   This block computes all derivative values.
        !
        df = 0
        do i=1, n
            do j =1,  grids%supportSPA(ixt)
                df(i)   = df(i)-((grids%posteriorSPA(ixt,ixtype,ixa, ixAIME,  ixy,numPointsSPA-grids%supportSPA(ixt)+j)*payoff(i,j))/(denom(j)))
            end do
        end do
        dg =  -1
        ifail = -1
        iter = 0
        do while (ifail  < 0 )
            if (iter==0)ifail= 0
            iter = iter + 1

            if ( ixy == 7 .and. ixAIME == 3 .and. ixt==14 .and. ixType == 1 .and. ixa == 1) then
                continue
                iprint = 4 
            end if
            call nlpqlp (    np,      m,     me,      m,      n, &
                n+1,   mnn2,      x,      f,      g, &
                df,     dg,      u,     xl,     xu, &
                c,      d,    acc,  accqp, stpmin, &
                maxfun,  maxit,  maxnm,    rho, iprint, &
                mode,   iout,  ifail,     wa,    lwa, &
                kwa,   lkwa, active, lactiv,    lql, &
                ql)

            if (ifail>0) then
                write(*,*) "failed"
            end if

            if (ifail==-1) then
                denom = 0.0
                do j =1,  grids%supportSPA(ixt)
                    do i=1, sizebuffer
                        denom(j) = denom(j) + x(i)*payoff(i,j)
                    end do
                end do
                f = 0
                do j =1,  grids%supportSPA(ixt)
                    f(1) = f(1) - grids%posteriorSPA(ixt,ixtype,ixa, ixAIME,  ixy,numPointsSPA-grids%supportSPA(ixt)+j)*log(denom(j))
                end do
                g(1) =  1.0- sum(x(1:n))
            elseif (ifail==-2) then
                df = 0
                do i=1, n
                    do j =1,  grids%supportSPA(ixt)
                        df(i)   = df(i)-((grids%posteriorSPA(ixt,ixtype,ixa, ixAIME,  ixy,numPointsSPA-grids%supportSPA(ixt)+j)*payoff(i,j))/(denom(j)))
                    end do
                end do
                dg =  -1
            end if
        end do


        if (sizebuffer>grids%supportSPA(ixt)) then
            continue

        end if
        if ( ixy == 7 .and. ixAIME == 2 .and. ixt==17) then
            continue
        end if
        if (maxval(x)<0.9) then
            continue

        end if

        ccp = 0.0
        do j =1,  grids%supportSPA(ixt)
            psum = 0.0
            do i=1,sizebuffer
                ccp((bufferInd(i)-1)*grids%supportSPA(ixt)+j) = x(i)*payoff(i,j)
                psum = psum + ccp((bufferInd(i)-1)*grids%supportSPA(ixt)+j)
            end do
            ccp((/((bufferInd(i)-1)*grids%supportSPA(ixt)+j, i=1,sizebuffer)/)) = ccp((/((bufferInd(i)-1)*grids%supportSPA(ixt)+j, i=1,sizebuffer)/))/ psum
        end do
        !Get defaults actiosn
        q= 0.0
        j = 1
        locl =1
        locA1 = 0
        do i=1,dimD
            locA1 = locA1 + 1
            if (locA1 > maxa(locl)) then
                locl = locl+1
                locA1=1
            end if
            if (i==bufferInd(j)) then
                q(locl,locA1) = x(j)
                j = j + 1
            end if
        end do
    else
        ccp = 0.0
        ccp((bufferInd(1)-1)*grids%supportSPA(ixt)+1:bufferInd(1)*grids%supportSPA(ixt)) = 1.0

        !Get defaults actiosn
        q= 0.0
        locl =1
        locA1 = 0
        do i=1,dimD
            locA1 = locA1 + 1
            if (locA1 > maxa(locl)) then
                locl = locl+1
                locA1=1
            end if
            if (i==bufferInd(1)) then
                q(locl,locA1) = 1.0
            end if
        end do
        allocate ( x(1))
        x(1) = 1.0
    end if

    !!Cache for use on next round
    if (ixy < numPointsY) then
        saveDim(unemp,:) = grids%supportSPA(ixt)*maxa
        !locinitialGuessRI(unemp,:)=0.0
        offset = 0
        do k=1,labourChoices
            !locinitialGuessRI(unemp,(k-1)*grids%supportSPA(ixt)*numPointsA+1:(k-1)*grids%supportSPA(ixt)*numPointsA+saveDim(unemp,k))= ccp((k-1)*offset+1:k*saveDim(unemp,k))
            offset = offset + saveDim(unemp,k)
        end do
    end if

    ccpOut=0.0
    ccpMat = 0.0
    !ccpMat=reshape(ccp, (/grids%supportSPA(ixt),labourChoices,ixa1-1/))
    offset = 0
    do i=1,labourChoices
        do j =1,maxmaxA
            if (j > maxa(i)) then
                ccpMat(:,i,j) = 0.0
                exit
            end if
            do k=1,grids%supportSPA(ixt)
                ccpMat(k,i,j) =  ccp(k+(j-1)*grids%supportSPA(ixt)+offset)
            end do
        end do
        offset = offset+saveDim(unemp,i)
    end do
    ccpOut(1:grids%supportSPA(ixt),1:labourChoices,1:maxmaxA)=ccpMat


    !do i=1,labourChoices
    !    do j =1,maxmaxA
    !        q(i,j) = dot_product(ccpOut(:,i,j),grids%posteriorSPA(ixt,ixtype,ixa, ixAIME,  ixy,numPointsSPA-grids%supportSPA(ixt)+1:numPointsSPA))
    !    end do
    !end do

    do i=1,grids%supportSPA(ixt)
        checkOther(i) = sum(ccpOut(i,:,:))
        !write (*,*) checkSave(i), checkOther(i)
        if (checkOther(i) > 1.05 .OR. checkOther(i)<0.95) then
            write (*,*) "CCP doesn't sum to 1", checkOther(i), i, converged,  ixy, ixA, ixAIME, ixt, unemp
            write (*,*) "converg check", checksave
            write (*,*) "converg check", checkother
            write (*,*) "ccp", ccp
            write (*,*) "ccpout", ccpOut
            stop
        end if
        !ccp((/(i+j*grids%supportSPA(ixt),j=0,Dimd-1)/))=ccp((/(i+j*grids%supportSPA(ixt),j=0,Dimd-1)/))*check**-1
    end do
    !if (max(checkOther(1),checkOther(2)) /= max(checkSave(1), checkSave(2)) .OR. min(checkOther(1),checkOther(2)) /= min(checkSave(1), checkSave(2)) )pause
    !Calculate continuation value
    do i=1,grids%supportSPA(ixt)
        test = 0.0
        do j=1,sizebuffer
            test(j)=x(j)* GrossUtil(bufferInd(j),i)
        end do
        v(i) = log(sum(test))+scale
        !v(i) /= v(i) .or.
        if (v(i)<-huge(v(i))) then
            v(i)=-huge(v(i))
            !continue
        end if
    end do


    if ( ixy == numPointsY .AND. ixA == numPointsA .AND. ixAIME == numAIME .AND. ixt==1 .AND. ixType == 4) then
        write (*,*) morethan
    end if
    contains
    !!function
    !function func(p) result(res)
    !use Header
    !implicit none
    !!input
    !real(kind=rk), intent(in) :: p(:)
    !!output
    !real(kind=rk), allocatable :: res(:)
    !!local
    !real(kind=hp) :: systemEq(dim),test(dim), loc!,grids%supportSPA(ixt)
    !integer:: d, y, locl, loca1, j
    !
    !allocate(res(dim))
    !locl =1
    !locA1 = 0
    !do d=1,dimD
    !    locA1 = locA1 + 1
    !    if (locA1 > maxa(locl)) then
    !        locl = locl+1
    !        locA1=1
    !    end if
    !    do y=1,grids%supportSPA(ixt)
    !        GrossUtil(d,y) = const(locl,locA1,y)*(dot_product(grids%posteriorSPA(ixt,numPointsSPA-grids%supportSPA(ixt)+1:numPointsSPA),p((d-1)*grids%supportSPA(ixt)+1:d*grids%supportSPA(ixt))))
    !    end do
    !end do
    !
    !do d=1,dimD
    !    do y=1,grids%supportSPA(ixt)
    !        res((d-1)*grids%supportSPA(ixt)+y) = GrossUtil(d,y) / sum((GrossUtil(:,y)))
    !    end do
    !end do
    !
    !end function
    !!Jacoben
    !FUNCTION Jacob(p)
    !USE header
    !IMPLICIT NONE
    !REAL(kind=rk), DIMENSION(:), INTENT(IN) :: p
    !REAL(kind=rk), DIMENSION(size(p),size(p)) :: Jacob
    !
    !real(kind=rk) :: l1, l2 ,l3, l4
    !integer :: row, column
    !
    !integer :: d1, y1, d2, y2, locl, locA1
    !!shouldn't do this twice
    !locl =1
    !locA1 = 0
    !do d1=1,dimD
    !    locA1 = locA1 + 1
    !    if (locA1 > maxa(locl)) then
    !        locl = locl+1
    !        locA1=1
    !    end if
    !    do y1=1,grids%supportSPA(ixt)
    !        GrossUtil(d1,y1) = const(locl,locA1,y1)*&
    !            (dot_product(grids%posteriorSPA(ixt,numPointsSPA-grids%supportSPA(ixt)+1:numPointsSPA),p((d1-1)*grids%supportSPA(ixt)+1:d1*grids%supportSPA(ixt))))
    !    end do
    !end do
    !locl =1
    !locA1 = 0
    !do d1=1,dimD
    !    locA1 = locA1 + 1
    !    if (locA1 > maxa(locl)) then
    !        locl = locl+1
    !        locA1=1
    !    end if
    !    do y1=1,grids%supportSPA(ixt)
    !        do d2=1,dimD
    !            do y2=1,grids%supportSPA(ixt)
    !                row = (d1-1)*grids%supportSPA(ixt)+y1
    !                column = (d2-1)*grids%supportSPA(ixt)+y2
    !                !locl=(d1-1)/(maxmaxA)+1
    !                !locA1=mod(d1-1,maxmaxA)+1
    !                l1= sum(const(locl,locA1,:))
    !                l2 = dot_product(grids%posteriorSPA(ixt,numPointsSPA-grids%supportSPA(ixt)+1:numPointsSPA),p((d1-1)*grids%supportSPA(ixt)+1:d1*grids%supportSPA(ixt)))
    !                l3 = sum(GrossUtil(d1,:))
    !                l4 = l3**2
    !                jacob((d1-1)*grids%supportSPA(ixt)+y1,(d2-1)*grids%supportSPA(ixt)+y2) = &
    !                    (-sum(const(locl,locA1,:))*dot_product(grids%posteriorSPA(ixt,numPointsSPA-grids%supportSPA(ixt)+1:numPointsSPA),p((d1-1)* &
    !                    grids%supportSPA(ixt)+1:d1*grids%supportSPA(ixt)))/(sum(GrossUtil(y1,:)))**2)
    !                if (d1==d2) then
    !                    jacob((d1-1)*grids%supportSPA(ixt)+y1,(d2-1)*grids%supportSPA(ixt)+y2) = jacob((d1-1)*grids%supportSPA(ixt)+y1,(d2-1)*grids%supportSPA(ixt)+y2) + &
    !                        + 1 /(sum(GrossUtil(y1,:)))
    !                end if
    !                jacob((d1-1)*grids%supportSPA(ixt)+y1,(d2-1)*grids%supportSPA(ixt)+y2)= const(locl,locA1,y1)*grids%posteriorSPA(ixt,numPointsSPA-y2+1)* &
    !                    jacob((d1-1)*grids%supportSPA(ixt)+y1,(d2-1)*grids%supportSPA(ixt)+y2)
    !            end do
    !        end do
    !    end do
    !end do
    !
    !END FUNCTION Jacob

    FUNCTION funcv(x)
    USE header
    IMPLICIT NONE
    REAL(kind=rk), DIMENSION(:), INTENT(IN) :: x
    REAL(kind=rk), DIMENSION(size(x)) :: funcv
    !funcv = x - func(x)
    END FUNCTION funcv
    end subroutine
    !! ---------------------------------------------------------------------------------------------------------!
    !!!solve period
    !!! when RI is imporant
    !! ---------------------------------------------------------------------------------------------------------!
    !subroutine solvePeriodRI(params, grids, ixy, ixA, ixAIME ,ixt,ixType,spa,  EV1, ccpOut, v, q, uConst)
    !implicit none
    !
    !!input
    !real(kind=rk), intent(in) ::  EV1(:,:,:)
    !type (structparamstype), intent(in) :: params
    !type (gridsType), intent(in) :: grids
    !integer, intent(in) :: ixt,ixType,ixa,ixAIME,ixy, spa
    !
    !!!output
    !real(kind=rk), intent(out) :: ccpOut(:,:,:), v(:), q(:,:), uConst(:,:,:)
    !!local
    !!integer, parameter :: hp = selected_real_kind(16)
    !integer :: ixl,  ixA1, i, j, iter, dim, dimD ,K, labourChoices, offset, testInt, maxmaxA !A1
    !integer, allocatable :: maxA(:)
    !real(kind=rk) :: Y,  ubA1, A,AIME, va1, vb1,  check,EVloc(numAIME), checkSave(grids%supportSPA(ixt)), checkOther(grids%supportSPA(ixt)) !cons
    !real(kind=rk) :: scale, toputil
    !real(kind=rk), allocatable :: values(:), ccpMat(:,:,:), eye(:,:)
    !real(kind=rk), allocatable :: const(:,:,:), GrossUtil(:,:)
    !real(kind=rk), allocatable  ::  parameters(:), ccp(:), test(:), utilVec(:)
    !integer, save :: saveDim(2,numPointsL), lastixT = 0, unemp, div, d, locl, locA1
    !integer, allocatable :: rang(:)
    !real(kind=rk), allocatable, save :: locinitialGuessRI(:,:)
    !logical :: converged, unchanged
    !character (len=2000) :: text(2)
    !
    !!temp
    !real(kind=rk), allocatable :: tempC1(:,:), tempC2(:,:)
    !real(kind=rk) :: diff, lowestUtil
    !
    !REAL(kind=rk) :: error
    !LOGICAL :: logcheck
    !ixa1 = 0
    !lowestUtil = 0.0
    !!If suffered unemployment shock then can't choose to work
    !labourChoices=mod(ixy,2)*numPointsL+(1-mod(ixy,2))*1
    !allocate(maxA(labourChoices), const(labourChoices,numPointsA,grids%supportSPA(ixt)))
    !maxA(1) = numPointsA !default to all asset level possible
    !!if (ixy == 4 .AND. ixt ==14 .AND. ixAIME == 2) then
    !!    write(*,*) 'stop'
    !!end if
    !do ixL = 0,(labourChoices-1),1           ! points of income choice
    !    ! Value of income and information for optimisation
    !    if (ixa1 > numPointsA) maxA(ixl) = numPointsA
    !    A    = grids%Agrid(ixType,ixt, ixA)            ! assets today
    !
    !    Y    = grids%Ygrid(ixType,ixt, ixY)
    !    Y    = ixL*Y
    !
    !    AIME = grids%AIMEgrid(ixType,ixt,ixAIME)
    !
    !    call nextAIME(grids,params,ixt,ixType,ixl,y,AIME)
    !
    !    Y    = ixL*Y+ abs(mod(ixy,2)==0 )*grids%benefit(ixt)
    !
    !    call gross2net(params,grids,Y,ixt,ixl,ixType,AIME,spa)
    !
    !    ubA1 = (A + Y - params%minCons)*(1+params%r);    ! upper bound: assets tomorrow
    !    do ixA1 = 1, numPointsA
    !        if (grids%Agrid(ixType,ixt+1, ixA1) > ubA1) then
    !            const(ixl+1,ixa1:numPointsA,:) = -huge(const(ixl+1,ixa1,1))
    !            uConst(:,ixl+1,ixa1:numPointsA) = -huge(const(ixl+1,ixa1,1))
    !            maxA(ixl+1) = ixa1 - 1
    !            exit
    !        end if
    !
    !        do i = 1, grids%supportSPA(ixt)-1
    !            EVloc =  (1-params%p)*EV1(ixA1,:,numPointsSPA-grids%supportSPA(ixt)+i)+params%p*EV1(ixA1,:,numPointsSPA-grids%supportSPA(ixt)+i+1)
    !            const(ixl+1,ixa1,i) = objectivefunc(params, grids,grids%Agrid(ixType,ixt+1, ixA1), A, Y,ixL,ixt,ixType, AIME,EVloc)
    !            !if (const(ixl+1,ixa1,i) <lowestUtil) lowestUtil = const(ixl+1,ixa1,i)
    !            uConst(numPointsSPA-grids%supportSPA(ixt)+i,ixl+1,ixa1) = const(ixl+1,ixa1,i)
    !        end do
    !        EVloc =  EV1(ixA1,:,numPointsSPA)
    !        const(ixl+1,ixa1,i) = objectivefunc(params, grids,grids%Agrid(ixType,ixt+1, ixA1), A, Y,ixL,ixt,ixType, AIME,EVloc)
    !        uConst(numPointsSPA-grids%supportSPA(ixt)+i,ixl+1,ixa1) = const(ixl+1,ixa1,i)
    !        !if (const(ixl+1,ixa1,i) < lowestUtil) lowestUtil = const(ixl+1,ixa1,i)
    !    end do
    !
    !end do
    !if (ixa1 > numPointsA) maxA(labourChoices) = numPointsA
    !maxmaxA = maxval(maxA)
    !
    !!allocate(tempC1(maxA(1),grids%supportSPA(ixt)))
    !!tempC1 = const(1,1:maxA(1),:)
    !!if (labourChoices==2) then
    !!    allocate(tempC2(maxA(2),grids%supportSPA(ixt)))
    !!    tempC2 = const(2,1:maxA(2),:)
    !!    if (rank==0) write(*,*) max(maxval(tempC1),maxval(tempC2)) - min(minval(tempC1),minval(tempC2))
    !!    diff = max(maxval(tempC1),maxval(tempC2)) - min(minval(tempC1),minval(tempC2))
    !!else
    !!    if (rank==0) write(*,*) maxval(tempC1) - minval(tempC1)
    !!    diff = maxval(tempC1) - minval(tempC1)
    !!end if
    !scale = maxval(const)
    !const = const-scale
    !const=exp(const)
    !
    !
    !dimD = sum(maxA)!labourChoices*(ixa1-1)
    !dim = grids%supportSPA(ixt)*dimD !*labourChoices*(ixa1-1)
    !allocate(values(dim+1),parameters(dim),ccp(dim),ccpMat(grids%supportSPA(ixt),labourChoices,maxmaxA), GrossUtil(dimD,grids%supportSPA(ixt)),test(dimD),utilVec(dim),rang(dimd))
    !!copy into util vec for comparision with edge solutions
    !ixl =1
    !ixa1 = 0
    !do d=1,dimD
    !    ixa1 = ixa1 + 1
    !    if ( ixa1 > maxa(ixl)) then
    !        ixl = ixl+1
    !        ixa1=1
    !    end if
    !    do i = 1, grids%supportSPA(ixt)
    !        utilVec((d-1)*grids%supportSPA(ixt)+i) = uConst(numPointsSPA-grids%supportSPA(ixt)+i,ixl,ixa1)
    !    end do
    !end do
    !
    !!Default ccp to zero for regions that are ignored
    !ccp = 0.0
    !!Want to save grid to reuse over employment states
    !unemp = mod(ixy-1,2)+1
    !if (ixt/=lastixT) then
    !    parameters = grids%initialGuessRI(1,1:dim)
    !    if (allocated(locinitialGuessRI)) deallocate(locinitialGuessRI)
    !    allocate(locinitialGuessRI(2,grids%supportSPA(ixt)*numPointsL*numPointsA))
    !else if (ixy<=2) then
    !    parameters = grids%initialGuessRI(unemp,1:dim)
    !else
    !    parameters = grids%initialGuessRI(unemp,1:dim)
    !    !parameters=0.0
    !    !offset = 0
    !    !do k=1,labourChoices
    !    !    if (saveDim(unemp,k) > grids%supportSPA(ixt)*(maxmaxA) ) saveDim(unemp,k) = grids%supportSPA(ixt)*maxa(k)
    !    !    parameters(offset+1:offset+saveDim(unemp,k))= &
    !    !        locinitialGuessRI(unemp,(k-1)*grids%supportSPA(ixt)*numPointsA+1:(k-1)*grids%supportSPA(ixt)*numPointsA+saveDim(unemp,k))
    !    !    offset = offset + saveDim(unemp,k)
    !    !end do
    !    !!Make sure we don't pass zeros across as gets stuck at them if give them as initial conds
    !    !parameters = parameters + 0.1/dim
    !end if
    !!store last calculated period to know when we change period
    !lastixT = ixt
    !
    !!parameters = grids%initialGuessRI(1,1:dim)
    !!if (mod(ixy,2)==0) then
    !!    parameters(grids%supportSPA(ixt)*(ixa1-1)+1:dim)=0.0
    !!end if
    !
    !!Normalise guess to be prob dist
    !do i=1,grids%supportspa(ixt)
    !    check = sum(parameters((/(i+j*grids%supportspa(ixt),j=0,dimd-1)/)))
    !    parameters((/(i+j*grids%supportspa(ixt),j=0,dimd-1)/))=parameters((/(i+j*grids%supportspa(ixt),j=0,dimd-1)/))*check**-1
    !end do
    !
    !check = 1.0
    !iter=0
    !converged = .FALSE.
    !div = 0
    !do while (.not. converged)
    !    !Convergence criteria should depend of cost of attention otherwise utility can get very small and it just assigns equiprobable dist
    !    do while (check >= (0.0001/10**div) ) !0.0001 !0.001!0.005!0.0001!0.00005 !0.00001 !/params%lambda *params%lambda
    !        iter=iter+1
    !        ccp = func(parameters)
    !        check = sum(abs(ccp-parameters))
    !        parameters=ccp
    !    end do
    !    !ccp = parameters
    !    !call newt(funcv,jacob,ccp,0.001_rk,logcheck,error)
    !
    !    !Check is good enough solution
    !    do i=1,grids%supportSPA(ixt)
    !        !write (text(i),*) (i+j*grids%supportSPA(ixt),j=0,Dimd-1)
    !        checkSave(i) = sum(ccp((/(i+j*grids%supportSPA(ixt),j=0,Dimd-1)/))) !grids%supportSPA(ixt)*maxa(i)
    !        if (checkSave(i) > 1.0001 .OR. checkSave(i)<0.9999) then
    !            write (*,*) "CCP doesn't sum to 1", checkSave(i), params%lambda, ixt
    !            stop
    !        else
    !            if (div > 0) write (*,*) "passed!", checkSave(i)
    !            converged = .TRUE.
    !        end if
    !    end do
    !    div = div+1
    !end do
    !
    !
    !!!compare to "edge" solutions 9should be corner but too hard for now
    !topUtil = 0
    !do i = 1,grids%supportSPA(ixt)! numPointsSPA-grids%supportSPA(ixt)+1,numPointsSPA
    !    rang = (/(i+j*grids%supportSPA(ixt),j=0,dimd-1)/)
    !    topUtil = topUtil + grids%posteriorSPA(ixt,numPointsSPA-grids%supportSPA(ixt)+i)*dot_product(ccp((/(i+j*grids%supportSPA(ixt),j=0,dimd-1)/)),utilvec((/(i+j*grids%supportSPA(ixt),j=0,dimd-1)/))) ! utilvec((i-1)*grids%supportSPA(ixt)+1:i*grids%supportSPA(ixt)))
    !end do
    !unchanged = .true.
    !locl =1
    !locA1 = 0
    !do i=1,dimD
    !    locA1 = locA1 + 1
    !    if (locA1 > maxa(locl)) then
    !        locl = locl+1
    !        locA1=1
    !    end if
    !    if (dot_product(grids%posteriorSPA(ixt,numPointsSPA-grids%supportSPA(ixt)+1:numPointsSPA), utilvec((i-1)*grids%supportSPA(ixt)+1:i*grids%supportSPA(ixt))) > toputil  ) then
    !        topUtil = dot_product(grids%posteriorSPA(ixt,numPointsSPA-grids%supportSPA(ixt)+1:numPointsSPA), utilvec((i-1)*grids%supportSPA(ixt)+1:i*grids%supportSPA(ixt)))
    !        ccp = 0.0
    !        ccp((i-1)*grids%supportSPA(ixt)+1:i*grids%supportSPA(ixt)) = 1.0
    !        unchanged = .false.
    !        GrossUtil = 0.0
    !        do j = 1,grids%supportSPA(ixt)
    !            GrossUtil(i,j) = const(locl,locA1,j)
    !        end do
    !    end if
    !end do
    !if (unchanged .AND. maxval(ccp)<0.999) then
    !    write (*,*) 'Unchanged!'
    !
    !    !stop
    !end if
    !
    !!!Cache for use on next round
    !if (ixy < numPointsY) then
    !    saveDim(unemp,:) = grids%supportSPA(ixt)*maxa
    !    locinitialGuessRI(unemp,:)=0.0
    !    offset = 0
    !    do k=1,labourChoices
    !        locinitialGuessRI(unemp,(k-1)*grids%supportSPA(ixt)*numPointsA+1:(k-1)*grids%supportSPA(ixt)*numPointsA+saveDim(unemp,k))= ccp((k-1)*offset+1:k*saveDim(unemp,k))
    !        offset = offset + saveDim(unemp,k)
    !    end do
    !end if
    !
    !ccpOut=0.0
    !ccpMat = 0.0
    !!ccpMat=reshape(ccp, (/grids%supportSPA(ixt),labourChoices,ixa1-1/))
    !offset = 0
    !do i=1,labourChoices
    !    do j =1,maxmaxA
    !        if (j > maxa(i)) then
    !            ccpMat(:,i,j) = 0.0
    !            exit
    !        end if
    !        do k=1,grids%supportSPA(ixt)
    !            ccpMat(k,i,j) =  ccp(k+(j-1)*grids%supportSPA(ixt)+offset)
    !        end do
    !    end do
    !    offset = offset+saveDim(unemp,i)
    !end do
    !ccpOut(1:grids%supportSPA(ixt),1:labourChoices,1:maxmaxA)=ccpMat
    !!Get defaults actiosn
    !do i=1,labourChoices
    !    do j =1,maxmaxA
    !        q(i,j) = dot_product(ccpOut(:,i,j),grids%posteriorSPA(ixt,numPointsSPA-grids%supportSPA(ixt)+1:numPointsSPA))
    !    end do
    !end do
    !
    !do i=1,grids%supportSPA(ixt)
    !    checkOther(i) = sum(ccpOut(i,:,:))
    !    !write (*,*) checkSave(i), checkOther(i)
    !    if (checkOther(i) > 1.05 .OR. checkOther(i)<0.95) then
    !        write (*,*) "CCP doesn't sum to 1", checkOther(i), i, converged,  ixy, ixA, ixAIME, ixt, unemp
    !        write (*,*) "converg check", checksave
    !        write (*,*) "converg check", checkother
    !        write (*,*) "ccp", ccp
    !        write (*,*) "ccpout", ccpOut
    !        stop
    !    end if
    !    !ccp((/(i+j*grids%supportSPA(ixt),j=0,Dimd-1)/))=ccp((/(i+j*grids%supportSPA(ixt),j=0,Dimd-1)/))*check**-1
    !end do
    !!if (max(checkOther(1),checkOther(2)) /= max(checkSave(1), checkSave(2)) .OR. min(checkOther(1),checkOther(2)) /= min(checkSave(1), checkSave(2)) )pause
    !!Calculate continuation value
    !do i=1,grids%supportSPA(ixt)
    !    test= GrossUtil(:,i)
    !    v(i) = log(sum(test))+scale
    !end do
    !
    !contains
    !!function
    !function func(p) result(res)
    !use Header
    !implicit none
    !!input
    !real(kind=rk), intent(in) :: p(:)
    !!output
    !real(kind=rk), allocatable :: res(:)
    !!local
    !real(kind=hp) :: systemEq(dim),test(dim), loc!,grids%supportSPA(ixt)
    !integer:: d, y, locl, loca1, j
    !
    !allocate(res(dim))
    !locl =1
    !locA1 = 0
    !do d=1,dimD
    !    locA1 = locA1 + 1
    !    if (locA1 > maxa(locl)) then
    !        locl = locl+1
    !        locA1=1
    !    end if
    !    do y=1,grids%supportSPA(ixt)
    !        GrossUtil(d,y) = const(locl,locA1,y)*(dot_product(grids%posteriorSPA(ixt,numPointsSPA-grids%supportSPA(ixt)+1:numPointsSPA),p((d-1)*grids%supportSPA(ixt)+1:d*grids%supportSPA(ixt))))
    !    end do
    !end do
    !
    !do d=1,dimD
    !    do y=1,grids%supportSPA(ixt)
    !        res((d-1)*grids%supportSPA(ixt)+y) = GrossUtil(d,y) / sum((GrossUtil(:,y)))
    !    end do
    !end do
    !
    !end function
    !!Jacoben
    !FUNCTION Jacob(p)
    !USE header
    !IMPLICIT NONE
    !REAL(kind=rk), DIMENSION(:), INTENT(IN) :: p
    !REAL(kind=rk), DIMENSION(size(p),size(p)) :: Jacob
    !
    !real(kind=rk) :: l1, l2 ,l3, l4
    !integer :: row, column
    !
    !integer :: d1, y1, d2, y2, locl, locA1
    !!shouldn't do this twice
    !locl =1
    !locA1 = 0
    !do d1=1,dimD
    !    locA1 = locA1 + 1
    !    if (locA1 > maxa(locl)) then
    !        locl = locl+1
    !        locA1=1
    !    end if
    !    do y1=1,grids%supportSPA(ixt)
    !        GrossUtil(d1,y1) = const(locl,locA1,y1)*&
    !            (dot_product(grids%posteriorSPA(ixt,numPointsSPA-grids%supportSPA(ixt)+1:numPointsSPA),p((d1-1)*grids%supportSPA(ixt)+1:d1*grids%supportSPA(ixt))))
    !    end do
    !end do
    !locl =1
    !locA1 = 0
    !do d1=1,dimD
    !    locA1 = locA1 + 1
    !    if (locA1 > maxa(locl)) then
    !        locl = locl+1
    !        locA1=1
    !    end if
    !    do y1=1,grids%supportSPA(ixt)
    !        do d2=1,dimD
    !            do y2=1,grids%supportSPA(ixt)
    !                row = (d1-1)*grids%supportSPA(ixt)+y1
    !                column = (d2-1)*grids%supportSPA(ixt)+y2
    !                !locl=(d1-1)/(maxmaxA)+1
    !                !locA1=mod(d1-1,maxmaxA)+1
    !                l1= sum(const(locl,locA1,:))
    !                l2 = dot_product(grids%posteriorSPA(ixt,numPointsSPA-grids%supportSPA(ixt)+1:numPointsSPA),p((d1-1)*grids%supportSPA(ixt)+1:d1*grids%supportSPA(ixt)))
    !                l3 = sum(GrossUtil(d1,:))
    !                l4 = l3**2
    !                jacob((d1-1)*grids%supportSPA(ixt)+y1,(d2-1)*grids%supportSPA(ixt)+y2) = &
    !                    (-sum(const(locl,locA1,:))*dot_product(grids%posteriorSPA(ixt,numPointsSPA-grids%supportSPA(ixt)+1:numPointsSPA),p((d1-1)* &
    !                    grids%supportSPA(ixt)+1:d1*grids%supportSPA(ixt)))/(sum(GrossUtil(y1,:)))**2)
    !                if (d1==d2) then
    !                    jacob((d1-1)*grids%supportSPA(ixt)+y1,(d2-1)*grids%supportSPA(ixt)+y2) = jacob((d1-1)*grids%supportSPA(ixt)+y1,(d2-1)*grids%supportSPA(ixt)+y2) + &
    !                        + 1 /(sum(GrossUtil(y1,:)))
    !                end if
    !                jacob((d1-1)*grids%supportSPA(ixt)+y1,(d2-1)*grids%supportSPA(ixt)+y2)= const(locl,locA1,y1)*grids%posteriorSPA(ixt,numPointsSPA-y2+1)* &
    !                    jacob((d1-1)*grids%supportSPA(ixt)+y1,(d2-1)*grids%supportSPA(ixt)+y2)
    !            end do
    !        end do
    !    end do
    !end do
    !
    !END FUNCTION Jacob
    !
    !FUNCTION funcv(x)
    !USE header
    !IMPLICIT NONE
    !REAL(kind=rk), DIMENSION(:), INTENT(IN) :: x
    !REAL(kind=rk), DIMENSION(size(x)) :: funcv
    !!funcv = x - func(x)
    !END FUNCTION funcv
    !end subroutine


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
    integer :: i !, seedIn, n
    real (kind=rk) ::  uniformRand(dimEstimation+1)
    !INTEGER, DIMENSION(:), ALLOCATABLE :: seed


    !!get uniform random number
    CALL RANDOM_NUMBER(uniformRand)
    uniformRand = -0.5+uniformRand

    if (rank==0) write (*,*) "Guess ", 1
    p(1,1) = params%nu
    p(1,2) = params%beta
    p(1,3) = params%gamma
    p(1,4) = params%thetab
    if (dimEstimation ==5) p(1,5) = params%lambda
    y(1) = gmm_criteria(p(1,:))

    do i = 2, dimEstimation+1
        if (rank==0) then
            write (*,*) "Guess ", i
        end if
        p(i,1) = params%nu*(1+uniformRand(i))
        p(i,2) = (0.985+0.1*uniformRand(i))!params%beta*(1+uniformRand(i))
        p(i,3) = params%gamma*(1+uniformRand(i))
        p(i,4) = params%thetab*(1+uniformRand(i))
        if (dimEstimation ==5) p(i,5) = params%lambda*(1+uniformRand(i))
        y(i) = gmm_criteria(p(i,:))
    end do

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
    if (minval(control(1:4)) < 0 )  then
        gmm_criteria = huge(gmm_criteria)
        return
    end if

    params%nu = control(1)
    params%beta = control(2)
    params%gamma = control(3)
    params%thetab = control(4)
    if (dimEstimation == 5) params%lambda = control(5)
    gmm_criteria = gmm(params,grids,moments,weights)

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
    real (kind=rk) :: workarea(4) !reemployRaw(4,10), unemployRaw(4,10),
    integer :: i, ios

    !    unemployRaw(1,:) = (/0.227, 0.049, 0.027, 0.014, 0.004, 0.007, 0.013, 0.017, 0.007, 0.01/)
    !    unemployRaw(2,:) = (/0.099, 0.026, 0.01, 0.013, 0.008, 0.012, 0.003, 0.003, 0.002, 0.0/)
    !    unemployRaw(3,:) = (/0.201, 0.024, 0.054, 0.032, 0.029, 0.038, 0.02, 0.005, 0.004, 0.01/)
    !    unemployRaw(4,:) = (/0.063, 0.016, 0.01, 0.006, 0.003, 0.002, 0.004, 0.001, 0.003, 0.0/)
    !
    !#if PROD_SIZE == 5
    !    !if (numPointsProd==5) then
    !    do i=1,numPointsType
    !        grids%unemploy(i,1) =unemployRaw(i,1)+unemployRaw(i,2)!(/0.1762, 0.0415,0.0257,0.0165,0.0090/)
    !        grids%unemploy(i,2) =unemployRaw(i,3)+unemployRaw(i,4)
    !        grids%unemploy(i,3) =unemployRaw(i,5)+unemployRaw(i,6)
    !        grids%unemploy(i,4) =unemployRaw(i,7)+unemployRaw(i,8)
    !        grids%unemploy(i,5) =unemployRaw(i,9)+unemployRaw(i,10)
    !    end do
    !#elif  PROD_SIZE == 10
    !    !else if (numPointsProd==10) then
    !    grids%unemploy = unemployRaw
    !    !end if
    !#endif
    !    reemployRaw(1,:) = (/0.13, 0.043, 0.022, 0.005, 0.002, 0.007, 0.01, 0.002, 0.0, 0.007 /)
    !    reemployRaw(2,:) = (/0.174, 0.005, 0.013, 0.017, 0.013, 0.017, 0.013, 0.003, 0.003, 0.0 /)
    !    reemployRaw(3,:) = (/0.118, 0.029, 0.037, 0.011, 0.008, 0.013, 0.003, 0.0, 0.005, 0.0 /)
    !    reemployRaw(4,:) = (/0.238, 0.032, 0.02, 0.008, 0.008, 0.012, 0.008, 0.004, 0.0, 0.004/)
    !
    !#if PROD_SIZE == 5
    !    !if (numPointsProd==5) then
    !    !grids%reemploy = (/ 0.1923, 0.0333, 0.0200, 0.0108,0.0047/)
    !    do i=1,numPointsType
    !        grids%reemploy(i,1) =reemployRaw(i,1)+reemployRaw(i,2)!(/0.1762, 0.0415,0.0257,0.0165,0.0090/)
    !        grids%reemploy(i,2) =reemployRaw(i,3)+reemployRaw(i,4)
    !        grids%reemploy(i,3) =reemployRaw(i,5)+reemployRaw(i,6)
    !        grids%reemploy(i,4) =reemployRaw(i,7)+reemployRaw(i,8)
    !        grids%reemploy(i,5) =reemployRaw(i,9)+reemployRaw(i,10)
    !    end do
    !#elif  PROD_SIZE == 10
    !    !else if (numPointsProd==10) then
    !    grids%reemploy = reemployRaw
    !    !end if
    !#endif


    open (unit = 1001,file=trim(pathInputs) // 'unemploymentTrans.txt', action='read', IOSTAT = ios)
    do i=1,numPointsType
        read (1001, *,  IOSTAT = ios) workarea
        grids%unemploy(i,:) = workarea(3) /100.0
        read (1001, *,  IOSTAT = ios) workarea
        grids%reemploy(i,:) = workarea(2) / 100.0
    end do
    close (unit = 1001)

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
    !REAL(kind=rk) :: temp(85),  temp2(86), temp3(numSims, Tperiods)
    integer :: ios, I
    real (kind=rk) ::   conditioningProb, sig_inc !sig_initial,
    real (kind=rk) :: e(Tperiods, numSims)   ! the innovations to log income
    real (kind=rk) :: logy1(1,numSims)        ! draws for initial income
    real (kind=rk) :: ly(Tperiods, numSims)           ! log income
    integer :: s, t,  workAge
    real (kind=rk) ::  uniformRand(Tperiods,numSims), workArea(8), spIncprm(3,5)
    !INTEGER, DIMENSION(:), ALLOCATABLE :: seed
    integer ::  n, productiv(1), typeSim, startN, endN, samples !requiredl,seedIn,finish,
    real (kind=rk), allocatable :: eTemp(:,:), lyTemp(:)
    logical :: unemployed
    character(len=1024) :: outFile
    character (100) :: tempC
    !integer, allocatable :: rangeSims(:)
    real (kind=rk) :: mortalmen(Tperiods)

    open (unit = 1001,file=trim(pathMoments)//'InitialAssets.txt', action='read', IOSTAT = ios)
    read (1001, *,  IOSTAT = ios) grids%initialAssets
    close (unit=1001)


    !write (*,*)'Mean emp int assets', sum(grids%initialassets(:,2))/real(obsInitDist,rk) !

    grids%fc = 0.0

    do n = 1, numPointsType
        !Get stochastic earning process
        write (tempC,'(I3)') n
        outFile = trim(pathInputs) // 'estparams_'// trim(adjustl(tempC)) // '.txt'
        outfile=ADJUSTL(outfile)
        open (unit = 1001,file=outfile, action='read', IOSTAT = ios)
        read (1001, *,  IOSTAT = ios) workArea
        close (unit=1001)
        params%rho(n) = workArea(1)
        params%sigma(n) = workArea(2)**0.5
        params%sigma_int(n) = workArea(4)**0.5

        !Get deterministic earning process
        outFile = trim(pathInputs) // 'wageReg'// trim(adjustl(tempC)) // '.txt'
        outfile=ADJUSTL(outfile)
        open (unit = 1001,file=outfile, action='read', IOSTAT = ios)
        read (1001, *,  IOSTAT = ios) workArea
        close (unit=1001)
        params%delta(n,2) = workArea(1)
        params%delta(n,1) = workArea(2)
        params%delta(n,3) = workArea(3)

        if (mod(n,2)==1) then
            !Get deterministic earning process
            outFile = trim(pathInputs) // 'SpousalIncome'// trim(adjustl(tempC)) // '.txt'
            outfile=ADJUSTL(outfile)
            open (unit = 1001,file=outfile, action='read', IOSTAT = ios)
            read (1001, *,  IOSTAT = ios) workArea
            close (unit=1001)
            spIncprm(n,4) = workArea(1)
            spIncprm(n,3) = workArea(2)
            spIncprm(n,2) = workArea(3)
            spIncprm(n,1) = workArea(4)
            spIncprm(n,5) = workArea(5)
        end if
    end do

    outFile = trim(pathInputs) // 'type.txt'
    outfile=ADJUSTL(outfile)
    open (unit = 1001,file=outfile, action='read', IOSTAT = ios)
    do n = 1, numPointsType
        read (1001, *,  IOSTAT = ios) params%fracType(n)
    end do
    close (unit=1001)

    !Get entitlement age for private pension
    outFile = trim(pathInputs) // 'ppAge.txt'
    outfile=ADJUSTL(outfile)
    open (unit = 1001,file=outfile, action='read', IOSTAT = ios)
    do n = 1, numPointsType
        read (1001, *,  IOSTAT = ios) params%ERA(n)
        params%ERA(n) = params%ERA(n) -startAge+1
    end do
    close (unit=1001)

    !get db pension level
    outFile = trim(pathInputs) // 'privatePen.txt'
    outfile=ADJUSTL(outfile)
    open (unit = 1001,file=outfile, action='read', IOSTAT = ios)
    read (1001, *,  IOSTAT = ios) workArea
    close (unit=1001)

    params%db(1,1) = workArea(1)
    params%db(2,1) = workArea(2)
    params%db(3,1) = workArea(3)
    params%db(4,1) = workArea(4)
    params%db(1,2) = workArea(5)
    params%db(2,2) = workArea(6)
    params%db(3,2) = workArea(7)
    params%db(4,2) = workArea(8)

    !get db pension level
    outFile = trim(pathInputs) // 'statePen.txt'
    outfile=ADJUSTL(outfile)
    open (unit = 1001,file=outfile, action='read', IOSTAT = ios)
    read (1001, *,  IOSTAT = ios) workArea
    close (unit=1001)

    params%pension(1,1) = workArea(1)
    params%pension(2,1) = workArea(2)
    params%pension(3,1) = workArea(3)
    params%pension(4,1) = workArea(4)
    params%pension(1,2) = workArea(5)
    params%pension(2,2) = workArea(6)
    params%pension(3,2) = workArea(7)
    params%pension(4,2) = workArea(8)

    params%mu = 0
    params%r= 0.04
    params%tol = 1e-10
    params%minCons =  50.0 !1e-5 !
    params%percentCons = 0.0 !0.1
    params%mu = 0
    params%hrsWrk = real(weeksWorking,rk)/52.0*(40.0/(16.0*7)) ! 0.3159
    params%p = 0.06

    call getIncomeGrid(params, grids)

    !need to redo surival
    !temp2 = (/0.0, 0.0, 0.0, 0.0, 0.0, &
    !    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
    !    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.011071, 0.011907, 0.012807, &
    !    0.013676, 0.01475, 0.015818, 0.016846, 0.018079, 0.019343, 0.020659, 0.0218, 0.023505, 0.025202, 0.02696, &
    !    0.028831, 0.031017, 0.033496, 0.036024, 0.038912, 0.042054, 0.045689, 0.049653, 0.054036, 0.05886, 0.064093, &
    !    0.069636, 0.07533, 0.081069, 0.086912, 0.093067, 0.099807, 0.107483, 0.116125, 0.125196, 0.134361, 0.143881, &
    !    0.1542, 0.165675, 0.17842, 0.192363, 0.2117, 0.1672,0.1565, 0.1485,0.1459, 1.0/);
    !grids%mortal = temp2(startAge-20+1:86)
    outFile = trim(pathInputs) // 'femaleCondSurive.txt'
    outfile=ADJUSTL(outfile)
    open (unit = 1001,file=outfile, action='read', IOSTAT = ios)
    read (1001, *,  IOSTAT = ios) grids%mortal
    close (unit=1001)
    grids%mortal(Tperiods+1) = 1.0

    startN = 1
    do n = 1, numPointsType
        endN =  floor(params%fracType(n)*numSims)
        samples = endN - StartN + 1
        allocate(lyTemp(samples),eTemp(Tperiods, samples))
        call genStdNorm(samples,lyTemp)
        lyTemp = params%sigma_int(n)*lyTemp
        call genStdNorm(samples*Tperiods,eTemp)
        eTemp = params%sigma(n)*eTemp
        logy1(1,startN:endN) = lyTemp
        e(:,startN:endN) = eTemp
        deallocate(lyTemp,eTemp)
    end do

    CALL RANDOM_NUMBER(uniformRand)

    outFile = trim(pathInputs) // 'maleCondSurive.txt'
    outfile=ADJUSTL(outfile)
    open (unit = 1001,file=outfile, action='read', IOSTAT = ios)
    read (1001, *,  IOSTAT = ios) mortalmen
    close (unit=1001)
    mortalmen = 1.0 - mortalmen
    do s = 1, numSims, 1                           ! loop through individuals
        if (numPointsType == 1) then
            typeSim = 1
        else
            if (real(s,rk)/numSims <= params%fracType(1)) then
                typeSim = 1
            else if (real(s,rk)/numSims <= params%fracType(2)) then
                typeSim = 2
            else if (real(s,rk)/numSims <= params%fracType(3)) then
                typeSim =3
            else
                typeSim =4
            end if
        end if
        sig_inc = params%sigma(typeSim)/((1-params%rho(typeSim)**2)**0.5)
        unemployed = .FALSE.
        ! Get all the incomes, recursively
        ly(1, s) = truncate(logy1(1, s), -normBnd*sig_inc,normBnd*sig_inc )
        grids%Simy(1, s) = weeksYear*exp(ly(1, s)+params%delta(typeSim,1)*(startAge)**2+params%delta(typeSim,2)*(startAge)+params%delta(typeSim,3)-grids%fc(1))
        do t = 1,Tperiods-1,1                              ! loop through time periods for a particular individual
            workAge = startAge + t
            !gen spouse inc
            if (mod(typeSim,2)==1) then
                params%spouseInc(typeSim,t) = weeksYear*exp(spIncprm(typeSim,1)*workage**4+spIncprm(typeSim,2)*workage**3+spIncprm(typeSim,3)*workage**2+spIncprm(typeSim,4)*workage+spIncprm(typeSim,5))
            else
                params%spouseInc(typeSim,t) = 0.0
            end if
            if (t>stopwrok) params%spouseInc(typeSim,t) = params%spouseInc(typeSim,t-1)
            params%spouseInc(typeSim,t) = product(mortalmen(1:t-1))*params%spouseInc(typeSim,t)
            if (unemployed) then
                if (uniformRand(t,s) < grids%reemploy(typeSim,productiv(1))) then
                    unemployed = .FALSE.
                    ly(t+1, s) = (1 -params%rho(typeSim)) * params%mu + params%rho(typeSim) * ly(t, s) + e(t + 1, s)
                    ly(t+1, s) = truncate(ly(t+1, s), -normBnd*sig_inc,normBnd*sig_inc )
                    grids%Simy(t+1, s) = weeksYear*exp( ly(t+1, s) + params%delta(typeSim,1)*(workAge+1)**2+params%delta(typeSim,2)*(workAge+1)+params%delta(typeSim,3)-grids%fc(t) )
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
                    ly(t+1, s) = (1 -params%rho(typeSim)) * params%mu + params%rho(typeSim) * ly(t, s) + e(t + 1, s)
                    ly(t+1, s) = truncate(ly(t+1, s), -normBnd*sig_inc,normBnd*sig_inc )
                    grids%Simy(t+1, s) = weeksYear*exp( ly(t+1, s) + params%delta(typeSim,1)*(workAge+1)**2+params%delta(typeSim,2)*(workAge+1)+params%delta(typeSim,3)-grids%fc(t) )
                end if
            end if
        end do ! t
        params%spouseInc(typeSim,t) = mortalmen(Tperiods)*params%spouseInc(typeSim,t-1)
    end do! s

    params%spouseInc = params%spouseInc +params%minCons
    if (TendRI>Tretire) then
        grids%supportSPA(1:Tretire-1) = numPointsSPA
        do s = Tretire, TendRI
            grids%supportSPA(s) = grids%supportSPA(s-1) - 1
        end do
        allocate(grids%posteriorSPA(TendRI,numPointsType,numPointsA, numAIME,  numPointsY,numPointsSPA))
        grids%posteriorSPA=0.0
        do t = 1, TendRI
            workAge = startAge - 20 + t
            if (t<Tretire) then
                do s=0,(numPointsSPA -2)
                    grids%posteriorSPA(t,:, :, :, :,s+1) =  choose(workAge,s)*(params%p**s)*((1-params%p)**(workAge-s))
                end do
                grids%posteriorSPA(t,:, :,:,  :,numPointsSPA)=1-sum(grids%posteriorSPA(t,:, :, :,  :,:),5)
            else
                grids%posteriorSPA(t,:,:, :,  :,1:(t-Tretire+1))=0.0
                conditioningProb=0.0
                do s=0,(t-Tretire)
                    conditioningProb =  conditioningProb + choose(workAge,s)*(params%p**s)*((1-params%p)**(workAge-s))
                end do
                do s=(t-Tretire+1),(numPointsSPA -2)
                    grids%posteriorSPA(t,:,:, :,  :,s+1) =  (choose(workAge,s)*(params%p**s)*((1-params%p)**(workAge-s)))/(1-conditioningProb)
                end do
                grids%posteriorSPA(t,:, :, :,  :,numPointsSPA)=1-sum(grids%posteriorSPA(t,:, :, :,  :,:),5)
            end if
        end do
    end if
    params%spaTransMat = 0.0

    do i=1,numPointsSPA  -1
        params%spaTransMat(i,i) = 1.0 - params%p
        params%spaTransMat(i,i+1) = params%p
    end do
    params%spaTransMat(numPointsSPA,numPointsSPA) = 1.0
    params%spaTransMat = transpose(params%spaTransMat)

    CALL RANDOM_NUMBER(grids%initialGuessRI)

    call RANDOM_NUMBER(grids%logitShock)

    end subroutine

    !!-----------------------------------------------------------------------------------------------------------!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Law motion AIME
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!-----------------------------------------------------------------------------------------------------------!
    subroutine nextAIME(grids,params,ixt,typeSim,ixl,yin, AIME)
    implicit none

    !inputs
    type (structparamstype), intent(in) :: params
    type (gridsType), intent(in) :: grids
    integer, intent(in) :: ixt, typeSIm, ixl
    real(kind=rk), intent(in):: yin

    !changing
    real (kind=rk), intent(inout) :: AIME

    !local
    integer :: workAge

    workAge = startAge - 20 + ixt
    ! Next peridos AIME
    if (ixL==1 .AND.  ixt < params%ERA(typeSim) ) then ! spouseretire
        AIME =  Yin/workAge + AIME * (workAge-1)/workAge
    else if ((ixL==0 .AND.  ixt < params%ERA(typeSim) )) then !spouseretire
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
    Y    = ixL*Y
    AIME = grids%AIMEgrid(typeSim,ixt,ixAIME)
    !Also doesn't matter if income is arbitrary as in that case they are too old for AIME to increase
    call nextAIME(grids,params,ixt,typeSim,ixl,y,AIME)
    !benefits don't count for  aime
    Y    = ixL*Y+ abs(mod(ixy,2)==0 )*grids%benefit(ixt)
    !set spa=60 as it shouldn't matter as this subroutine only called if above SPA
    call gross2net(params,grids,Y,ixt,ixl,typeSim,AIME,spa)

    ubA1 = (A + Y - params%minCons)*(1+params%r)    ! upper bound: assets tomorrow
    policyA1temp = 0
    do ixA1 = 1, numPointsA
        if (grids%Agrid(typeSim,ixt+1, ixA1) > ubA1) then
            !if (policyA1temp ==0) then
            !    policyA1temp = ixa1
            !    EV1  = EV(ixA1,:)
            !    negV = objectivefunc(params, grids,grids%Agrid(typeSim,ixt+1, ixA1), A, Y,ixL,ixt,typeSim, AIME,EV1)
            !end if
            exit
        end if
        EV1  = EV(ixA1,:)  ! relevant section of EV matrix (in assets tomorrow) because income is irrelevant
        negVtemp = objectivefunc(params, grids,grids%Agrid(typeSim,ixt+1, ixA1), A, Y,ixL,ixt,typeSim, AIME,EV1)
        if (negVtemp > negV) then
            negV = negVtemp
            policyA1temp = ixa1
        end if
    end do
    !!-----------------------------------------------------------------------------------------------------------!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Calculate Standard errors
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!-----------------------------------------------------------------------------------------------------------!
    end subroutine

    subroutine getStandardErrors(numParams, params, grids, empricMoments, weights, sterrs)
    implicit none

#ifdef mpiBuild
    include 'mpif.h'
#endif
    type (structparamstype), intent(inout) :: params
    type (gridsType), intent(inout) :: grids

    real(kind=rk), intent(in) :: empricMoments(2,24)
    real(kind=rk), intent(in) :: weights(2,24)


    ! The NAG variables
    integer, intent(in) ::  numParams
    !real (kind=rk), intent(in) :: Params(numParams)
    real (kind=rk), intent(out) :: sterrs(numParams)

    ! For programme
    real (kind=rk) :: vcv(numParams, numParams)
    real (kind=rk) ::  paramsPert(numParams), paramsLoc(numParams) !transParams(numParams),
    integer :: ixParams, ixMom1, ixMom2
    real (kind=rk), allocatable :: momentsPertUp(:), momentsPertDown(:), momentsForChiSq(:,:)
    real (kind=rk), allocatable :: deriv(:), allDerivs(:,:), allDerivsT(:,:), W(:, :), DW(:,:), DWD(:,:), DWDinv(:,:)
    real (kind=rk), allocatable :: gW_opt(:,:), gWg_opt(:, :), gRinv(:,:), gRinvg(:, :), G_GSGinv(:,:), G_GSGinv_Gprime(:,:), G_GSGinv_GprimeW(:,:), I(:, :), P(:,:)
    real (kind=rk), allocatable :: PW(:,:), R(:,:),Rinv(:,:)
    real (kind=rk) ::  unWeightedUnSquaredMoments(nummoments), f
    character (100) :: out_sterrs, out_vcv
    !character (10) :: charixtype
    !integer :: ierrorm, startpoint, endpoint
    real (kind=rk) :: W_InvertedSquaredSDs(nummoments)
    !integer :: outErr

    ! Redundant but required by getMomentsNAG_NM
    !real (kind = rk) :: sumofsquares, ruser
    !integer :: iuser
    integer ::  ix2 !

    allocate (momentsForChiSq(numMoments, 1))
    allocate (momentsPertUp(numMoments))
    allocate (momentsPertDown(numMoments))
    allocate (deriv(numMoments))
    allocate (allDerivs(numParams, numMoments))
    allocate (allDerivsT(numMoments, numParams))
    allocate (W(numMoments, numMoments))
    allocate (DW(numParams, numMoments))
    allocate (DWD(numParams, numParams))
    allocate (DWDinv(numParams, numParams))
    allocate (gRinv(1, numMoments), gRinvg(1, 1))
    allocate (gW_opt(1, numMoments), gWg_opt(1, 1))
    allocate (G_GSGinv(numMoments, numParams))
    allocate (G_GSGinv_Gprime(numMoments, numMoments))
    allocate (G_GSGinv_GprimeW(numMoments, numMoments))
    allocate (I(numMoments, numMoments))
    allocate (P(numMoments,numMoments))

    allocate (PW(numMoments, numMoments))
    allocate (R(numMoments, numMoments))
    allocate (Rinv(numMoments, numMoments))

#ifdef mpiBuild
    call mpi_barrier(mpi_comm_world, ierror)
    if (ierror.ne.0) stop 'mpi problem173'
#endif  


    paramsLoc(1) =  params%nu
    paramsLoc(2) =  params%beta
    paramsLoc(3) =  params%gamma
    paramsloc(4) =  params%thetab
    if (numParams >= 5) paramsLoc(5) =    params%lambda

    ! Allocate


    f= gmm(params,grids,empricMoments,weights, .true., .true.,unWeightedUnSquaredMoments)
    momentsForChiSq(:, 1) = unWeightedUnSquaredMoments

    W_InvertedSquaredSDs(1:24) = weights(1,:)
    W_InvertedSquaredSDs(25:48) = weights(2,:)
    do ixParams = 1, numParams
        if (rank==0) write (*,*) "Do parameter ", ixParams
        paramsPert = paramsLoc

        paramsPert(ixParams) = paramsLoc(ixParams) + 0.1_rk * paramsLoc(ixParams)
        call copytoParams(paramsPert, params)
        !if (ixParams==1) params%nu=0.95
        !call gmm(params,grids,numerical,earnProc,empricMoments,weights,f,outErr,unWeightedUnSquaredMoments)
        f= gmm(params,grids,empricMoments,weights, .true., .true.,unWeightedUnSquaredMoments)
        momentsPertUp = unWeightedUnSquaredMoments

        paramsPert(ixParams) = paramsLoc(ixParams) - 0.1_rk * paramsLoc(ixParams)
        call copytoParams(paramsPert, params)
        !if (ixParams==1) params%nu=0.05
        f= gmm(params,grids,empricMoments,weights, .true., .true.,unWeightedUnSquaredMoments)
        !call gmm(params,grids,numerical,earnProc,empricMoments,weights,f,outErr,unWeightedUnSquaredMoments)

        momentsPertDown = unWeightedUnSquaredMoments

        deriv = (momentsPertUp - momentsPertDown) / (0.2_rk*paramsLoc(ixParams))

        allDerivs(ixParams, :) = deriv
        allDerivsT(:, ixParams) = deriv

        !if ((ixParams.eq.1)) then !(rank.eq.0).and.
        !write (*,*) rank, "Derivative", maxval(deriv)
        !write (*,*) rank, momentsPertUp(1), momentsPertDown(1)
        !write (*,*) rank, momentsPertUp(2), momentsPertDown(2)
        !end if
        if (rank == 0 ) write (*,*) deriv
    end do

#ifdef mpiBuild
    call mpi_barrier(mpi_comm_world, ierror)
    if (ierror.ne.0) stop 'mpi problem'
#endif
    if (rank == 0 )  write (*,*) 'Done loop'
    !map weights to W_InvertedSquaredSDs
    do ixMom1 = 1, numMoments
        do ixMom2 = 1, numMoments
            if (ixMom1.eq.ixMom2) then
                W(ixMom1, ixMom2) = W_InvertedSquaredSDs(ixMom1)**2
                I(ixMom1, ixMom2) = 1.0_rk
            else
                W(ixMom1, ixMom2) = 0.0_rk
                I(ixMom1, ixMom2) = 0.0_rk
            end if
        end do
    end do

    ! Get standard errors
    if (rank==0) write(*,*) '-----------------'
    if (rank==0) write(*,*) W
    if (rank==0) write(*,*) '-----------------'
    if (rank==0) write(*,*) allDerivs
    if (rank==0) write(*,*) '-----------------'
    DW = matmul(allDerivs, W)
    DWD = matmul(DW,  allDerivsT)
    !write (*,*) "Matrix ", DWD
    call getmatrixinverse(DWD, numParams, numParams + 1, DWDinv)
    vcv = DWDinv


    do ixParams = 1, numParams
        sterrs(ixParams) = vcv(ixParams,ixParams)**(0.5)
    end do

    ! Get chi2 - using Lockwood's notation

    ! This isn't working - when I try to take the inverse of R, I'm getting that the matrix isn't PD
    ! I don't know what is causing this. While G_GSGinv_Gprime has whole rows of zeros (it will as long as one derivative is zero)
    ! This doesn't mean that P has rows of zeros - P is I - G_GSGinv_GprimeW
    ! So I don't know what the story is - so I'm leaving it
    ! Notation follows Lockwood (T:\newassettest\readings) page 51


    !G_GSGinv = matmul(allDerivsT, DWDinv)
    !G_GSGinv_Gprime = matmul(G_GSGinv, allDerivs)
    !G_GSGinv_GprimeW = matmul(G_GSGinv_Gprime, W)

    !P = I - G_GSGinv_GprimeW
    !PW = matmul(P, W)
    !R = matmul(PW, P)


    !call getmatrixinverse(R, numMoments, numMoments + 1, Rinv)
    !gRinv = matmul( transpose(momentsForChiSq),Rinv)
    !gRinvg =  matmul( gRinv, momentsForChiSq)


    ! Over identification test errroneously assuming optimal weighting matrix
    gW_opt = matmul( transpose(momentsForChiSq), W)
    gWg_opt =  matmul( gW_opt, momentsForChiSq)

    if (rank.eq.0) write (*,*) 'Overidentification statistics is', gWg_opt


    !write(charixtype, '(I10)') globaltype
    out_sterrs = trim(path) // 'sterrs'  // '.txt'
    out_vcv = trim(path) // 'vcv'  // '.txt'


    if (rank.eq.0) then

        write(*,*) '-----------------'
        do ix2 = 1, numMoments
            write(*,*) allDerivs(1, ix2), allDerivs(2, ix2)
        end do

        write(*,*)
        write(*,*) '-----------------'
        do ix2 = 1, numMoments
            write(*,*) W_InvertedSquaredSDs(ix2)
        end do

        write(*,*)


        write(*,*) DWD

        write(*,*)

        write(*,*) '-----------------'
        write(*,*) 'Standard errors'
        write(*,*) '-----------------'


        write(*,*) sterrs
    end if

    if (rank.eq.0) then
        open (unit=846, file=out_sterrs, action='write', status='unknown')
        write(846,*) sterrs
        close(846)
        open (unit=846, file=out_vcv, action='write', status='unknown')
        write(846,*) vcv
        close(846)
    end if

#ifdef mpiBuild
    call mpi_barrier(mpi_comm_world, ierror)
    if (ierror.ne.0) stop 'mpi problem'
#endif


    deallocate (momentsForChiSq)
    deallocate (momentsPertUp)
    deallocate (momentsPertDown)
    deallocate (deriv)
    deallocate (allDerivs)
    deallocate (allDerivsT)
    deallocate (W)
    deallocate (DW)
    deallocate (DWD)
    deallocate (DWDinv)
    deallocate (gW_opt)
    deallocate (gWg_opt)
    deallocate (gRinv)
    deallocate (gRinvg)
    deallocate (G_GSGinv)
    deallocate (G_GSGinv_Gprime)
    deallocate (G_GSGinv_GprimeW)
    deallocate (I)
    deallocate (P)
    deallocate (PW)
    deallocate (R)
    deallocate (Rinv)

    contains
    !copy pertubation to params
    subroutine copytoParams(pert, params)
    implicit none
    type (structparamstype), intent(inout) :: params
    real(kind=rk) :: pert(numparams)

    params%nu = pert(1)
    params%beta = pert(2)
    params%gamma = pert(3)
    params%thetab = pert(4)
    if (numParams >= 5) params%lambda = pert(5)

    end subroutine
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
    subroutine solveOuterLoop(params, grids )
    implicit none
#ifdef mpiBuild
    include 'mpif.h'
#endif     
    type (structparamstype) :: params
    type (gridsType) :: grids


    !local
    logical :: delCache
    real (kind=rk) :: ypath(Tperiods, numSims) !income
    real (kind=rk) :: cpath(Tperiods, numSims)  !consumption
    integer :: lpath(Tperiods, numSims) !labour supply
    real (kind=rk) :: vpath(Tperiods, numSims)  !value
    integer :: apath(Tperiods + 1,numSims) !this is the path at the start of each period, so we include the 'start' of death
    integer :: combinedData(3,24,numSims)
    real (kind=rk) :: AIME(Tperiods + 1,numSims), error
    integer ::  ios, requiredl !,typeSim,action,
    integer ::SPA !,i, lamb
    logical:: recal = .true.
    integer :: pos(7,numSims,1), prio(7,numSims,1)
    real (kind=rk) ::post(7,numSims), prior(7,numSims,2), beta(2,1)
    logical :: mask(7,numSims)



    if (recal) call solveValueFunction( params, grids, .TRUE., .TRUE. )
    !simulate
#ifdef mpiBuild
    call mpi_barrier(mpi_comm_world, ierror)
    if (ierror.ne.0) stop 'mpi problem180'
#endif         
    if (rank == 0) then
        !delCache = .FALSE.
        !if (modelChoice>1) then
        !if (counterFact) then
        call simWithUncer(params, grids, grids%Simy, cpath, apath,  lpath, ypath, AIME, 1, .FALSE. , error, pos, prio) !5 vpath,
        !post = pos(:,:,1) + 59.0
        !prior(:,:,2) = 1.0
        !prior(:,:,1) = prio(:,:,1) + 59.0
        !mask = .true.
        !call doReg(post, prior, numSims, 7, 2, mask, .true., beta)
        !write (*,*) "beta is ", beta

        write(*,*) "error = ", error
        call writetofile(grids, ypath, cpath, apath, vpath, lpath, grids%Simy,AIME)
        call writetofileByType(grids,params,ypath, cpath, apath, vpath, lpath, grids%Simy,AIME)
        !else
        !    delCache = .FALSE.
        !    do SPA =1,3
        !        write (*,*) "Sim True SPA", 59+SPA
        !        if (SPA==3) delCache = .TRUE.
        !        call simWithUncer(params, grids, grids%Simy, cpath, apath, lpath, ypath, AIME, SPA, .FALSE. ) !vpath,
        !        combinedData(SPA,:,:) = apath(52-startAge+1:75-startAge+1,:)*(2*lpath(52-startAge+1:75-startAge+1,:)-1)
        !        call writetofileAge(grids,params,ypath, cpath, apath, vpath, lpath, grids%Simy,AIME, 59+SPA)
        !    end do
        !    inquire (iolength=requiredl) combinedData
        !    !print *, path
        !    !print *, trim(path) // 'CombinedData'
        !    open (unit=201, form="unformatted", file=trim(path) // 'CombinedData', ACCESS="STREAM", action='write', IOSTAT = ios)
        !    write (201)  combinedData
        !    close( unit=201)
        !end if
        !else
        !    call simWithUncer(params, grids, grids%Simy, cpath, apath,  lpath, ypath, AIME, TrueSPA, .FALSE. ) !vpath,
        !    !do i = 1, Tperiods
        !    !    amax(i) =  maxval( grids%Agrid(1,1, apath(i,:) ))
        !    !end do
        !    !write (*,*) maxval(amax)
        !    call writetofile(grids,ypath, cpath, apath, vpath, lpath, grids%Simy,AIME)
        !    call writetofileByType(grids,params,ypath, cpath, apath, vpath, lpath, grids%Simy,AIME)
        !end if
    end if
    do while (error > 1.0d-9)
        if (recal) call solveValueFunction( params, grids, .false., .TRUE. )
        call simWithUncer(params, grids, grids%Simy, cpath, apath,  lpath, ypath, AIME, 1, .FALSE. , error ) !5 vpath,
        write(*,*) "error = ", error
    end do

    end subroutine
    end module routines

