    module Header
    use Policy
    implicit none
    !#include "../globalMacros.txt"
#include "globalMacros.txt"

    integer, parameter :: sp = 4
    integer, parameter :: rk = selected_real_kind(15)
    integer, parameter :: hp = selected_real_kind(20)
    integer, parameter :: li = 8

    integer, parameter :: numPointsType = TYPES_SIZE !1!
#ifdef _WIN64
    integer, parameter :: numPointsA =   30 !90 !  120!  !140 !110 !120 !30 !30 !120!10!30!100!50!
#else
    integer, parameter :: numPointsA = 120!120!10!30!100!50!
#endif 

    integer, parameter :: numPointsProd = PROD_SIZE !5!10 !
    integer, parameter :: numPointsY = 2*numPointsProd !20
    integer, parameter :: numAIME = 8!numPointsA  !10 !5
    integer, parameter :: numPointsL = 2
    integer, parameter :: numPointsSPA = 8 !11!8 !9!
    integer, parameter :: numPointsPost = 5
    real (kind=rk) :: probBlock = 1.0_rk/real(numPointsPost,rk)
    integer, parameter :: numSims =  1000!50000 !1016!  5000!250! 10000!
    integer, parameter :: startAge =  52 !20!
    integer, parameter :: weeksYear = 52
    integer, parameter :: weeksWorking = 48
    integer, parameter :: endAge = 101
    integer, parameter :: Tperiods = endAge -startAge
    integer, parameter :: Tretire =60 -startAge +1
    integer, protected :: TrueSPA = 1
    integer, parameter :: TendRI = Tretire +  numPointsSPA - 1!CHANGE THIS ONE FOR RE
    integer, parameter :: normBnd = 3 !4 is to big as the lowest point will be at the bound
    integer, protected  :: dimEstimation = 4
    !integer, parameter :: spouseretire = 65 -startAge+1
    integer, parameter :: stopwrok = 80 -startAge+1
    integer, parameter :: obsInitDist = 838 !874
    integer, parameter :: numMoments = 48
    integer :: modelChoice
    integer, protected :: EndPeriodRI
    logical, protected :: uncerRE
    logical, parameter :: counterFact = .TRUE.!.FALSE. !
    integer, protected :: sizeDimSPA
    integer, protected :: model
    integer, parameter :: numStates = 5
    integer, parameter :: numChoices = 2
    integer, protected :: numRIflags
    integer, protected :: dimVal(Tperiods,numStates)
    integer, protected :: dimPol(Tperiods,numStates)
    integer, parameter :: flgRcvd = 1
    integer, parameter :: flgNotRcvd = 2
    !integer, parameter :: knn = 2 !5 !

#ifdef _WIN64          
    character(len=250), parameter :: pathMoments = 'C:\Users\Uctphen\Dropbox\SourceCode\upgradeProject\moments\'
    !character(len=250), parameter :: pathErrors = 'C:\Users\Uctphen\Dropbox\SourceCode\upgradeProject\VSProj - Copy\'
    character(len=250), parameter :: pathInputs = 'C:\Users\Uctphen\Dropbox\SourceCode\upgradeProject\input\'
    character(len=250), parameter :: path_bobyqa = 'C:\Users\Uctphen\bobyqa\'
#else
    character(len=250), parameter :: pathMoments = '/home/jamie/Dropbox/SourceCode/upgradeProject/moments/'
    !character(len=250), parameter :: pathErrors = 'C:\Users\Uctphen\Dropbox\SourceCode\upgradeProject\VSProj - Copy\'
    character(len=250), parameter :: pathInputs = '/home/jamie/Dropbox/SourceCode/upgradeProject/input/'
    character(len=250), parameter :: path_bobyqa = 'C:\Users\Uctphen\bobyqa\'

    !character(len=250), parameter :: pathMoments = '~/FortrCodeRI/moments/'
    !!character(len=250), parameter :: pathErrors = '~/FortrCodeRI/Errors/'
    !character(len=250), parameter :: pathInputs = '~/FortrCodeRI/input/'
    !character(len=250), parameter :: path_bobyqa = '/scratch/scratch/uctphen/store/'
#endif 
    character(len=250), protected :: pathDataStore
    character(len=250), protected :: path
    !Holds the structural parameters that will eventually be estimated (at least in some cases)
    type structparamstype
        !Personal
        real (kind=rk) :: gamma
        real (kind=rk) :: r
        real (kind=rk) :: beta
        real (kind=rk) :: mu
        real (kind=rk) :: rho(numPointsType)
        real (kind=rk) :: sigma(numPointsType)
        real (kind=rk) :: sigma_int(numPointsType)
        real (kind=rk) :: fracType(numPointsType)
        real (kind=rk) :: ERA(numPointsType)
        real (kind=rk) :: nu
        !Instutional
        real (kind=rk) :: delta(numPointsType,3)
        real (kind=rk) :: pension(numPointsType,2)
        real (kind=rk) :: hrsWrk
        real (kind=rk) :: spouseInc(numPointsType,Tperiods)
        real (kind=rk) :: minCons
        real (kind=rk) :: db(numPointsType,2)
        real (kind=rk) :: startA
        real (kind=rk) :: thetab
        real (kind=rk) :: k
        real (kind=rk) :: tol
        real (kind=rk) :: p
        real (kind=rk) :: lambda
        real (kind=rk) :: percentCons
        real (kind=rk) :: utlityShifter
        real (kind=rk) :: spaTransMat(numPointsSPA,numPointsSPA)
        !bound
        real (kind=rk) :: BnDnu(2)
        real (kind=rk) :: BnDbeta(2)
        real (kind=rk) :: BnDgamma(2)
        real (kind=rk) :: BnDdb1(2)
        real (kind=rk) :: BnDdb2(2)
        real (kind=rk) :: BnDthetab(2)
    end type structparamstype

    type beliefStoreType
        integer :: sizePost
        real (kind=rk), allocatable :: storePost(:,:)
        real (kind=rk), allocatable :: distMean(:)
        integer(kind=li) :: storePostBin(184756) !storePostBin(:)
    end type beliefStoreType

    type gridsType
        real (kind=rk) :: Agrid(numPointsType,Tperiods+1,numPointsA)
        real (kind=rk) :: YgridExtrap(numPointsType,Tperiods,numPointsProd)
        real (kind=rk) :: Ygrid(numPointsType,Tperiods,numPointsY)
        real (kind=rk) :: incTransitionMrx(numPointsY,numPointsY)
        real (kind=rk) :: AIMEgrid(numPointsType,Tperiods+1,numAIME)
        real (kind=rk) :: benefit(Tperiods)
        real (kind=rk) :: fc(Tperiods)
        real (kind=rk) :: maxInc(numPointsType,Tperiods)
        real (kind=rk) :: initialAssets(obsInitDist,2)
        real (kind=rk) :: mortal(Tperiods+1)
        real (kind=rk) :: Simy(Tperiods, numSims)
        real (kind=rk) :: unemploy(numPointsType,numPointsProd)
        real (kind=rk) :: reemploy(numPointsType,numPointsProd)
        real (kind=rk), allocatable :: posteriorSPA(:,:,:, :,  :,:)
        !real (kind=rk) :: posteriorSPA(TendRI,1,1,1,1,numPointsSPA) Tretire
        integer :: supportSPA(TendRI)
        !real (kind=rk) :: initialGuessRI(numPointsSPA+1,numPointsSPA)
        real (kind=rk) :: initialGuessRI(numPointsSPA*numPointsL*numPointsA+1,numPointsSPA*numPointsL*numPointsA)
        real (kind=rk) :: logitShock(numSims,Tperiods)
        type (beliefStoreType) :: beliefStore(1:TendRI-1)
    end type gridsType

    !type sparseCOOType
    !    integer :: col
    !    integer :: row
    !    real(kind=rk) :: val
    !end type sparseCOOType
    !
    !type policyType
    !    type(sparseCOOType), allocatable :: COO(:)
    !end type policyType

    type modelObjectsType
        real (kind=sp), allocatable :: V(:, :, :,:,:)
        !real (kind=sp), allocatable :: policy(:, :, :, :,:,:,:)
        type(policyType), allocatable :: policy(:, :, :,:,:)
        !real (kind=rk), allocatable :: u(:, :, :, :,:,:)
        !real (kind=rk), allocatable :: q(:, :, :, :,:)
        real (kind=sp), allocatable :: EV(:, :, :,:,:)
    end type modelObjectsType

    !allocate(modelObjects%EV(numPointsA, numAIME, numPointsSPA,numPointsY))
    !allocate(modelObjects%policy(numPointsA, numAIME, numPointsSPA, numPointsY,numPointsL,numPointsA))
    !allocate(modelObjects%V(numPointsA, numAIME, numPointsSPA,numPointsY))

    !! For mpi
    integer :: rank, ierror, procsize
    integer, parameter :: mpiDim =  numPointsA * numAIME

    real (kind=rk), protected :: timeHack
    REAL (kind=rk), protected :: rate
    INTEGER (KIND=8) , protected  :: c1,c2,cr,cm

    contains
    ! ---------------------------------------------------------------------------------------------------------!
    !!Set global proctied variables
    !----------------------------------------------------------------------------------------------------------!
    subroutine setModel
    implicit none

#ifdef mpiBuild
    include 'mpif.h'
    CHARACTER(len=32) :: arg
#endif 
    integer :: globalCounter, globalSize
    ! First initialize the system_clock
    CALL system_clock(count_rate=cr)
    CALL system_clock(count_max=cm)
    rate = REAL(cr)
    CALL SYSTEM_CLOCK(c1)

#ifdef _WIN64
    if (rank == 0) then
        write (*,*) "Press 1 for RE non uncert SPA, 2 for RE+uncert SPA, 3 for RI"

        read (*,*) modelChoice
    end if
#else
    CALL get_command_argument(1, arg)
    read(arg(1:1),'(i)') modelChoice
#endif 

#ifdef mpiBuild
    call MPI_Bcast( modelChoice, 1, MPI_INTEGER ,0, mpi_comm_world, ierror)
    if (ierror.ne.0) stop 'mpi problem171'
#endif     


    model = modelChoice
    select case(modelChoice)
    case(1)
        dimVal(:,:) = reshape(((/(numPointsA,globalCounter=1,Tperiods), (numAIME,globalCounter=1,Tperiods), &
            (1,globalCounter=1,Tperiods),(1,globalCounter=1,Tperiods),(numPointsY,globalCounter=1,Tperiods)/)),(/Tperiods,numStates/))
        dimPol= dimVal 
        !dimPol(:,:) = reshape(((/(numPointsA,globalCounter=1,Tperiods), (numAIME,globalCounter=1,Tperiods), &
        !    (1,globalCounter=1,Tperiods),(1,globalCounter=1,Tperiods),(numPointsY,globalCounter=1,Tperiods), &
        !    (numPointsL,globalCounter=1,Tperiods),(numPointsA,globalCounter=1,Tperiods)/)),(/Tperiods,numStates+ numChoices/))
        !(numpointsA,globalCounter=1,Tperiods)
        sizeDimSPA = 1
        EndPeriodRI = 1
        numRIflags = 1
#ifdef _WIN64
        pathDataStore = "C:\Users\Uctphen\DataStore\PolicyFuncsBaseline\"
        path = "C:\Users\Uctphen\Dropbox\SourceCode\upgradeProject\VSProj - Copy\outBaseline\"
#else
        pathDataStore = "/scratch/scratch/uctphen/PolicyFuncsBaseline/"
        path = "/scratch/scratch/uctphen/outBaseline/"
#endif         
    case(2)
        sizeDimSPA = numPointsSPA
        EndPeriodRI = TendRI
        uncerRE = .TRUE.
        numRIflags = 1
#ifdef _WIN64
        pathDataStore = "C:\Users\Uctphen\DataStore\PolicyFuncsRE\"
        path = "C:\Users\Uctphen\Dropbox\SourceCode\upgradeProject\VSProj - Copy\outRE\"
#else
        pathDataStore = "/scratch/scratch/uctphen/policyFuncsRE/"
        path = "/scratch/scratch/uctphen/outRE/"
#endif 
        dimVal(:,:) = reshape(((/(numPointsA,globalCounter=1,Tperiods), (numAIME,globalCounter=1,Tperiods), &
            (1,globalCounter=1,Tperiods),(numPointsSPA,globalCounter=1,TendRI-1),(1,globalCounter=EndPeriodRI,Tperiods), &
            (numPointsY,globalCounter=1,Tperiods)/)),(/Tperiods,numStates/))
        dimpol = dimval
        !dimPol(:,:) = reshape(((/(numPointsA,globalCounter=1,Tperiods), (numAIME,globalCounter=1,Tperiods), &
        !    (1,globalCounter=1,Tperiods),(numPointsSPA,globalCounter=1,TendRI-1),(1,globalCounter=EndPeriodRI,Tperiods),&
        !    (numPointsY,globalCounter=1,Tperiods), &
        !    (numPointsL,globalCounter=1,Tperiods),(numPointsA,globalCounter=1,Tperiods)/)),(/Tperiods,numStates+ numChoices/))
    case default
        numRIflags = 2
        sizeDimSPA = 999
        EndPeriodRI = TendRI
        uncerRE = .FALSE.
        globalSize =  choose(numPointsSPA - 1 + numPointsPost,numPointsPost)
        write(*,*) numPointsProd, numPointsY
        do globalCounter = 1, Tperiods
            if (globalCounter >= Tretire) then
                globalSize = choose(numPointsSPA - 1 + numPointsPost   - (globalCounter-Tretire+1),numPointsPost)
            end if
            if (globalCounter >= EndPeriodRI -1) globalSize = 1
            if (globalCounter == EndPeriodRI) numRIflags = 1
            if (globalCounter < EndPeriodRI) then
                dimVal(globalCounter,:) =  (/numPointsA,numAIME,globalSize,numpointsspa,numPointsY/)
               ! dimPol(globalCounter,:) =  (/numPointsA,numAIME,globalSize,numpointsspa,numPointsY,numPointsL,numPointsA/)
            else
                dimVal(globalCounter,:) =  (/numPointsA,numAIME,globalSize,1,numPointsY/)
                !dimPol(globalCounter,:) =  (/numPointsA,numAIME,globalSize,1,numPointsY,numPointsL,numPointsA/)
            end if
            dimpol=dimval
        end do
        numRIflags = 2
#ifdef _WIN64
        pathDataStore = "C:\Users\Uctphen\DataStore\PolicyFuncs\"
        path = "C:\Users\Uctphen\Dropbox\SourceCode\upgradeProject\VSProj - Copy\out\"
#else
        !pathDataStore = "/scratch/scratch/uctphen/PolicyFuncs/"
        !path = "/scratch/scratch/uctphen/out/"
        path = "/home/jamie/Dropbox/SourceCode/upgradeProject/VSProj - Copy/out/"
        pathDataStore = "/home/jamie/tempFolder/"
#endif 
        dimEstimation = dimEstimation + 1
    end select
    contains
    function choose (n, k) result (res)

    implicit none
    integer, intent (in) :: n
    integer, intent (in) :: k
    real(kind=rk):: res

    res = factorial (n) / (factorial (k) * factorial (n - k))

    end function choose
    function factorial (num) result (res)

    implicit none
    integer, intent (in) :: num
    real(KIND=rk) :: res
    integer :: i

    !Why doesn't this work?
    !res = product ((/(i, i = 1, num)/))
    res=1
    if (num > 0) then
        do i =1,num
            res=res*i
        end do
    end if
    end function factorial
    end subroutine
    subroutine setSPA(newSPA)
    implicit none
    integer, intent(in) :: newSPA
    TrueSPA=newSPA
    end subroutine
    end module Header
