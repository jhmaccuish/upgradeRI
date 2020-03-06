    module Header

    !#include "../globalMacros.txt"
#include "globalMacros.txt"

    integer, parameter :: rk = selected_real_kind(15)
    integer, parameter :: hp = selected_real_kind(20)

    integer, parameter :: numPointsType = TYPES_SIZE !1!
#ifdef _WIN64      
    integer, parameter :: numPointsA =  30 !90 !  !140 !110 !120 !30 !30 !120!10!30!100!50!
#else
    integer, parameter :: numPointsA = 120!120!10!30!100!50!
#endif 

    integer, parameter :: numPointsProd = PROD_SIZE !5!10 !
    integer, parameter :: numPointsY = 2*numPointsProd !20
    integer, parameter :: numAIME = 8!numPointsA  !10 !5
    integer, parameter :: numPointsL = 2
    integer, parameter :: numPointsSPA = 11
    integer, parameter :: numSims =  50000 !1016!  5000!250! 10000!
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
#ifdef _WIN64          
    character(len=250), parameter :: pathMoments = 'C:\Users\Uctphen\Dropbox\SourceCode\upgradeProject\moments\'
    !character(len=250), parameter :: pathErrors = 'C:\Users\Uctphen\Dropbox\SourceCode\upgradeProject\VSProj - Copy\'
    character(len=250), parameter :: pathInputs = 'C:\Users\Uctphen\Dropbox\SourceCode\upgradeProject\input\'   
    character(len=250), parameter :: path_bobyqa = 'C:\Users\Uctphen\bobyqa\'
#else
    character(len=250), parameter :: pathMoments = '~/FortrCodeRI/moments/'
    !character(len=250), parameter :: pathErrors = '~/FortrCodeRI/Errors/'
    character(len=250), parameter :: pathInputs = '~/FortrCodeRI/input/'  
    character(len=250), parameter :: path_bobyqa = '/scratch/scratch/uctphen/store/'
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
        real (kind=rk) :: posteriorSPA(TendRI,numPointsSPA)
        integer :: supportSPA(TendRI)
        !real (kind=rk) :: initialGuessRI(numPointsSPA+1,numPointsSPA)
        real (kind=rk) :: initialGuessRI(numPointsSPA*numPointsL*numPointsA+1,numPointsSPA*numPointsL*numPointsA)
        real (kind=rk) :: logitShock(numSims,Tperiods)
    end type gridsType


    type modelObjectsType
        real (kind=rk), allocatable :: V(:, :, :,:)
        real (kind=rk), allocatable :: policy(:, :, :, :,:,:)
        real (kind=rk), allocatable :: u(:, :, :, :,:,:)
        real (kind=rk), allocatable :: q(:, :, :, :,:)
        real (kind=rk), allocatable :: EV(:, :, :,:)
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
    CHARACTER(len=32) :: arg
#ifdef mpiBuild
    include 'mpif.h'
#endif 
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



    select case(modelChoice)
    case(1)
        EndPeriodRI = 1
#ifdef _WIN64
        pathDataStore = "C:\Users\Uctphen\DataStore\PolicyFuncsBaseline\"
        path = "C:\Users\Uctphen\Dropbox\SourceCode\upgradeProject\VSProj - Copy\outBaseline\"
#else
        pathDataStore = "/scratch/scratch/uctphen/PolicyFuncsBaseline/"
        path = "/scratch/scratch/uctphen/outBaseline/"
#endif         
    case(2)
        EndPeriodRI = TendRI
        uncerRE = .TRUE.
#ifdef _WIN64
        pathDataStore = "C:\Users\Uctphen\DataStore\PolicyFuncsRE\"
        path = "C:\Users\Uctphen\Dropbox\SourceCode\upgradeProject\VSProj - Copy\outRE\"
#else
        pathDataStore = "/scratch/scratch/uctphen/policyFuncsRE/"
        path = "/scratch/scratch/uctphen/outRE/"
#endif 
        case default
        EndPeriodRI = TendRI
        uncerRE = .FALSE.
#ifdef _WIN64
        pathDataStore = "C:\Users\Uctphen\DataStore\PolicyFuncs\"
        path = "C:\Users\Uctphen\Dropbox\SourceCode\upgradeProject\VSProj - Copy\out\"
#else
        pathDataStore = "/scratch/scratch/uctphen/PolicyFuncs/"
        path = "/scratch/scratch/uctphen/out/"
#endif 
        dimEstimation = dimEstimation + 1
    end select

    end subroutine
    subroutine setSPA(newSPA)
    implicit none
    integer, intent(in) :: newSPA
    TrueSPA=newSPA
    end subroutine
    end module Header