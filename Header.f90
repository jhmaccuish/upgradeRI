    module Header

    !#include "../globalMacros.txt"
#include "globalMacros.txt"

    integer, parameter :: rk = selected_real_kind(15)
    integer, parameter :: hp = selected_real_kind(20)

    integer, parameter :: numPointsType = TYPES_SIZE !1!
    integer, parameter :: numPointsA = 30!10!30!100!50!
    
    integer, parameter :: numPointsProd = PROD_SIZE !5!10 !
    integer, parameter :: numPointsY = 2*numPointsProd !20
    integer, parameter :: numAIME = 8!numPointsA  !10 !5
    integer, parameter :: numPointsL = 2
    integer, parameter :: numPointsSPA = 11
    integer, parameter :: numSims =  1016!  5000!250! 10000!
    integer, parameter :: startAge =  52 !20!
    integer, parameter :: endAge = 105
    integer, parameter :: Tperiods = endAge -startAge
    integer, parameter :: Tretire =60 -startAge +1
    integer, parameter :: TrueSPA = 1
    integer, parameter :: TendRI = Tretire +  numPointsSPA - 1!CHANGE THIS ONE FOR RE
    integer, parameter :: normBnd = 4
    integer, protected  :: dimEstimation = 6
    integer, parameter :: spouseretire = 65 -startAge+1
    integer, parameter :: stopwrok = 80 -startAge+1
    integer :: modelChoice
    integer, protected :: EndPeriodRI
    logical, protected :: uncerRE
    logical, parameter :: counterFact = .TRUE.!.FALSE. !
    character(len=250), parameter :: pathMoments = 'C:\Users\Uctphen\Dropbox\SourceCode\upgradeProject\moments\'
    character(len=250), parameter :: pathErrors = 'C:\Users\Uctphen\Dropbox\SourceCode\upgradeProject\VSProj - Copy\'
    character(len=250), protected :: pathDataStore
    character(len=250), protected :: path
    !Holds the structural parameters that will eventually be estimated (at least in some cases)
    type structparamstype
        !Personal
        real (kind=rk) :: gamma
        real (kind=rk) :: r
        real (kind=rk) :: beta
        real (kind=rk) :: sigma
        real (kind=rk) :: mu
        real (kind=rk) :: rho
        real (kind=rk) :: nu
        !Instutional
        real (kind=rk) :: delta(numPointsType,3)
        real (kind=rk) :: pension
        real (kind=rk) :: hrsWrk
        real (kind=rk) :: spouseInc(numPointsType)
        real (kind=rk) :: minCons
        real (kind=rk) :: db(2)
        real (kind=rk) :: startA
        real (kind=rk) :: thetab
        real (kind=rk) :: k
        real (kind=rk) :: tol
        real (kind=rk) :: p
        real (kind=rk) :: lambda
        real (kind=rk) :: percentCons
        real (kind=rk) :: utlityShifter
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
        real (kind=rk) :: initialAssets(1016)
        real (kind=rk) :: mortal(Tperiods+1)
        real (kind=rk) :: Simy(Tperiods, numSims)
        real (kind=rk) :: unemploy(numPointsType,numPointsProd)
        real (kind=rk) :: reemploy(numPointsType,numPointsProd)
        real (kind=rk) :: posteriorSPA(TendRI,numPointsSPA)
        integer :: supportSPA(TendRI)
        !real (kind=rk) :: initialGuessRI(numPointsSPA+1,numPointsSPA)
        real (kind=rk) :: initialGuessRI(numPointsSPA*numPointsL*numPointsA+1,numPointsSPA*numPointsL*numPointsA)
    end type gridsType

    !Made allocatable to allow for very large arrays
    !type modelObjectsType
    !    real (kind=rk):: V(numPointsA, numAIME, numPointsSPA,numPointsY)
    !    real (kind=rk) :: policy(numPointsA, numAIME, numPointsSPA, numPointsY,numPointsL,numPointsA)
    !    real (kind=rk) :: EV(numPointsA, numAIME, numPointsSPA,numPointsY)
    !end type modelObjectsType

    !type locmodelObjectsType
    !    real (kind=rk):: V(numPointsA, numAIME, numPointsSPA,numPointsY/2)
    !    real (kind=rk) :: policy(numPointsA, numAIME, numPointsSPA, numPointsY/2,numPointsL,numPointsA)
    !    real (kind=rk) :: EV(numPointsA, numAIME, numPointsSPA,numPointsY/2)
    !end type locmodelObjectsType

    type modelObjectsType
        real (kind=rk), allocatable :: V(:, :, :,:)
        real (kind=rk), allocatable :: policy(:, :, :, :,:,:)
        real (kind=rk), allocatable :: EV(:, :, :,:)
    end type modelObjectsType

    !allocate(modelObjects%EV(numPointsA, numAIME, numPointsSPA,numPointsY))
    !allocate(modelObjects%policy(numPointsA, numAIME, numPointsSPA, numPointsY,numPointsL,numPointsA))
    !allocate(modelObjects%V(numPointsA, numAIME, numPointsSPA,numPointsY))

    !! For mpi
    integer :: rank, ierror, procsize
    integer, parameter :: mpiDim =  numPointsA * numAIME

    !! Test controls
    logical, parameter :: fullLifeCycle = .FALSE.
    !logical, parameter :: intermediateToFile = .FALSE.

    contains
    subroutine setModel
    implicit none

#ifdef mpi
    include 'mpif.h'
#endif 

    if (rank == 0) then
        write (*,*) "Press 1 for RE non uncert SPA, 2 for RE+uncert SPA, 3 for RI"

        read (*,*) modelChoice
    end if
#ifdef mpi
    call MPI_Bcast( modelChoice, 1, MPI_INTEGER ,0, mpi_comm_world, ierror)
    if (ierror.ne.0) stop 'mpi problem171'
#endif     
    select case(modelChoice)
    case(1)
        EndPeriodRI = 1
        pathDataStore = "C:\Users\Uctphen\DataStore\PolicyFuncsBaseline\"
        path = "C:\Users\Uctphen\Dropbox\SourceCode\upgradeProject\VSProj - Copy\outBaseline\"
    case(2)
        EndPeriodRI = TendRI
        uncerRE = .TRUE.
        pathDataStore = "C:\Users\Uctphen\DataStore\PolicyFuncsRE\"
        path = "C:\Users\Uctphen\Dropbox\SourceCode\upgradeProject\VSProj - Copy\outRE\"
        case default
        EndPeriodRI = TendRI
        uncerRE = .FALSE.
        pathDataStore = "C:\Users\Uctphen\DataStore\PolicyFuncs\"
        path = "C:\Users\Uctphen\Dropbox\SourceCode\upgradeProject\VSProj - Copy\out\"
        dimEstimation = dimEstimation + 1
    end select

    end subroutine

    end module Header
