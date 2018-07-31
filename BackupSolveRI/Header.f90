    module Header

    integer, parameter :: rk = selected_real_kind(15)

    integer, parameter :: numPointsType = 4
    integer, parameter :: numPointsA = 15 !30 !20 !30 !45
    integer, parameter :: numPointsProd = 5
    integer, parameter :: numPointsY = 2*numPointsProd !20
    integer, parameter :: numAIME = numPointsA  !10 !5
    integer, parameter :: numPointsL = 2
    integer, parameter :: numPointsSPA = 11
    integer, parameter :: numSims =  5000!1016!  250! 10000!
    integer, parameter :: startAge =  50 !52 !20!
    integer, parameter :: endAge = 105
    integer, parameter :: Tperiods = endAge -startAge
    integer, parameter :: Tretire =60 -startAge
    integer, parameter :: TendRI = Tretire +  numPointsSPA - 1
    integer, parameter :: normBnd = 4
    integer, parameter :: dimEstimation = 6
    integer, parameter :: spouseretire = 65 -startAge
    integer, parameter :: stopwrok = 80 -startAge

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
        real (kind=rk) :: initialGuessRI(numPointsSPA*numPointsL*numPointsA+1,numPointsSPA*numPointsL*numPointsA)
    end type gridsType

    type modelObjectsType
        real (kind=rk) :: V(Tperiods+1, numPointsType,numPointsA, numAIME, numPointsSPA,numPointsY)
        real (kind=rk) :: policyA1(Tperiods, numPointsType, numPointsA, numAIME, numPointsSPA, numPointsY,numPointsA)
        !real (kind=rk) :: policyC(Tperiods, numPointsType, numPointsA, numAIME,numPointsSPA, numPointsY)
        integer :: policyL(Tperiods, numPointsType, numPointsA, numAIME,numPointsSPA, numPointsY,numPointsL)
        real (kind=rk) :: EV(Tperiods+1, numPointsType, numPointsA, numAIME, numPointsSPA,numPointsY)
    end type modelObjectsType

    !! For mpi
    integer :: rank, ierror, procsize
    integer, parameter :: mpiDim =  numPointsA * numAIME !*numPointsSPA

    !! Test controls
    logical, parameter :: fullLifeCycle = .FALSE.
    logical, parameter :: intermediateToFile = .TRUE.

    end module Header
