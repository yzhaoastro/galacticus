!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
!!    Andrew Benson <abenson@carnegiescience.edu>
!!
!! This file is part of Galacticus.
!!
!!    Galacticus is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!
!!    Galacticus is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with Galacticus.  If not, see <http://www.gnu.org/licenses/>.

  !!{
  Implements a node operator class that applies granule effect to orbiting satellite halos.
  !!}

  use :: Dark_Matter_Particles   , only : darkMatterParticleClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass
  
  !![
  <nodeOperator name="nodeOperatorSatelliteGranule">
   <description>A node operator class that applies granule effect to orbiting satellite halos.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorSatelliteGranule
     private
     double precision                                     :: massParticle
     class           (darkMatterParticleClass  ), pointer :: darkMatterParticle_   => null()
     class           (darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_ => null()
     integer                                              :: axnCoreID      , aynCoreID , aznCoreID   , &
                     &                                       axnOuterID     , aynOuterID, aznOuterID  , &
                     &                                       densityCoreID  , massCoreID, radiusCoreID, &
                     &                                       radiusSolitonID, tauID     , maxTimeID   , &
                     &                                       frequencyCoreID, stepNumID
     double precision                                     :: alphaCore      , alphaOuter
   contains
     final     :: satelliteGranuleDestructor
     procedure :: acceleration                => granuleAcceleration
     procedure :: nodeTreeInitialize          => satelliteGranuleNodeTreeInitialize
     procedure :: differentialEvolution       => satelliteGranuleDifferentialEvolution
     procedure :: differentialEvolutionScales => satelliteGranuleDifferentialEvolutionScales
  end type nodeOperatorSatelliteGranule

  interface nodeOperatorSatelliteGranule
     module procedure satelliteGranuleConstructorParameters
     module procedure satelliteGranuleConstructorInternal
  end interface nodeOperatorSatelliteGranule

contains

  function satelliteGranuleConstructorParameters(parameters) result(self)
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorSatelliteGranule)                       :: self
    type            (inputParameters             ), intent(inout)        :: parameters
    class           (darkMatterParticleClass     ),              pointer :: darkMatterParticle_
    class           (darkMatterProfileDMOClass   ),              pointer :: darkMatterProfileDMO_
    double precision                                                     :: alphaCore            , alphaOuter
  
    !![
    <inputParameter>
      <name>alphaCore</name>
      <defaultValue>5.0d0</defaultValue>
      <source>parameters</source>
      <description>Dimensionless scaling factor for the granule-induced acceleration in the solitonic core region..</description>
    </inputParameter>
    <inputParameter>
      <name>alphaOuter</name>
      <defaultValue>0.3d0</defaultValue>
      <source>parameters</source>
      <description>Dimensionless scaling factor for the granule-induced acceleration in the outer halo region.</description>
    </inputParameter>
    <objectBuilder class="darkMatterParticle"   name="darkMatterParticle_"   source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    !!]
    self =nodeOperatorSatelliteGranule(darkMatterParticle_,darkMatterProfileDMO_,alphaCore,alphaOuter)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterParticle_"   />
    <objectDestructor name="darkMatterProfileDMO_"/>
    !!]
  end function satelliteGranuleConstructorParameters

  function satelliteGranuleConstructorInternal(darkMatterParticle_,darkMatterProfileDMO_,alphaCore,alphaOuter) result(self)
    use :: Dark_Matter_Particles         , only : darkMatterParticleFuzzyDarkMatter
    use :: Numerical_Constants_Prefixes  , only : kilo
    implicit none
    type            (nodeOperatorSatelliteGranule)                     :: self
    class           (darkMatterParticleClass     ), intent(in), target :: darkMatterParticle_
    class           (darkMatterProfileDMOClass   ), intent(in), target :: darkMatterProfileDMO_
    double precision                              , intent(in)         :: alphaCore            , alphaOuter

    !![
    <constructorAssign variables="*darkMatterParticle_,*darkMatterProfileDMO_,alphaCore,alphaOuter"/>
    <addMetaProperty component="basic"             name="tau"                  id="self%tauID"            isEvolvable="yes" isCreator="yes"/>
    <addMetaProperty component="basic"             name="maxTime"              id="self%maxTimeID"        isEvolvable="no"  isCreator="yes"/>
    <addMetaProperty component="basic"             name="frequencyCore"        id="self%frequencyCoreID"  isEvolvable="no"  isCreator="yes"/>
    <addMetaProperty component="basic"             name="axnCore"              id="self%axnCoreID"        rank="1"          isCreator="yes"/>
    <addMetaProperty component="basic"             name="aynCore"              id="self%aynCoreID"        rank="1"          isCreator="yes"/>
    <addMetaProperty component="basic"             name="aznCore"              id="self%aznCoreID"        rank="1"          isCreator="yes"/>
    <addMetaProperty component="basic"             name="axnOuter"             id="self%axnOuterID"       rank="1"          isCreator="yes"/>
    <addMetaProperty component="basic"             name="aynOuter"             id="self%aynOuterID"       rank="1"          isCreator="yes"/>
    <addMetaProperty component="basic"             name="aznOuter"             id="self%aznOuterID"       rank="1"          isCreator="yes"/>
    <addMetaProperty component="basic"             name="stepNum"              type="integer"         id="self%stepNumID"   isCreator="yes"/>
    <addMetaProperty component="darkMatterProfile" name="solitonDensityCore"   id="self%densityCoreID"    isEvolvable="no"  isCreator="no" />
    <addMetaProperty component="darkMatterProfile" name="solitonMassCore"      id="self%massCoreID"       isEvolvable="no"  isCreator="no" />
    <addMetaProperty component="darkMatterProfile" name="solitonRadiusCore"    id="self%radiusCoreID"     isEvolvable="no"  isCreator="no" />
    <addMetaProperty component="darkMatterProfile" name="solitonRadiusSoliton" id="self%radiusSolitonID"  isEvolvable="no"  isCreator="no" />
    !!]

    select type (darkMatterParticle__ => self%darkMatterParticle_)
    class is (darkMatterParticleFuzzyDarkMatter)
       self%massParticle=+darkMatterParticle__%mass()*kilo
    class default
       call Error_Report('expected a `darkMatterParticleFuzzyDarkMatter` dark matter particle object'//{introspection:location})
    end select
  end function satelliteGranuleConstructorInternal


  subroutine satelliteGranuleDestructor(self)
    implicit none
    type(nodeOperatorSatelliteGranule), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterParticle_"   />
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    !!]
  end subroutine satelliteGranuleDestructor

  subroutine satelliteGranuleDifferentialEvolutionScales(self,node)
    !!{
    Set the absolute ODE solver scale for the solitonic core mass evolution,
    using a fraction of the minimum core mass as reference, following \cite{chan_diversity_2022}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentSatellite, treeNode
    implicit none
    class           (nodeOperatorSatelliteGranule), intent(inout)          :: self
    type            (treeNode                    ), intent(inout)          :: node
    class           (nodeComponentBasic          )               , pointer :: basic
    class           (nodeComponentSatellite      )               , pointer :: satellite

    if (.not. node%isSatellite()) return
    
    satellite => node%satellite()
    basic     => node%basic    ()
    call basic%floatRank0MetaPropertyScale(            &
         &                                 self%tauID, &
         &                                 1.0d0       &
         &                                )
    return
  end subroutine satelliteGranuleDifferentialEvolutionScales

  subroutine satelliteGranuleDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    use :: Galacticus_Nodes                , only : nodeComponentBasic, nodeComponentSatellite, treeNode
    use :: Numerical_Constants_Math        , only : Pi
    implicit none
    class(nodeOperatorSatelliteGranule), intent(inout), target :: self
    type(treeNode), intent(inout), target                      :: node
    class           (nodeComponentBasic               ),                pointer :: basic
    logical, intent(inout)                                     :: interrupt
    procedure(interruptTask), intent(inout), pointer           :: functionInterrupt
    integer, intent(in)                                        :: propertyType
    class (nodeComponentSatellite), pointer                     :: satellite
    double precision :: tauRate, maxTime
    double precision, dimension(3) :: acceleration
    !$GLC attributes unused :: interrupt, functionInterrupt, propertyType

    if (.not. node%isSatellite()) return
    
    satellite => node%satellite()
    basic => node%basic()
    maxTime = basic%floatRank0MetaPropertyGet(self%maxTimeID)

    tauRate = 2*Pi/maxTime

    call basic%floatRank0MetaPropertyRate(          &
         &                               self%tauID, &
         &                               tauRate           &
         &                               )

    acceleration = self%acceleration(node)
    call satellite%velocityRate(acceleration)
    return
  end subroutine satelliteGranuleDifferentialEvolution

  subroutine satelliteGranuleNodeTreeInitialize(self,node)
    !!{
    Initialize solitonic core properties for a tree node. Computes the initial core mass using the analytic core–halo relation,
    records it as the minimum core mass for tolerance scaling, and stores the value in the meta-property database of the node.    
    !!}
    use :: Coordinates                     , only : coordinateSpherical      , assignment(=)
    use :: Galacticus_Nodes                , only : nodeComponentBasic            , treeNode
    use :: Galactic_Structure_Options      , only : componentTypeDarkHalo         , massTypeDark
    use :: Statistics_Distributions         , only : distributionFunction1DNormal
    use :: Mass_Distributions              , only : massDistributionClass    , kinematicsDistributionClass
    implicit none
    class(nodeOperatorSatelliteGranule), intent(inout), target :: self
    type            (treeNode                         ), intent(inout), target  :: node
    type            (treeNode                         ),                pointer :: nodeHost
    class           (nodeComponentBasic               ),                pointer :: basic, basicHost
    class           (massDistributionClass             ), pointer       :: massDistributionHost_
    type            (coordinateSpherical                )               :: coordinates
    type (distributionFunction1DNormal)          :: normalCore, normalOuter
    integer :: stepNum
    integer :: i
    double precision    , dimension(:  ), allocatable :: axnCore, aynCore, aznCore, axnOuter, aynOuter, aznOuter
    double precision :: maxTime, densityCentral, frequencyCore

    if (.not. node%isSatellite()) return

    nodeHost           => node%parent
    do while (associated(nodeHost%parent))
       nodeHost => nodeHost%parent
    end do
    basicHost => nodeHost%basic()
    maxTime = basicHost%time ()

    nodeHost           => node%parent
    !massDistributionHost_        => self             %darkMatterProfileDMO_%get         (                      nodeHost)
    massDistributionHost_ =>  nodeHost%massDistribution()
    
    !![
    <objectDestructor name="massDistributionHost_"/>
    !!]
    
    !coordinates        = [0.0d0,0.0d0,0.0d0]
    !densityCentral = +massDistributionHost_%density(coordinates)
    densityCentral = 10.0d10
    frequencyCore = 10.94d0*(sqrt(densityCentral)) !Gyr^-1
    stepNum  = ceiling(4.0d0*frequencyCore*maxTime)
    maxTime = stepNum/(4.0d0*frequencyCore)

    stepNum = ceiling(100.0d0)
    maxTime =14.0d0
    
    basic              => node%basic()
    call basic%floatRank0MetaPropertySet(self%frequencyCoreID ,frequencyCore)
    call basic%floatRank0MetaPropertySet(self%maxTimeID ,maxTime)
    call basic%integerRank0MetaPropertySet(self%stepNumID ,stepNum)

    normalCore = distributionFunction1DNormal(0.0d0, self%alphaCore**2)
    normalOuter = distributionFunction1DNormal(0.0d0, self%alphaOuter**2)

    allocate(axnCore(0:stepNum))
    allocate(aynCore(0:stepNum))
    allocate(aznCore(0:stepNum))
    allocate(axnOuter(0:stepNum))
    allocate(aynOuter(0:stepNum))
    allocate(aznOuter(0:stepNum))
    
    do i = 0, stepNum
       axnCore(i) = normalCore%sample(randomNumberGenerator_=node%hostTree%randomNumberGenerator_)
       aynCore(i) = normalCore%sample(randomNumberGenerator_=node%hostTree%randomNumberGenerator_)
       aznCore(i) = normalCore%sample(randomNumberGenerator_=node%hostTree%randomNumberGenerator_)

       axnOuter(i) = normalOuter%sample(randomNumberGenerator_=node%hostTree%randomNumberGenerator_)
       aynOuter(i) = normalOuter%sample(randomNumberGenerator_=node%hostTree%randomNumberGenerator_)
       aznOuter(i) = normalOuter%sample(randomNumberGenerator_=node%hostTree%randomNumberGenerator_)
    end do
    
    call basic%floatRank1MetaPropertySet(self%axnCoreID ,axnCore)
    call basic%floatRank1MetaPropertySet(self%aynCoreID ,aynCore)
    call basic%floatRank1MetaPropertySet(self%aznCoreID ,aznCore)
    call basic%floatRank1MetaPropertySet(self%axnOuterID,axnOuter)
    call basic%floatRank1MetaPropertySet(self%aynOuterID,aynOuter)
    call basic%floatRank1MetaPropertySet(self%aznOuterID,aznOuter)
    return
  end subroutine satelliteGranuleNodeTreeInitialize

  function granuleAcceleration(self,node)
    !!{
    Return a test granule acceleration for satellites.
    Units: km/s/Gyr (same as other dynamical friction acceleration models).
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentSatellite, nodeComponentDarkMatterProfile, treeNode
    use :: Coordinates                     , only : coordinateSpherical      , coordinateCartesian        , assignment(=)
    use :: Error_Functions                 , only : Error_Function
    use :: Galactic_Structure_Options      , only : componentTypeDarkHalo         , massTypeDark, componentTypeAll, massTypeAll
    use :: Mass_Distributions              , only : massDistributionClass    , kinematicsDistributionClass
    use :: Numerical_Constants_Astronomical, only : gigaYear                 , gravitationalConstant_internal, megaParsec
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Numerical_Constants_Physical    , only : plancksConstant
    use :: Vectors                         , only : Vector_Magnitude
    
    implicit none
    double precision                                     , dimension(3)          :: granuleAcceleration
    class(nodeOperatorSatelliteGranule), intent(inout), target :: self
    type            (treeNode                                   ), intent(inout)         :: node
    type            (treeNode                  )               , pointer :: nodeHost
    class           (nodeComponentBasic        )               , pointer :: basicHost,  basic
    class           (nodeComponentSatellite    )               , pointer :: satellite
    class           (massDistributionClass             ), pointer       :: massDistribution_     , massDistributionTotal_
    class           (kinematicsDistributionClass       ), pointer       :: kinematics_           , kinematicsTotal_
    class           (nodeComponentDarkMatterProfile), pointer               :: darkMatterProfileHost
    type            (coordinateSpherical                )               :: coordinates, coordinatesSoliton
    !type            (coordinateCartesian                 )              :: velocity

    double precision :: velocityDispersion, velocityDispersionSoliton, sigmaJeans
    double precision :: densityCore, massCore, radiusSoliton, radiusCore
    double precision :: amplCore, amplOuter
    double precision :: frequencyCore, frequencyOuter
    double precision :: radius, velocityNorm, tau
    double precision, dimension(3) :: positionSatellite, velocitySatellite
    double precision, dimension(:), allocatable :: arrayN
    double precision, dimension(:), allocatable :: axnCore, aynCore, aznCore, axnOuter, aynOuter, aznOuter
    double precision :: axCore, ayCore, azCore, axOuter, ayOuter, azOuter
    double precision , dimension(3) :: aRawCore, aRawOuter, velocityDirection
    double precision :: Atotal, Aperp, Apar
    double precision, dimension(3) :: aPar, aPerp
    integer :: stepNum
    
    !$GLC attributes initialized :: axnCore, aynCore, aznCore, axnOuter, aynOuter, aznOuter
    !$GLC attributes unused :: self

    if (.not. node%isSatellite()) then
      granuleAcceleration = [0.0d0, 0.0d0, 0.0d0]
      return
    end if

    satellite => node     %satellite(        )
    basic => node%basic      (                 )
    positionSatellite  =  satellite     %position   (                 )
    radius             =  Vector_Magnitude          (positionSatellite)
    velocitySatellite  =  satellite     %velocity   (                 )
    velocityNorm        =  Vector_Magnitude          (velocitySatellite)
    !massBound          =  satellite     %boundMass  (                 )
    
    nodeHost           => node%parent
    basicHost          => nodeHost      %basic      (                 )
    darkMatterProfileHost => nodeHost%darkMatterProfile()

    ! massDistribution_           =>  self%darkMatterProfileDMO_%get(node  )
    ! kinematics_                 =>  massDistribution_ %kinematicsDistribution(   )
    
    massDistribution_ => node%massDistribution(componentType=componentTypeDarkHalo, massType=massTypeDark)
    massDistributionTotal_   => node%massDistribution(componentType=componentTypeAll,       massType=massTypeAll)
    kinematics_                 =>  massDistribution_ %kinematicsDistribution(   )
    kinematicsTotal_                 =>  massDistributionTotal_ %kinematicsDistribution(   )
    coordinates        = [radius,0.0d0,0.0d0]
    velocityDispersion =  +kinematics_    %velocityDispersion1D  (coordinates    ,massDistribution_    ,massDistributionTotal_    )

    !![
    <objectDestructor name="massDistribution_"    />
    <objectDestructor name="massDistributionTotal_"/>
    <objectDestructor name="kinematics_"          />
    <objectDestructor name="kinematicsTotal_"      />
    !!]

    axnCore = basic%floatRank1MetaPropertyGet(self%axnCoreID)
    aynCore = basic%floatRank1MetaPropertyGet(self%aynCoreID)
    aznCore = basic%floatRank1MetaPropertyGet(self%aznCoreID)
    axnOuter = basic%floatRank1MetaPropertyGet(self%axnOuterID)
    aynOuter = basic%floatRank1MetaPropertyGet(self%aynOuterID)
    aznOuter = basic%floatRank1MetaPropertyGet(self%aznOuterID)

    densityCore = darkMatterProfileHost%floatRank0MetaPropertyGet(self%densityCoreID)
    massCore    = darkMatterProfileHost%floatRank0MetaPropertyGet(self%massCoreID)
    radiusCore  = darkMatterProfileHost%floatRank0MetaPropertyGet(self%radiusCoreID    )
    radiusSoliton  = darkMatterProfileHost%floatRank0MetaPropertyGet(self%radiusSolitonID )
    frequencyCore = basic%floatRank0MetaPropertyGet(self%frequencyCoreID)
    tau = basic%floatRank0MetaPropertyGet(self%tauID)
    stepNum   = basic%integerRank0MetaPropertyGet(self%stepNumID)

    arrayN = arrayCos(tau, stepNum)

    coordinates        = [radius,0.0d0,0.0d0]
    coordinatesSoliton = [radiusSoliton,0.0d0,0.0d0]

    velocityDispersion =  +kinematics_    %velocityDispersion1D  (coordinates    ,massDistribution_    ,massDistribution_    )
    velocityDispersionSoliton =  +kinematics_    %velocityDispersion1D  (coordinatesSoliton    ,massDistribution_    ,massDistribution_    )
    sigmaJeans  = sqrt(2.0d0)*velocityDispersion

    if (radius<=radiusCore) then
        amplCore  = gravitationalConstant_internal*massCore/(radiusCore**2.0d0)
        !velocityDispersion =  +kinematics_    %velocityDispersion1D  (coordinates    ,massDistribution_    ,massDistribution_    )
        amplOuter = +amplitudeOuter(self,velocityDispersion,densityCore)
        !velocityDispersion =  +kinematics_    %velocityDispersion1D  (coordinatesSoliton    ,massDistribution_    ,massDistribution_    )
        frequencyOuter = self%massParticle*sqrt(2.0d0)*(velocityDispersionSoliton**2.0d0)/0.35d0/plancksConstant
    else if (radius<=radiusSoliton) then
        amplCore = 0.0d0
        !coordinates        =  [radius,0.0d0,0.0d0]
        !velocityDispersion =  +kinematics_    %velocityDispersion1D  (coordinates    ,massDistribution_    ,massDistribution_    )
        amplOuter = +amplitudeOuter(self,velocityDispersion,densityCore)
        !velocityDispersion =  +kinematics_    %velocityDispersion1D  (coordinatesSoliton    ,massDistribution_    ,massDistribution_    )
        frequencyOuter = self%massParticle*sqrt(2.0d0)*(velocityDispersionSoliton**2.0d0)/0.35d0/plancksConstant
    else
        amplCore = 0.0d0
        !velocityDispersion =  +kinematics_    %velocityDispersion1D  (coordinatesSoliton    ,massDistribution_    ,massDistribution_    )
        amplOuter = +amplitudeOuter(self,velocityDispersionSoliton,densityCore)
        !velocityDispersion =  +kinematics_    %velocityDispersion1D  (coordinates    ,massDistribution_    ,massDistribution_    )
        frequencyOuter = self%massParticle*sqrt(2.0d0)*(velocityDispersion**2.0d0)/0.35d0/plancksConstant
    end if

    axCore = dot_product(axnCore, arrayN) / sqrt(frequencyCore * tau)
    ayCore = dot_product(aynCore, arrayN) / sqrt(frequencyCore * tau)
    azCore = dot_product(aznCore, arrayN) / sqrt(frequencyCore * tau)

    axOuter = dot_product(axnOuter, arrayN) / sqrt(frequencyOuter * tau)
    ayOuter = dot_product(aynOuter, arrayN) / sqrt(frequencyOuter * tau)
    azOuter = dot_product(aznOuter, arrayN) / sqrt(frequencyOuter * tau)

    aRawCore = [axCore, ayCore, azCore]
    aRawOuter = [axOuter, ayOuter, azOuter]

    if (velocityNorm > 0.0d0) then
        velocityDirection = velocitySatellite/velocityNorm
    else
        velocityDirection = [0.0d0, 0.0d0, 0.0d0]
    end if

    Atotal = sqrt( 1.0d0 + (velocityNorm / sigmaJeans)**2 )
    Aperp  = sqrt( 1.0d0 + 0.5d0 * (velocityNorm / sigmaJeans)**2 )
    Apar   = 1.0d0 / (3.0d0 / Atotal - 2.0d0 / Aperp)

    ! The velocity-dependent factor is applied only to a_raw_outer.
    ! The acceleration from the core is set to 0 at d > r_c.
    aPar  = dot_product(aRawOuter, velocityDirection)*velocityDirection/sqrt(Apar)      ! parallel
    aPerp = (aRawOuter - aPar)/sqrt(Aperp)                    ! perpendicular

    granuleAcceleration = amplCore*aRawCore + amplOuter*(aPar+aPerp)
    return
  end function granuleAcceleration

  function arrayCos(tau, m) result(arr)
    !!{
    The first element is 1/sqrt(2). The next n elements are cos(nt/2) for n = 1, 2, ..., N.
    !!}
    implicit none
    integer, intent(in) :: m
    double precision, intent(in) :: tau
    double precision, allocatable, dimension(:) :: arr
    integer :: n

    allocate(arr(0:m))
    arr(0) = 1.0d0 / sqrt(2.0d0)
    do n = 1, m
        arr(n) = cos(dble(n) * tau / 2.0d0)
    end do
    return
  end function arrayCos

  function amplitudeOuter(self,velocityDispersion,densityCore)
    use :: Numerical_Constants_Astronomical, only : gigaYear                 , gravitationalConstant_internal, megaParsec
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Numerical_Constants_Physical    , only : plancksConstant
    use :: Vectors                         , only : Vector_Magnitude
    implicit none
    double precision :: amplitudeOuter
    class(nodeOperatorSatelliteGranule), intent(inout), target :: self
    double precision, intent(in) :: velocityDispersion,densityCore
    double precision :: sigmaJeans, LambdaDB, deff, meff
    
    sigmaJeans  = sqrt(2.0d0)*velocityDispersion
    LambdaDB    = plancksConstant*sigmaJeans/self%massParticle
    deff        = 0.35d0*LambdaDB
    meff        = densityCore*(0.282d0*LambdaDB)**3.0d0

    amplitudeOuter = gravitationalConstant_internal*meff/(deff**2.0d0)
  return
  end function amplitudeOuter
