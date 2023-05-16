
abstract type StellarType end
abstract type Star             <: StellarType end
abstract type CompactObject    <: StellarType end
abstract type SubStellarObject <: StellarType end
abstract type Other            <: StellarType end

struct DeeplyOrFullyConvectiveLowMassMainSequence <: StellarType
    key::Int
    DeeplyOrFullyConvectiveLowMassMainSequence() = new(0)
end

struct MainSequence <: StellarType
    key::Int
    MainSequence() = new(1)
end

struct HertzsprungGap <: StellarType
    key::Int
    HertzsprungGap() = new(2)
end

struct FirstGiantBranch <: StellarType
    key::Int
    FirstGiantBranch() = new(3)
end

struct CoreHeliumBurning <: StellarType
    key::Int
    CoreHeliumBurning() = new(4)
end

struct FirstAsymptoticGiantBranch <: StellarType
    key::Int
    FirstAsymptoticGiantBranch() = new(5)
end

struct SecondAsymptoticGiantBranch <: StellarType
    key::Int
    SecondAsymptoticGiantBranch() = new(6)
end

struct MainSequenceNakedHelium <: StellarType
    key::Int
    MainSequenceNakedHelium() = new(7)
end

struct HertzsprungGapNakedHelium <: StellarType
    key::Int
    HertzsprungGapNakedHelium() = new(8)
end

struct GiantBranchNakedHelium <: StellarType
    key::Int
    GiantBranchNakedHelium() = new(9)
end

struct HeliumWhiteDwarf <: CompactObject
    key::Int
    HeliumWhiteDwarf() = new(10)
end

struct CarbonOxygenWhiteDwarf <: CompactObject
    key::Int
    CarbonOxygenWhiteDwarf() = new(11)
end

struct OxygenNeonWhiteDwarf <: CompactObject
    key::Int
    OxygenNeonWhiteDwarf() = new(12)
end

struct NeutronStar <: CompactObject
    key::Int
    NeutronStar() = new(13)
end

struct BlackHole <: CompactObject
    key::Int
    BlackHole() = new(14)<: Star
end

struct MasslessSupernova <: Other
    key::Int
    MasslessSupernova() = new(15)
end

struct UnknownStellarType <: Other
    key::Int
    UnknownStellarType() = new(16)
end

struct PreMainSequence <: Star
    key::Int
    PreMainSequence() = new(17)
end

struct Planet <: SubStellarObject
    key::Int
    Planet() = new(18)
end

struct BrownDwarf <: SubStellarObject
    key::Int
    BrownDwarf() = new(19)
end


















