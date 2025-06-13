
abstract type AbstractStellarType end

"""
Deeply or fully convective low main sequence star with index 0. 
"""
struct DeeplyOrFullyConvectiveLowMassMainSequence <: AbstractStellarType
    number::Int
    DeeplyOrFullyConvectiveLowMassMainSequence() = new(0)
end

"""
Main sequence star with index 1. 
"""
struct MainSequence <: AbstractStellarType
    number::Int
    MainSequence() = new(1)
end

"""
Hertzsprung gap star with index 2. 
"""
struct HertzsprungGap <: AbstractStellarType# Giant
    number::Int
    HertzsprungGap() = new(2)
end

"""
First giant branch star with index 3. 
"""
struct FirstGiantBranch <: AbstractStellarType# Giant
    number::Int
    FirstGiantBranch() = new(3)
end

"""
Core helium burning star with index 4. 
"""
struct CoreHeliumBurning <: AbstractStellarType# Giant
    number::Int
    CoreHeliumBurning() = new(4)
end

"""
First asymptotic gant branch star with index 5.
"""
struct FirstAsymptoticGiantBranch <: AbstractStellarType# Giant
    number::Int
    FirstAsymptoticGiantBranch() = new(5)
end

"""
Second asymptotic giant branch star with index 6.
"""
struct SecondAsymptoticGiantBranch <: AbstractStellarType# Giant
    number::Int
    SecondAsymptoticGiantBranch() = new(6)
end

"""
Main sequence naked helium star with index 7.
"""
struct MainSequenceNakedHelium <: AbstractStellarType
    number::Int
    MainSequenceNakedHelium() = new(7)
end

"""
Hertzsprung gap naked helium star with index 8.
"""
struct HertzsprungGapNakedHelium <: AbstractStellarType# Giant
    number::Int
    HertzsprungGapNakedHelium() = new(8)
end

"""
Giant branch naked helium star with index 9.
"""
struct GiantBranchNakedHelium <: AbstractStellarType# Giant
    number::Int
    GiantBranchNakedHelium() = new(9)
end

"""
Helium white dwarf with index 10.
"""
struct HeliumWhiteDwarf <: AbstractStellarType
    number::Int
    HeliumWhiteDwarf() = new(10)
end

"""
Carbon/oxygen white dwarf with index 11.
"""
struct CarbonOxygenWhiteDwarf <: AbstractStellarType
    number::Int
    CarbonOxygenWhiteDwarf() = new(11)
end

"""
Oxygen/neon white dwarf with index 12.
"""
struct OxygenNeonWhiteDwarf <: AbstractStellarType
    number::Int
    OxygenNeonWhiteDwarf() = new(12)
end

"""
Neutron star with index 13.
"""
struct NeutronStar <: AbstractStellarType
    number::Int
    NeutronStar() = new(13)
end

"""
Black hole with index 14.
"""
struct BlackHole <: AbstractStellarType
    number::Int
    BlackHole() = new(14)
end

"""
Massless supernova with index 15.
"""
struct MasslessSupernova <: AbstractStellarType
    number::Int
    MasslessSupernova() = new(15)
end

"""
Unknown stellar type with index 16.
"""
struct UnknownStellarType <: AbstractStellarType
    number::Int
    UnknownStellarType() = new(16)
end

"""
Pre main sequence star with index 17.
"""
struct PreMainSequence <: AbstractStellarType
    number::Int
    PreMainSequence() = new(17)
end

"""
Planet with index 18.
"""
struct Planet <: AbstractStellarType
    number::Int
    Planet() = new(18)
end

"""
Brown dwarf with index 19.
"""
struct BrownDwarf <: AbstractStellarType
    number::Int
    BrownDwarf() = new(19)
end

# abstract type Star                 <: StellarType end
# abstract type Giant                <: Star        end
# abstract type StellarRemnant       <: StellarType end
# abstract type SubStellarObject     <: StellarType end
# abstract type Other                <: StellarType end
# abstract type WhiteDwarf           <: StellarRemnant  end
# abstract type CompactObject        <: StellarRemnant  end

const Star = Union{DeeplyOrFullyConvectiveLowMassMainSequence, MainSequence, HertzsprungGap, FirstGiantBranch, CoreHeliumBurning, FirstAsymptoticGiantBranch, SecondAsymptoticGiantBranch, MainSequenceNakedHelium, HertzsprungGapNakedHelium, GiantBranchNakedHelium}
const Giant = Union{HertzsprungGap, FirstGiantBranch, CoreHeliumBurning, FirstAsymptoticGiantBranch, SecondAsymptoticGiantBranch, HertzsprungGapNakedHelium, GiantBranchNakedHelium}
const StellarRemnant = Union{HeliumWhiteDwarf, CarbonOxygenWhiteDwarf, OxygenNeonWhiteDwarf, NeutronStar, BlackHole, MasslessSupernova}
const SubStellarObject = Union{Planet, BrownDwarf}
const WhiteDwarf = Union{HeliumWhiteDwarf, CarbonOxygenWhiteDwarf, OxygenNeonWhiteDwarf}
const CompactObject = Union{NeutronStar, BlackHole}

const  stellar_types = Dict{Int, AbstractStellarType}(
                                0  => DeeplyOrFullyConvectiveLowMassMainSequence(),
                                1  => MainSequence(),
                                2  => HertzsprungGap(),
                                3  => FirstGiantBranch(),
                                4  => CoreHeliumBurning(),
                                5  => FirstAsymptoticGiantBranch(),
                                6  => SecondAsymptoticGiantBranch(),
                                7  => MainSequenceNakedHelium(),
                                8  => HertzsprungGapNakedHelium(),
                                9  => GiantBranchNakedHelium(),
                                10 => HeliumWhiteDwarf(),
                                11 => CarbonOxygenWhiteDwarf(),
                                12 => OxygenNeonWhiteDwarf(),
                                13 => NeutronStar(),
                                14 => BlackHole(),
                                15 => MasslessSupernova(),
                                16 => UnknownStellarType(),
                                17 => PreMainSequence(),
                                18 => Planet(),
                                19 => BrownDwarf()
                            )

function stellar_type_from_index(index)
    @assert (0 <= index <= 19) "Stellar index must be in the range {0, 19}."
    stellar_types[index]
end



