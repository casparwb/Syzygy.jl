
abstract type StellarType end
abstract type Star             <: StellarType end
abstract type CompactObject    <: StellarType end
abstract type SubStellarObject <: StellarType end
abstract type Other            <: StellarType end
abstract type WhiteDwarf       <: CompactObject end

"""
Deeply or fully convective low main sequence star with index 0. 
"""
struct DeeplyOrFullyConvectiveLowMassMainSequence <: Star
    index::Int
    DeeplyOrFullyConvectiveLowMassMainSequence() = new(0)
end

"""
Main sequence star with index 1. 
"""
struct MainSequence <: Star
    index::Int
    MainSequence() = new(1)
end

"""
Hertzsprung gap star with index 2. 
"""
struct HertzsprungGap <: Star
    index::Int
    HertzsprungGap() = new(2)
end

"""
First giant branch star with index 3. 
"""
struct FirstGiantBranch <: Star
    index::Int
    FirstGiantBranch() = new(3)
end

"""
Core helium burning star with index 4. 
"""
struct CoreHeliumBurning <: Star
    index::Int
    CoreHeliumBurning() = new(4)
end

"""
First asymptotic gant branch star with index 5.
"""
struct FirstAsymptoticGiantBranch <: Star
    index::Int
    FirstAsymptoticGiantBranch() = new(5)
end

"""
Second asymptotic giant branch star with index 6.
"""
struct SecondAsymptoticGiantBranch <: Star
    index::Int
    SecondAsymptoticGiantBranch() = new(6)
end

"""
Main sequence naked helium star with index 7.
"""
struct MainSequenceNakedHelium <: Star
    index::Int
    MainSequenceNakedHelium() = new(7)
end

"""
Hertzsprung gap naked helium star with index 8.
"""
struct HertzsprungGapNakedHelium <: Star
    index::Int
    HertzsprungGapNakedHelium() = new(8)
end

"""
Giant branch naked helium star with index 9.
"""
struct GiantBranchNakedHelium <: Star
    index::Int
    GiantBranchNakedHelium() = new(9)
end

"""
Helium white dwarf with index 10.
"""
struct HeliumWhiteDwarf <: WhiteDwarf
    index::Int
    HeliumWhiteDwarf() = new(10)
end

"""
Carbon/oxygen white dwarf with index 11.
"""
struct CarbonOxygenWhiteDwarf <: WhiteDwarf
    index::Int
    CarbonOxygenWhiteDwarf() = new(11)
end

"""
Oxygen/neon white dwarf with index 12.
"""
struct OxygenNeonWhiteDwarf <: WhiteDwarf
    index::Int
    OxygenNeonWhiteDwarf() = new(12)
end

"""
Neutron star with index 13.
"""
struct NeutronStar <: CompactObject
    index::Int
    NeutronStar() = new(13)
end

"""
Black hole with index 14.
"""
struct BlackHole <: CompactObject
    index::Int
    BlackHole() = new(14)
end

"""
Massless supernova with index 15.
"""
struct MasslessSupernova <: Other
    index::Int
    MasslessSupernova() = new(15)
end

"""
Unknown stellar type with index 16.
"""
struct UnknownStellarType <: Other
    index::Int
    UnknownStellarType() = new(16)
end

"""
Pre main sequence star with index 17.
"""
struct PreMainSequence <: Star
    index::Int
    PreMainSequence() = new(17)
end

"""
Planet with index 18.
"""
struct Planet <: SubStellarObject
    index::Int
    Planet() = new(18)
end

"""
Brown dwarf with index 19.
"""
struct BrownDwarf <: SubStellarObject
    index::Int
    BrownDwarf() = new(19)
end

const  stellar_types = Dict(
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










