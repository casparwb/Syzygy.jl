
abstract type StellarType end
abstract type Star                 <: StellarType end
abstract type Giant                <: Star        end
abstract type StellarRemnant       <: StellarType end
abstract type SubStellarObject     <: StellarType end
abstract type Other                <: StellarType end
abstract type WhiteDwarf           <: StellarRemnant  end
abstract type CompactObject        <: StellarRemnant  end

"""
Deeply or fully convective low main sequence star with index 0. 
"""
struct DeeplyOrFullyConvectiveLowMassMainSequence <: Star
    number::Int
    DeeplyOrFullyConvectiveLowMassMainSequence() = new(0)
end

"""
Main sequence star with index 1. 
"""
struct MainSequence <: Star
    number::Int
    MainSequence() = new(1)
end

"""
Hertzsprung gap star with index 2. 
"""
struct HertzsprungGap <: Star# Giant
    number::Int
    HertzsprungGap() = new(2)
end

"""
First giant branch star with index 3. 
"""
struct FirstGiantBranch <: Star# Giant
    number::Int
    FirstGiantBranch() = new(3)
end

"""
Core helium burning star with index 4. 
"""
struct CoreHeliumBurning <: Star# Giant
    number::Int
    CoreHeliumBurning() = new(4)
end

"""
First asymptotic gant branch star with index 5.
"""
struct FirstAsymptoticGiantBranch <: Star# Giant
    number::Int
    FirstAsymptoticGiantBranch() = new(5)
end

"""
Second asymptotic giant branch star with index 6.
"""
struct SecondAsymptoticGiantBranch <: Star# Giant
    number::Int
    SecondAsymptoticGiantBranch() = new(6)
end

"""
Main sequence naked helium star with index 7.
"""
struct MainSequenceNakedHelium <: Star
    number::Int
    MainSequenceNakedHelium() = new(7)
end

"""
Hertzsprung gap naked helium star with index 8.
"""
struct HertzsprungGapNakedHelium <: Star# Giant
    number::Int
    HertzsprungGapNakedHelium() = new(8)
end

"""
Giant branch naked helium star with index 9.
"""
struct GiantBranchNakedHelium <: Star# Giant
    number::Int
    GiantBranchNakedHelium() = new(9)
end

"""
Helium white dwarf with index 10.
"""
struct HeliumWhiteDwarf <: WhiteDwarf
    number::Int
    HeliumWhiteDwarf() = new(10)
end

"""
Carbon/oxygen white dwarf with index 11.
"""
struct CarbonOxygenWhiteDwarf <: WhiteDwarf
    number::Int
    CarbonOxygenWhiteDwarf() = new(11)
end

"""
Oxygen/neon white dwarf with index 12.
"""
struct OxygenNeonWhiteDwarf <: WhiteDwarf
    number::Int
    OxygenNeonWhiteDwarf() = new(12)
end

"""
Neutron star with index 13.
"""
struct NeutronStar <: CompactObject
    number::Int
    NeutronStar() = new(13)
end

"""
Black hole with index 14.
"""
struct BlackHole <: CompactObject
    number::Int
    BlackHole() = new(14)
end

"""
Massless supernova with index 15.
"""
struct MasslessSupernova <: Other
    number::Int
    MasslessSupernova() = new(15)
end

"""
Unknown stellar type with index 16.
"""
struct UnknownStellarType <: Other
    number::Int
    UnknownStellarType() = new(16)
end

"""
Pre main sequence star with index 17.
"""
struct PreMainSequence <: Star
    number::Int
    PreMainSequence() = new(17)
end

"""
Planet with index 18.
"""
struct Planet <: SubStellarObject
    number::Int
    Planet() = new(18)
end

"""
Brown dwarf with index 19.
"""
struct BrownDwarf <: SubStellarObject
    number::Int
    BrownDwarf() = new(19)
end

const  stellar_types = Dict{Int, StellarType}(
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



