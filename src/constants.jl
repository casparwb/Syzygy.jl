using StaticArrays

const GRAVCONST = 6.6743015e-11u"m^3/kg/s^2" |> upreferred
const UNITLESS_G = GRAVCONST.val
const G² = UNITLESS_G^2
const G³ = UNITLESS_G^3
const G⁴ = UNITLESS_G^4

const speed_of_light = 299_792_458u"m/s" |> upreferred
const UNITLESS_c = ustrip(unit_length/unit_time, speed_of_light)
const c² = UNITLESS_c*UNITLESS_c
const c³ = c²*UNITLESS_c
const c⁴ = c³*UNITLESS_c
const c⁻² = 1/c²
const c⁻³ = 1/c³
const c⁻⁴ = c⁻²/c²
const c⁻⁵ = c⁻⁴/UNITLESS_c
const c⁻⁶ = c⁻⁵/UNITLESS_c
const c⁻⁷ = c⁻⁶/UNITLESS_c

const π² = π^2

const surface_gravity_unit_conversion_factor = ustrip(u"cm/s^2", 1.0*unit_length/unit_time^2)

const stellar_type = Dict(  "deeply or fully convective low mass MS star" => 0,
                            "Main Sequence Star" => 1,
                            "Hertzsprung Gap" => 2,
                            "First Giant Branch" => 3,
                            "Core Helium Burning" => 4,
                            "First Asymptotic Giant Branch" => 5,
                            "Second Asymptotic Giant Branch" => 6,
                            "Main Sequence Naked Helium Star" => 7,
                            "Hertzsprung Gap Naked Helium Star" => 8,
                            "Giant Branch Naked Helium Star" => 9,
                            "Helium White Dwarf" => 10,
                            "Carbon/Oxygen White Dwarf" => 11,
                            "Oxygen/Neon White Dwarf" => 12,
                            "Neutron Star" => 13,
                            "Black Hole" => 14,
                            "Massless Supernova" => 15,
                            "Unknown stellar type" => 16,
                            "Pre-main-sequence Star" => 17,
                            "Planet" => 18,
                            "Brown Dwarf" => 19)


const stellar_type_index = Dict(5  => "First Asymptotic Giant Branch",
                                16 => "Unknown stellar type",
                                12 => "Oxygen/Neon White Dwarf",
                                8  => "Hertzsprung Gap Naked Helium star",
                                17 => "Pre-main-sequence Star",
                                1  => "Main Sequence star",
                                19 => "Brown Dwarf",
                                0  => "deeply or fully convective low mass MS star",
                                6  => "Second Asymptotic Giant Branch",
                                11 => "Carbon/Oxygen White Dwarf",
                                9  => "Giant Branch Naked Helium star",
                                14 => "Black Hole",
                                3  => "First Giant Branch",
                                7  => "Main Sequence Naked Helium star",
                                13 => "Neutron Star",
                                4  => "Core Helium Burning",
                                15 => "Massless Supernova",
                                2  => "Hertzsprung Gap",
                                10 => "Helium White Dwarf",
                                18 => "Planet")


const tidal_structure_coefficients = Dict{Tuple{Float64, Int64}, SVector{6, Float64}}(
                                     (1.5, 2) => SA[-0.397, 1.678, 1.277, -12.42, 9.446, -5.550],
                                     (1.5, 3) => SA[-0.909, 1.574, 12.37, -57.40, 80.10, -46.43],
                                     (2.0, 2) => SA[-0.517, -0.906, 23.88, -93.49, 112.3, -44.15],
                                     (2.0, 3) => SA[-1.040, -1.354, 37.64, -139.9, 168.2, -66.53],
                                     (3.0, 2) => SA[-1.124, 0.877, -13.37, 21.55, -16.48, 4.124],
                                     (3.0, 3) => SA[-1.703, 2.653, -14.34, 12.85, -0.492, -3.60],
                                     (0.0, 2) => SA[0.0,    0.0,   0.0,    0.0,   0.0,    0.0],
                                     (0.0, 3) => SA[0.0,    0.0,   0.0,    0.0,   0.0,    0.0],

)

const orbital_elements = (:a, :P, :e, :ω, :i, :Ω, :ν)

const hierarchy_labels = ["Primary", "Secondary", "Tertiary", 
                          "Quaternary", "Quinary", "Senary", 
                          "Septenary", "Octonary", "Oonary", "Denary"]

