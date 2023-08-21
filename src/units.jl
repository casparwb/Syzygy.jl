using Unitful, UnitfulAstro
import DynamicQuantities

CGS_units() = Unitful.preferunits(u"cm, g, s")

@unit stp "stp" StellarType 1.0u"1" false
Unitful.register(Syzygy)