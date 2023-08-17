using Unitful, UnitfulAstro

CGS_units() = Unitful.preferunits(u"cm, g, s")

@unit stp "stp" StellarType 1.0u"1" false
Unitful.register(Syzygy)