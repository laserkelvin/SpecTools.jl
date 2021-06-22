
using PhysicalConstants.CODATA2018: c_0, h, k_B

# Boltzmann constant in cm/K
const k_cm = Float64((k_B * (1 / h) * (1 / (c_0 * 100))).val)

# grabbed from https://physics.nist.gov/cgi-bin/cuu/Value?kshhz
const k_mhz = 2.083661912e7
# c, h, k_B = c_0.val, h.val, k_B.val

export
    c_0,
    k_B,
    h,
    k_cm,
    MHz2wavenumber

MHz2wavenumber(ν::AbstractFloat) = ν / (c_0.val / 1e4)
