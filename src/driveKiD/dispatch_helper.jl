"""
    Helper functions to arrange CloudMicrophysics.jl
    autoconversion arguments for easier dispatch.
"""
conv_args(::CMP.KK2000, params) = (N_d = params.prescribed_Nd,)
conv_args(::CMP.B1994, params) = (N_d = params.prescribed_Nd, smooth_transition = true)
conv_args(::CMP.TC1980, params) = (N_d = params.prescribed_Nd, smooth_transition = true)
conv_args(::CMP.LD2004, params) = (N_d = params.prescribed_Nd, smooth_transition = true)
conv_args(::CMP.VarTimescaleAcnv, params) = (N_d = params.prescribed_Nd,)

"""
    Helper functions to arrange CloudMicrophysics.jl
    accretion arguments for easier dispatch.
"""
accr_args(rf::CMP.KK2000, ql, qr, ρ) = (rf, ql, qr, ρ)
accr_args(rf::CMP.B1994, ql, qr, ρ) = (rf, ql, qr, ρ)
accr_args(rf::CMP.TC1980, ql, qr, ρ) = (rf, ql, qr)
