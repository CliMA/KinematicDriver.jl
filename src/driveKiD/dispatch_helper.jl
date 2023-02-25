"""
    Helper functions to arrange CloudMicrophysics.jl
    autoconversion arguments for easier dispatch.
"""
conv_args(::CMT.KK2000Type, params) = (N_d = params.prescribed_Nd,)
conv_args(::CMT.B1994Type, params) = (N_d = params.prescribed_Nd, smooth_transition = true)
conv_args(::CMT.TC1980Type, params) = (N_d = params.prescribed_Nd, smooth_transition = true)
conv_args(::CMT.LD2004Type, params) = (N_d = params.prescribed_Nd, smooth_transition = true)

"""
    Helper functions to arrange CloudMicrophysics.jl
    accretion arguments for easier dispatch.
"""
accr_args(rf::CMT.KK2000Type, ql, qr, ρ) = (rf, ql, qr, ρ)
accr_args(rf::CMT.B1994Type, ql, qr, ρ) = (rf, ql, qr, ρ)
accr_args(rf::CMT.TC1980Type, ql, qr, ρ) = (rf, ql, qr)
