#================================
FILE DEFINING THE FOLLOWING:
 - ODE MODEL
 - FUNCTION TO ESTIMATE 
   HEMOGLOBIN FROM RBC LEVELS
================================#

function model_iron_homeostasis!(du, u, p, t)

    HMO, FOS, Inulin, VitC, LfFe, LumenFe3, LumenFe2, LIP, Ft, FtFe, B, P, RBC, TFe, H, MFe = u; # Variables

    #=====
    Parameters
    =====#

    # HMO
    D_H = p[1]; κ_H = p[2]; K_H = p[3]; ϕ_H = p[4];

    # FOS
    D_F = p[5]; κ_F = p[6]; K_F = p[7]; ϕ_F = p[8];

    # Inulin
    D_I = p[9]; κ_I = p[10]; K_I = p[11]; ϕ_I = p[12];

    # VitC
    D_V = p[13]; ϕ_V = p[14];

    # LfFe
    D_L = p[15]; k_L = p[16]; ϕ_L = p[17];

    # LumenFe3
    D_Fe3 = p[18]; c_max = p[19]; K_pH = p[20]; K_Vit = p[21]; ϕ_Fe3 = p[22];

    # LumenFe2
    D_Fe2 = p[23]; κ_Fe = p[24]; K_Fe = p[25]; k_enter = p[26]; β = p[27]; ϕ_Fe2 = p[28];

    # LIP
    k_1 = p[29]; k_2 = p[30]; k_utilize = p[31]; k_fpn = p[32]; θ = p[33];

    # B
    μ_B = p[34]; η_B0 = p[35]; a_H = p[36]; a_F = p[37]; a_I = p[38]; α_BP = p[39];

    # P
    μ_P = p[40]; η_P0 = p[41]; a_Fe = p[42];

    # RBC
    p_max = p[43]; ϵ = p[44]; ζ = p[45]; δ_R = p[46];

    # TFe
    f_p = p[47]; γ = p[48]; δ_TFe = p[49];

    # H
    h = p[50]; K_T = p[51]; K_p = p[52]; δ_H = p[53];

    # MFe
    δ_M = p[54];

    # pH
    ω = p[55]; ν = p[56];

    #=====
    Equations
    =====#

    a_HB  = 1/κ_H;
    a_HP  = 6/(7*κ_H);
    a_F   = 1.5/κ_F;
    a_I   = 1/κ_I;
    a_Fe  = 2/κ_Fe;

    η_B     = η_B0*(1 + a_HB*κ_H*HMO/(K_H + HMO))*(1 + a_F*κ_F*FOS/(K_F + FOS))*(1 + a_I*κ_I*Inulin/(K_I + Inulin));
    η_P     = η_P0*(1 - a_HP*κ_H*HMO/(K_H + HMO))*(1 + a_Fe*κ_Fe*LumenFe2/(K_Fe + LumenFe2));
    pRBC    = p_max*(ϵ/(ϵ + RBC))*(TFe/(ζ + TFe));
    pH      = ω - ν*log10(B);

    du[1]  = dHMO        = D_H - κ_H*HMO*(B + P)/(K_H + HMO) - ϕ_H*HMO;
    du[2]  = dFOS        = D_F - κ_F*FOS*B/(K_F + FOS) - ϕ_F*FOS;
    du[3]  = dInulin     = D_I - κ_I*Inulin*B/(K_I + Inulin) - ϕ_I*Inulin;
    du[4]  = dVitC       = D_V - ϕ_V*VitC;
    du[5]  = dLfFe       = D_L - k_L*LfFe - ϕ_L*LfFe;
    du[6]  = dLumenFe3   = D_Fe3 - c_max*LumenFe3*(K_pH/(K_pH + pH))*(1 + (VitC/(K_Vit + VitC))) - ϕ_Fe3*LumenFe3;
    du[7]  = dLumenFe2   = D_Fe2 + c_max*LumenFe3*(K_pH/(K_pH + pH))*(1 + (VitC/(K_Vit + VitC))) - κ_Fe*P*LumenFe2/(K_Fe + LumenFe2) - k_enter*LumenFe2*β/(β + FtFe) - ϕ_Fe2*LumenFe2;
    du[8]  = dLIP        = k_enter*LumenFe2*β/(β + FtFe) + k_L*LfFe + k_1*FtFe - k_2*Ft*LIP - k_utilize*LIP - k_fpn*LIP*θ/(θ + H);
    du[9]  = dFt         = k_1*FtFe - k_2*LIP*Ft;
    du[10] = dFtFe       = k_2*LIP*Ft - k_1*FtFe;
    du[11] = dB          = μ_B*B*(1 - B/η_B) - α_BP*B*P;
    du[12] = dP          = μ_P*P*(1 - P/η_P) - α_BP*B*P;
    du[13] = dRBC        = pRBC - δ_R*RBC;
    du[14] = dTFe        = k_fpn*LIP*θ/(θ + H) - pRBC*f_p*TFe + δ_R*RBC*f_p*TFe*γ/(γ + H) - δ_TFe*TFe;
    du[15] = dH          = h*(TFe/(K_T + TFe))*(1 - (pRBC)/(K_p + pRBC)) - δ_H*H;
    du[16] = dMFe        = δ_R*RBC*f_p*TFe*(1 - (γ)/(γ + H)) - δ_M*MFe;

end

function Hb_estimator!(x)
    # x is the number of RBCs per kgbw

    η₁ = 8.2E-14;       # Volume of an RBC (L)
    η₂ = 0.560/7.0;     # Blood volume (L) per kgbw
    Hb = (100/3)*(η₁/η₂)*x;
    return Hb; # Hemoglobin levels in g/dL
end

