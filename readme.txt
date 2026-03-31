================================================================================
kOmegaDynamic RAS Turbulence Model: Implementation and Configuration
================================================================================

This document reports the implemented equations and the user-configurable 
parameters for the kOmegaDynamic turbulence model.

--------------------------------------------------------------------------------
1. MATHEMATICAL FORMULATION
--------------------------------------------------------------------------------

1.1 Transport Equations
The transport equations for the dynamic k-omega model with the zero-decay 
adjustment applied to the filtered variables (subscript t) are:

k-equation:
dk_t/dt + d(k_t * <u_j>_t)/dx_j = 
    2 * nu_t * <S_ij>_t * <S_ij>_t 
    - c_mu* * beta* * k_t * omega_t 
    + d/dx_j [ (nu + nu_t/sigma_k) * dk_t/dx_j ]

omega-equation:
d_omega_t/dt + d(omega_t * <u_j>_t)/dx_j = 
    2 * gamma * (omega_t/k_t) * nu_t * <S_ij>_t * <S_ij>_t 
    - c_mu* * beta * omega_t^2 
    + d/dx_j [ (nu + nu_t/sigma_w) * d_omega_t/dx_j ]

Note: In OpenFOAM, alphaK = 1/sigma_k and alphaOmega = 1/sigma_w.

1.2 Dynamic Procedure and Eddy Viscosity
The turbulent eddy viscosity is defined as:
    nu_t = c_mu * k_t / (beta* * omega_t)

The dynamic coefficient c_mu (bounded 0.0 <= c_mu <= 0.2) is evaluated using 
the Germano-like procedure:
    c_mu = (M_ij * L_ij) / (M_ij * M_ij)

Where:
    L_ij = <<u_i>_t * <u_j>_t>_T - <<u_i>_t>_T * <<u_j>_t>_T
    M_ij = 2 * <(k_t / (beta* * omega_t)) * <S_ij>_t>_T 
           - 2 * (k_T_t / (beta* * omega_T_t)) * <<S_ij>_t>_T

1.3 Zero-Decay Modification
The normalized dynamic coefficient c_mu* is used in the destruction terms:
    c_mu* = c_mu / c_mu_0 (where c_mu_0 = 0.09)

If "Zero-Decay" is set to "no", c_mu* defaults to 1.0.

--------------------------------------------------------------------------------
2. OPENFOAM CONFIGURATION (constant/momentumTransport)
--------------------------------------------------------------------------------

simulationType RAS;
RAS
{
    RASModel        kOmegaDynamic;
    turbulence      on;
    printCoeffs     on; 
    kOmegaDynamicCoeffs
    {
        // Standard k-omega parameters
        betaStar        0.09;
        beta            0.075;       
        gamma           0.55;
        alphaK          0.5;        
        alphaOmega      0.5;        

        // Dynamic procedure parameters
        Cmu_0           0.09;       // Base reference value
        dynamicCmuMin   0.0;        // Lower clipping limit 
        dynamicCmuMax   0.2;        // Upper clipping limit 

        // Averaging and stabilization
        nWindow         2;          // Buffer size for temporal averaging
        nStart          10;         // Time steps before dynamic procedure
        
        // Model variant switch
        Zero-Decay      yes;        // Activate Zero-Decay modification (or no)
    }
}

--------------------------------------------------------------------------------
3. PARAMETERS  
--------------------------------------------------------------------------------

| Parameter    | Description                                      | Default |
|--------------|--------------------------------------------------|---------|
| betaStar     | Constant for turbulence dissipation rate         | 0.09    |
| beta         | Constant for specific dissipation destruction    | 0.075   |
| gamma        | Constant for specific dissipation production     | 0.55    |
| alphaK       | Inverse Prandtl number for k (1/sigma_k)         | 0.5     |
| alphaOmega   | Inverse Prandtl number for omega (1/sigma_w)     | 0.5     |
| Cmu_0        | Base reference value for the Cmu coefficient     | 0.09    |
| nWindow      | Time steps for temporal averaging window         | 2       |
| nStart       | Warm-up steps before dynamic procedure starts    | 10      |
| Zero-Decay   | Toggle for Cmu* normalization (yes/no)           | yes     |

================================================================================
