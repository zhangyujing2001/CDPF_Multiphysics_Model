/*-----------------------*- Fluent UDF -*-----------------------*\
  @File name        chemicalReactions.c
  @Author           Yujing Zhang
  @Email            zhangyujing@tongji.edu.cn
  @Created on       2025-10-24
  @Last modified    2025-10-28
\*--------------- Copyright (c) 2025 Yujing Zhang --------------*/

//- Header files
//  "udf.h" header file contains:
//  1. Definitions for DEFINE macros;
//  2. #include compiler directives for C and C++ library function
//  header files (e.g., "stdio.h", "math.h");
//  3. Header files (e.g., "mem.h") for other ANSYS Fluent-supplied
//  macros and functions.
#include "udf.h"

/* Constants */
#define GAS_CONST           8.314         // Molar gas constant (J.mol-1.K-1)
#define MOLAR_MASS_O2       0.032         // Molar mass of O2 (kg.mol-1)
#define MOLAR_MASS_NO2      0.046         // Molar mass of NO2 (kg.mol-1)
#define MOLAR_MASS_NO       0.030         // Molar mass of NO (kg.mol-1)
#define MOLAR_MASS_CO2      0.044         // Molar mass of CO2 (kg.mol-1)
#define MOLAR_MASS_N2       0.028         // Molar mass of N2 (kg.mol-1)
#define DELTA_H_R1         -393.5e3       // Standard enthalpy change of reaction R1 (J.mol-1)
#define DELTA_H_R2         -279.3e3       // Standard enthalpy change of reaction R2 (J.mol-1)
#define DELTA_H_R3         -336.4e3       // Standard enthalpy change of reaction R3 (J.mol-1)
#define DELTA_H_R4         -57.1e3        // Standard enthalpy change of reaction R4 (J.mol-1)
#define REF_TEMP            298.15        // Reference temperature (K)
#define PORE_SIZE           1.0e-5        // Pore size (m)
#define POROSITY            0.6           // Porosity (-)
#define OP_PRESSURE         101325        // Operating pressure (Pa)
#define SMALL_VALUE_1       1.0e-12       // Small value 1
#define SMALL_VALUE_2       1.0e-30       // Small value 2

/* Set species index */
enum speciesIndex
{
    OXYGEN,
    NITROGEN_DIOXIDE,
    NITRIC_OXIDE,
    CARBON_DIOXIDE,
    NITROGEN
};

/* Define a 'struct' to hold species info */
typedef struct
{
    int speciesID;           // Species index
    real molarMass;          // Molar mass (kg.mol-1)
    real switchTemp;         // Temperature at which coefficients switch (K)
    real lowTempCoeff[5];    // coefficients of Shomate equation at low temperature
    real highTempCoeff[5];   // coefficients of Shomate equation at high temperature
}
SpeciesInfo;

/* Define a 'struct' to hold reaction info */
typedef struct
{
    int reactantNum;        // Number of reactants (Max = 2)
    int reactants[2];       // Array of reactant species indices
    real reactOrders[2];    // Array of reaction orders for each reactant
    real freqFactor;        // Frequency factor (s-1)
    real actEnergy;         // Activation energy (J.mol-1)
}
ReactionInfo;

/* Basic species info from NIST Chemistry WebBook */
const SpeciesInfo O2info = {OXYGEN, MOLAR_MASS_O2, 700.0, 
                            {31.32234, -20.23531, 57.86644, -36.50624, -0.007374}, 
                            {30.03235, 8.772972, -3.988133, 0.788313, -0.741599}};
const SpeciesInfo NO2info = {NITROGEN_DIOXIDE, MOLAR_MASS_NO2, 1200.0, 
                            {16.10857, 75.89525, -54.38740, 14.30777, 0.239423}, 
                            {56.82541, 0.738053, -0.144721, 0.009777, -5.459911}};
const SpeciesInfo NOinfo = {NITRIC_OXIDE, MOLAR_MASS_NO, 1200.0, 
                            {23.83491, 12.58878, -1.139011, -1.497459, 0.214194}, 
                            {35.99169, 0.957170, -0.148032, 0.009974, -3.004088}};
const SpeciesInfo CO2info = {CARBON_DIOXIDE, MOLAR_MASS_CO2, 1200.0, 
                            {24.99735, 55.18696, -33.69137, 7.948387, -0.136638}, 
                            {58.16639, 2.720074, -0.492289, 0.038844, -6.447293}};
const SpeciesInfo N2info = {NITROGEN, MOLAR_MASS_N2, 500.0, 
                            {28.98641, 1.853978, -9.647459, 16.63537, 0.000117}, 
                            {19.50583, 19.88705, -8.598535, 1.369784, 0.527601}};

/* Reaction information */
const ReactionInfo reactionR11 = {1, {OXYGEN, 0}, {0.9, 0.0}, 8.50e7, 1.64e5};
const ReactionInfo reactionR12 = {1, {OXYGEN, 0}, {0.3, 0.0}, 1.19e5, 1.14e5};
const ReactionInfo reactionR21 = {1, {NITROGEN_DIOXIDE, 0}, {0.6, 0.0}, 0.48, 2.67e4};
const ReactionInfo reactionR22 = {1, {NITROGEN_DIOXIDE, 0}, {0.6, 0.0}, 0.51, 2.68e4};
const ReactionInfo reactionR31 = {2, {OXYGEN, NITROGEN_DIOXIDE}, {0.3, 1.0}, 1.395e3, 4.78e4};
const ReactionInfo reactionR32 = {2, {OXYGEN, NITROGEN_DIOXIDE}, {0.3, 0.4}, 51.4, 5.22e4};
const ReactionInfo reactionR4f = {2, {OXYGEN, NITRIC_OXIDE}, {0.5, 1.0}, 3.61e6, 8.18e4};
const ReactionInfo reactionR4r = {1, {NITROGEN_DIOXIDE, 0}, {1.0, 0.0}, 3.61e6, 8.18e4};

/**
 * Function to calculate the SSA (specific surface area).
 * @param poreSize     Pore size (m).
 * @param porosity     Porosity (-).
 * @return specific surface area (m-1).
 */
real getSSA(real poreSize, real porosity)
{
    real ssa = 4.0 * (1.0 - porosity) / poreSize;
    return ssa;
}

/**
 * Function to calculate the mole ratio n_i = w_i / M_i.
 * @param massFrac      Mass fraction of species i (-).
 * @param s             Pointer to 'SpeciesInfo' struct.
 * @return mole ratio (mol.kg-1).
 */
real getMoleRatio(real massFrac, const SpeciesInfo *s)
{   
    if (massFrac < 0.0)
    {
        massFrac = 0.0;
    }
    return massFrac / s->molarMass;
}

/** 
 * Function to calculate the partial pressure of species i.
 * @param moleFrac      Mole fraction of species i (-).
 * @param totalPress    Total pressure (Pa).
 * @return partial pressure of species i (Pa).
 */
real getPartialPressure(real moleFrac, real totalPressure)
{   
    // Small value
    real smallValue = SMALL_VALUE_1 * OP_PRESSURE;

    // Protection against invalid mole fraction
    if (!isfinite(moleFrac)) moleFrac = 0.0;
    if (moleFrac < 0.0) moleFrac = 0.0;
    if (moleFrac > 1.0) moleFrac = 1.0;

    // Protection against invalid total pressure
    if (!isfinite(totalPressure) || totalPressure <= 0.0)
    {
        totalPressure = OP_PRESSURE;
    }

    // Calculate partial pressure
    real p = moleFrac * totalPressure;
    if (!isfinite(p) || p < smallValue)
    {
        p = smallValue;
    }

    return p;
}

/**
 * Function to calculate the equilibrium constant.
 * @param nSpecies       Number of species.
 * @param pressure       Array that stores partial pressures of species (Pa).
 * @param stoichCoeffs   Array that stores stoichiometric coefficients of species.
 * @return equilibrium constant K_eq.
 */
real getEquilConst(int nSpecies, real pressure[], real stoichCoeffs[])
{   
    real smallValue = SMALL_VALUE_1 * OP_PRESSURE;
    real logEquilConst = 0.0;
    for (int i = 0; i < nSpecies; i++)
    {   
        real p = pressure[i];
        if (p < smallValue)
        {
            p = smallValue;
        }
        logEquilConst += stoichCoeffs[i] * log(p / OP_PRESSURE);
    }
    return exp(logEquilConst);
}

/**
 * Function to calculate the rate constant based on Arrhenius equation.
 * @param temperature    Temperature (K).
 * @param freqFactor     Frequency factor (s-1).
 * @param actEnergy      Activation energy (J.mol-1).
 * @return rate constant (s-1).
 */
real getRateConst(real temperature, real freqFactor, real actEnergy)
{
    real rateConst = freqFactor * exp(-actEnergy / (GAS_CONST * temperature));
    return rateConst;
}

/** 
 * Function to calculate the rate of reaction.
 * @param c                 UDF data-structure, cell index.
 * @param t                 UDF data-structure, pointer to cell thread.
 * @param r                 Pointer to 'ReactionInfo' struct.
 * @return reaction rate (mol.m-3.s-1).
 */
real getReactRate(cell_t c, Thread *t, const ReactionInfo *r)
{
    // Variables
    real rhoMixture = C_R(c, t);
    real tempMixture = C_T(c, t);
    real rateConst = getRateConst(tempMixture, r->freqFactor, r->actEnergy);
    real totMolarMass = MOLAR_MASS_O2 + MOLAR_MASS_NO2 + MOLAR_MASS_NO + 
                        MOLAR_MASS_CO2 + MOLAR_MASS_N2;
    real concMixture = rhoMixture / totMolarMass;

    real massFracProduct = 1.0;
    for (int i = 0; i < r->reactantNum; i++)
    {
        real massFrac = C_YI(c, t, r->reactants[i]);
        massFracProduct *= pow(massFrac, r->reactOrders[i]);
    }

    return rateConst * massFracProduct * concMixture;
}

/**
 * Calculate the specific heat capacity at standard condition using Shomate equation.
 * @param temperature        Temperature (K).
 * @param s                  Pointer to 'SpeciesInfo' struct.
 * @return specific heat capacity (J.mol-1.K-1).
 */
real getSpecificHeat(real temperature, const SpeciesInfo *s)
{
    // Variables
    real t = temperature / 1000.0;
    real coeffs[5];
    real cp;

    // Calculate Cp based on Shomate equation
    if (temperature <= s->switchTemp)
    {
        for (int i = 0; i < 5; i++)
        {
            coeffs[i] = s->lowTempCoeff[i];
        }
        cp = coeffs[0] + 
             coeffs[1] * t + 
             coeffs[2] * pow(t, 2) +
             coeffs[3] * pow(t, 3) + 
             coeffs[4] / pow(t, 2);
    }
    else
    {
        for (int i = 0; i < 5; i++)
        {
            coeffs[i] = s->highTempCoeff[i];
        }
        cp = coeffs[0] + 
             coeffs[1] * t + 
             coeffs[2] * pow(t, 2) +
             coeffs[3] * pow(t, 3) + 
             coeffs[4] / pow(t, 2);
    }
    return cp;
}

/**
 * Calculate the sensible enthaply at standard condition using numerical integration.
 * @param temperature         Temperature (K).
 * @param refTemp             Reference temperature (K).
 * @param s                   Pointer to 'SpeciesInfo' struct.
 * @return sensible enthalpy (J.mol-1).
 */
real getSensibleEnthalpy(real temperature, real refTemp, const SpeciesInfo *s)
{
    // Variables
    real steps = 100.0;
    real delta = (temperature - refTemp) / steps;
    real sum = 0.5 * (getSpecificHeat(refTemp, s) + getSpecificHeat(temperature, s));
    for (int i = 1; i < steps; i++)
    {
        real temp = refTemp + i * delta;
        sum += getSpecificHeat(temp, s);
    }
    return sum * delta;
}

/**
 * Calculate the partial pressure of O2, NO2, NO. This function is used in 
 * the source term calculation of O2, NO2, NO transport equations (R4 only).
 * @param c                 UDF data-structure, cell index.
 * @param t                 UDF data-structure, pointer to cell thread.
 * @param pressureArray     Array to store partial pressures of O2, NO2, NO (Pa).
 * @return void.
 */
void getPressureArray(cell_t c, Thread *t, real pressureArray[3])
{   
    // Mole ratio of each species in mixture.
    real moleRatioO2 = getMoleRatio(C_YI(c, t, OXYGEN), &O2info);
    real moleRatioNO2 = getMoleRatio(C_YI(c, t, NITROGEN_DIOXIDE), &NO2info);
    real moleRatioNO = getMoleRatio(C_YI(c, t, NITRIC_OXIDE), &NOinfo);
    real moleRatioCO2 = getMoleRatio(C_YI(c, t, CARBON_DIOXIDE), &CO2info);
    real moleRatioN2 = getMoleRatio(C_YI(c, t, NITROGEN), &N2info);
    
    // Protection against negative mole ratios
    moleRatioO2 = (moleRatioO2 > 0.0) ? moleRatioO2 : 0.0;
    moleRatioNO2 = (moleRatioNO2 > 0.0) ? moleRatioNO2 : 0.0;
    moleRatioNO = (moleRatioNO > 0.0) ? moleRatioNO : 0.0;
    moleRatioCO2 = (moleRatioCO2 > 0.0) ? moleRatioCO2 : 0.0;
    moleRatioN2 = (moleRatioN2 > 0.0) ? moleRatioN2 : 0.0;

    // Molar fraction of each species in R4.
    real totMoleRatio = moleRatioO2 + moleRatioNO2 +
                        moleRatioNO + moleRatioCO2 +
                        moleRatioN2;
    if (totMoleRatio < SMALL_VALUE_2)
    {
        totMoleRatio = SMALL_VALUE_2;
    }
    real moleFracO2 = moleRatioO2 / totMoleRatio;
    real moleFracNO2 = moleRatioNO2 / totMoleRatio;
    real moleFracNO = moleRatioNO / totMoleRatio;

    // Partial pressure of each species in R4 (Pa).
    real totPressure = C_P(c, t) + OP_PRESSURE;
    pressureArray[0] = getPartialPressure(moleFracO2, totPressure);
    pressureArray[1] = getPartialPressure(moleFracNO2, totPressure);
    pressureArray[2] = getPartialPressure(moleFracNO, totPressure);
}

/**
 * Mass flux of O2 due to reaction R11 + R12 + R31 + R32.
 * @param massFluxOxygen          UDF name.
 * @param t                       Pointer to thread on which boundary condition is to be applied.
 * @param i                       Index that identifies the variable that is to be defined.
 * @return void.
 */
DEFINE_PROFILE(massFluxOxygen, t, i)
{
#if !RP_HOST
    // Variables
    face_t f;
    real ssa = getSSA(PORE_SIZE, POROSITY);

    // Loop over faces in the thread,
    // compute the face value for the boundary variable,
    // and store the value in memory.
    begin_f_loop(f, t)
    {
        if (PRINCIPAL_FACE_P(f, t))
        {
            // Adjacent cell thread and index.
            Thread *t0 = THREAD_T0(t);
            cell_t c0 = F_C0(f, t);

            // Rate of reaction R11 & R12.
            real reactRateR11 = getReactRate(c0, t0, &reactionR11);
            real reactRateR12 = getReactRate(c0, t0, &reactionR12);

            // Rate of reaction R31.
            real reactRateR31 = getReactRate(c0, t0, &reactionR31);
            real reactRateR32 = getReactRate(c0, t0, &reactionR32);

            // Mass flux of O2 due to reactions R11, R12, R31 & R32. (kg.m-2.s-1).
            real massFluxO2 = MOLAR_MASS_O2 * (-reactRateR11 - reactRateR12 - 
                              0.5 * reactRateR31 - 0.5 * reactRateR32) / ssa;

            // 'F_PROFILE'
            F_PROFILE(f, t, i) = massFluxO2;
        }
    }
    end_f_loop(f, t)
#endif
}

/**
 * Mass flux of NO2 due to reaction R21 + R22 + R31 + R32.
 * @param massFluxNitrogenDioxide          UDF name.
 * @param t                                Pointer to thread on which boundary condition is to be applied.
 * @param i                                Index that identifies the variable that is to be defined.
 * @return void.
 */
DEFINE_PROFILE(massFluxNitrogenDioxide, t, i)
{
#if !RP_HOST
    // Variables
    face_t f;
    real ssa = getSSA(PORE_SIZE, POROSITY);

    // Loop over faces in the thread,
    // compute the face value for the boundary variable,
    // and store the value in memory.
    begin_f_loop(f, t)
    {
        if (PRINCIPAL_FACE_P(f, t))
        {
            // Adjacent cell thread and index.
            Thread *t0 = THREAD_T0(t);
            cell_t c0 = F_C0(f, t);

            // Rate of reaction R21 & R22.
            real reactRateR21 = getReactRate(c0, t0, &reactionR21);
            real reactRateR22 = getReactRate(c0, t0, &reactionR22);

            // Rate of reaction R31 & R32.
            real reactRateR31 = getReactRate(c0, t0, &reactionR31);
            real reactRateR32 = getReactRate(c0, t0, &reactionR32);

            // Mass flux of NO2 due to reactions R21, R22, R31 & R32 (kg.m-2.s-1).
            real massFluxNO2 = MOLAR_MASS_NO2 * (-2.0 * reactRateR21 -2.0 * reactRateR22 - 
                               reactRateR31 - reactRateR32) / ssa;

            // 'F_PROFILE'
            F_PROFILE(f, t, i) = massFluxNO2;
        }
    }
    end_f_loop(f, t)
#endif
}

/**
 * Mass flux of NO due to reaction R21 + R22 + R31 + R32.
 * @param massFluxNitricOxide          UDF name.
 * @param t                            Pointer to thread on which boundary condition is to be applied.
 * @param i                            Index that identifies the variable that is to be defined.
 * @return void.
 */
DEFINE_PROFILE(massFluxNitricOxide, t, i)
{
#if !RP_HOST
    // Variables
    face_t f;
    real ssa = getSSA(PORE_SIZE, POROSITY);

    // Loop over faces in the thread,
    // compute the face value for the boundary variable,
    // and store the value in memory.
    begin_f_loop(f, t)
    {
        if (PRINCIPAL_FACE_P(f, t))
        {
            // Adjacent cell thread and index.
            Thread *t0 = THREAD_T0(t);
            cell_t c0 = F_C0(f, t);

            // Rate of reaction R21 & R22.
            real reactRateR21 = getReactRate(c0, t0, &reactionR21);
            real reactRateR22 = getReactRate(c0, t0, &reactionR22);

            // Rate of reaction R31 & R32.
            real reactRateR31 = getReactRate(c0, t0, &reactionR31);
            real reactRateR32 = getReactRate(c0, t0, &reactionR32);

            // Mass flux of NO due to reactions R21 and R31 (kg.m-2.s-1).
            real massFluxNO = MOLAR_MASS_NO * (2.0 * reactRateR21 + 2.0 * reactRateR22 + 
                              reactRateR31 + reactRateR32) / ssa;

            // 'F_PROFILE'
            F_PROFILE(f, t, i) = massFluxNO;
        }
    }
    end_f_loop(f, t)
#endif
}

/**
 * Mass flux of CO2 due to reaction R11 + R12 + R21 + R22 + R31 + R32.
 * @param massFluxCarbonDioxide          UDF name.
 * @param t                              Pointer to thread on which boundary condition is to be applied.
 * @param i                              Index that identifies the variable that is to be defined.            
 * @return void.
 */
DEFINE_PROFILE(massFluxCarbonDioxide, t, i)
{
#if !RP_HOST
    // Variables
    face_t f;
    real ssa = getSSA(PORE_SIZE, POROSITY);

    // Loop over faces in the thread,
    // compute the face value for the boundary variable,
    // and store the value in memory.
    begin_f_loop(f, t)
    {
        if (PRINCIPAL_FACE_P(f, t))
        {
            // Adjacent cell thread and index.
            Thread *t0 = THREAD_T0(t);
            cell_t c0 = F_C0(f, t);

            // Rate of reaction R11 & R12.
            real reactRateR11 = getReactRate(c0, t0, &reactionR11);
            real reactRateR12 = getReactRate(c0, t0, &reactionR12);

            // Rate of reaction R21 & R22.
            real reactRateR21 = getReactRate(c0, t0, &reactionR21);
            real reactRateR22 = getReactRate(c0, t0, &reactionR22);

            // Rate of reaction R31 & R32.
            real reactRateR31 = getReactRate(c0, t0, &reactionR31);
            real reactRateR32 = getReactRate(c0, t0, &reactionR32);

            // Mass flux of CO2 due to reactions R11, R12, R21, R22, R31 & R32. (kg.m-2.s-1).
            real massFluxCO2 = MOLAR_MASS_CO2 * (reactRateR11 + reactRateR12 + 
                               reactRateR21 + reactRateR22 + reactRateR31 + reactRateR32) / ssa;

            // 'F_PROFILE'
            F_PROFILE(f, t, i) = massFluxCO2;
        }
    }
    end_f_loop(f, t)
#endif
}

/**
 * Heat flux due to chemical reactions.
 * @param heatFlux          UDF name.
 * @param t                 Pointer to thread on which boundary condition is to be applied.
 * @param i                 Index that identifies the variable that is to be defined.
 */
DEFINE_PROFILE(heatFlux, t, i)
{
#if !RP_HOST
    /* Variables */
    face_t f;
    real A[ND_ND], Amag;                     // Geometric variables for boundaries
    real gradT[ND_ND];                       // Temperature gradient vector
    real ssa = getSSA(PORE_SIZE, POROSITY);  // Specific surface area
    Thread *t0 = t->t0;                      // Adjacent cell thread

    /* Do nothing if areas are not computed yet or not next to fluid */
    if (!Data_Valid_P() || !FLUID_THREAD_P(t0))
    {
        return;
    }

    /* Loop over faces in the thread,
     * compute the face value for the boundary variable,
     * and store the value in memory.
     */
    begin_f_loop(f, t)
    {
        if (PRINCIPAL_FACE_P(f, t))
        {
            // Adjacent cell index
            cell_t c0 = F_C0(f, t);

            // Temperature
            real tempMixture = C_T(c0, t0);

            // Rate of reaction R11, R12, R21, R22, R31 & R32.
            real reactRateR11 = getReactRate(c0, t0, &reactionR11) / ssa;
            real reactRateR12 = getReactRate(c0, t0, &reactionR12) / ssa;
            real reactRateR21 = getReactRate(c0, t0, &reactionR21) / ssa;
            real reactRateR22 = getReactRate(c0, t0, &reactionR22) / ssa;
            real reactRateR31 = getReactRate(c0, t0, &reactionR31) / ssa;
            real reactRateR32 = getReactRate(c0, t0, &reactionR32) / ssa;

            // Sensible enthalpy changes of reaction R1, R2 and R3 (J.mol-1).
            real enthalpyChangeR1 = getSensibleEnthalpy(tempMixture, REF_TEMP, &CO2info) - 
                                     getSensibleEnthalpy(tempMixture, REF_TEMP, &O2info);
            real enthalpyChangeR2 = getSensibleEnthalpy(tempMixture, REF_TEMP, &CO2info) + 
                                     2.0 * getSensibleEnthalpy(tempMixture, REF_TEMP, &NOinfo) - 
                                     2.0 * getSensibleEnthalpy(tempMixture, REF_TEMP, &NO2info);
            real enthalpyChangeR3 = getSensibleEnthalpy(tempMixture, REF_TEMP, &CO2info) + 
                                     getSensibleEnthalpy(tempMixture, REF_TEMP, &NOinfo) - 
                                     getSensibleEnthalpy(tempMixture, REF_TEMP, &NO2info) - 
                                     0.5 * getSensibleEnthalpy(tempMixture, REF_TEMP, &O2info);
            
            // Standard enthalpy changes of reaction R1, R2 and R3 (J.mol-1).
            real stdEnthalpyR1 = DELTA_H_R1 + enthalpyChangeR1;
            real stdEnthalpyR2 = DELTA_H_R2 + enthalpyChangeR2;
            real stdEnthalpyR3 = DELTA_H_R3 + enthalpyChangeR3;

            // Heat flux due to chemical reactions (W.m-2).
            real reactHeatFlux = -((reactRateR11 + reactRateR12) * stdEnthalpyR1 + 
                                   (reactRateR21 + reactRateR22) * stdEnthalpyR2 + 
                                   (reactRateR31 + reactRateR32) * stdEnthalpyR3);
            
            // Info of boundary face geometry
            F_AREA(A, f, t);
            Amag = NV_MAG(A);

            // Reconstruction gradient of temperature at the boundary face
            if (NULL != T_STORAGE_R_NV(t0, SV_T_RG))
            {
                NV_V(gradT, =, C_T_RG(c0, t0));
                real lambdaMixture = C_K_L(c0, t0);
                real conductHeatFlux = -lambdaMixture * NV_DOT(gradT, A) / Amag;
                F_PROFILE(f, t, i) = reactHeatFlux + conductHeatFlux;
            }
            else
            {
                F_PROFILE(f, t, i) = reactHeatFlux;
            }
        }
    }
    end_f_loop(f, t)
#endif
}

/**
 * Source term for species (O2) transport equation.
 * @param oxygenSource          UDF name.
 * @param c                     Index that identifies cell on which the source term is to be applied.
 * @param t                     Pointer to cell thread.
 * @param dS                    Not used.
 * @param eqn                   Equation index.
 * @return source term for O2 transport equation (kg.m-3.s-1).
 */
DEFINE_SOURCE(oxygenSource, c, t, dS, eqn)
{
    // Rate of reaction R4f (forward).
    real reactRateR4f = getReactRate(c, t, &reactionR4f);

    // Partial pressures of O2, NO2, NO in R4 (Pa).
    real pressure[3];
    getPressureArray(c, t, pressure);

    // Stoichiometric coefficients of each species in R4.
    real stoichCoeffs[] = {-0.5, 1.0, -1.0};

    // Equilibrium constant of reaction R4.
    real equilConst = getEquilConst(3, pressure, stoichCoeffs);

    // Rate of reaction R4r (reverse).
    real reactRateR4r = getReactRate(c, t, &reactionR4r) / equilConst;

    // Net rate of production/consumption of O2 due to reaction R4.
    real massRateO2 = -0.5 * (reactRateR4f - reactRateR4r) * MOLAR_MASS_O2;

    // Return source term for O2 transport equation (kg.m-3.s-1).
    real source = massRateO2;
    return source;
}

/**
 * Source term for species (NO2) transport equation.
 * @param nitrogenDioxideSource    UDF name.
 * @param c                        Index that identifies cell on which the source term is to be applied.
 * @param t                        Pointer to cell thread.
 * @param dS                       Not used.
 * @param eqn                      Equation index.
 * @return source term for NO2 transport equation (kg.m-3.s-1).
 */
DEFINE_SOURCE(nitrogenDioxideSource, c, t, dS, eqn)
{
    // Rate of reaction R4f (forward).
    real reactRateR4f = getReactRate(c, t, &reactionR4f);

    // Partial pressures of O2, NO2, NO in R4 (Pa).
    real pressure[3];
    getPressureArray(c, t, pressure);

    // Stoichiometric coefficients of each species in R4.
    real stoichCoeffs[] = {-0.5, 1.0, -1.0};

    // Equilibrium constant of reaction R4.
    real equilConst = getEquilConst(3, pressure, stoichCoeffs);

    // Rate of reaction R4r (reverse).
    real reactRateR4r = getReactRate(c, t, &reactionR4r) / equilConst;

    // Net rate of production/consumption of NO2 due to reaction R4.
    real massRateNO2 = (reactRateR4f - reactRateR4r) * MOLAR_MASS_NO2;

    // Return source term for NO2 transport equation (kg.m-3.s-1).
    real source = massRateNO2;
    return source;
}

/**
 * Source term for species (NO) transport equation.
 * @param nitricOxideSource        UDF name.
 * @param c                        Index that identifies cell on which the source term is to be applied.
 * @param t                        Pointer to cell thread.
 * @param dS                       Not used.
 * @param eqn                      Equation index.
 * @return source term for NO transport equation (kg.m-3.s-1).
 */
DEFINE_SOURCE(nitricOxideSource, c, t, dS, eqn)
{
    // Rate of reaction R4f (forward).
    real reactRateR4f = getReactRate(c, t, &reactionR4f);

    // Partial pressures of O2, NO2, NO in R4 (Pa).
    real pressure[3];
    getPressureArray(c, t, pressure);

    // Stoichiometric coefficients of each species in R4.
    real stoichCoeffs[] = {-0.5, 1.0, -1.0};

    // Equilibrium constant of reaction R4.
    real equilConst = getEquilConst(3, pressure, stoichCoeffs);

    // Rate of reaction R4r (reverse).
    real reactRateR4r = getReactRate(c, t, &reactionR4r) / equilConst;

    // Net rate of production/consumption of NO due to reaction R4.
    real massRateNO = -1.0 * (reactRateR4f - reactRateR4r) * MOLAR_MASS_NO;

    // Return source term for NO transport equation (kg.m-3.s-1).
    real source = massRateNO;
    return source;
}

/**
 * Source term for energy equation due to reaction R4.
 * @param heatReleaseFromR4        UDF name.
 * @param c                        Index that identifies cell on which the source term is to be applied.
 * @param t                        Pointer to cell thread.
 * @param dS                       Not used.
 * @param eqn                      Equation index.
 */
DEFINE_SOURCE(heatReleaseFromR4, c, t, dS, eqn)
{
     // Rate of reaction R4f (forward).
    real reactRateR4f = getReactRate(c, t, &reactionR4f);

    // Partial pressures of O2, NO2, NO in R4 (Pa).
    real pressure[3];
    getPressureArray(c, t, pressure);

    // Stoichiometric coefficients of each species in R4.
    real stoichCoeffs[] = {-0.5, 1.0, -1.0};

    // Equilibrium constant of reaction R4.
    real equilConst = getEquilConst(3, pressure, stoichCoeffs);

    // Rate of reaction R4r (reverse).
    real reactRateR4r = getReactRate(c, t, &reactionR4r) / equilConst;

    // Overall rate of reaction R4.
    real netReactRateR4 = reactRateR4f - reactRateR4r;

    // Sensible enthalpy change of reaction R4 (J.mol-1).
    real tempMixture = C_T(c, t);
    real enthalpyChangeR4 = getSensibleEnthalpy(tempMixture, REF_TEMP, &NO2info) -
                            getSensibleEnthalpy(tempMixture, REF_TEMP, &NOinfo) - 
                            0.5 * getSensibleEnthalpy(tempMixture, REF_TEMP, &O2info);

    // Standard enthalpy change of reaction R4 (J.mol-1).
    real stdEnthalpyR4 = DELTA_H_R4 + enthalpyChangeR4;

    // Heat release due to reaction R4 (W.m-3).
    real heatReleaseR4 = -netReactRateR4 * stdEnthalpyR4;

    // Return source term for energy equation due to reaction R4 (W.m-3).
    real source = heatReleaseR4;
    return source;
}