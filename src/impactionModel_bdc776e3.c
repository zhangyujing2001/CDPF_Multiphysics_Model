/*-----------------------*- Fluent UDF -*-----------------------*\
  @File name        particleTransport.c
  @Author           Yujing Zhang
  @Email            zhangyujing@tongji.edu.cn
  @Created on       2025-11-13
  @Last modified    2025-11-18
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
#define PI          3.141592       // Approximate value of Pi
#define RHO_1       1770           // Particle density, kg.m-3
#define D           0.4e-9         // Minimum separation, 0.4nm
#define D_0         0.165e-9       // "Cut-off" separation, 0.165nm
#define GAMMA_1     0.095          // Surface energy of particle (graphite), J.m-2
#define GAMMA_2     2.561          // Surface energy of substrate (gamma-alumina), J.m-2
#define UDM_IMPACT  0              // 'UDM_0' for recording the number of impacts
#define UDM_DEPOSIT 1              // 'UDM_1' for recording the number of deposits
#define UDM_NUMBER  2              // Number of UDMs used

/**
 * Function to calculate the Hamaker constant.
 * @param gamma1   Surface energy of particle.
 * @param gamma2   Surface energy of substrate.
 * @return Hamaker constant.
 */
real getHamaker(real gamma1, real gamma2)
{
    // Effective surface energy
    real gammaEff = sqrt(gamma1 * gamma2);

    // Hamaker constant
    real A = 24 * PI * gammaEff * pow(D_0, 2);
    return A;
}

/**
 * Function to calculate the critical velocity based on B-H theory.
 * @param dp       Particle diameter.
 * @return Critical velocity.
 */
real getCritVelocity(real dp)
{
    // Hamaker constant
    real A = getHamaker(GAMMA_1, GAMMA_2);

    // Critical velocity
    real critVel = sqrt(A / (PI * RHO_1 * D * pow(dp, 2)));
    return critVel;
}

/**
 * Function to compute the restitution coefficient.
 * @param vi       Incident velocity of the particle.
 * @param vcr      Critical velocity.
 * @return Restitution coefficient.
 */
real getRestitutionCoeff(real vi, real vcr)
{
    real e = sqrt(1 - pow((vcr / vi), 2));
    return e;
}

/**
 * `DEFINE_ON_DEMAND` macro to initialize UDMs for particle impact and deposition counts.
 * @param initUDMs    UDF name.
 * @return void.
 */
DEFINE_ON_DEMAND(initUDMs)
{
    // Variables.
    Domain *d = Get_Domain(1);
    Thread *ct;
    cell_t c;

    // Loop over all cell threads in the domain, and initialize UDMs to zero.
    thread_loop_c(ct, d)
    {
        begin_c_loop(c, ct)
        {
            for (int i = 0; i < UDM_NUMBER; i++)
            {
                C_UDMI(c, ct, i) = 0.0;
            }
        }
        end_c_loop(c, ct)
    }
}

/**
 * `DEFINE_DPM_BC` macro to define a boundary condition for particle-wall interactions.
 * The function is executed every time a particle touches a boundary of the domain, 
 * except for symmetric or periodic boundaries.
 * @param impactDepositRebound    UDF name.
 * @param tp                      Pointer to the 'Tracked_Particle' data structure.
 * @param t                       Pointer to the face thread the particle is hitting.
 * @param f                       Index of the face that the particle is hitting.
 * @param f_normal                Array containing the unit vector normal to the face.
 * @param dim                     Dimension of the simulation (2D or 3D).
 * @return the status of the tracked particle.
 */
DEFINE_DPM_BC(impactDepositRebound, tp, t, f, f_normal, dim)
{
    // Variables
    cell_t c = TP_CELL(tp);              // cell index where the particle is located
    Thread *tcell = TP_CELL_THREAD(tp);  // thread pointer for the cell
    real vn = 0.0;                       // normal velocity
    real normal[ND_ND];                  // array to store 'f_normal' vector
    real NV_VEC(x);                      // particle position vector
    int idim = dim;                      // dimension of the simulation
    int thread_id = THREAD_ID(t);        // id of the face thread

    // Store the normal vector
    for (int i = 0; i < idim; i++)
    {
        normal[i] = f_normal[i];
    }

    // Main
    if (TP_TYPE(tp) == DPM_TYPE_INERT)
    {
        // Detect particle impact
        if ((NNULLP(t)) && (THREAD_TYPE(t) == THREAD_F_WALL))
        {
            // Record the coordinates of impact and increment impact count
            F_CENTROID(x, f, t);
            C_UDMI(c, tcell, UDM_IMPACT) += 1.0;
        }
        // Write impact information to file
        FILE *fp0;
        char impactFile[128];
        sprintf(impactFile, "impact_info_%d.csv", thread_id);
        fp0 = fopen(impactFile, "a");
        if (fp0 != NULL)
        {
            fprintf(fp0,
                    "%lld, %g, %g, %g, %g, %g, %g\n",
                    TP_ID(tp),
                    TP_POS(tp)[0], TP_POS(tp)[1], 
                    TP_VEL(tp)[0], TP_VEL(tp)[1], 
                    TP_DIAM(tp), TP_TIME(tp));
            fclose(fp0);
        }
        else
        {
            Message("Warning: Unable to open file %s\n", impactFile);
        }

        // Normal velocity calculation
        for (int i = 0; i < idim; i++)
        {
            vn += TP_VEL(tp)[i] * normal[i];
        }

        // Critical velocity calculation
        real dp = TP_DIAM(tp);
        real vcr = getCritVelocity(dp);

        // Determine particle behavior based on normal velocity and critical velocity
        if (fabs(vn) <= vcr)
        {
            if (NNULLP(t) && (THREAD_TYPE(t) == THREAD_F_WALL))
            {
                C_UDMI(c, tcell, UDM_DEPOSIT) += 1;
            }
            FILE *fp1;
            char depositFile[128];
            sprintf(depositFile, "deposit_info_%d.csv", thread_id);
            fp1 = fopen(depositFile, "a");
            if (fp1 != NULL)
            {
                fprintf(fp1,
                        "%lld, %g, %g, %g, %g, %g, %g\n",
                        TP_ID(tp),
                        TP_POS(tp)[0], TP_POS(tp)[1], 
                        TP_VEL(tp)[0], TP_VEL(tp)[1], 
                        TP_DIAM(tp), TP_TIME(tp));
                fclose(fp1);
            }
            else
            {
                Message("Warning: Unable to open file %s\n", depositFile);
            }
            return PATH_END;  // remove particle from domain after recording its info
        }
        else
        {
            // Compute restitution coefficient
            real vi = fabs(vn);
            real restitution = getRestitutionCoeff(vi, vcr);

            // Compute tangential velocity before rebound
            for (int i = 0; i < idim; i++)
            {
                TP_VEL(tp)[i] -= vn * normal[i];
            }

            // Compute tangential velocity after rebound
            for (int i = 0; i < idim; i++)
            {
                TP_VEL(tp)[i] *= restitution;
            }

            // Compute normal velocity after rebound
            for (int i = 0; i < idim; i++)
            {
                TP_VEL(tp)[i] -= restitution * vn * normal[i];
            }

            // Update particle velocity
            for (int i = 0; i < idim; i++)
            {
                TP_VEL0(tp)[i] = TP_VEL(tp)[i];
            }
            return PATH_ACTIVE;  // continue tracking the particle
        }
    }
    return PATH_ABORT;  // Kill the process if not inert particle
}