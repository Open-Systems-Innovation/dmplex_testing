#include "petscdmlabel.h"
#include "petscsys.h"
#include "petscsystypes.h"
static char help[] =
    "2D FVM elastic wave problem with custom Riemann Solver and Slope Limiter\n";

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Include statements 
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#include <petscdmplex.h>
#include <petscoptions.h>
#include <petscviewerhdf5.h>
#include "petscds.h"
#include "petscts.h"
#include "petscdmplex.h"   

#define DIM 2 /* Geometric dimension */
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Physics Context 
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
typedef struct {
  PetscReal sigma_xx;
  PetscReal sigma_yy;
  PetscReal sigma_xy;
  PetscReal velocity[DIM];
  PetscReal lambda;
  PetscReal mu;
  PetscReal density;
} ElasticNode;

typedef struct {
  ElasticNode elasticnode;
  PetscReal vals[DIM +1];
} ElasticNodeUnion;

typedef struct _n_Physics *Physics;

struct _n_Physics {
  void (*riemann)(PetscInt, PetscInt, const PetscReal[], const PetscReal[], const PetscScalar[], const PetscScalar[], PetscInt, const PetscScalar[], PetscScalar[], void *);
  PetscReal maxspeed; /* kludge to pick initial time step, need to add
                         monitoring and step control */
  PetscInt dim; // dim
  PetscInt Nf;  // number of fields
  PetscReal blob_lambda;
  PetscReal blob_mu;
  PetscReal blob_density;
  PetscReal silicone_lambda;
  PetscReal silicone_mu;
  PetscReal silicone_density;
};

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   User Context 
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
typedef struct _n_User *User;

struct _n_User {
  Physics   physics;
  PetscReal maxspeed;
  PetscInt  monitorStepOffset;
  char      outputBasename[PETSC_MAX_PATH_LEN]; /* Basename for output files */
  PetscInt  vtkInterval;                        /* For monitor */
  PetscBool vtkmon;
};

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Riemann Solver 
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
//static void PhysicsRiemann_Elastic_HLL(PetscInt dim, PetscInt Nf, const PetscReal *qp, const PetscReal *n, const PetscScalar *xL, const PetscScalar *xR, PetscInt numConstants, const PetscScalar constants[], PetscScalar *flux, Physics phys)
//{
//  Physics_Elastic *elastic = (Physics_Elastic *)phys->data;
//  PetscReal        cL, cR;
//  PetscReal        nn[DIM];
//  const ElasticNode *uL = (const ElasticNode *)xL, *uR = (const ElasticNode *)xR;
//  ElasticNodeUnion fL, fR;
//  PetscInt         i;
//  PetscReal        zero = 0.;
//  PetscErrorCode   ierr;
//
//  /* Check for valid density to prevent division by zero */
//  if (uL->rho <= 0 || uR->rho <= 0) {
//    for (i = 0; i < dim + 2; i++) flux[i] = zero;
//    return;
//  }
//
//  /* Normalize the normal vector */
//  nn[0] = n[0];
//  nn[1] = n[1];
//  Normalize2Real(nn);
//
//  /* Calculate fluxes for each side */
//  ElasticFlux(phys, nn, uL, &fL.elasticnode);
//  ElasticFlux(phys, nn, uR, &fR.elasticnode);
//
//  /* Compute wave speeds based on Lamé parameters */
//  PetscReal cp_L = PetscSqrtReal((elastic->lambda + 2 * elastic->mu) / uL->rho);
//  PetscReal cs_L = PetscSqrtReal(elastic->mu / uL->rho);
//  PetscReal cp_R = PetscSqrtReal((elastic->lambda + 2 * elastic->mu) / uR->rho);
//  PetscReal cs_R = PetscSqrtReal(elastic->mu / uR->rho);
//
//  /* Max wave speeds on each side */
//  cL = PetscMax(cp_L, cs_L);
//  cR = PetscMax(cp_R, cs_R);
//
//  /* Projected velocities along the normal */
//  PetscReal v_L = Dot2Real(uL->velocity, nn);
//  PetscReal v_R = Dot2Real(uR->velocity, nn);
//
//  /* Calculate sL and sR for HLL */
//  PetscReal sL = PetscMin(v_L - cL, v_R - cR);
//  PetscReal sR = PetscMax(v_L + cL, v_R + cR);
//
//  /* HLL Riemann solver logic */
//  if (sL > zero) {
//    for (i = 0; i < dim + 2; i++) flux[i] = fL.vals[i] * Norm2Real(n);
//  } else if (sR < zero) {
//    for (i = 0; i < dim + 2; i++) flux[i] = fR.vals[i] * Norm2Real(n);
//  } else {
//    for (i = 0; i < dim + 2; i++) {
//      flux[i] = ((sR * fL.vals[i] - sL * fR.vals[i] + sR * sL * (xR[i] - xL[i])) / (sR - sL)) * Norm2Real(n);
//    }
//  }
//}

static inline PetscReal Dot2Real(const PetscReal *x, const PetscReal *y)
{
  return x[0] * y[0] + x[1] * y[1];
}

static inline PetscReal Norm2Real(const PetscReal *x)
{
  return PetscSqrtReal(PetscAbsReal(Dot2Real(x, x)));
}

static inline void Normalize2Real(PetscReal *x)
{
  PetscReal a = 1. / Norm2Real(x);
  x[0] *= a;
  x[1] *= a;
}

static inline PetscReal DotDIMReal(const PetscReal *x, const PetscReal *y)
{
  PetscInt  i;
  PetscReal prod = 0.0;

  for (i = 0; i < DIM; i++)
    prod += x[i] * y[i];
  return prod;
}

static PetscErrorCode Flux(Physics phys, const PetscReal *n, const ElasticNode *x, ElasticNode *f)
{
  PetscFunctionBeginUser;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* Compute elastic wave speeds: P-wave speed (cp), S-wave speed (cs), and max speed (c_max) */
PetscErrorCode ElasticWaveSpeed(PetscReal lambda, PetscReal mu, PetscReal rho, PetscReal *cp, PetscReal *cs)
{
  PetscFunctionBeginUser;

  if (rho <= 0) {
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Density (rho) must be positive.");
  }
  if (mu < 0 || lambda < 0) {
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Lamé parameters (lambda and mu) must be non-negative.");
  }

  /* Calculate P-wave (c_p) and S-wave (c_s) speeds */
  PetscReal bulk = lambda + 2.0 * mu;
  *cp = PetscSqrtReal(bulk / rho);
  *cs = PetscSqrtReal(mu / rho);

  PetscFunctionReturn(PETSC_SUCCESS);
}

//static void NewElasticRiemann(
//    PetscInt dim, PetscInt Nf, const PetscReal *qp, const PetscReal *n,
//    const PetscScalar *xL, const PetscScalar *xR, PetscInt numConstants,
//    const PetscScalar constants[], PetscScalar *flux, void *ctx)
//{
//  /* Input Parameters:
//       dim          - The spatial dimension
//       Nf           - The number of fields
//       x            - The coordinates at a point on the interface
//       n            - The normal vector to the interface
//       uL           - The state vector to the left of the interface
//       uR           - The state vector to the right of the interface
//       numConstants - number of constant parameters
//       constants    - constant parameters
//       ctx          - optional user context
//       flux         - output array of flux through the interface
//  */
//  //PetscPrintf(PETSC_COMM_WORLD, "dim = %d, Nf = % d, qp = (%f, %f), n = (%f, %f), \n uL =  (%f, %f, %f, %f, %f), \n uR = (%f, %f, %f, %f, %f), numConstants = %d\n", dim, Nf, qp[0], qp[1], n[0], n[1], uL[0], uL[1], uL[2], uL[3], uL[4], uR[0], uR[1], uR[2], uR[3], uR[4], numConstants );
//
//  Physics phys = (Physics) ctx;  // Cast context to your struct
//  
//  PetscReal       cpL, cpR, csL, csR, nn[DIM], s2;
//  PetscReal       zero = 0.0;
//  PetscInt        i;
//  PetscErrorCode  ierr;
//
//  PetscReal du[5];
//  PetscReal detP, detS;
//  PetscReal a1, a2, a3, a4; // alpha values
//  
//  // Compute jumps in the states (uL is left, uR is right)
//  // and store in du
//
//  /* Normalize the normal vector */
//  for (i = 0, s2 = 0.; i < DIM; i++) {
//    nn[i] = n[i];
//    s2 += nn[i] * nn[i];
//  }
//  s2 = PetscSqrtReal(s2); /* |n|_2 = sum(n^2)^1/2 */
//  for (i = 0.; i < DIM; i++)
//    nn[i] /= s2;
//
//  // find the difference of each component of u
//  for (PetscInt i = 0; i < 5; ++i) {
//      du[i] = xR[i] - xL[i];
//  }
//
//  const ElasticNode *uL = (const ElasticNode *)xL;
//  const ElasticNode *uR = (const ElasticNode *)xR;
//  ElasticNodeUnion   fL, fR;
//
//  PetscReal lambdaR = 4;
//  PetscReal muR = 1;
//  PetscReal rhoR = 1;
//  PetscReal lambdaL= 4;
//  PetscReal muL= 1;
//  PetscReal rhoL = 1;
//
//  ierr = Flux(phys, nn, uL, &fL.elasticnode);
//  if (ierr) {
//    PetscCallVoid(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
//    for (i = 0; i < 2 + dim; i++) fL.vals[i] = zero / zero;
//    PetscCallVoid(PetscFPTrapPop());
//  }
//
//  ierr = Flux(phys, nn, uR, &fR.elasticnode);
//  if (ierr) {
//    PetscCallVoid(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
//    for (i = 0; i < 2 + dim; i++) fR.vals[i] = zero / zero;
//    PetscCallVoid(PetscFPTrapPop());
//  }
//  /* Calculate wave speeds */
//  ierr = ElasticWaveSpeed(lambdaL, muL, rhoL, &cpL, &csL);
//  if (ierr) {
//    PetscCallVoid(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
//    cpL = zero / zero;
//    csL = zero / zero;
//    PetscCallVoid(PetscFPTrapPop());
//  }
//  ierr = ElasticWaveSpeed(lambdaR, muR, rhoR, &cpR, &csR);
//  if (ierr) {
//    PetscCallVoid(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
//    cpR = zero / zero;
//    csR = zero / zero;
//    PetscCallVoid(PetscFPTrapPop());
//  }
//  PetscReal velL  = DotDIMReal(uL->velocity, nn);
//  PetscReal velR  = DotDIMReal(uR->velocity, nn);
//  //PetscReal speed = PetscMax(velR + cR, velL + cL);
//  PetscReal speed = PetscMax(velR, velL);
//
//  for (i = 0; i < 2 + dim; i++)
//    flux[i] =
//      0.5 * ((fL.vals[i] + fR.vals[i]) + speed * (xL[i] - xR[i])) * s2;
//
//}
static void ElasticRiemann(
    PetscInt dim, PetscInt Nf, const PetscReal *qp, const PetscReal *n,
    const PetscScalar *uL, const PetscScalar *uR, PetscInt numConstants,
    const PetscScalar constants[], PetscScalar *flux, void *ctx)
{
    // Normalize the normal vector
    PetscReal nn[dim];
    nn[0] = n[0];
    nn[1] = n[1];
    Normalize2Real(nn);

    // Obtain material parameters
    PetscReal lambdaL = uL[5], muL = uL[6], densityL = uL[7];
    PetscReal lambdaR = uR[5], muR = uR[6], densityR = uR[7];

    // Compute P-wave and S-wave speeds
    PetscReal bulkL = lambdaL + 2.0 * muL;
    PetscReal bulkR = lambdaR + 2.0 * muR;
    PetscReal cpL = PetscSqrtReal(bulkL / densityL);
    PetscReal csL = PetscSqrtReal(muL / densityL);
    PetscReal cpR = PetscSqrtReal(bulkR / densityR);
    PetscReal csR = PetscSqrtReal(muR / densityR);

    // Estimate wave speeds
    PetscReal SL = PetscMin(-cpL, -csL);
    PetscReal SR = PetscMax(cpR, csR);

    // Calculate fluxes for left and right states
    PetscReal FL[5], FR[5];
    FL[0] = uL[3] * nn[0] + uL[4] * nn[1]; // flux component example
    FL[1] = uL[4] * nn[0] + uL[3] * nn[1];
    FL[2] = 0.0; // update this based on your problem setup
    FL[3] = lambdaL * (uL[0] * nn[0] + uL[1] * nn[1]);
    FL[4] = muL * (uL[0] * nn[1] - uL[1] * nn[0]);

    FR[0] = uR[3] * nn[0] + uR[4] * nn[1];
    FR[1] = uR[4] * nn[0] + uR[3] * nn[1];
    FR[2] = 0.0;
    FR[3] = lambdaR * (uR[0] * nn[0] + uR[1] * nn[1]);
    FR[4] = muR * (uR[0] * nn[1] - uR[1] * nn[0]);

    // Compute the HLL flux
    for (int i = 0; i < 5; ++i) {
        if (SL >= 0) {
            flux[i] = FL[i]; // Left state flux
        } else if (SR <= 0) {
            flux[i] = FR[i]; // Right state flux
        } else {
            // HLL flux
            flux[i] = (SR * FL[i] - SL * FR[i] + SL * SR * (uR[i] - uL[i])) / (SR - SL);
        }
    }
    flux[5] = 0.0;
    flux[6] = 0.0;
    flux[7] = 0.0;
}


//static void ElasticRiemann(
//    PetscInt dim, PetscInt Nf, const PetscReal *qp, const PetscReal *n,
//    const PetscScalar *uL, const PetscScalar *uR, PetscInt numConstants,
//    const PetscScalar constants[], PetscScalar *flux, void *ctx)
//{
//  /* Input Parameters:
//       dim          - The spatial dimension
//       Nf           - The number of fields
//       x            - The coordinates at a point on the interface
//       n            - The normal vector to the interface
//       uL           - The state vector to the left of the interface
//       uR           - The state vector to the right of the interface
//       numConstants - number of constant parameters
//       constants    - constant parameters
//       ctx          - optional user context
//       flux         - output array of flux through the interface
//  */
//  //PetscPrintf(PETSC_COMM_WORLD, "dim = %d, Nf = % d, qp = (%f, %f), n = (%f, %f), \n uL =  (%f, %f, %f, %f, %f), \n uR = (%f, %f, %f, %f, %f), numConstants = %d\n", dim, Nf, qp[0], qp[1], n[0], n[1], uL[0], uL[1], uL[2], uL[3], uL[4], uR[0], uR[1], uR[2], uR[3], uR[4], numConstants );
//  //PetscPrintf(PETSC_COMM_WORLD, "dim = %d, Nf = %d\n", dim, Nf);
//  //PetscPrintf(PETSC_COMM_WORLD, "qp = (%f, %f), n = (%f, %f)\n", qp[0], qp[1], n[0], n[1]);
//  //PetscPrintf(PETSC_COMM_WORLD, "uL = (%f, %f, %f, %f, %f)\n", uL[0], uL[1], uL[2], uL[3], uL[4]);
//  //PetscPrintf(PETSC_COMM_WORLD, "uR = (%f, %f, %f, %f, %f)\n", uR[0], uR[1], uR[2], uR[3], uR[4]);
//  //PetscPrintf(PETSC_COMM_WORLD, "numConstants = %d\n", numConstants);
//  Physics phys = (Physics) ctx;  // Cast context to your struct
//  
//  PetscReal du[5];
//  PetscReal detP, detS;
//  PetscReal nn[dim];
//  PetscReal a1, a2, a3, a4; // alpha values
//  
//  // Compute jumps in the states (uL is left, uR is right)
//  // and store in du
//
//  /* Normalize the normal vector */
//  nn[0] = n[0];
//  nn[1] = n[1];
//  Normalize2Real(nn);
//
//  for (PetscInt i = 0; i < 5; ++i) {
//      du[i] = uR[i] - uL[i];
//  }
//
//  // obtain material paramter values for L and R sides
//  PetscReal lambdaL = uL[5];
//  PetscReal muL = uL[6];
//  PetscReal densityL = uL[7];
//  PetscReal lambdaR = uR[5];
//  PetscReal muR = uR[6];
//  PetscReal densityR = uR[7];
//
//  //PetscPrintf(PETSC_COMM_WORLD, "LambdaR = %f\n", lambdaR);
//  //PetscPrintf(PETSC_COMM_WORLD, "qp[0]= %f, qp[1] = %f\n", qp[0], qp[1]);
//  //if (lambdaR < 0.0) {
//    //PetscPrintf(PETSC_COMM_WORLD, "LambdaR = %f\n", lambdaR);
//  //  lambdaR = 2;
//  //  muR = 1;
//  //  densityR = 1;
//  //PetscPrintf(PETSC_COMM_WORLD, "qp[0]= %f, qp[1] = %f\n", qp[0], qp[1]);
//  //PetscPrintf(PETSC_COMM_WORLD, "n[0]= %f, n[1] = %f\n", n[0], n[1]);
// // }
// // if (lambdaR < 2.0) {
// //   PetscPrintf(PETSC_COMM_WORLD, "qp[0]= %f, qp[1] = %f\n", qp[0], qp[1]);
// // }
//  //if (lambdaR != lambdaL) {
//  //  PetscPrintf(PETSC_COMM_WORLD, "lambdaR = %f,  muR = %f,  densityR = %f\n", lambdaR, muR, densityR);
//  //  PetscPrintf(PETSC_COMM_WORLD, "lambdaL = %f,  muL = %f,  densityL = %f\n", lambdaL, muL, densityL);
//  //}
//  // Calculate cp (P-wave speed) and cs (S-wave speed)
//  PetscReal bulkL = lambdaL + 2.0 * muL;
//  PetscReal bulkR = lambdaR + 2.0 * muR;
//  PetscReal cpL = PetscSqrtReal(bulkL / densityL);
//  PetscReal csL = PetscSqrtReal(muL / densityL);
//  PetscReal cpR = PetscSqrtReal(bulkR / densityR);
//  PetscReal csR = PetscSqrtReal(muR / densityR);
//  /* Max wave speeds on each side */
//  PetscReal cL = PetscMax(cpL, csL);
//  PetscReal cR = PetscMax(cpR, csR);
//  PetscReal maxspeed = (cL > cR) ? cL : cR;
//  //if (maxspeed > 2.0) {
//  // PetscPrintf(PETSC_COMM_WORLD, "Max wavespeed: %f\n", (cL > cR) ? cL : cR);
//  //}
//
//  //PetscPrintf(PETSC_COMM_WORLD, "cpL = %f,  csL = %f,  cpR = %f,  csR = %f\n", cpL, csL, cpR, csR);
//  // Calculate useful multiples of the norm vector components
//  PetscReal nx = nn[0];
//  PetscReal ny = nn[1];
//  PetscReal nx2 = nx * nx;
//  PetscReal ny2 = ny * ny;
//  PetscReal nxy = nx * ny;
//  
////  PetscPrintf(PETSC_COMM_WORLD, "nx = %f   ny = %f  nx2 = %f  ny2 = %f   nxy = %f \n", nx, ny, nx2, ny2, nxy);
//  // Define Eigenvalues (wave speeds)
//  PetscReal s[4] = {-cpL, cpR, -csL, csR};
//
//  // Define the 4 eigenvectors (from columns of Matrix R)
//  PetscReal r1[5] = {lambdaL + 2 * muL * nx2, lambdaL + 2 * muL * ny2,
//                     2 * muL * nxy, nx * cpL, ny * cpL};
//  PetscReal r2[5] = {lambdaR + 2 * muR * nx2, lambdaR + 2 * muR * ny2,
//                     2 * muR * nxy, -nx * cpR, -ny * cpR};
//  PetscReal r3[5] = {-2 * muL * nxy, 2 * muL * nxy, muL * (nx2 - ny2),
//                     -ny * csL, nx * csL};
//  PetscReal r4[5] = {-2 * muR * nxy, 2 * muR * nxy, muR * (nx2 - ny2), ny * csR,
//                     -nx * csR};
//
//  // Compute the 4 alphas
//  detP = cpR * bulkL + cpL * bulkR;
//  detS = csR * muL + csL * muR;
//  //PetscPrintf(PETSC_COMM_WORLD, "detP = %f  detS = %f \n", detP, detS);
//  // P wave strengths
//  a1 = (cpR * (du[0] * nx2 + du[1] * ny2 + 2 * nxy * du[2]) +
//         bulkR * (nx * du[3] + ny * du[4])) / detP;
//  a2 = (cpL * (du[0] * nx2 + du[1] * ny2 + 2 * nxy * du[2]) -
//         bulkL * (nx * du[3] + ny * du[4])) / detP;
//  // S wave strengths
//  a3 = (csR * (du[2] * (nx2 - ny2) + nxy * (du[1] - du[0])) +
//         muR * (nx * du[4] - ny * du[3])) / detS;
//  a4 = (csL * (du[2] * (nx2 - ny2) + nxy * (du[1] - du[0])) -
//         muL * (nx * du[4] - ny * du[3])) / detS;
//  if (a1 > 0.0) {
//    //PetscPrintf(PETSC_COMM_WORLD, "a1 = %f,  a2 = %f,  a3 = %f,  a4 = %f\n", a1,a2,a3,a4); 
//    //PetscPrintf(PETSC_COMM_WORLD, "du[0] = %f,  du[1] = %f,  du[2] = %f,  du[3] = %f,  du[4] = %f\n", du[0],du[1],du[2],du[3],du[4]); 
//  }
//  // Compute the waves
//// Compute the waves
//  PetscReal W1[5], W2[5], W3[5], W4[5];
//  for (int i = 0; i < 5; ++i) {
//      W1[i] = a1 * r1[i];
//      W2[i] = a2 * r2[i];
//      W3[i] = a3 * r3[i];
//      W4[i] = a4 * r4[i];
//  }
//  // First wave (P-wave, left-going)
//  // Second wave (P-wave, right-going)
//  // Third wave (S-wave, left-going)
//  // Fourth wave (S-wave, right-going)
//  // Use the waves to update the flux
//  flux[0] = s[1] * W2[0] + s[3] * W4[0] + s[0] * W1[0] + s[2] * W3[0];
//  flux[1] = s[1] * W2[1] + s[3] * W4[1] + s[0] * W1[1] + s[2] * W3[1];
//  flux[2] = s[1] * W2[2] + s[3] * W4[2] + s[0] * W1[2] + s[2] * W3[2];
//  flux[3] = s[1] * W2[3] + s[3] * W4[3] + s[0] * W1[3] + s[2] * W3[3];
//  flux[4] = s[1] * W2[4] + s[3] * W4[4] + s[0] * W1[4] + s[2] * W3[4];
//  flux[5] = 0.0;
//  flux[6] = 0.0;
//  flux[7] = 0.0;
//  
////  PetscPrintf(PETSC_COMM_WORLD, "flux = (%f, %f, %f, %f, %f, %f, %f, %f)\n", flux[0], flux[1], flux[2], flux[3], flux[4], flux[5], flux[6], flux[7]);
//}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Boundary Conditions
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static PetscErrorCode LeftBoundary(PetscReal time, const PetscReal *c, const PetscReal *n, const PetscScalar *xI, PetscScalar *xG, void *ctx)
{
  PetscFunctionBeginUser;
  Physics         phys   = (Physics)ctx;

 // PetscPrintf(PETSC_COMM_WORLD, "time= %f, c = (%f, %f), n = (%f, %f) \n xI = (%f, %f, %f, %f, %f, %f, %f, %f),  xG = (%f, %f, %f, %f, %f, %f, %f, %f)\n", time, c[0], c[1], n[0], n[1], xI[0], xI[1], xI[2], xI[3], xI[4], xI[5], xI[6], xI[7], xG[0], xG[1], xG[2], xG[3], xG[4], xG[5], xG[6], xG[7]);
  xG[0] =  xI[0];
  xG[1] =  xI[1];
  xG[2] =  xI[2];
  xG[3] =  xI[3];
  xG[4] =  xI[4];
  xG[5] =  phys->silicone_lambda; 
  xG[6] =  phys->silicone_mu;
  xG[7] =  phys->silicone_density;
//  if (time < 0.025) {
//    xG[0] =  PetscSignReal((0.1 * PETSC_PI * time)/0.025); 
// }

//  PetscPrintf(PETSC_COMM_WORLD, "time= %f, c = (%f, %f), n = (%f, %f) \n  xG = (%f, %f, %f, %f, %f, %f, %f, %f)\n", time, c[0], c[1], n[0], n[1], xG[0], xG[1], xG[2], xG[3], xG[4], xG[5], xG[6], xG[7]);
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode BoundaryOutflow(PetscReal time, const PetscReal *c, const PetscReal *n, const PetscScalar *xI, PetscScalar *xG, void *ctx)
{
  PetscFunctionBeginUser;
  Physics         phys   = (Physics)ctx;

 // PetscPrintf(PETSC_COMM_WORLD, "time= %f, c = (%f, %f), n = (%f, %f) \n xI = (%f, %f, %f, %f, %f, %f, %f, %f),  xG = (%f, %f, %f, %f, %f, %f, %f, %f)\n", time, c[0], c[1], n[0], n[1], xI[0], xI[1], xI[2], xI[3], xI[4], xI[5], xI[6], xI[7], xG[0], xG[1], xG[2], xG[3], xG[4], xG[5], xG[6], xG[7]);
  xG[0] =  xI[0];
  xG[1] =  xI[1];
  xG[2] =  xI[2];
  xG[3] =  xI[3];
  xG[4] =  xI[4];
  xG[5] =  phys->silicone_lambda; 
  xG[6] =  phys->silicone_mu;
  xG[7] =  phys->silicone_density;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SetUpBC(DM dm, PetscDS ds, Physics phys)
{
  DMLabel        label;
  PetscInt       boundaryid;  // Physical group label for boundary
  
  PetscFunctionBeginUser;

  // Add absorbing boundary conditions
  // Check if the label exists
  PetscCall(DMGetLabel(dm, "left_boundary", &label));
  // if it doesn't exist, then throw error
  if (!label) {
      PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Error: Label 'left_boundary' not found\n"));
      PetscFunctionReturn(PETSC_ERR_ARG_WRONG);
  }

  // add the boundary condition to the left 
  boundaryid = 4;
  PetscCall(PetscDSAddBoundary(ds, DM_BC_NATURAL_RIEMANN, "left_boundary", label, 1, &boundaryid, 0, 0, NULL, (void (*)(void))LeftBoundary, NULL, phys, NULL));

  // add the boundary condition to the bottom 
  boundaryid = 1;
  PetscCall(DMGetLabel(dm, "bottom_boundary", &label));
  PetscCall(PetscDSAddBoundary(ds, DM_BC_NATURAL_RIEMANN, "bottom_boundary", label, 1, &boundaryid, 0, 0, NULL, (void (*)(void))BoundaryOutflow, NULL, phys, NULL));

  // add the boundary condition to the right
  boundaryid = 2;
  PetscCall(DMGetLabel(dm, "right_boundary", &label));
  PetscCall(PetscDSAddBoundary(ds, DM_BC_NATURAL_RIEMANN, "right_boundary", label, 1, &boundaryid, 0, 0, NULL, (void (*)(void))BoundaryOutflow, NULL, phys, NULL));

  // add the boundary condition to the top
  boundaryid = 3;
  PetscCall(DMGetLabel(dm, "top_boundary", &label));
  PetscCall(PetscDSAddBoundary(ds, DM_BC_NATURAL_RIEMANN, "top_boundary", label, 1, &boundaryid, 0, 0, NULL, (void (*)(void))BoundaryOutflow, NULL, phys, NULL));

  //PetscCall(DMGetLabel(dm, "Face Sets", &label));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Initial Conditions
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static PetscErrorCode blobMaterial (PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u, void *ctx)
{
  Physics         phys   = (Physics)ctx;

  PetscFunctionBeginUser;

  u[0] = 0;
  u[1] = 0;
  u[2] = 0;
  u[3] = 0;
  u[4] = 0;
  u[5] = phys->blob_lambda;
  u[6] = phys->blob_mu;
  u[7] = phys->blob_density;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode siliconeMaterial (PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u, void *ctx)
{
  Physics         phys   = (Physics)ctx;
  PetscFunctionBeginUser;

  u[0] = 0;
  u[1] = 0;
  u[2] = 0;
  u[3] = 0;
  u[4] = 0;
  u[5] = phys->silicone_lambda;
  u[6] = phys->silicone_mu;
  u[7] = phys->silicone_density;
  PetscFunctionReturn(PETSC_SUCCESS);
}
static PetscErrorCode recieverMaterial (PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u, void *ctx)
{
  Physics         phys   = (Physics)ctx;

  PetscFunctionBeginUser;
  u[0] = 0;
  u[1] = 0;
  u[2] = 0;
  u[3] = 1;
  u[4] = 0;
  u[5] = phys->silicone_lambda;
  u[6] = phys->silicone_mu;
  u[7] = phys->silicone_density;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode SetUpMaterialProperties(DM dm, Vec X, User user)
{
  DM                 plex, dmCell;
  //Model              mod = user->model;
  Vec                cellgeom;
  const PetscScalar *cgeom;
  PetscScalar       *x;
  PetscInt           cStart, cEnd, c;
  DMLabel            siliconeLabel, outerLayerLabel, recieverLabel;

  PetscFunctionBeginUser;
  PetscCall(DMConvert(dm, DMPLEX, &plex));
  PetscCall(DMPlexGetGeometryFVM(plex, NULL, &cellgeom, NULL));
  PetscCall(VecGetDM(cellgeom, &dmCell));
  PetscCall(DMPlexGetSimplexOrBoxCells(dm, 0, &cStart, &cEnd));
  PetscCall(VecGetArrayRead(cellgeom, &cgeom));
  PetscCall(VecGetArray(X, &x));

  // Get the labels for each material
  PetscCall(DMGetLabel(plex, "silicone", &siliconeLabel));
  PetscCall(DMGetLabel(plex, "outer_layer", &outerLayerLabel));
  PetscCall(DMGetLabel(plex, "reciever", &recieverLabel));

  for (c = cStart; c < cEnd; ++c) {
    const PetscFVCellGeom *cg;
    PetscScalar           *xc;
    PetscInt               siliconeValue, outerLayerValue, recieverValue;

    // Read the geometry data for the current cell
    PetscCall(DMPlexPointLocalRead(dmCell, c, cgeom, &cg));
    PetscCall(DMPlexPointGlobalRef(dm, c, x, &xc));
    if (!xc) continue; // Skip if the cell is not owned

    // Check the material label for the current cell
//    PetscCall(DMGetLabelValue(siliconeLabel, c, &siliconeValue));
    PetscCall(DMGetLabelValue(dm, "silicone", c, &siliconeValue));
    PetscCall(DMGetLabelValue(dm, "outer_layer", c, &outerLayerValue));
    PetscCall(DMGetLabelValue(dm, "reciever", c, &recieverValue));

   // PetscPrintf(PETSC_COMM_WORLD, "Xc = %f, %f, %f, %f\n", xc[0], xc[1], xc[2], xc[3]);
 //   // Set initial condition based on the material type

      xc[0] = 0.0;
      xc[1] = 0.0;
      xc[2] = 0.0;
      xc[3] = 0.0;
    if (siliconeValue == 9) {
 //     PetscCall(siliconeMaterial(dim, 0.0, cg->centroid, xc, mod->solutionctx));
      xc[5] = 2;
      xc[6] = 1;
      xc[7] = 1;
    } else if (outerLayerValue == 10) {
      xc[5] = 2;
      xc[6] = 1;
      xc[7] = 1;
    } else if (recieverValue == 11) {
      xc[0] = 1;
      xc[5] = 2;
      xc[6] = 1;
      xc[7] = 1;
    } else {
      xc[5] = 2;
      xc[6] = 1;
      xc[7] = 1;
    }
  }

  PetscCall(VecRestoreArrayRead(cellgeom, &cgeom));
  PetscCall(VecRestoreArray(X, &x));
  PetscCall(DMDestroy(&plex));
  PetscFunctionReturn(PETSC_SUCCESS);
}


//PetscErrorCode SetUpMaterialProperties(DM dm, Vec X, Physics phys)
//{
//  DMLabel label;
//  PetscInt materialid; 
//
//  // declare an array of function pointers called func
//  // one function per field
//  // we have to do it like this becaues DMProjectFunction expects an
//  // array of function pointers because C is like that
//  PetscErrorCode (*func[1])(PetscInt dim, PetscReal time, const PetscReal x[],
//                            PetscInt Nf, PetscScalar *u, void *ctx);
//  // now func[0] is a pointer to a function
//  void *ctx[1];  // array of pointers for DMProjectFunctionLabel
//  ctx[0] = phys;
//
//  PetscFunctionBeginUser;
//
//  // blob material
////  materialid = 6;
////  func[0] = blobMaterial;
////  ctx[0] = (void *)phys;
////  PetscCall(DMGetLabel(dm, "blob", &label));
////  PetscCall(DMProjectFunctionLabel(dm, 0.0, label, 1, &materialid, 0, NULL, func, ctx, ADD_VALUES, X));
//
//  // silicone material
//  materialid = 9;
//  func[0] = siliconeMaterial;
//  PetscCall(DMGetLabel(dm, "silicone", &label));
//  PetscCall(DMProjectFunctionLabel(dm, 0.0, label, 1, &materialid, 0, NULL, func, ctx, INSERT_VALUES, X));
//
//  materialid = 10;
//  func[0] = siliconeMaterial;
//  PetscCall(DMGetLabel(dm, "outer_layer", &label));
//  PetscCall(DMProjectFunctionLabel(dm, 0.0, label, 1, &materialid, 0, NULL, func, ctx, INSERT_VALUES, X));
//  
//  materialid = 11;
//  func[0] = recieverMaterial;
//  PetscCall(DMGetLabel(dm, "reciever", &label));
//  PetscCall(DMProjectFunctionLabel(dm, 0.0, label, 1, &materialid, 0, NULL, func, ctx, INSERT_VALUES, X));
//
//
//
//  // boundary
// // materialid = 1;
// // PetscCall(DMGetLabel(dm, "Face Sets", &label));
// // DMLabelView(label, PETSC_VIEWER_STDOUT_WORLD);
// // PetscCall(DMProjectFunctionLabel(dm, 0.0, label, 1, &materialid, 0, NULL, func, ctx, ADD_VALUES, X));
// // PetscCall(VecViewFromOptions(X, NULL, "-X_view"));
//
//  PetscFunctionReturn(PETSC_SUCCESS);
//}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Initial Conditions 
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
PetscErrorCode zero_vector(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u, void *ctx)
{
  // Initialize all values to zero initially
  for (PetscInt d = 0; d < 5; ++d) {
    u[d] = 0.0;
  }
  return 0;
}

PetscErrorCode SetInitialConditions(DM dm, Vec X, User user)
{
  PetscFunctionBeginUser;
  PetscErrorCode (*funcs[1])(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx);

  funcs[0] = zero_vector;
  PetscCall(DMProjectFunction(dm, 0.0, funcs, NULL, INSERT_ALL_VALUES, X));
  
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Output to VTK to view in Paraview 
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static PetscErrorCode OutputVTK(DM dm, const char *filename, PetscViewer *viewer)
{
  PetscFunctionBeginUser;
  PetscCall(PetscViewerCreate(PetscObjectComm((PetscObject)dm), viewer));
  PetscCall(PetscViewerSetType(*viewer, PETSCVIEWERVTK));
  PetscCall(PetscViewerFileSetName(*viewer, filename));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MonitorVTK(TS ts, PetscInt stepnum, PetscReal time, Vec X, void *ctx)
{
  User        user = (User)ctx;
  DM          dm;
  PetscViewer viewer;
  char        filename[PETSC_MAX_PATH_LEN];
  PetscReal   xnorm;
  PetscBool   rollback;

  PetscFunctionBeginUser;

  // Check for rollback
  PetscCall(TSGetStepRollBack(ts, &rollback));
  if (rollback)
   PetscFunctionReturn(PETSC_SUCCESS);

  // Get the current solution
  PetscCall(PetscObjectSetName((PetscObject)X, "u"));
  PetscCall(VecGetDM(X, &dm));

  // Find the norm of the solution vector for summary printing
  PetscCall(VecNorm(X, NORM_INFINITY, &xnorm));

  // Adjust step iteration number by user offset
  if (stepnum >= 0) stepnum += user->monitorStepOffset;

  // Output to VTK at regular intervals or at the final time
  if (user->vtkInterval < 1) PetscFunctionReturn(PETSC_SUCCESS);
  // the final time is represented by stepnum = -1
  if ((stepnum == -1) ^ (stepnum % user->vtkInterval == 0)) {
    if (stepnum == -1) {
      PetscCall(TSGetStepNumber(ts, &stepnum));  // Adjust for final time
    }

    // Generate the VTK filename and open the VTK viewer
    PetscCall(PetscSNPrintf(filename, sizeof(filename), "%s-%03" PetscInt_FMT ".vtu", user->outputBasename, stepnum));
    PetscCall(OutputVTK(dm, filename, &viewer));

    // Write the solution vector to VTK using VecView
    PetscCall(VecView(X, viewer));

    // Clean up
    PetscCall(PetscViewerDestroy(&viewer));
  }

  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MonitorHDF5(TS ts, PetscInt stepnum, PetscReal time,
                                  Vec X, void *ctx) {
  PetscViewer viewer;
  User user = (User)ctx;
  DM                dm;

  PetscFunctionBeginUser;

  // open mesh.h5 for HDF5 input/output
  PetscCall(PetscViewerHDF5Open(PETSC_COMM_WORLD, "mesh.h5", FILE_MODE_WRITE, &viewer));

  //PetscCall(PetscViewerHDF5PushGroup(viewer, "/big_dm"));
  PetscCall(TSGetDM(ts, &dm));
  PetscCall(DMSetOutputSequenceNumber(dm, stepnum, time));
  PetscCall(DMView(dm, viewer));
  PetscCall(VecView(X, viewer));
  //PetscCall(PetscViewerHDF5PopGroup(viewer));

  //PetscCall(PetscViewerHDF5PopTimestepping(viewer));
  //PetscCall(PetscViewerHDF5PopGroup(viewer));

  PetscCall(PetscViewerDestroy(&viewer));
  
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Initialize Time stepping (TS) object
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static PetscErrorCode InitializeTS(DM dm, User user, TS *ts)
{
  PetscFunctionBeginUser;
  PetscCall(TSCreate(PetscObjectComm((PetscObject)dm), ts));
  PetscCall(TSSetType(*ts, TSSSP)); // use Runge-Kutta, -ts_ssp_type {rks2,rks3,rk104}
  PetscCall(TSSetDM(*ts, dm));
  if (user->vtkmon) PetscCall(TSMonitorSet(*ts, MonitorVTK, user, NULL));
  PetscCall(DMTSSetBoundaryLocal(dm, DMPlexTSComputeBoundary, user));
  PetscCall(DMTSSetRHSFunctionLocal(dm, DMPlexTSComputeRHSFunctionFVM, user));
  PetscCall(TSSetMaxTime(*ts, 0.4));
  PetscCall(TSSetExactFinalTime(*ts, TS_EXACTFINALTIME_STEPOVER));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Process User Options 
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Main program
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
int main(int argc, char **argv)
{
  // Declare variables
  MPI_Comm          comm;
  PetscDS           ds;
  PetscFV           fv;
  User              user;
  Physics           phys;
  DM                dm, plex;
  PetscReal         ftime, cfl, dt, minRadius, maxspeed;
  PetscInt          nsteps;
  PetscInt          dim = 2;
  PetscInt          numComponents = 8;
  TS                ts;
  TSConvergedReason reason;
  Vec               X;
  PetscViewer       viewer;
  //PetscViewer       viewer;

  // Initialize PETSc code
  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &argv, NULL, help));

  comm = PETSC_COMM_WORLD;
      
  PetscCall(PetscNew(&user));
  PetscCall(PetscNew(&user->physics));
  phys              = user->physics;
  phys->riemann     = (PetscRiemannFunc)ElasticRiemann;
  phys->dim = 2;  // dim
  phys->Nf = 1;   // number of fields
  user->vtkmon      = PETSC_TRUE;
  user->vtkInterval = 1;
  PetscStrcpy(user->outputBasename, "paraview/paraview");
  phys->blob_lambda = 200;
  phys->blob_mu = 100;
  phys->blob_density = 1;
  phys->silicone_lambda = 2;
  phys->silicone_mu = 1;
  phys->silicone_density = 1;

  
  // Read in 3D mesh from file
  //PetscCall(DMPlexCreateFromFile(comm, user->infile, NULL, PETSC_TRUE,&dm));
  PetscCall(DMPlexCreateGmshFromFile(comm, "mesh.msh", PETSC_TRUE, &dm));
  PetscCall(DMSetFromOptions(dm));
  PetscCall(DMViewFromOptions(dm, NULL, "-orig_dm_view"));
  PetscCall(DMGetDimension(dm, &dim));

  // add ghost points
  {
    DM gdm;
    PetscInt numGhost;
    PetscCall(DMPlexConstructGhostCells(dm, NULL, &numGhost, &gdm));

    PetscPrintf(comm, "number of ghost points  = %d\n", numGhost);
    PetscCall(DMDestroy(&dm));
    dm = gdm;
  }

  // Create finite volume
  PetscCall(PetscFVCreate(comm, &fv));
  PetscCall(PetscFVSetFromOptions(fv));
  PetscCall(PetscFVSetNumComponents(fv, numComponents));
  PetscCall(PetscFVSetSpatialDimension(fv, dim));
  PetscCall(PetscObjectSetName((PetscObject)fv, "finite_volume"));
  
  // Define component names for pressure and velocity
  {
    // Set names for the components of the field 
    PetscCall(PetscFVSetComponentName(fv, 0, "Stress_11")); // 11 component of stress 
    PetscCall(PetscFVSetComponentName(fv, 1, "Stress_22")); // 22 component of stress 
    PetscCall(PetscFVSetComponentName(fv, 2, "Stress_12")); // 12 component of stress 
    PetscCall(PetscFVSetComponentName(fv, 3, "Velocity_x")); // X component of velocity
    PetscCall(PetscFVSetComponentName(fv, 4, "Velocity_y")); // Y component of velocity
    PetscCall(PetscFVSetComponentName(fv, 5, "Lambda")); // Y component of velocity
    PetscCall(PetscFVSetComponentName(fv, 6, "Mu")); // Y component of velocity
    PetscCall(PetscFVSetComponentName(fv, 7, "Density")); // Y component of velocity
  }
  PetscCall(PetscFVViewFromOptions(fv, NULL, "-fv_view"));
  
  // add FV to DM
  PetscCall(DMAddField(dm, NULL, (PetscObject)fv));

  // Create the Discrete Systems (DS)
  // each DS object has a set of fields with a PetscVW discretization
  PetscCall(DMCreateDS(dm));
  PetscCall(DMGetDS(dm, &ds));

  PetscCall(DMViewFromOptions(dm, NULL, "-dm_view"));

  // set the pointwise functions. The actual equations to be enforced
  // on each region
  PetscCall(PetscDSSetRiemannSolver(ds, 0, user->physics->riemann));
  PetscCall(PetscDSSetContext(ds, 0, user->physics));
  PetscCall(SetUpBC(dm, ds, user->physics));
  PetscCall(PetscDSSetFromOptions(ds));

  // initialize TS object
  PetscCall(InitializeTS(dm, user, &ts));

  // create solution vector
  PetscCall(DMCreateGlobalVector(dm, &X));
  PetscCall(PetscObjectSetName((PetscObject)X, "solution"));
  PetscCall(SetInitialConditions(dm, X, user));
  PetscCall(SetUpMaterialProperties(dm, X, phys));

 // MAYBE DO THIS NEXT TO VIEW THE DM IN PARAVIEW  
//  PetscCall(DMConvert(dm, DMPLEX, &plex));
//  if (user->vtkmon) {
//    DM  dmCell;
//    Vec cellgeom, partition;
// 
//    PetscCall(DMPlexGetGeometryFVM(plex, NULL, &cellgeom, NULL));
//    PetscCall(OutputVTK(dm, "ex11-cellgeom.vtk", &viewer));
//    PetscCall(VecView(cellgeom, viewer));
//    PetscCall(PetscViewerDestroy(&viewer));
//    
// //   PetscCall(CreatePartitionVec(dm, &dmCell, &partition));
//    //PetscCall(OutputVTK(dmCell, "ex11-partition.vtk", &viewer));
//    //PetscCall(VecView(partition, viewer));
//    //PetscCall(PetscViewerDestroy(&viewer));
//    //PetscCall(VecDestroy(&partition));
//    //PetscCall(DMDestroy(&dmCell));
//   }
//PetscCall(MPIU_Allreduce(&phys->maxspeed, &mod->maxspeed, 1, MPIU_REAL, MPIU_MAX, PetscObjectComm((PetscObject)ts)));

  // CFL condition for time stepping
  //cfl = 0.9*4;
  PetscCall(DMPlexGetGeometryFVM(dm, NULL, NULL, &minRadius));
  PetscReal silicone_speed = PetscSqrtReal((phys->silicone_lambda + 2.0 * phys->silicone_mu) / phys->silicone_density);
  PetscReal blob_speed = PetscSqrtReal((phys->blob_lambda + 2.0 * phys->blob_mu) / phys->blob_density);
  maxspeed = (silicone_speed > blob_speed) ? silicone_speed : blob_speed;
  dt = (minRadius / maxspeed)/5;
  PetscPrintf(comm, "dt = %f\n", dt);
  PetscPrintf(comm, "maxspeed = %f\n", maxspeed);
  PetscPrintf(comm, "minRadius = %f\n", minRadius);
  PetscCall(TSSetTimeStep(ts, dt));
  PetscCall(TSSetFromOptions(ts));

  PetscCall(TSSetSolution(ts, X));
  PetscCall(VecViewFromOptions(X, NULL, "-X_view"));
  PetscCall(VecDestroy(&X));
  PetscCall(TSViewFromOptions(ts, NULL, "-ts_view"));
  PetscCall(PetscDSViewFromOptions(ds, NULL, "-ds_view"));
  PetscCall(TSSolve(ts, NULL));

  PetscCall(TSGetSolveTime(ts, &ftime));
  PetscCall(TSGetStepNumber(ts, &nsteps));
  PetscCall(TSGetConvergedReason(ts, &reason));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "%s at time %g after %" PetscInt_FMT " steps\n", TSConvergedReasons[reason], (double)ftime, nsteps));

  // Free objects from memory
  PetscCall(TSDestroy(&ts));
  PetscCall(DMDestroy(&dm));
  PetscCall(PetscFVDestroy(&fv));
  
  // End main program
  PetscFinalize();
  return 0;
}

