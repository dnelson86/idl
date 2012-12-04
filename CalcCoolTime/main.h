/*
 * CalcCoolTime Routine for IDL
 * dnelson nov.2012
 * based on cooling/cooling.c (Arepo r21365)
 */
 
#define NCOOLTAB          2000
#define SMALLNUM          1.0e-60
#define MAX_TABLESIZE     250          /* Max # of lines in TREECOOL */
#define HYDROGEN_MASSFRAC 0.76
#define PROTONMASS        1.67262178e-24
#define BOLTZMANN         1.38065e-16
#define GAMMA             (5.0/3)  /**< adiabatic index of simulated gas */
#define GAMMA_MINUS1      (GAMMA-1)
#define HUBBLE_PARAM      0.7

/* proto */
double CoolingRate(double logT, double rho, double nelec, double scalefac);
void   find_abundances_and_rates(double logT, double rho, double ne_input);
void   IonizeParamsUVB(double scalefac);
void   MakeRateTable(void);
int    ReadIonizeParams(char *fname, int which);

/* data for gas state */
typedef struct
{
 double ne, necgs, nHcgs;
 double bH0, bHep, bff, aHp, aHep, aHepp, ad, geH0, geHe0, geHep;
 double gJH0ne, gJHe0ne, gJHepne;
 double nH0, nHp, nHep, nHe0, nHepp;
 double XH, yhelium;
 double mhboltz;
 double ethmin;   /* minimum internal energy for neutral gas */
} GasState;

/* tabulated rates */
typedef struct
{
 double BetaH0, BetaHep, Betaff;
 double AlphaHp, AlphaHep, Alphad, AlphaHepp;
 double GammaeH0, GammaeHe0, GammaeHep;
} RateTable;

/* photo-ionization/heating rate table */
typedef struct
{
  float variable;       /* logz for UVB */
  float gH0, gHe, gHep; /* photo-ionization rates */
  float eH0, eHe, eHep; /* photo-heating rates */
} PhotoTable;

/* current interpolated photo-ionization/heating rates */
typedef struct
{
 char J_UV;
 double gJH0, gJHep, gJHe0, epsH0, epsHep, epsHe0;
} PhotoCurrent;
