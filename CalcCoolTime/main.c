/*
 * CalcCoolTime Routine for IDL
 * dnelson nov.2012
 * based on cooling/cooling.c (Arepo r21365)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
 
#include "idl_export.h"
#include "main.h"

static double Tmin = 0.0;     /**< min temperature in log10 */ 
static double Tmax = 9.0;     /**< max temperature in log10 */ 
static double deltaT;         /**< log10 of temperature spacing in the interpolation tables */ 

static GasState gs;           /**< gas state */
static RateTable *RateT;      /**< tabulated rates */
static PhotoTable *PhotoTUVB; /**< photo-ionization/heating rate table for UV background*/
static PhotoCurrent pc;       /**< current interpolated photo rates */
static int NheattabUVB;       /**< length of UVB photo table */

/* IDL usage:
 *
 *  npts = long(n_elements(u))
 *  u = float(u) ; code units
 *  rho = float(rho) ; code units
 *  ne  = float(ne)
 *  scalefac = float(time) ; Header.Time
 *  coolrate_out = fltarr(npts)
 *
 *  ret = Call_External('/n/home07/dnelson/idl/CalcCoolTime/CalcCoolTime.so', 'CalcCoolTime', $
                        npts,u,rho,ne,cooltime_out,scalefac,/CDECL)
 */

int CalcCoolTime(int argc, void* argv[])
{
  int npts;
  float *u, *rho, *ne;
  float *cooltime_out;
  float scalefac;
  char buf[128];
  
  /* --- IDL glue --- */
  
  if (argc != 6)
  {
    sprintf(buf,"Wrong number of arguments (%d)!\n",argc);
    IDL_Message(IDL_M_GENERIC,IDL_MSG_RET,buf);
    return 0;
  } else {
    //IDL_Message(IDL_M_GENERIC,IDL_MSG_RET,"CalcCoolTime Loaded.");
  }
  
  npts = *(int *)argv[0];
  u    = (float *)argv[1];
  rho  = (float *)argv[2];
  ne   = (float *)argv[3];
  
  cooltime_out = (float *)argv[4];
  scalefac     = *(float *)argv[5];
  
  /* --- init --- */
  
  // if scalefac>1, this is a flag that u actually contains log10(temp) already
  int flag = 0;
  if(scalefac > 1.0) {
    flag = 1;
    scalefac -= 1; // true scalefac is scalefac-1
    //IDL_Message(IDL_M_GENERIC,IDL_MSG_RET,"Using direct log10(temp).");
  }
  
  gs.XH = HYDROGEN_MASSFRAC; /* set default hydrogen mass fraction */

  memset(&pc, 0, sizeof(PhotoCurrent)); /* zero photo-ionization/heating rates */

  RateT = (RateTable*) malloc( (NCOOLTAB + 1) * sizeof(RateTable) );
  MakeRateTable(); /* allocate and construct rate table */

  if (!ReadIonizeParams("/n/home07/dnelson/idl/CalcCoolTime/TREECOOL_fg_dec11", 0)) /* read photo tables */
    return 0;

  IonizeParamsUVB(scalefac); /* set UV background ionization state for this redshift */
  
  /* --- cooling rate --- */
  
  double temp = 0.0;
  double ratefact, LambdaNet;
  int i;
  
  for ( i = 0; i < npts; i++)
  {
    //rho *= All.UnitDensity_in_cgs * HUBBLE_PARAM * HUBBLE_PARAM; // handle in IDL
    //u_old *= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs; // handle in IDL

    gs.nHcgs = gs.XH * rho[i] / PROTONMASS;	/* hydrogen number dens in cgs units */
    ratefact = gs.nHcgs * gs.nHcgs / rho[i];

    if(flag)
    {
      // convert u to temp (if scalefac>1 this a flag that u is actually already in log10(temp))
      temp = u[i];
      // calculate thermal energy for u/udot
      u[i] = pow(10.0,temp) * (1 + ne[i] + gs.yhelium) / (1 + 4 * gs.yhelium) / gs.mhboltz / GAMMA_MINUS1;
    }
    else
      temp = log10(GAMMA_MINUS1 * u[i] * gs.mhboltz * (1 + 4 * gs.yhelium) / (1 + ne[i] + gs.yhelium));
  
    // calculate lambda
    LambdaNet = CoolingRate(temp, rho[i], ne[i], scalefac);

    if(LambdaNet >= 0)		/* ups, we have actually heating due to UV background */
      cooltime_out[i] = 0;
    else 
      cooltime_out[i] = u[i] / (-ratefact * LambdaNet);
    
    //cooltime_out[i] *= HUBBLE_PARAM / All.UnitTime_in_s; // handle in IDL
  }
  
  return 0;
}

/** \brief Computes the actual abundance ratios.
 * 
 *  The chemical composition of the gas is primordial (no metals are present)
 * 
 *  \param logT     log10 of gas temperature
 *  \param rho      gas density
 *  \param ne_guess electron number density relative to hydrogen number density
 */
void find_abundances_and_rates(double logT, double rho, double ne_input)
{
  double neold;
  int j, niter;
  double flow, fhi, t;

  double logT_input, rho_input;

  logT_input = logT;
  rho_input = rho;

  if(logT <= Tmin)		/* everything neutral */
    {
      gs.nH0 = 1.0;
      gs.nHe0 = gs.yhelium;
      gs.nHp = 0;
      gs.nHep = 0;
      gs.nHepp = 0;
      gs.ne = 0;
      return;
    }

  if(logT >= Tmax)		/* everything is ionized */
    {
      gs.nH0 = 0;
      gs.nHe0 = 0;
      gs.nHp = 1.0;
      gs.nHep = 0;
      gs.nHepp = gs.yhelium;
      gs.ne = gs.nHp + 2.0 * gs.nHepp;
      return;
    }

  t = (logT - Tmin) / deltaT;
  j = (int) t;
  fhi = t - j;
  flow = 1 - fhi;

  gs.nHcgs = gs.XH * rho / PROTONMASS;	/* hydrogen number dens in cgs units */

  gs.ne = ne_input;
  neold = gs.ne;
  niter = 0;
  gs.necgs = gs.ne * gs.nHcgs;

  /* evaluate number densities (cf KWH eqns 33-38) in units of nH */
  gs.aHp = flow * RateT[j].AlphaHp + fhi * RateT[j + 1].AlphaHp;
  gs.aHep = flow * RateT[j].AlphaHep + fhi * RateT[j + 1].AlphaHep;
  gs.aHepp = flow * RateT[j].AlphaHepp + fhi * RateT[j + 1].AlphaHepp;
  gs.ad = flow * RateT[j].Alphad + fhi * RateT[j + 1].Alphad;
  gs.geH0 = flow * RateT[j].GammaeH0 + fhi * RateT[j + 1].GammaeH0;
  gs.geHe0 = flow * RateT[j].GammaeHe0 + fhi * RateT[j + 1].GammaeHe0;
  gs.geHep = flow * RateT[j].GammaeHep + fhi * RateT[j + 1].GammaeHep;

  if(gs.necgs <= 1.e-25 || pc.J_UV == 0) 
  {
    gs.gJH0ne = gs.gJHe0ne = gs.gJHepne = 0;
  }
  else
  {
    gs.gJH0ne = pc.gJH0 / gs.necgs;
    gs.gJHe0ne = pc.gJHe0 / gs.necgs;
    gs.gJHepne = pc.gJHep / gs.necgs;
  }

  gs.nH0 = gs.aHp / (gs.aHp + gs.geH0 + gs.gJH0ne);	/* eqn (33) */
  gs.nHp = 1.0 - gs.nH0;		/* eqn (34) */

  if((gs.gJHe0ne + gs.geHe0) <= SMALLNUM)	/* no ionization at all */
  {
    gs.nHep = 0.0;
    gs.nHepp = 0.0;
    gs.nHe0 = gs.yhelium;
  }
  else
  {
    gs.nHep = gs.yhelium / (1.0 + (gs.aHep + gs.ad) / (gs.geHe0 + gs.gJHe0ne) + (gs.geHep + gs.gJHepne) / gs.aHepp);	/* eqn (35) */
    gs.nHe0 = gs.nHep * (gs.aHep + gs.ad) / (gs.geHe0 + gs.gJHe0ne);	/* eqn (36) */
    gs.nHepp = gs.nHep * (gs.geHep + gs.gJHepne) / gs.aHepp;	/* eqn (37) */
  }

  gs.bH0 = flow * RateT[j].BetaH0 + fhi * RateT[j + 1].BetaH0;
  gs.bHep = flow * RateT[j].BetaHep + fhi * RateT[j + 1].BetaHep;
  gs.bff = flow * RateT[j].Betaff + fhi * RateT[j + 1].Betaff;
}



/** \brief  Calculate (heating rate-cooling rate)/n_h^2 in cgs units.
 * 
 *  \param logT     log10 of gas temperature
 *  \param rho      gas density
 *  \param nelec    electron number density relative to hydrogen number density
 *  \return         (heating rate-cooling rate)/n_h^2
 */
double CoolingRate(double logT, double rho, double nelec, double scalefac)
{
  double Lambda, Heat;
  double LambdaExc, LambdaIon, LambdaRec, LambdaFF, LambdaCmptn = 0.0;
  double LambdaExcH0, LambdaExcHep, LambdaIonH0, LambdaIonHe0, LambdaIonHep;
  double LambdaRecHp, LambdaRecHep, LambdaRecHepp, LambdaRecHepd;
  double redshift;
  double T;

  if(logT <= Tmin)
    logT = Tmin + 0.5 * deltaT;	/* floor at Tmin */

  gs.nHcgs = gs.XH * rho / PROTONMASS;	/* hydrogen number dens in cgs units */

  if(logT < Tmax)
    {
      find_abundances_and_rates(logT, rho, nelec);

      /* Compute cooling and heating rate (cf KWH Table 1) in units of nH**2 */
      T = pow(10.0, logT);

      LambdaExcH0 = gs.bH0 * gs.ne * gs.nH0;
      LambdaExcHep = gs.bHep * gs.ne * gs.nHep;
      LambdaExc = LambdaExcH0 + LambdaExcHep;	/* excitation */
      LambdaIonH0 = 2.18e-11 * gs.geH0 * gs.ne * gs.nH0;
      LambdaIonHe0 = 3.94e-11 * gs.geHe0 * gs.ne * gs.nHe0;
      LambdaIonHep = 8.72e-11 * gs.geHep * gs.ne * gs.nHep;
      LambdaIon = LambdaIonH0 + LambdaIonHe0 + LambdaIonHep;	/* ionization */
      LambdaRecHp = 1.036e-16 * T * gs.ne * (gs.aHp * gs.nHp);
      LambdaRecHep = 1.036e-16 * T * gs.ne * (gs.aHep * gs.nHep);
      LambdaRecHepp = 1.036e-16 * T * gs.ne * (gs.aHepp * gs.nHepp);
      LambdaRecHepd = 6.526e-11 * gs.ad * gs.ne * gs.nHep;
      LambdaRec = LambdaRecHp + LambdaRecHep + LambdaRecHepp + LambdaRecHepd;
      LambdaFF = gs.bff * (gs.nHp + gs.nHep + 4 * gs.nHepp) * gs.ne;
      Lambda = LambdaExc + LambdaIon + LambdaRec + LambdaFF;

      if(scalefac)
	{
	  redshift = 1 / scalefac - 1;
	  LambdaCmptn = 5.65e-36 * gs.ne * (T - 2.73 * (1. + redshift)) * pow(1. + redshift, 4.) / gs.nHcgs;

	  Lambda += LambdaCmptn;
	}
      else
	  LambdaCmptn = 0;

      Heat = 0;
      if(pc.J_UV != 0)
	  Heat += (gs.nH0 * pc.epsH0 + gs.nHe0 * pc.epsHe0 + gs.nHep * pc.epsHep) / gs.nHcgs;

    }
  else				/* here we're outside of tabulated rates, T>Tmax K */
    {
      /* at high T (fully ionized); only free-free and Compton cooling are present. Assumes no heating. */

      Heat = 0;

      LambdaExcH0 = LambdaExcHep = LambdaIonH0 = LambdaIonHe0 = LambdaIonHep = LambdaRecHp = LambdaRecHep = LambdaRecHepp = LambdaRecHepd = 0;

      /* very hot: H and He both fully ionized */
      gs.nHp = 1.0;
      gs.nHep = 0;
      gs.nHepp = gs.yhelium;
      gs.ne = gs.nHp + 2.0 * gs.nHepp;

      T = pow(10.0, logT);
      LambdaFF = 1.42e-27 * sqrt(T) * (1.1 + 0.34 * exp(-(5.5 - logT) * (5.5 - logT) / 3)) * (gs.nHp + 4 * gs.nHepp) * gs.ne;

      if(scalefac)
	{
	  redshift = 1 / scalefac - 1;
	  /* add inverse Compton cooling off the microwave background */
	  LambdaCmptn = 5.65e-36 * gs.ne * (T - 2.73 * (1. + redshift)) * pow(1. + redshift, 4.) / gs.nHcgs;
	}
      else
	  LambdaCmptn = 0;

      Lambda = LambdaFF + LambdaCmptn;
    }

  return (Heat - Lambda);
}


/** \brief Make cooling rates interpolation table.
 *
 *  Set up interpolation tables in T for cooling rates given in KWH, ApJS, 105, 19
 *  Some coefficients are updated for GFM_PRIMORDIAL_RATES
 */
void MakeRateTable(void)
{
  int i;
  double T;
  double Tfact;

  gs.yhelium = (1 - gs.XH) / (4 * gs.XH);
  gs.mhboltz = PROTONMASS / BOLTZMANN;
  Tmin = 1.0;
  deltaT = (Tmax - Tmin) / NCOOLTAB;
  gs.ethmin = pow(10.0, Tmin) * (1. + gs.yhelium) / ((1. + 4. * gs.yhelium) * gs.mhboltz * GAMMA_MINUS1);
  /* minimum internal energy for neutral gas */

  for(i = 0; i <= NCOOLTAB; i++)
    {
      RateT[i].BetaH0 = RateT[i].BetaHep = RateT[i].Betaff = RateT[i].AlphaHp = RateT[i].AlphaHep = RateT[i].AlphaHepp = RateT[i].Alphad = RateT[i].GammaeH0 = RateT[i].GammaeHe0 = RateT[i].GammaeHep = 0;

      T = pow(10.0, Tmin + deltaT * i);
      Tfact = 1.0 / (1 + sqrt(T / 1.0e5));

      /* collisional excitation */
      /* Cen 1992 */
      if(118348 / T < 70)
	RateT[i].BetaH0 = 7.5e-19 * exp(-118348 / T) * Tfact;
      if(473638 / T < 70)
	RateT[i].BetaHep = 5.54e-17 * pow(T, -0.397) * exp(-473638 / T) * Tfact;

      /* free-free */
      RateT[i].Betaff = 1.43e-27 * sqrt(T) * (1.1 + 0.34 * exp(-(5.5 - log10(T)) * (5.5 - log10(T)) / 3));

      /* recombination */

      /* Cen 1992 */
      /* Hydrogen II */
      RateT[i].AlphaHp = 8.4e-11 * pow(T / 1000, -0.2) / (1. + pow(T / 1.0e6, 0.7)) / sqrt(T);
      /* Helium II */
      RateT[i].AlphaHep = 1.5e-10 * pow(T, -0.6353);
      /* Helium III */
      RateT[i].AlphaHepp = 4. * RateT[i].AlphaHp;

      /* Cen 1992 */
      /* dielectric recombination */
      if(470000 / T < 70)
	RateT[i].Alphad = 1.9e-3 * pow(T, -1.5) * exp(-470000 / T) * (1. + 0.3 * exp(-94000 / T));

      /* collisional ionization */

      /* Cen 1992 */
      /* Hydrogen */
      if(157809.1 / T < 70)
	RateT[i].GammaeH0 = 5.85e-11 * sqrt(T) * exp(-157809.1 / T) * Tfact;
      /* Helium */
      if(285335.4 / T < 70)
	RateT[i].GammaeHe0 = 2.38e-11 * sqrt(T) * exp(-285335.4 / T) * Tfact;
      /* Hellium II */
      if(631515.0 / T < 70)
	RateT[i].GammaeHep = 5.68e-12 * sqrt(T) * exp(-631515.0 / T) * Tfact;

    }
}

/** \brief Read table input for ionizing parameters.
 * 
 *  \param file that contains the tabulated parameters
 *  \param which flag used to identify the type of the ionizing backgroung (0 = UV background,
 *               1 = AGN background)
 */
int ReadIonizeParams(char *fname, int which)
{
  int iter, i;
  FILE *fdcool;
  float dummy;

  if (which == 0)
  {
  NheattabUVB = 0;

  for (iter = 0, i = 0; iter < 2; iter++)
    {
      if(!(fdcool = fopen(fname, "r"))) {
        printf("ERROR: Cannot read ionization table in file `%s'\n", fname);
        return 0;
      }
      if (iter == 0)
        while (fscanf(fdcool, "%g %g %g %g %g %g %g", &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy) != EOF)
          NheattabUVB++;
      if (iter == 1)
        while (fscanf(fdcool, "%g %g %g %g %g %g %g",&PhotoTUVB[i].variable,&PhotoTUVB[i].gH0,&PhotoTUVB[i].gHe,&PhotoTUVB[i].gHep,&PhotoTUVB[i].eH0,&PhotoTUVB[i].eHe,&PhotoTUVB[i].eHep) != EOF)
          i++;
      fclose(fdcool);

      if (iter == 0)
        {
          PhotoTUVB = (PhotoTable *) malloc( NheattabUVB * sizeof(PhotoTable) );
        }
    }
    /* ignore zeros at end of treecool file */
    for (i = 0; i < NheattabUVB; ++i)
      if(PhotoTUVB[i].gH0 == 0.0)
        break;

    NheattabUVB = i;
    //printf("COOLING: using %d ionization table entries from file `%s'.\n", NheattabUVB, fname);
    
    return 1;
  }

  return 0;
}


/** \brief Set the ionization parameters for the UV background.
 */
void IonizeParamsUVB(double scalefac)
{
  int i, ilow;
  double logz, dzlow, dzhi;
  double redshift;

  if(scalefac)
    redshift = 1 / scalefac - 1;
  else
    {
      memset(&pc, 0, sizeof(PhotoCurrent)); // SetZeroIonization
      return;
    }

  logz = log10(redshift + 1.0);
  ilow = 0;
  for(i = 0; i < NheattabUVB; i++)
    {
      if(PhotoTUVB[i].variable < logz)
	ilow = i;
      else
	break;
    }

  dzlow = logz - PhotoTUVB[ilow].variable;
  dzhi = PhotoTUVB[ilow +1].variable - logz;

  if(logz > PhotoTUVB[NheattabUVB -1].variable || PhotoTUVB[ilow].gH0 == 0 || PhotoTUVB[ilow + 1].gH0 == 0 || NheattabUVB == 0)
    {
      memset(&pc, 0, sizeof(PhotoCurrent)); // SetZeroIonization
      return;
    }
  else
    pc.J_UV = 1; 

  pc.gJH0 = pow(10., (dzhi * log10(PhotoTUVB[ilow].gH0) + dzlow * log10(PhotoTUVB[ilow + 1].gH0)) / (dzlow + dzhi));
  pc.gJHe0 = pow(10., (dzhi * log10(PhotoTUVB[ilow].gHe) + dzlow * log10(PhotoTUVB[ilow + 1].gHe)) / (dzlow + dzhi));
  pc.gJHep = pow(10., (dzhi * log10(PhotoTUVB[ilow].gHep) + dzlow * log10(PhotoTUVB[ilow + 1].gHep)) / (dzlow + dzhi));
  pc.epsH0 = pow(10., (dzhi * log10(PhotoTUVB[ilow].eH0) + dzlow * log10(PhotoTUVB[ilow + 1].eH0)) / (dzlow + dzhi));
  pc.epsHe0 = pow(10., (dzhi * log10(PhotoTUVB[ilow].eHe) + dzlow * log10(PhotoTUVB[ilow + 1].eHe)) / (dzlow + dzhi));
  pc.epsHep = pow(10., (dzhi * log10(PhotoTUVB[ilow].eHep) + dzlow * log10(PhotoTUVB[ilow + 1].eHep)) / (dzlow + dzhi));

  return;
}
