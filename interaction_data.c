// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
/** \file interaction_data.c
    Implementation of \ref interaction_data.h "interaction_data.h"
 */
#include <string.h>
#include <stdlib.h>
#include "config.h"
#include "debug.h"
#include "interaction_data.h"
#include "communication.h"
#include "grid.h"
#include "p3m.h"
#include "debye_hueckel.h"
#include "mmm1d.h"
#include "lj.h"
#include "ljcos.h"
#include "gb.h"
#include "parser.h"

/****************************************
 * variables
 *****************************************/
int n_particle_types = 0;
int n_interaction_types = 0;
IA_parameters *ia_params = NULL;

#ifdef ELECTROSTATICS
Coulomb_parameters coulomb = { 0.0, COULOMB_NONE };
Debye_hueckel_params dh_params = { 0.0, 0.0, 0.0, 0.0 };
#endif

int n_bonded_ia = 0;
Bonded_ia_parameters *bonded_ia_params = NULL;

double max_cut;

double lj_force_cap = 0.0;

#ifdef CONSTRAINTS
int n_constraints       = 0;
Constraint *constraints = NULL;
#endif 

/*****************************************
 * functions
 *****************************************/

/** Initialize interaction parameters. */
void initialize_ia_params(IA_parameters *params) {
  params->LJ_eps =
    params->LJ_sig =
    params->LJ_cut =
    params->LJ_shift =
    params->LJ_offset = 
    params->LJ_capradius = 0;
  
  params->LJCOS_eps =
    params->LJCOS_sig =
    params->LJCOS_cut = 
    params->LJCOS_offset = 
    params->LJCOS_alfa = 
    params->LJCOS_beta = 
    params->LJCOS_rmin = 0 ;
    
  params->GB_eps =
    params->GB_sig =
    params->GB_cut =
    params->GB_k1 =
    params->GB_k2 =   
    params->GB_mu = 
    params->GB_nu =
    params->GB_chi1 = 
    params->GB_chi2 = 0 ;

}

/** Copy interaction parameters. */
void copy_ia_params(IA_parameters *dst, IA_parameters *src) {
  dst->LJ_eps = src->LJ_eps;
  dst->LJ_sig = src->LJ_sig;
  dst->LJ_cut = src->LJ_cut;
  dst->LJ_shift = src->LJ_shift;
  dst->LJ_offset = src->LJ_offset;
  dst->LJ_capradius = src->LJ_capradius;

  dst->LJCOS_eps = src->LJCOS_eps;
  dst->LJCOS_sig = src->LJCOS_sig;
  dst->LJCOS_cut = src->LJCOS_cut;
  dst->LJCOS_offset = src->LJCOS_offset;
  dst->LJCOS_alfa = src->LJCOS_alfa;
  dst->LJCOS_beta = src->LJCOS_beta;
  dst->LJCOS_rmin = src->LJCOS_rmin;
  
  dst->GB_eps = src->GB_eps;
  dst->GB_sig = src->GB_sig;
  dst->GB_cut = src->GB_cut;
  dst->GB_k1 = src->GB_k1;
  dst->GB_k2 = src->GB_k2; 
  dst->GB_mu = src->GB_mu;
  dst->GB_nu = src->GB_nu;
  dst->GB_chi1 = src->GB_chi1;
  dst->GB_chi2 = src->GB_chi2; 

}

/** returns non-zero if particles of type i and j have a nonbonded interaction */
int checkIfParticlesInteract(int i, int j) {
  IA_parameters *data = get_ia_param(i, j);

  if (data->LJ_cut != 0)
    return 1;
  
  if (data->LJCOS_cut != 0)
    return 1;

  if (data->GB_cut != 0)
    return 1;

  return 0;
}

char *get_name_of_bonded_ia(int i) {
  switch (i) {
  case BONDED_IA_FENE:
    return "FENE";
  case BONDED_IA_ANGLE:
    return "angle ";
  case BONDED_IA_DIHEDRAL:
    return "dihedral";
  case BONDED_IA_HARMONIC:
    return "HARMONIC";
  default:
    fprintf(stderr, "%d: internal error: name of unknown interaction %d requested\n",
	    this_node, i);
    errexit();
  }
  /* just to keep the compiler happy */
  return "";
}

/** This function increases the LOCAL ia_params field
    to the given size. This function is not exported
    since it does not do this on all nodes. Use
    make_particle_type_exist for that.
*/
void realloc_ia_params(int nsize)
{
  int i, j;
  IA_parameters *new_params;

  if (nsize <= n_particle_types)
    return;

  new_params = (IA_parameters *) malloc(nsize*nsize*sizeof(IA_parameters));
  if (ia_params) {
    /* if there is an old field, copy entries and delete */
    for (i = 0; i < nsize; i++)
      for (j = 0; j < nsize; j++) {
	if ((i < n_particle_types) && (j < n_particle_types))
	  copy_ia_params(&new_params[i*nsize + j],
			 &ia_params[i*n_particle_types + j]);
	else
	  initialize_ia_params(&new_params[i*nsize + j]);
      }
    free(ia_params);
  }
  else {
    /* new field, just init */
    for (i = 0; i < nsize; i++)
      for (j = 0; j < nsize; j++)
	initialize_ia_params(&new_params[i*nsize + j]);
  }

  n_particle_types = nsize;
  ia_params = new_params;
}

int printBondedIAToResult(Tcl_Interp *interp, int i)
{
  Bonded_ia_parameters *params = &bonded_ia_params[i];
  char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];

  sprintf(buffer, "%d ", i);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  
  switch (params->type) {
  case BONDED_IA_FENE:
    Tcl_PrintDouble(interp, params->p.fene.k, buffer);
    Tcl_AppendResult(interp, "FENE ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, params->p.fene.r, buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    return (TCL_OK);
  case BONDED_IA_ANGLE:
    Tcl_PrintDouble(interp, params->p.angle.bend, buffer);
    Tcl_AppendResult(interp, "angle ", buffer, (char *) NULL);
    return (TCL_OK);
  case BONDED_IA_DIHEDRAL:
    Tcl_AppendResult(interp, "dihedral",(char *) NULL);
    return (TCL_OK);
  case BONDED_IA_HARMONIC:
    Tcl_PrintDouble(interp, params->p.harmonic.k, buffer);
    Tcl_AppendResult(interp, "HARMONIC ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, params->p.harmonic.r, buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    return (TCL_OK);
  case BONDED_IA_NONE:
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "unknown bonded interaction number ",buffer,
		     (char *) NULL);
    return (TCL_ERROR);
  default:
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "unknown bonded interaction type",(char *) NULL);
    return (TCL_ERROR);
  }
  return (TCL_ERROR);
}

int printNonbondedIAToResult(Tcl_Interp *interp, int i, int j)
{
  char buffer[TCL_DOUBLE_SPACE + 2*TCL_INTEGER_SPACE];
  IA_parameters *data = get_ia_param(i, j);

  if (!data) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "interaction does not exist",
		     (char *) NULL);
    return (TCL_ERROR);
  }

  sprintf(buffer, "%d %d ", i, j);
  Tcl_AppendResult(interp, buffer, (char *) NULL);
  if (data->LJ_cut != 0) {
    Tcl_PrintDouble(interp, data->LJ_eps, buffer);
    Tcl_AppendResult(interp, "lennard-jones ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->LJ_sig, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->LJ_cut, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->LJ_shift, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->LJ_offset, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->LJ_capradius, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);  
  }
  if (data->LJCOS_cut != 0) {
    Tcl_PrintDouble(interp, data->LJCOS_eps, buffer);
    Tcl_AppendResult(interp, "lj-cos ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->LJCOS_sig, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->LJCOS_cut, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->LJCOS_offset, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->LJCOS_alfa, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->LJCOS_beta, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->LJCOS_rmin, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);  
  }
  if (data->GB_cut != 0) {
    Tcl_PrintDouble(interp, data->GB_eps, buffer);
    Tcl_AppendResult(interp, "gay-berne ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->GB_sig, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->GB_cut, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->GB_k1, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->GB_k2, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->GB_mu, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->GB_nu, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->GB_chi1, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->GB_chi2, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
 
  }
  return (TCL_OK);
}

int printCoulombIAToResult(Tcl_Interp *interp) 
{
#ifdef ELECTROSTATICS
  char buffer[TCL_DOUBLE_SPACE + 2*TCL_INTEGER_SPACE];
  if (coulomb.bjerrum == 0.0) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "coulomb 0.0",
		     (char *) NULL);
    return (TCL_OK);
  }
  Tcl_PrintDouble(interp, coulomb.bjerrum, buffer);
  Tcl_AppendResult(interp, "coulomb ", buffer, " ", (char *) NULL);
  if (coulomb.method == COULOMB_P3M) {
    Tcl_PrintDouble(interp, p3m.r_cut, buffer);
    Tcl_AppendResult(interp, "p3m ", buffer, " ", (char *) NULL);
    sprintf(buffer,"%d",p3m.mesh[0]);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    sprintf(buffer,"%d",p3m.cao);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, p3m.alpha, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, p3m.accuracy, buffer);
    Tcl_AppendResult(interp, buffer, "} ", (char *) NULL);

    Tcl_PrintDouble(interp, p3m.epsilon, buffer);
    Tcl_AppendResult(interp, "{coulomb epsilon ", buffer, " ", (char *) NULL);
    sprintf(buffer,"%d",p3m.inter);
    Tcl_AppendResult(interp, "n_interpol ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, p3m.mesh_off[0], buffer);
    Tcl_AppendResult(interp, "mesh_off ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, p3m.mesh_off[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, p3m.mesh_off[2], buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
  }
  else if (coulomb.method == COULOMB_DH) {
    Tcl_PrintDouble(interp, dh_params.kappa, buffer);
    Tcl_AppendResult(interp, "dh ", buffer, " ",(char *) NULL);
    Tcl_PrintDouble(interp, dh_params.r_cut, buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
  }
  else if (coulomb.method == COULOMB_MMM1D) {
    Tcl_PrintDouble(interp, sqrt(mmm1d_params.far_switch_radius_2), buffer);
    Tcl_AppendResult(interp, "mmm1d ", buffer, " ",(char *) NULL);
    sprintf(buffer, "%d", mmm1d_params.bessel_cutoff);
    Tcl_AppendResult(interp, buffer, " ",(char *) NULL);
    Tcl_PrintDouble(interp, mmm1d_params.maxPWerror, buffer);
    Tcl_AppendResult(interp, buffer,(char *) NULL);
  }
#else
  Tcl_AppendResult(interp, "ELECTROSTATICS not compiled (see config.h)",(char *) NULL);
#endif
  return (TCL_OK);
}

void calc_maximal_cutoff()
{
  int i, j;
  max_cut = -1.0;
  /* bonded */
  for (i = 0; i < n_bonded_ia; i++) {
    switch (bonded_ia_params[i].type) {
    case BONDED_IA_FENE:
      if(max_cut < bonded_ia_params[i].p.fene.r)
	max_cut = bonded_ia_params[i].p.fene.r;
      break;
    case BONDED_IA_HARMONIC:
      if(max_cut < bonded_ia_params[i].p.harmonic.r)
	max_cut = bonded_ia_params[i].p.harmonic.r;
      break;
    default:
      break;
    }
  }
  /* non bonded */
  for (i = 0; i < n_particle_types; i++)
     for (j = i; j < n_particle_types; j++) {
       	if (checkIfParticlesInteract(i, j)) {
	  IA_parameters *data = get_ia_param(i, j);
    if (data->LJ_cut != 0) {
	    if(max_cut < (data->LJ_cut+data->LJ_offset) ) 
	      max_cut = (data->LJ_cut+data->LJ_offset);
	  }
	  if (data->LJCOS_cut != 0) {
	    if(max_cut < (data->LJCOS_cut+data->LJCOS_offset) ) 
	      max_cut = (data->LJCOS_cut+data->LJCOS_offset);
	  }
	  if (data->GB_cut != 0) {
	    if(max_cut < (data->GB_cut) ) 
	      max_cut = (data->GB_cut);
	  }	  
	}
     }
#ifdef ELECTROSTATICS
  /* real space electrostatic */
  switch (coulomb.method) {
  case COULOMB_P3M:
    if (max_cut < p3m.r_cut) 
      max_cut = p3m.r_cut;
    break;
  case COULOMB_DH:
    if (max_cut < dh_params.r_cut) 
      max_cut = dh_params.r_cut;
    break;
  case COULOMB_MMM1D:
    /* uses mi interactions anyways */
    if (max_cut < 0)
      max_cut = 0;
    break;
  }
#endif
}

int inter_print_all(Tcl_Interp *interp)
{
  int i, j, start = 1;

  for (i = 0; i < n_bonded_ia; i++) {
    if (bonded_ia_params[i].type != BONDED_IA_NONE) {
      if (start) {
	Tcl_AppendResult(interp, "{", (char *)NULL);
	start = 0;
      }
      else
	Tcl_AppendResult(interp, " {", (char *)NULL);

      printBondedIAToResult(interp, i);
      Tcl_AppendResult(interp, "}", (char *)NULL);
    }
  }
  for (i = 0; i < n_particle_types; i++)
    for (j = i; j < n_particle_types; j++) {
      if (checkIfParticlesInteract(i, j)) {
	if (start) {
	  Tcl_AppendResult(interp, "{", (char *)NULL);
	  start = 0;
	}
	else
	  Tcl_AppendResult(interp, " {", (char *)NULL);
	printNonbondedIAToResult(interp, i, j);
	Tcl_AppendResult(interp, "}", (char *)NULL);
      }
    }
#ifdef ELECTROSTATICS
  if(coulomb.bjerrum > 0.0) {
    if (start) {
      Tcl_AppendResult(interp, "{", (char *)NULL);
      start = 0;
    }
    else
      Tcl_AppendResult(interp, " {", (char *)NULL);
    printCoulombIAToResult(interp);
    Tcl_AppendResult(interp, "}", (char *)NULL);
  }
#endif
  if(lj_force_cap != 0.0) {
    char buffer[TCL_DOUBLE_SPACE];
    
    if (start) {
      Tcl_AppendResult(interp, "{", (char *)NULL);
      start = 0;
    }
    else
      Tcl_AppendResult(interp, " {", (char *)NULL);
    if (lj_force_cap == -1.0)
      Tcl_AppendResult(interp, "ljforcecap individual");
    else {
      Tcl_PrintDouble(interp, lj_force_cap, buffer);
      Tcl_AppendResult(interp, "ljforcecap ", buffer, (char *) NULL);
    }
    Tcl_AppendResult(interp, "}", (char *)NULL);
  }
  return (TCL_OK);
}

int inter_print_bonded(Tcl_Interp *interp, int i)
{
  char buffer[TCL_INTEGER_SPACE];

  Tcl_ResetResult(interp);

  if(i < 0) {
    Tcl_AppendResult(interp, "interaction type must be nonnegative",
		     (char *) NULL);
    return (TCL_ERROR);
  }
  
  make_bond_type_exist(i);
      
  /* print specific interaction information */
  if(i<n_bonded_ia) {
    printBondedIAToResult(interp, i);
    return TCL_OK;
  }

  sprintf(buffer, "%d", i);
  Tcl_AppendResult(interp, "unknown bonded interaction number ", buffer,
		   (char *) NULL);
  return TCL_ERROR;
}

int lj_cos_set_params(int part_type_a, int part_type_b,
		      double eps, double sig, double cut,
		      double offset)
{
  IA_parameters *data, *data_sym;

  double facsq;

  make_particle_type_exist(part_type_a);
  make_particle_type_exist(part_type_b);
    
  data     = get_ia_param(part_type_a, part_type_b);
  data_sym = get_ia_param(part_type_b, part_type_a);
  
  if (!data || !data_sym) {
    return TCL_ERROR;
  }

  /* LJCOS should be symmetrically */
  data_sym->LJCOS_eps    = data->LJCOS_eps    = eps;
  data_sym->LJCOS_sig    = data->LJCOS_sig    = sig;
  data_sym->LJCOS_cut    = data->LJCOS_cut    = cut;
  data_sym->LJCOS_offset = data->LJCOS_offset = offset;

  /* Calculate dependent parameters */
  facsq = driwu2*SQR(sig);
  data_sym->LJCOS_rmin = data->LJCOS_rmin = sqrt(driwu2)*sig;
  data_sym->LJCOS_alfa = data->LJCOS_alfa = PI/(SQR(data->LJCOS_cut)-facsq);
  data_sym->LJCOS_beta = data->LJCOS_beta = PI*(1.-(1./(SQR(data->LJCOS_cut)/facsq-1.)));

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);
  mpi_bcast_ia_params(part_type_b, part_type_a);
  
  return TCL_OK;
}

int lennard_jones_set_params(int part_type_a, int part_type_b,
			     double eps, double sig, double cut,
			     double shift, double offset,
			     double cap_radius)
{
  IA_parameters *data, *data_sym;

  make_particle_type_exist(part_type_a);
  make_particle_type_exist(part_type_b);
    
  data     = get_ia_param(part_type_a, part_type_b);
  data_sym = get_ia_param(part_type_b, part_type_a);

  if (!data || !data_sym) {
    return TCL_ERROR;
  }

  /* LJ should be symmetrically */
  data->LJ_eps    = data_sym->LJ_eps    = eps;
  data->LJ_sig    = data_sym->LJ_sig    = sig;
  data->LJ_cut    = data_sym->LJ_cut    = cut;
  data->LJ_shift  = data_sym->LJ_shift  = shift;
  data->LJ_offset = data_sym->LJ_offset = offset;
 
  if (cap_radius > 0) {
    data->LJ_capradius = cap_radius;
    data_sym->LJ_capradius = cap_radius;
  }

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);
  mpi_bcast_ia_params(part_type_b, part_type_a);

  if (lj_force_cap != -1.0)
    mpi_lj_cap_forces(lj_force_cap);

  return TCL_OK;
}

int gay_berne_set_params(int part_type_a, int part_type_b,
			     double eps, double sig, double cut,
			     double k1, double k2,
			     double mu, double nu)
{
  IA_parameters *data, *data_sym;

  make_particle_type_exist(part_type_a);
  make_particle_type_exist(part_type_b);
    
  data     = get_ia_param(part_type_a, part_type_b);
  data_sym = get_ia_param(part_type_b, part_type_a);

  if (!data || !data_sym) {
    return TCL_ERROR;
  }

  /* GB should be symmetrically */
  data->GB_eps    = data_sym->GB_eps    = eps;
  data->GB_sig    = data_sym->GB_sig    = sig;
  data->GB_cut    = data_sym->GB_cut    = cut;
  data->GB_k1     = data_sym->GB_k1     = k1;
  data->GB_k2     = data_sym->GB_k2     = k2;
  data->GB_mu     = data_sym->GB_mu     = mu;
  data->GB_nu     = data_sym->GB_nu     = nu;
 
   /* Calculate dependent parameters */

  data->GB_chi1 = data_sym->GB_chi1 = ((data->GB_k1*data->GB_k1) - 1) / ((data->GB_k1*data->GB_k1) + 1);
  data->GB_chi2 = data_sym->GB_chi2 = (pow(data->GB_k2,(1/data->GB_mu))-1)/(pow(data->GB_k2,(1/data->GB_mu))+1);

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);
  mpi_bcast_ia_params(part_type_b, part_type_a);

  return TCL_OK;
}

int dihedral_set_params(int bond_type, double bend)
{
  if(bond_type < 0)
    return TCL_ERROR;

  make_bond_type_exist(bond_type);

  bonded_ia_params[bond_type].type = BONDED_IA_DIHEDRAL;
  bonded_ia_params[bond_type].num = 0;

  bend = 0;

  mpi_bcast_ia_params(bond_type, -1); 

  return TCL_OK;
}

int angle_set_params(int bond_type, double bend)
{
  if(bond_type < 0)
    return TCL_ERROR;

  make_bond_type_exist(bond_type);

  bonded_ia_params[bond_type].p.angle.bend = bend;
  bonded_ia_params[bond_type].type = BONDED_IA_ANGLE;
  bonded_ia_params[bond_type].num = 2;
 
  /* broadcast interaction parameters */
  mpi_bcast_ia_params(bond_type, -1); 

  return TCL_OK;
}

int harmonic_set_params(int bond_type, double k, double r)
{
  if(bond_type < 0)
    return TCL_ERROR;

  make_bond_type_exist(bond_type);

  bonded_ia_params[bond_type].p.harmonic.k = k;
  bonded_ia_params[bond_type].p.harmonic.r = r;
  bonded_ia_params[bond_type].type = BONDED_IA_HARMONIC;
  bonded_ia_params[bond_type].p.harmonic.r2 = SQR(bonded_ia_params[bond_type].p.harmonic.r);
  bonded_ia_params[bond_type].num  = 1;

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(bond_type, -1); 

  return TCL_OK;
}

int fene_set_params(int bond_type, double k, double r)
{
  if(bond_type < 0)
    return TCL_ERROR;

  make_bond_type_exist(bond_type);

  bonded_ia_params[bond_type].p.fene.k = k;
  bonded_ia_params[bond_type].p.fene.r = r;

  bonded_ia_params[bond_type].type = BONDED_IA_FENE;
  bonded_ia_params[bond_type].p.fene.r2 = SQR(bonded_ia_params[bond_type].p.fene.r);
  bonded_ia_params[bond_type].num  = 1;

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(bond_type, -1); 
  
  return TCL_OK;
}

int inter_print_non_bonded(Tcl_Interp * interp,
			   int part_type_a, int part_type_b)
{
  IA_parameters *data, *data_sym;

  Tcl_ResetResult(interp);

  make_particle_type_exist(part_type_a);
  make_particle_type_exist(part_type_b);
    
  data     = get_ia_param(part_type_a, part_type_b);
  data_sym = get_ia_param(part_type_b, part_type_a);

  if (!data || !data_sym) {
    Tcl_AppendResult(interp, "particle types must be nonnegative",
		     (char *) NULL);
    return TCL_ERROR;
  }

  return printNonbondedIAToResult(interp, part_type_a, part_type_b);
}


int inter_parse_non_bonded(Tcl_Interp * interp,
			   int part_type_a, int part_type_b,
			   int argc, char ** argv)
{
  double eps, sig, cut, shift, offset, cap_radius;
  double k1, k2, mu, nu; /* parameters needed for Gay-Berne*/

  Tcl_ResetResult(interp);

  if (argc <= 0) {
    Tcl_AppendResult(interp, "wrong # args:  should be \"",
		     "inter <type 1> <type 2> ?interaction? ?values?\"",
		     (char *) NULL);
    return TCL_ERROR;
  }

  /* get interaction parameters */

  /* parse 
   *                        lennard-jones
   * interaction
   */
  if (ARG0_IS_S("lennard-jones")) {

    /* get lennard-jones interaction type */
    if (argc < 6 || argc > 7) {
      Tcl_AppendResult(interp, "lennard-jones needs 5 parameters: "
		       "<lj_eps> <lj_sig> <lj_cut> <lj_shift> <lj_offset>",
		       (char *) NULL);
      return TCL_ERROR;
    }

    /* copy lennard-jones parameters */
    if ((! ARG_IS_D(1, eps))   ||
	(! ARG_IS_D(2, sig))   ||
	(! ARG_IS_D(3, cut))   ||
	(! ARG_IS_D(4, shift)) ||
	(! ARG_IS_D(5, offset)    ))
      {
	Tcl_AppendResult(interp, "lennard-jones needs 5 DOUBLE parameters: "
			 "<lj_eps> <lj_sig> <lj_cut> <lj_shift> <lj_offset>",
			 (char *) NULL);
	return TCL_ERROR;
      }
	
    cap_radius = -1.0;

    if (argc == 7 && (! ARG_IS_D(6, cap_radius))) {
      Tcl_AppendResult(interp, "invalid value for lennard-jones capradius",
		       (char *) NULL);
      return TCL_ERROR;
    }

    CHECK_VALUE(lennard_jones_set_params(part_type_a, part_type_b,
					 eps, sig, cut, shift, offset,
					 cap_radius),
		"particle types must be nonnegative");
  }


  /* parse 
   *                        lj-cos
   * interaction
   */
  if (ARG0_IS_S("lj-cos")) {
  	
    /* this is a quick fix for the inconsistency in the ljcos parameters 
       there are 7 parameters for ljcos, but you read in only four of them.
       The rest is calculated in lj_cos_set_params.
       This is a problem with the blockfile format (Mehmet) 
    */

    if (argc < 5 || argc > 8) {
      Tcl_AppendResult(interp, "lj-cos needs 4 parameters: "
		       "<ljcos_eps> <ljcos_sig> <ljcos_cut> <ljcos_offset>",
		       (char *) NULL);
      return TCL_ERROR;
    }

    /* copy lj-cos parameters */
    if ((! ARG_IS_D(1, eps))   ||
	(! ARG_IS_D(2, sig))   ||
	(! ARG_IS_D(3, cut))   ||
	(! ARG_IS_D(4, offset)    ))
      {
	Tcl_AppendResult(interp, "lj-cos needs 4 DOUBLE parameters: "
			 "<ljcos_eps> <ljcos_sig> <ljcos_cut> <ljcos_offset>",
			 (char *) NULL);
	return TCL_ERROR;
      }

    CHECK_VALUE(lj_cos_set_params(part_type_a, part_type_b, eps, sig, cut, offset),
		"particle types must be nonnegative");
  }

  /* parse 
   *                        gay-berne
   * interaction
   */
  if (ARG0_IS_S("gay-berne")) {
  	
    /* there are 9 parameters for gay-berne, but you read in only 7 of them.
       The rest is calculated in gay_berne_set_params.
    */

    if (argc < 8 || argc > 9) {
      Tcl_AppendResult(interp, "gay-berne needs 7 parameters: "
		       "<gb_eps> <gb_sig> <gb_cut> <gb_k1> <gb_k2> <gb_mu> <gb_nu>",
		       (char *) NULL);
      return TCL_ERROR;
    }

    /* copy gay-berne parameters */
    if ((! ARG_IS_D(1, eps))   ||
	(! ARG_IS_D(2, sig))   ||
	(! ARG_IS_D(3, cut))   ||
	(! ARG_IS_D(4, k1 ))   ||
	(! ARG_IS_D(5, k2 ))   ||
	(! ARG_IS_D(6, mu ))   ||	
	(! ARG_IS_D(7, nu )    ))
      {
	Tcl_AppendResult(interp, "gay-berne needs 7 DOUBLE parameters: "
			 "<gb_eps> <gb_sig> <gb_cut> <gb_k1> <gb_k2> <gb_mu> <gb_nu>",
			 (char *) NULL);
	return TCL_ERROR;
      }

    CHECK_VALUE(gay_berne_set_params(part_type_a, part_type_b, eps, sig, cut, k1, k2, mu, nu),
		"particle types must be nonnegative");
  }
      
  Tcl_AppendResult(interp, "unknown interaction type \"", argv[0],
		   "\"", (char *)NULL);
  return TCL_ERROR;
}

int inter_print_partner_num(Tcl_Interp *interp, int bond_type)
{
  Bonded_ia_parameters * params;
  char buffer[TCL_INTEGER_SPACE];

  if(bond_type < 0) {
    Tcl_AppendResult(interp, "interaction type must be nonnegative",
		     (char *) NULL);
    return TCL_ERROR;
  }
  make_bond_type_exist(bond_type);

  if(bond_type < n_bonded_ia) {
    params = &bonded_ia_params[bond_type];
    sprintf(buffer, "%d", params->num);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    return TCL_OK;
  }
 
  sprintf(buffer, "%d", bond_type);
  Tcl_AppendResult(interp, "unknown bonded interaction number ", buffer,
		   (char *) NULL);
  return TCL_ERROR;
}

int inter_parse_bonded(Tcl_Interp *interp,
		       int bond_type,
		       int argc, char ** argv)
{
  double k, r, bend;

  if (ARG0_IS_S("num")) {
    if (argc == 1)
      return inter_print_partner_num(interp, bond_type);
    else {
	Tcl_AppendResult(interp, "to manny parameters",
			 (char *) NULL);
	return TCL_ERROR;
    }
  }
  
  if (ARG0_IS_S("fene")) {

      if (argc != 3) {
	Tcl_AppendResult(interp, "fene needs 2 parameters: "
			 "<k_fene> <r_fene>", (char *) NULL);
	return TCL_ERROR;
      }

      if ((! ARG_IS_D(1, k)) || (! ARG_IS_D(2, r))) 
	{
	  Tcl_AppendResult(interp, "fene needs 2 DOUBLE parameters: "
			   "<k_fene> <r_fene>", (char *) NULL);
	  return TCL_ERROR;
	}

      CHECK_VALUE(fene_set_params(bond_type, k, r), "bond type must be nonnegative");
  }
    
  if (ARG0_IS_S("harmonic")) {

      if (argc != 3) {
	Tcl_AppendResult(interp, "harmonic needs 2 parameters: "
			 "<k_harmonic> <r_harmonic>", (char *) NULL);
	return TCL_ERROR;
      }

      if ((! ARG_IS_D(1, k)) || (! ARG_IS_D(2, r))) {
	Tcl_AppendResult(interp, "harmonic needs 2 DOUBLE parameters: "
			 "<k_harmonic> <r_harmonic>", (char *) NULL);
	return TCL_ERROR;
      }

      CHECK_VALUE(harmonic_set_params(bond_type, k, r), "bond type must be nonnegative");
  }

  if (ARG0_IS_S("angle")) {

    if (argc != 2) {
      Tcl_AppendResult(interp, "angle needs 1 parameter: "
		       "<bend> ", (char *) NULL);
      return (TCL_ERROR);
    }

    if (! ARG_IS_D(1, bend)) {
      Tcl_AppendResult(interp, "angle needs 1 DOUBLE parameter: "
		       "<bend> ", (char *) NULL);
      return TCL_ERROR;
    }
    
    CHECK_VALUE(angle_set_params(bond_type, bend), "bond type must be nonnegative");
  }
    
  if (ARG0_IS_S("dihedral")) {
 
      /* remove this lines after the implementation of dihedral */
      Tcl_AppendResult(interp, "Do not use interaction type \"", argv[0],
		       "\"", (char *) NULL);
      return TCL_ERROR;

      CHECK_VALUE(dihedral_set_params(bond_type, bend), "bond type must be nonnegative");
  }

  Tcl_AppendResult(interp, "unknown interaction type \"", argv[0],
		   "\"", (char *) NULL);
  return TCL_ERROR;
}

int ljforcecap_set_params(double ljforcecap)
{
  if (lj_force_cap != -1.0)
    mpi_lj_cap_forces(lj_force_cap);

  return TCL_OK;
}

int inter_parse_ljforcecap(Tcl_Interp * interp, int argc, char ** argv)
{
  char buffer[TCL_DOUBLE_SPACE];

  if (argc == 0) {
    if (lj_force_cap == -1.0)
      Tcl_AppendResult(interp, "ljforcecap individual", (char *) NULL);
    else {
      Tcl_PrintDouble(interp, lj_force_cap, buffer);
      Tcl_AppendResult(interp, "ljforcecap ", buffer, (char *) NULL);
    }
    return TCL_OK;
  }

  if (argc > 1) {
    Tcl_AppendResult(interp, "inter ljforcecap takes at most 1 parameter",
		     (char *) NULL);      
    return TCL_ERROR;
  }
  
  if (ARG0_IS_S("individual"))
      lj_force_cap = -1.0;
  else if (! ARG0_IS_D(lj_force_cap) || lj_force_cap < 0) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "force cap must be a nonnegative double value or \"individual\"",
		     (char *) NULL);
    return TCL_ERROR;
  }

  CHECK_VALUE(ljforcecap_set_params(lj_force_cap),
	      "If you can read this, you should change it. (Use the source Luke!)");
  return TCL_ERROR;
}

#ifdef ELECTROSTATICS
void p3m_set_tune_params(double r_cut, int mesh, int cao,
			 double alpha, double accuracy)
{
  if (r_cut >= 0)
    p3m.r_cut = r_cut;

  if (mesh >= 0)
    p3m.mesh[2] = p3m.mesh[1] = p3m.mesh[0] = mesh;

  if (cao >= 0)
    p3m.cao = cao;

  if (alpha >= 0)
    p3m.alpha = alpha;

  if (accuracy >= 0)
    p3m.accuracy = accuracy;

  mpi_bcast_coulomb_params();
}

int p3m_set_params(double r_cut, int mesh, int cao,
		   double alpha, double accuracy)
{
  if(r_cut < 0)
    return -1;

  if(mesh < 0)
    return -2;

  if(cao < 1 || cao > 7 || cao > mesh)
    return -3;

  p3m.r_cut = r_cut;
  p3m.mesh[2] = p3m.mesh[1] = p3m.mesh[0] = mesh;
  p3m.cao = cao;

  if (alpha > 0)
    p3m.alpha = alpha;
  else
    if (alpha != -1.0)
      return -4;

  if (accuracy > 0)
    p3m.accuracy = accuracy;
  else
    if (accuracy != -1.0)
      return -5;

  mpi_bcast_coulomb_params();

  return 0;
}

int p3m_set_mesh_offset(double x, double y, double z)
{
  if(x < 0.0 || x > 1.0 ||
     y < 0.0 || y > 1.0 ||
     z < 0.0 || z > 1.0 )
    return TCL_ERROR;

  p3m.mesh_off[0] = x;
  p3m.mesh_off[1] = y;
  p3m.mesh_off[2] = z;

  mpi_bcast_coulomb_params();

  return TCL_OK;
}

int p3m_set_eps(double eps)
{
  p3m.epsilon = eps;

  mpi_bcast_coulomb_params();

  return TCL_OK;
}

int p3m_set_ninterpol(int n)
{
  if (n < 0)
    return TCL_ERROR;

  p3m.inter = n;

  mpi_bcast_coulomb_params();

  return TCL_OK;
}

int coulomb_set_bjerrum(double bjerrum)
{
  if (bjerrum < 0.0)
    return TCL_ERROR;
  
  coulomb.bjerrum = bjerrum;

  if (coulomb.bjerrum == 0.0) {

    if (coulomb.method == COULOMB_P3M) {

      p3m.bjerrum = 0.0;
      p3m.alpha   = 0.0;
      p3m.r_cut   = 0.0;
      p3m.mesh[0] = 0;
      p3m.mesh[1] = 0;
      p3m.mesh[2] = 0;
      p3m.cao     = 0;

    } else if (coulomb.method == COULOMB_DH) {

      dh_params.bjerrum = 0.0;
      dh_params.r_cut   = 0.0;
      dh_params.kappa   = 0.0;

    } else if (coulomb.method == COULOMB_MMM1D) {

      mmm1d_params.maxPWerror = 1e40;
      mmm1d_params.bjerrum = 0.0;
      mmm1d_params.bessel_cutoff = 0;

    }
 
    mpi_bcast_coulomb_params();
    coulomb.method = COULOMB_NONE;
    mpi_bcast_coulomb_params();

  }

  return TCL_OK;
}

int inter_parse_p3m_tune_params(Tcl_Interp * interp, int argc, char ** argv)
{
  int mesh, cao;
  double r_cut, accuracy;

  mesh = cao = -1;
  r_cut = accuracy = -1.0;

  while(argc > 0) {

    if(ARG0_IS_S("r_cut")) {
      if (! (argc > 1 && ARG1_IS_D(r_cut))) {
	  Tcl_AppendResult(interp, "r_cut expects double",
			   (char *) NULL);
	  return TCL_ERROR;
      }
 
    } else if(ARG0_IS_S("mesh")) {
      if(! (argc > 1 && ARG1_IS_I(mesh))) {
	Tcl_AppendResult(interp, "mesh expects integer",
			   (char *) NULL);
	return TCL_ERROR;
      }
      
    } else if(ARG0_IS_S("cao")) {
      if(! (argc > 1 && ARG1_IS_I(cao))) {
	Tcl_AppendResult(interp, "cao expects integer",
			   (char *) NULL);
	return TCL_ERROR;
      } 

    } else if(ARG0_IS_S("accuracy")) {
      if(! (argc > 1 && ARG1_IS_D(accuracy))) {
	Tcl_AppendResult(interp, "accuracy expects double",
			 (char *) NULL);
	  return TCL_ERROR;
      }

    } else {
      Tcl_AppendResult(interp, "Unkwon p3m tune parameter \"",argv[0],"\"",
		       (char *) NULL);
      return TCL_ERROR;  
    }
    
    argc -= 2;
    argv += 2;
  }

  p3m_set_tune_params(r_cut, mesh, cao, -1.0, accuracy);

  if(P3M_tune_parameters(interp) == TCL_ERROR) 
    return TCL_ERROR;

  return TCL_OK;
}

int inter_parse_p3m(Tcl_Interp * interp, int argc, char ** argv)
{
  double r_cut, alpha, accuracy = -1.0;
  int mesh, cao, i;

  coulomb.method = COULOMB_P3M;
  p3m.bjerrum    = coulomb.bjerrum;
    
#ifdef PARTIAL_PERIODIC
  if(periodic[0] == 0 ||
     periodic[1] == 0 ||
     periodic[2] == 0)
    {
      Tcl_AppendResult(interp, "Need periodicity (1,1,1) with Coulomb P3M",
		       (char *) NULL);
      return TCL_ERROR;  
    }
#endif

  if(node_grid[0] < node_grid[1] || node_grid[1] < node_grid[2]) {
    Tcl_AppendResult(interp, "Node grid not suited for Coulomb P3M. Node grid must be sorted, largest first.", (char *) NULL);
    return TCL_ERROR;  
  }

  if (ARG0_IS_S("tune"))
    return inter_parse_p3m_tune_params(interp, argc-1, argv+1);
      
  if(! ARG0_IS_D(r_cut)) {
    Tcl_AppendResult(interp, "Unknown p3m parameter \"", argv[0],"\"",
		     (char *) NULL);
    return TCL_ERROR;  
  }

  if(argc < 3) {
    Tcl_AppendResult(interp, "p3m needs at least 3 parameters: <r_cut> <mesh> <cao> [<alpha> [accuracy]]",
		     (char *) NULL);
    return TCL_ERROR;  
  }

  if((! ARG_IS_I(1, mesh)) || (! ARG_IS_I(2, cao))) {
    Tcl_AppendResult(interp, "integer expected", (char *) NULL);
    return TCL_ERROR;
  }
	
  if(argc > 3) {
    if(! ARG_IS_D(3, alpha))
      return TCL_ERROR;
  }
  else {
    Tcl_AppendResult(interp, "Automatic p3m tuning not implemented.",
		     (char *) NULL);
    return TCL_ERROR;  
  }

  if(argc > 4) {
    if(! ARG_IS_D(4, accuracy)) {
      Tcl_AppendResult(interp, "double expected", (char *) NULL);
      return TCL_ERROR;
    }
  }

  if ((i = p3m_set_params(r_cut, mesh, cao, alpha, accuracy)) < 0) {
    switch (i) {
    case -1:
      Tcl_AppendResult(interp, "r_cut must be positive", (char *) NULL);
      break;
    case -2:
      Tcl_AppendResult(interp, "mesh must be positive", (char *) NULL);
      break;
    case -3:
      Tcl_AppendResult(interp, "cao must be between 1 and 7 and less than mesh",
		       (char *) NULL);
      break;
    case -4:
      Tcl_AppendResult(interp, "alpha must be positive", (char *) NULL);
      break;
    case -5:
      Tcl_AppendResult(interp, "accuracy must be positive", (char *) NULL);
      break;
    default:;
      Tcl_AppendResult(interp, "unspecified error", (char *) NULL);
    }

    return TCL_ERROR;

  } else {
    return TCL_OK;
  }
  return TCL_ERROR;
}

int inter_parse_p3m_opt_params(Tcl_Interp * interp, int argc, char ** argv)
{
  int i; double d1, d2, d3;

  Tcl_ResetResult(interp);

  while (argc > 0) {
    /* p3m parameter: inter */
    if (ARG0_IS_S("n_interpol")) {
      
      if(argc < 2) {
	Tcl_AppendResult(interp, argv[0], " needs 1 parameter",
			 (char *) NULL);
	return TCL_ERROR;
      }
      
      if (! ARG1_IS_I(i)) {
	Tcl_AppendResult(interp, argv[0], " needs 1 INTEGER parameter",
			 (char *) NULL);
	return TCL_ERROR;
      }
      
      if (p3m_set_ninterpol(i) == TCL_ERROR) {
	Tcl_AppendResult(interp, argv[0], " argument must be positive",
			 (char *) NULL);
	return TCL_ERROR;
      }

      argc -= 2;
      argv += 2;
    }
    
    /* p3m parameter: mesh_off */
    else if (ARG0_IS_S("mesh_off")) {
      
      if(argc < 4) {
	Tcl_AppendResult(interp, argv[0], " needs 3 parameters",
			 (char *) NULL);
	return TCL_ERROR;
      }
	
      if ((! ARG_IS_D(1, d1)) ||
	  (! ARG_IS_D(2, d2)) ||
	  (! ARG_IS_D(3, d3)))
	{
	  Tcl_AppendResult(interp, argv[0], " needs 3 DOUBLE parameters",
			   (char *) NULL);
	  return TCL_ERROR;
	}

      if (p3m_set_mesh_offset(d1, d2 ,d3) == TCL_ERROR)
	{
	  Tcl_AppendResult(interp, argv[0], " parameters have to be between 0.0 an 1.0",
			   (char *) NULL);
	  return TCL_ERROR;
	}

      argc -= 4;
      argv += 4;
    }
    
    /* p3m parameter: epsilon */
    else if(ARG0_IS_S( "epsilon")) {

      if(argc < 2) {
	Tcl_AppendResult(interp, argv[0], " needs 1 parameter",
			 (char *) NULL);
	return TCL_ERROR;
      }

      if (! ARG1_IS_D(d1)) {
	Tcl_AppendResult(interp, argv[0], " needs 1 DOUBLE parameter",
			 (char *) NULL);
	return TCL_ERROR;
      }
	
      if (p3m_set_eps(d1) == TCL_ERROR) {
	Tcl_AppendResult(interp, argv[0], " There is no error msg yet!",
			 (char *) NULL);
	return TCL_ERROR;
      }

      argc -= 2;
      argv += 2;	    
    }
    else {
      Tcl_AppendResult(interp, "Unknown coulomb p3m parameter: \"",argv[0],"\"",(char *) NULL);
      return TCL_ERROR;
    }
  }

  return TCL_OK;
}

int dh_set_params(double kappa, double r_cut)
{
  if(dh_params.kappa < 0.0)
    return -1;

  if(dh_params.r_cut < 0.0)
    return -2;

  dh_params.kappa = kappa;
  dh_params.r_cut = r_cut;

  mpi_bcast_coulomb_params();

  return 1;
}

int inter_parse_dh(Tcl_Interp * interp, int argc, char ** argv)
{
  double kappa, r_cut;
  int i;

  if(argc < 2) {
    Tcl_AppendResult(interp, "Not enough parameters: inter coulomb dh <kappa> <r_cut>", (char *) NULL);
    return TCL_ERROR;
  }
  
  coulomb.method = COULOMB_DH;
  dh_params.bjerrum = coulomb.bjerrum;

  if(! ARG0_IS_D(kappa))
    return TCL_ERROR;
  if(! ARG1_IS_D(r_cut))
    return TCL_ERROR;

  if ( (i = dh_set_params(kappa, r_cut)) < 0) {
    switch (i) {
    case -1:
      Tcl_AppendResult(interp, "dh kappa must be positiv.",(char *) NULL);
      break;
    case -2:
      Tcl_AppendResult(interp, "dh r_cut must be positiv.",(char *) NULL);
      break;
    default:
      Tcl_AppendResult(interp, "unspecified error",(char *) NULL);
    }
    
    return TCL_ERROR;
  }

  return TCL_OK;
}

int inter_parse_mmm1d(Tcl_Interp * interp, int argc, char ** argv)
{
  double switch_rad, maxPWerror;
  int bessel_cutoff;

  if (argc < 2) {
    Tcl_AppendResult(interp, "Not enough parameters: inter coulomb mmm1d <switch radius> {<bessel cutoff>} <maximal error for near formula> | tune  <maximal pairwise error>", (char *) NULL);
    return TCL_ERROR;
  }

  coulomb.method = COULOMB_MMM1D;

  if (ARG0_IS_S("tune")) {
    /* autodetermine bessel cutoff AND switching radius */
    if (! ARG_IS_D(1, maxPWerror))
      return TCL_ERROR;
    bessel_cutoff = -1;
    switch_rad = -1;
  }
  else {
    if (argc == 2) {
      /* autodetermine bessel cutoff */
      if ((! ARG_IS_D(0, switch_rad)) ||
	  (! ARG_IS_D(1, maxPWerror))) 
	return TCL_ERROR;
      bessel_cutoff = -1;
    }
    else if (argc == 3) {
      if((! ARG_IS_D(0, switch_rad)) ||
	 (! ARG_IS_I(1, bessel_cutoff)) ||
	 (! ARG_IS_D(2, maxPWerror))) 
	return TCL_ERROR;
    }
    else {
      Tcl_AppendResult(interp, "Too many parameters: inter coulomb mmm1d <switch radius> {<bessel cutoff>} <maximal error for near formula> | tune  <maximal pairwise error>", (char *) NULL);
      return TCL_ERROR;
    }
  }

  return set_mmm1d_params(interp, coulomb.bjerrum, switch_rad,
			  bessel_cutoff, maxPWerror);
}

int inter_parse_coulomb(Tcl_Interp * interp, int argc, char ** argv)
{
  double d1;

  Tcl_ResetResult(interp);

  if(argc == 0) {
    Tcl_AppendResult(interp, "{", (char *)NULL);
    printCoulombIAToResult(interp);
    Tcl_AppendResult(interp, "}", (char *)NULL);
    return TCL_OK;
  }
  
  if (! ARG0_IS_D(d1)) {
    if (coulomb.method == COULOMB_P3M)
      return inter_parse_p3m_opt_params(interp, argc, argv);
    else {
      Tcl_AppendResult(interp, "(P3M not enabled) Expect: inter coulomb <bjerrum>",
		       (char *) NULL);
      return TCL_ERROR;
    }
  }

  if (coulomb_set_bjerrum(d1) == TCL_ERROR) {
    Tcl_AppendResult(interp, argv[0], "bjerrum length must be positive",
		     (char *) NULL);
    return TCL_ERROR;
  }
    
  argc -= 1;
  argv += 1;

  if (d1 == 0.0 && argc == 0)
    return TCL_OK;

  if(argc < 1) {
    Tcl_AppendResult(interp, "wrong # args for inter coulomb.",
		     (char *) NULL);
    return TCL_ERROR;
  }

  /* check method */
  if(ARG0_IS_S("p3m"))    
    return inter_parse_p3m(interp, argc-1, argv+1);

  if (ARG0_IS_S("dh"))
    return inter_parse_dh(interp, argc-1, argv+1);    
    
  if (ARG0_IS_S("mmm1d"))
    return inter_parse_mmm1d(interp, argc-1, argv+1);

  coulomb.bjerrum = 0.0;
  coulomb.method  = COULOMB_NONE;

  mpi_bcast_coulomb_params();

  Tcl_AppendResult(interp, "Do not know coulomb method \"",argv[1],
		   "\": coulomb switched off", (char *) NULL);

  return TCL_ERROR;
}
#endif

int inter_parse_rest(Tcl_Interp * interp, int argc, char ** argv)
{

  if(ARG0_IS_S("ljforcecap"))
    return inter_parse_ljforcecap(interp, argc-1, argv+1);
  
  if(ARG0_IS_S("coulomb")) {
#ifdef ELECTROSTATICS
    return inter_parse_coulomb(interp, argc-1, argv+1);
#else
    Tcl_AppendResult(interp, "ELECTROSTTICS not compiled (see config.h)", (char *) NULL);
#endif
  }

  Tcl_AppendResult(interp, "unknown interaction type \"", argv[0],
		   "\"", (char *) NULL);

  return TCL_ERROR;
}

int inter(ClientData _data, Tcl_Interp *interp,
	  int argc, char **argv)
{
  int i, j;
  double dummy;

  Tcl_ResetResult(interp);
  
  /* first we handle the special cases 

     1. no parameters
     2. one parameter
     3. two parameters

     then the bonded interactions
     then the non bonded interactions
     then the rest
   */

  /* no argument -> print all interaction informations. */
  if (argc == 1)
    return inter_print_all(interp);

  /* There is only 1 parameter */
  if (argc == 2) {

    if (ARG1_IS_I(i))
      return inter_print_bonded(interp, i);
    
    if (ARG1_IS_D(dummy)) {
      Tcl_AppendResult(interp, "integer or string expected",
		       (char *) NULL);
      return TCL_ERROR;
    }
    
    return inter_parse_rest(interp, argc-1, argv+1);
  }
  
  /* There are only 2 parameters */
  if (argc == 3) {
    
    if (ARG_IS_I(1, i) && ARG_IS_I(2, j))
      return inter_print_non_bonded(interp, i, j);
    
    if (ARG_IS_I(1, i)) {
      Tcl_AppendResult(interp, "not enough arguments",
		       (char *) NULL);
      return TCL_ERROR;
    }
    
    if (ARG1_IS_D(dummy)) {
      Tcl_AppendResult(interp, "integer or string expected",
		       (char *) NULL);
      return TCL_ERROR;
    }
      
    return inter_parse_rest(interp, argc-1, argv+1);
  }

  /****************************************************
   * Here we have more than 2 parameters
   ****************************************************/

  // non bonded interactions
  if (ARG_IS_I(1, i) && ARG_IS_I(2, j))
    return inter_parse_non_bonded(interp, i, j, argc-3, argv+3);

  Tcl_ResetResult(interp);

  // bonded interactions
  if (ARG_IS_I(1, i) && ! ARG_IS_D(2, dummy)) {
    Tcl_ResetResult(interp);
    return inter_parse_bonded(interp, i, argc-2, argv+2);
  }

  Tcl_ResetResult(interp);

  if (ARG_IS_D(1, dummy)) {
      Tcl_AppendResult(interp, "integer or string expected",
		       (char *) NULL);
      return TCL_ERROR;
  }
  
  // named interactions
  return inter_parse_rest(interp, argc-1, argv+1);
}

void make_particle_type_exist(int type)
{
  int ns = type + 1;
  if (ns <= n_particle_types)
    return;

  mpi_bcast_n_particle_types(ns);
}

void make_bond_type_exist(int type)
{
  int i, ns = type + 1;
  if(ns <= n_bonded_ia)
    return;
  bonded_ia_params = (Bonded_ia_parameters *)realloc(bonded_ia_params,
						     ns*sizeof(Bonded_ia_parameters));
  for (i = n_bonded_ia; i < ns; i++)
    bonded_ia_params[i].type = BONDED_IA_NONE;

  n_bonded_ia = ns;
}

#ifdef CONSTRAINTS
int printConstraintToResult(Tcl_Interp *interp, int i)
{
  Constraint *con = &constraints[i];
  char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
  sprintf(buffer, "%d ", i);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  
  switch (con->type) {
  case CONSTRAINT_WAL:
    Tcl_PrintDouble(interp, con->c.wal.n[0], buffer);
    Tcl_AppendResult(interp, "wall normal ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.wal.n[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.wal.n[2], buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    Tcl_PrintDouble(interp, con->c.wal.d, buffer);
    Tcl_AppendResult(interp, " dist ", buffer, (char *) NULL);
    sprintf(buffer, "%d", con->part_rep.r.type);
    Tcl_AppendResult(interp, " type ", buffer, (char *) NULL);
    break;
  case CONSTRAINT_SPH:
    Tcl_PrintDouble(interp, con->c.sph.pos[0], buffer);
    Tcl_AppendResult(interp, "sphere center ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.sph.pos[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.sph.pos[2], buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    Tcl_PrintDouble(interp, con->c.sph.rad, buffer);
    Tcl_AppendResult(interp, " radius ", buffer, (char *) NULL);
    Tcl_PrintDouble(interp, con->c.sph.direction, buffer);
    Tcl_AppendResult(interp, " direction ", buffer, (char *) NULL);
    sprintf(buffer, "%d", con->part_rep.r.type);
    Tcl_AppendResult(interp, " type ", buffer, (char *) NULL);
    break;
  case CONSTRAINT_CYL:
    Tcl_PrintDouble(interp, con->c.cyl.pos[0], buffer);
    Tcl_AppendResult(interp, "cylinder center ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.cyl.pos[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.cyl.pos[2], buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    Tcl_PrintDouble(interp, con->c.cyl.axis[0], buffer);
    Tcl_AppendResult(interp, " axis ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.cyl.axis[1], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.cyl.axis[2], buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    Tcl_PrintDouble(interp, con->c.cyl.rad, buffer);
    Tcl_AppendResult(interp, " radius ", buffer, (char *) NULL);
    Tcl_PrintDouble(interp, con->c.cyl.length, buffer);
    Tcl_AppendResult(interp, " length ", buffer, (char *) NULL);
    Tcl_PrintDouble(interp, con->c.cyl.direction, buffer);
    Tcl_AppendResult(interp, " direction ", buffer, (char *) NULL);
    sprintf(buffer, "%d", con->part_rep.r.type);
    Tcl_AppendResult(interp, " type ", buffer, (char *) NULL);
    break;
  case CONSTRAINT_ROD:
    Tcl_PrintDouble(interp, con->c.rod.pos[0], buffer);
    Tcl_AppendResult(interp, "rod center ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, con->c.rod.pos[1], buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    Tcl_PrintDouble(interp, con->c.rod.rad, buffer);
    Tcl_AppendResult(interp, " radius ", buffer, (char *) NULL);
    Tcl_PrintDouble(interp, con->c.rod.lambda, buffer);
    Tcl_AppendResult(interp, " lambda ", buffer, (char *) NULL);
    sprintf(buffer, "%d", con->part_rep.r.type);
    Tcl_AppendResult(interp, " type ", buffer, (char *) NULL);
    break;
  default:
    sprintf(buffer, "%d", con->type);
    Tcl_AppendResult(interp, "unknown constraint type ", buffer, ".", (char *) NULL);
    return (TCL_OK);
  }

  return (TCL_OK);
}

int constraint_print_all(Tcl_Interp *interp)
{
  int i;
  if(n_constraints>0) Tcl_AppendResult(interp, "{", (char *)NULL);
  for (i = 0; i < n_constraints; i++) {
    if(i>0) Tcl_AppendResult(interp, " {", (char *)NULL);
    printConstraintToResult(interp, i);
    Tcl_AppendResult(interp, "}", (char *)NULL);
  }
  return (TCL_OK);
}

void printConstraintForceToResult(Tcl_Interp *interp, int con)
{
  double f[3];
  char buffer[TCL_DOUBLE_SPACE];

  mpi_get_constraint_force(con, f);

  Tcl_PrintDouble(interp, f[0], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, f[1], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, f[2], buffer);
  Tcl_AppendResult(interp, buffer, (char *) NULL);
}

Constraint *generate_constraint()
{
  n_constraints++;
  constraints = realloc(constraints,n_constraints*sizeof(Constraint));
  constraints[n_constraints-1].type = CONSTRAINT_NONE;
  constraints[n_constraints-1].part_rep.r.identity = -n_constraints;
  
  return &constraints[n_constraints-1];
}

int constraint_wall(Constraint *con, Tcl_Interp *interp,
		    int argc, char **argv)
{
  int i;
  double norm;
  con->type = CONSTRAINT_WAL;
  /* invalid entries to start of */
  con->c.wal.n[0] = 
    con->c.wal.n[1] = 
    con->c.wal.n[2] = 0;
  con->c.wal.d = 0;
  con->part_rep.r.type = -1;
  while (argc > 0) {
    if(!strncmp(argv[0], "normal", strlen(argv[0]))) {
      if(argc < 4) {
	Tcl_AppendResult(interp, "constraint wall normal <nx> <ny> <nz> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(con->c.wal.n[0])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[2], &(con->c.wal.n[1])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[3], &(con->c.wal.n[2])) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 4; argv += 4;
    }
    else if(!strncmp(argv[0], "dist", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint wall dist <d> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(con->c.wal.d)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "type", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint wall type <t> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetInt(interp, argv[1], &(con->part_rep.r.type)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else
      break;
  }
  norm = SQR(con->c.wal.n[0])+SQR(con->c.wal.n[1])+SQR(con->c.wal.n[2]);
  if (norm < 1e-10 || con->part_rep.r.type < 0) {
    Tcl_AppendResult(interp, "usage: constraint wall normal <nx> <ny> <nz> dist <d> type <t>",
		     (char *) NULL);
    return (TCL_ERROR);    
  }
  for(i=0;i<3;i++) con->c.wal.n[i] /= sqrt(norm);

  make_particle_type_exist(con->type);

  return (TCL_OK);
}

int constraint_sphere(Constraint *con, Tcl_Interp *interp,
		      int argc, char **argv)
{
  con->type = CONSTRAINT_SPH;

  /* invalid entries to start of */
  con->c.sph.pos[0] = 
    con->c.sph.pos[1] = 
    con->c.sph.pos[2] = 0;
  con->c.sph.rad = 0;
  con->c.sph.direction = -1;
  con->part_rep.r.type = -1;

  while (argc > 0) {
    if(!strncmp(argv[0], "center", strlen(argv[0]))) {
      if(argc < 4) {
	Tcl_AppendResult(interp, "constraint sphere center <x> <y> <z> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(con->c.sph.pos[0])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[2], &(con->c.sph.pos[1])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[3], &(con->c.sph.pos[2])) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 4; argv += 4;
    }
    else if(!strncmp(argv[0], "radius", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint sphere radius <r> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(con->c.sph.rad)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "direction", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "1/-1 or inside/outside is expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (!strncmp(argv[1], "inside", strlen(argv[1])))
	con->c.sph.direction = -1;
      else if (!strncmp(argv[1], "outside", strlen(argv[1])))
	con->c.sph.direction = 1;
      else if (Tcl_GetDouble(interp, argv[1], &(con->c.sph.direction)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "type", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint sphere type <t> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetInt(interp, argv[1], &(con->part_rep.r.type)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else
      break;
  }

  if (con->c.sph.rad < 0. || con->part_rep.r.type < 0) {
    Tcl_AppendResult(interp, "usage: constraint sphere center <x> <y> <z> radius <d> direction <direction> type <t>",
		     (char *) NULL);
    return (TCL_ERROR);    
  }

  make_particle_type_exist(con->type);

  return (TCL_OK);
}

int constraint_cylinder(Constraint *con, Tcl_Interp *interp,
		    int argc, char **argv)
{
  double axis_len;
  int i;

  con->type = CONSTRAINT_CYL;
  /* invalid entries to start of */
  con->c.cyl.pos[0] = 
    con->c.cyl.pos[1] = 
    con->c.cyl.pos[2] = 0;
  con->c.cyl.axis[0] = 
    con->c.cyl.axis[1] = 
    con->c.cyl.axis[2] = 0;
  con->c.cyl.rad = 0;
  con->c.cyl.length = 0;
  con->c.cyl.direction = 0;
  con->part_rep.r.type = -1;
  while (argc > 0) {
    if(!strncmp(argv[0], "center", strlen(argv[0]))) {
      if(argc < 4) {
	Tcl_AppendResult(interp, "constraint cylinder center <x> <y> <z> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(con->c.cyl.pos[0])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[2], &(con->c.cyl.pos[1])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[3], &(con->c.cyl.pos[2])) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 4; argv += 4;
    }
    else if(!strncmp(argv[0], "axis", strlen(argv[0]))) {
      if(argc < 4) {
	Tcl_AppendResult(interp, "constraint cylinder axis <rx> <ry> <rz> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(con->c.cyl.axis[0])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[2], &(con->c.cyl.axis[1])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[3], &(con->c.cyl.axis[2])) == TCL_ERROR)
	return (TCL_ERROR);

      argc -= 4; argv += 4;    
    }
    else if(!strncmp(argv[0], "radius", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint cylinder radius <rad> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(con->c.cyl.rad)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "length", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint cylinder length <len> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(con->c.cyl.length)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "direction", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint cylinder direction <dir> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (!strncmp(argv[1], "inside", strlen(argv[1])))
	con->c.cyl.direction = -1;
      else if (!strncmp(argv[1], "outside", strlen(argv[1])))
	con->c.cyl.direction = 1;
      else if (Tcl_GetDouble(interp, argv[1], &(con->c.cyl.direction)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "type", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint cylinder type <t> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetInt(interp, argv[1], &(con->part_rep.r.type)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else
      break;
  }

  axis_len=0.;
  for (i=0;i<3;i++)
    axis_len += SQR(con->c.cyl.axis[i]);

  if (con->c.cyl.rad < 0. || con->part_rep.r.type < 0 || axis_len < 1e-30 ||
      con->c.cyl.direction == 0 || con->c.cyl.length <= 0) {
    Tcl_AppendResult(interp, "usage: constraint cylinder center <x> <y> <z> axis <rx> <ry> <rz> radius <rad> length <length> direction <direction> type <t>",
		     (char *) NULL);
    return (TCL_ERROR);    
  }

  /*normalize the axis vector */
      axis_len = sqrt (axis_len);
      for (i=0;i<3;i++) {
	con->c.cyl.axis[i] /= axis_len;
      }
      

  make_particle_type_exist(con->type);

  return (TCL_OK);
}

int constraint_rod(Constraint *con, Tcl_Interp *interp,
		   int argc, char **argv)
{
  con->type = CONSTRAINT_ROD;
  /* invalid entries to start of */
  con->c.rod.pos[0] = con->c.rod.pos[1] = 0;
  con->c.rod.rad = -1;
  con->c.rod.lambda = 0;
  con->part_rep.r.type = -1;
  while (argc > 0) {
    if(!strncmp(argv[0], "center", strlen(argv[0]))) {
      if(argc < 3) {
	Tcl_AppendResult(interp, "constraint rod center <px> <py> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(con->c.rod.pos[0])) == TCL_ERROR ||
	  Tcl_GetDouble(interp, argv[2], &(con->c.rod.pos[1])) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 3; argv += 3;
    }
    else if(!strncmp(argv[0], "radius", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint rod radius <r> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(con->c.rod.rad)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "type", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint rod type <t> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetInt(interp, argv[1], &(con->part_rep.r.type)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else if(!strncmp(argv[0], "lambda", strlen(argv[0]))) {
      if (argc < 1) {
	Tcl_AppendResult(interp, "constraint rod lambda <l> expected", (char *) NULL);
	return (TCL_ERROR);
      }
      if (Tcl_GetDouble(interp, argv[1], &(con->c.rod.lambda)) == TCL_ERROR)
	return (TCL_ERROR);
      argc -= 2; argv += 2;
    }
    else
      break;
  }
  if (con->c.rod.rad < 0 || con->part_rep.r.type < 0) {
    Tcl_AppendResult(interp, "usage: constraint rod center <px> <py> rad <r> type <t>",
		     (char *) NULL);
    return (TCL_ERROR);
  }

  make_particle_type_exist(con->type);

  return (TCL_OK);
}

#endif

int constraint(ClientData _data, Tcl_Interp *interp,
	       int argc, char **argv)
{
#ifdef CONSTRAINTS
  int status, c_num;

  if (argc < 2) return constraint_print_all(interp);
  
  if(!strncmp(argv[1], "wall", strlen(argv[1]))) {
    status = constraint_wall(generate_constraint(),interp, argc - 2, argv + 2);
    mpi_bcast_constraint(-1);
    return status;
  }
  else if(!strncmp(argv[1], "sphere", strlen(argv[1]))) {
    status = constraint_sphere(generate_constraint(),interp, argc - 2, argv + 2);
    mpi_bcast_constraint(-1);
    return status;
  }
  else if(!strncmp(argv[1], "cylinder", strlen(argv[1]))) {
    status = constraint_cylinder(generate_constraint(),interp, argc - 2, argv + 2);
    mpi_bcast_constraint(-1);
    return status;
  }
  else if(!strncmp(argv[1], "rod", strlen(argv[1]))) {
    status = constraint_rod(generate_constraint(),interp, argc - 2, argv + 2);
    mpi_bcast_constraint(-1);
    return status;
  }
  else if(!strncmp(argv[1], "force", strlen(argv[1]))) {
    if(argc < 3) {
      Tcl_AppendResult(interp, "which particles force?",(char *) NULL);
      return (TCL_ERROR);
    }
    if(Tcl_GetInt(interp, argv[2], &(c_num)) == TCL_ERROR) return (TCL_ERROR);
    if(c_num < 0 || c_num >= n_constraints) {
      Tcl_AppendResult(interp, "constraint does not exist",(char *) NULL);
      return (TCL_ERROR);
    }
    printConstraintForceToResult(interp, c_num);
    return (TCL_OK);
  }
  else if(!strncmp(argv[1], "delete", strlen(argv[1]))) {
    if(argc < 3) {
      /* delete all */
      mpi_bcast_constraint(-2);
      return (TCL_OK);
    }
    if(Tcl_GetInt(interp, argv[2], &(c_num)) == TCL_ERROR) return (TCL_ERROR);
    if(c_num < 0 || c_num >= n_constraints) {
      Tcl_AppendResult(interp, "Can not delete non existing constraint",(char *) NULL);
      return (TCL_ERROR);
    }
    mpi_bcast_constraint(c_num);
    return (TCL_OK);    
  }
  else if (argc == 2 && Tcl_GetInt(interp, argv[1], &c_num) == TCL_OK) {
    printConstraintToResult(interp, c_num);
    return (TCL_OK);
  }

  Tcl_AppendResult(interp, "possible constraints: wall sphere cylinder or constraint delete {c} to delete constraint(s)",(char *) NULL);
  return (TCL_ERROR);
#else
  Tcl_AppendResult(interp, "Constraints not compiled in!" ,(char *) NULL);
  return (TCL_ERROR);
#endif
}
