/*
  Copyright (C) 2012,2013,2014 The ESPResSo project
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
#ifndef _OBJECT_IN_FLUID_AREA_FORCE_GLOBAL_H
#define _OBJECT_IN_FLUID_AREA_FORCE_GLOBAL_H
/** \file area_force_global.hpp
 *  Routines to calculate the AREA_FORCE_GLOBAL energy or/and and force 
 *  for a particle triple (triangle from mesh). (Dupin2007)
 *  \ref forces.cpp
*/

#include "utils.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "cells.hpp"
#include "lb.hpp"
#include "grid.hpp"
#include "errorhandling.hpp"

/** set parameters for the AREA_FORCE_GLOBAL potential. 
*/
int area_force_global_set_params(int bond_type, double A0_g, double ka_g);

/************************************************************/

/** called in force_calc() from within forces.cpp
 *  calculates the global area for a cell before the forces are handled
 *  sums up parts for area with mpi_reduce from local triangles
 *  synchronization with allreduce
 *
 *  !!! loop over particles from domain_decomposition !!!
 */  

inline void calc_area_global(double *area, int molType){ //first-fold-then-the-same approach
	double partArea=0.;

	/** loop over particles */
	int c, np, i ,j;
	Cell *cell;
	Particle *p, *p1, *p2, *p3;
	double p11[3],p22[3],p33[3];
	int img[3];
	double AA[3],BB[3];
	Bonded_ia_parameters *iaparams;
    int type_num, n_partners,id;
    BondedInteraction type;

	int test=0;

	/* Loop local cells */
	for (c = 0; c < local_cells.n; c++) {
		cell = local_cells.cell[c];
		p   = cell->part;
		np  = cell->n;
		
		/* Loop cell particles */
		for(i=0; i < np; i++) {				
			j = 0;
			p1=&p[i];
			while(j<p1->bl.n){
				/* bond type */
				type_num = p1->bl.e[j++];
				iaparams = &bonded_ia_params[type_num];			
				type = iaparams->type;
				n_partners = iaparams->num;
				id=p1->p.mol_id;
				//printf("neigh=%d, type=%d type_num=%d\n", p1->bl.n-1, type, type_num);
				if(type == BONDED_IA_AREA_FORCE_GLOBAL && id == molType){ // BONDED_IA_AREA_FORCE_GLOBAL with correct molType  
					test++;
					/* fetch particle 2 */
					p2 = local_particles[p1->bl.e[j++]];
                    if (!p2) {
                        ostringstream msg;
                        msg <<"area calc: bond broken between particles " << p1->p.identity << " and " << p1->bl.e[j-1] << " (particles not stored on the same node - area_force_global1); n " << p1->bl.n << " max " << p1->bl.max ;
                        runtimeError(msg);
						return;
					}
					/* fetch particle 3 */
					//if(n_partners>2){
					p3 = local_particles[p1->bl.e[j++]];
                    if (!p3) {
                        ostringstream msg;
                        msg <<"area calc: bond broken between particles " << p1->p.identity << ", " << p1->bl.e[j-2] << " and " << p1->bl.e[j-1] << " (particles not stored on the same node - area_force_global1); n " << p1->bl.n << " max " << p1->bl.max ;
                        runtimeError(msg);
						return;
					}
					// remaining neighbors fetched
					
					// getting unfolded positions of all particles
					#ifdef GHOST_FLAG
					// first find out which particle out of p1, p2 (possibly p3, p4) is not a ghost particle. In almost all cases it is p1, however, it might be other one. we call this particle reference particle.
					if (p1->l.ghost != 1) {
						//unfold non-ghost particle using image, because for physical particles, the structure p->l.i is correctly set
						memcpy(p11, p1->r.p, 3*sizeof(double));
						memcpy(img, p1->l.i, 3*sizeof(int));
						unfold_position(p11,img);
						// other coordinates are obtained from its relative positions to the reference particle
						get_mi_vector(AA, p2->r.p, p11);
						get_mi_vector(BB, p3->r.p, p11);
						for (int i=0; i < 3; i++) { p22[i] = p11[i] + AA[i]; p33[i] = p11[i] + BB[i]; }
					} else {
						// in case the first particle is a ghost particle
						if (p2->l.ghost != 1) {
							memcpy(p22, p2->r.p, 3*sizeof(double));
							memcpy(img, p2->l.i, 3*sizeof(int));
							unfold_position(p22,img);
							get_mi_vector(AA, p1->r.p, p22);
							get_mi_vector(BB, p3->r.p, p22);
							for (int i=0; i < 3; i++) { p11[i] = p22[i] + AA[i]; p33[i] = p22[i] + BB[i]; }
						} else {
							// in case the first and the second particle are ghost particles
							if (p3->l.ghost != 1) {
								memcpy(p33, p3->r.p, 3*sizeof(double));
								memcpy(img, p3->l.i, 3*sizeof(int));
								unfold_position(p33,img);
								get_mi_vector(AA, p1->r.p, p33);
								get_mi_vector(BB, p2->r.p, p33);
								for (int i=0; i < 3; i++) { p11[i] = p33[i] + AA[i]; p22[i] = p33[i] + BB[i]; }
							} else {
								printf("Something wrong in area_force_local.hpp: All particles in a bond are ghost particles, impossible to unfold the positions...");
								return;
							}
						}
					}
					#endif
					#ifndef GHOST_FLAG
						// if ghost flag was not defined we have no other option than to assume the first particle is a physical one.
						memcpy(p11, p1->r.p, 3*sizeof(double));
						memcpy(img, p1->l.i, 3*sizeof(int));
						unfold_position(p11,img);
						// other coordinates are obtained from its relative positions to the reference particle
						get_mi_vector(AA, p2->r.p, p11);
						get_mi_vector(BB, p3->r.p, p11);
						for (int i=0; i < 3; i++) { p22[i] = p11[i] + AA[i]; p33[i] = p11[i] + BB[i]; }
					#endif
					// unfolded positions correct
					partArea += area_triangle(p11,p22,p33);
				}
				else{
					j+=n_partners;
				}	
			}
		}
    }

	MPI_Allreduce(&partArea, area, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}


inline void add_area_global_force(double area, int molType){  //first-fold-then-the-same approach
	double aa, force1[3], force2[3], force3[3], rh[3], hn, h[3];
	int k;
	
	/** loop over particles */
	int c, np, i ,j;
	Cell *cell;
	Particle *p, *p1, *p2, *p3;
	double p11[3],p22[3],p33[3];
	double AA[3],BB[3];
	int img[3];

	Bonded_ia_parameters *iaparams;
    int type_num, n_partners,id;
    BondedInteraction type;

	int test=0;
	
	/* Loop local cells */
	for (c = 0; c < local_cells.n; c++) {
		cell = local_cells.cell[c];
		p   = cell->part;
		np  = cell->n;
		
		/* Loop cell particles */
		for(i=0; i < np; i++) {				
			j = 0;
			p1=&p[i];
			//printf("i=%d neigh=%d\n", i, p1->bl.n);
			while(j<p1->bl.n){
				/* bond type */
				type_num = p1->bl.e[j++];
				iaparams = &bonded_ia_params[type_num];
				type = iaparams->type;
				n_partners = iaparams->num;
				id=p1->p.mol_id;
				//printf("neigh=%d, type=%d type_num=%d\n", p1->bl.n-1, type, type_num);
				//printf("id %d molType %d\n", id, molType); 
				if(type == BONDED_IA_AREA_FORCE_GLOBAL && id == molType){ // BONDED_IA_VOLUME_FORCE with correct molType
					test++;
					/* fetch particle 2 */
					p2 = local_particles[p1->bl.e[j++]];
                    if (!p2) {
                        ostringstream msg;
                        msg <<"add area: bond broken between particles " << p1->p.identity << " and " << p1->bl.e[j-1] << " (particles not stored on the same node - area_force_global2); n " << p1->bl.n << " max " << p1->bl.max ;
                        runtimeError(msg);
						return;
					}
					/* fetch particle 3 */
					//if(n_partners>2){
					p3 = local_particles[p1->bl.e[j++]];
                    if (!p3) {
                        ostringstream msg;
                        msg <<"add area: bond broken between particles " << p1->p.identity << ", " << p1->bl.e[j-2] << " and " << p1->bl.e[j-1] << " (particles not stored on the same node); n " << p1->bl.n << " max " << p1->bl.max;
                        runtimeError(msg);
						return;
					}
					
					// getting unfolded positions of all particles
					#ifdef GHOST_FLAG
					// first find out which particle out of p1, p2 (possibly p3, p4) is not a ghost particle. In almost all cases it is p1, however, it might be other one. we call this particle reference particle.
					if (p1->l.ghost != 1) {
						//unfold non-ghost particle using image, because for physical particles, the structure p->l.i is correctly set
						memcpy(p11, p1->r.p, 3*sizeof(double));
						memcpy(img, p1->l.i, 3*sizeof(int));
						unfold_position(p11,img);
						// other coordinates are obtained from its relative positions to the reference particle
						get_mi_vector(AA, p2->r.p, p11);
						get_mi_vector(BB, p3->r.p, p11);
						for (int i=0; i < 3; i++) { p22[i] = p11[i] + AA[i]; p33[i] = p11[i] + BB[i]; }
					} else {
						// in case the first particle is a ghost particle
						if (p2->l.ghost != 1) {
							memcpy(p22, p2->r.p, 3*sizeof(double));
							memcpy(img, p2->l.i, 3*sizeof(int));
							unfold_position(p22,img);
							get_mi_vector(AA, p1->r.p, p22);
							get_mi_vector(BB, p3->r.p, p22);
							for (int i=0; i < 3; i++) { p11[i] = p22[i] + AA[i]; p33[i] = p22[i] + BB[i]; }
						} else {
							// in case the first and the second particle are ghost particles
							if (p3->l.ghost != 1) {
								memcpy(p33, p3->r.p, 3*sizeof(double));
								memcpy(img, p3->l.i, 3*sizeof(int));
								unfold_position(p33,img);
								get_mi_vector(AA, p1->r.p, p33);
								get_mi_vector(BB, p2->r.p, p33);
								for (int i=0; i < 3; i++) { p11[i] = p33[i] + AA[i]; p22[i] = p33[i] + BB[i]; }
							} else {
								printf("Something wrong in area_force_local.hpp: All particles in a bond are ghost particles, impossible to unfold the positions...");
								return;
							}
						}
					}
					#endif
					#ifndef GHOST_FLAG
						// if ghost flag was not defined we have no other option than to assume the first particle is a physical one.
						memcpy(p11, p1->r.p, 3*sizeof(double));
						memcpy(img, p1->l.i, 3*sizeof(int));
						unfold_position(p11,img);
						// other coordinates are obtained from its relative positions to the reference particle
						get_mi_vector(AA, p2->r.p, p11);
						get_mi_vector(BB, p3->r.p, p11);
						for (int i=0; i < 3; i++) { p22[i] = p11[i] + AA[i]; p33[i] = p11[i] + BB[i]; }
					#endif
					// unfolded positions correct
					
					for(k=0;k<3;k++){
						h[k]=1.0/3.0 *(p11[k]+p22[k]+p33[k]);
					}
					
					aa=( area - iaparams->p.area_force_global.A0_g) / iaparams->p.area_force_global.A0_g;

					//aminusb(3,h,p11,rh);				// area_forces for each triangle node
					vecsub(h,p11,rh);				// area_forces for each triangle node
					hn=normr(rh);
					for(k=0;k<3;k++) {
						force1[k] =  iaparams->p.area_force_global.ka_g * aa * rh[k]/hn;
						//(&part1)->f.f[k]+=force[k];
					}
					//aminusb(3,h,p22,rh);				// area_forces for each triangle node
					vecsub(h,p22,rh);				// area_forces for each triangle node
					hn=normr(rh);
					for(k=0;k<3;k++) {
						force2[k] =  iaparams->p.area_force_global.ka_g * aa * rh[k]/hn;
						//(&part2)->f.f[k]+=force[k];
					}
					//aminusb(3,h,p33,rh);				// area_forces for each triangle node
					vecsub(h,p33,rh);				// area_forces for each triangle node
					hn=normr(rh);
					for(k=0;k<3;k++) {
						force3[k] =  iaparams->p.area_force_global.ka_g * aa * rh[k]/hn;
						//(&part3)->f.f[k]+=force[k];
					}
	
					for(k=0;k<3;k++) {
						p1->f.f[k] += force1[k]; 
						p2->f.f[k] +=  force2[k];
						p3->f.f[k] +=  force3[k];
					}

				}
				else{
					j+=n_partners;
				}
			}
		}
    }
	
}

#endif 
