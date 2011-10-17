#include "particle_data.h"
#include "interaction_data.h"
#include "virtual_sites_relative.h"
#include "collision.h"
#include "virtual_sites.h"
#include "integrate.h"
#include "cells.h"


#ifdef COLLISION_DETECTION


// During force calculation, colliding particles are recorded in thequeue
// The queue is processed after force calculation, when it is save to add
// particles
collision_struct * collision_queue;

// Number of collisions recoreded in the queue
int number_of_collisions;

// Bond type used for marker bonds
int collision_detection_bond_marker_real=0;


// bond type used between virtual sites 
int collision_detection_bond_marker_virtual=1;


// Detect a collision between the given particles.
// In case of a collision, a bond is added between them as marker
// and the collision is recorded in the queue
void detect_collision(Particle* p1, Particle* p2)
{
  //printf("in collsiion_detction"); 
  double dist_betw_part, vec21[3], collisioncriter=1.15;
  int part1, part2, the_bond_type_added_on_collision=0, size;

  // Obtain distance between particles
  dist_betw_part = distance2vec(p1->r.p, p2->r.p, vec21);
  if (dist_betw_part > collisioncriter)
    return;

  part1 = p1->p.identity;
  part2 = p2->p.identity;
      
  // Retrieving the particles from local_particles is necessary, because the particle might be a
  // ghost, and those don't contain bonding info
  p1 = local_particles[part1];
  p2 = local_particles[part2];

  // Ignore virtual particles
  if ((p1->p.isVirtual) || (p2->p.isVirtual))
    return;


  // Check, if there's already a bond between the particles
  // First check the bonds of p1 
  int i = 0;
  int found = 0;
  while(i < p1->bl.n) {
    size = bonded_ia_params[p1->bl.e[i]].num;

    if (p1->bl.e[i] == the_bond_type_added_on_collision &&
        p1->bl.e[i + 1] == part2) {
      // There's a bond, already. Nothing to do for these particles
      return;
    }
    i += size + 1;
  }
      
  // Check, if a bond is already stored in p2
  i = 0;
  while(i < p2->bl.n) {
    size = bonded_ia_params[p2->bl.e[i]].num;

    /* COMPARE P2 WITH P1'S BONDED PARTICLES*/

    if (p2->bl.e[i] == the_bond_type_added_on_collision &&
        p2->bl.e[i + 1] == part1) {
      return;
    }
    i += size + 1;
  }


  // If we're still here, there is no previous bond between the particles

  // Create the marker bond between the particles
  int bondG[2];
  bondG[0]=collision_detection_bond_marker_real;
  bondG[1]=part2;
  local_change_bond(part1, bondG, 0);
  //printf("real bond was created\n");   
  // Insert collision info into the queue
      
  // Point of collision
  double new_position[3];
  for (i=0;i<3;i++) {
    new_position[i] =p1->r.p[i] - vec21[i] * 0.50;
    //printf("connection point is %f\n",new_position[i]);
  }
       
  number_of_collisions = number_of_collisions+1;
  //	printf("number of collisions each time are %d\n",number_of_collisions);       
  // Allocate mem for the new collision info
  collision_queue = (collision_struct *) realloc (collision_queue, (number_of_collisions) * sizeof(collision_struct));
      
  // Save the collision      
  collision_queue[number_of_collisions-1].pp1 = part1;
  collision_queue[number_of_collisions-1].pp2 = part2;
  for (i=0;i<3;i++) {
    collision_queue[number_of_collisions-1].point_of_collision[i] = new_position[i]; 
  }

}

void prepare_collision_queue()
{
  
  number_of_collisions=0;

  collision_queue = (collision_struct *) malloc (sizeof(collision_struct));

}


// Handle the collisions stored in the queue
void handle_collisions ()
{
  //	printf("number of collisions in handle collision are %d\n",number_of_collisions);  
  int delete =0, bondG[2], i;

  if (number_of_collisions > 0) {
    // Go through the queue
    for (i=0;i<number_of_collisions;i++) {
      //  printf("Handling collision of particles %d %d\n", collision_queue[i].pp1, collision_queue[i].pp2);
      //  fflush(stdout);
   
      int j;

      // The following lines will remove the relative velocity from
      // colliding particles
      //   double v[3];
      //   for (j=0;j<3;j++)
      //   {
      //    v[j] =0.5 *((local_particles[collision_queue[i].pp1])->m.v[j] +(local_particles[collision_queue[i].pp2])->m.v[j]);
      //    (local_particles[collision_queue[i].pp1])->m.v[j] =v[j];
      //    (local_particles[collision_queue[i].pp2])->m.v[j] =v[j];
      //   }

      // Create virtual sites and bind them gotether
  
      // Virtual site related to first particle in the collision
      place_particle(max_seen_particle+1,collision_queue[i].point_of_collision);
      vs_relate_to(max_seen_particle,collision_queue[i].pp1);
      (local_particles[max_seen_particle])->p.isVirtual=1;
  
      //printf("virtual1 was created: %d rel_to %d\n", max_seen_particle, collision_queue[i].pp1); 
      // Virtual particle related to 2nd particle of the collision
      place_particle(max_seen_particle+1,collision_queue[i].point_of_collision);
      vs_relate_to(max_seen_particle,collision_queue[i].pp2);
      (local_particles[max_seen_particle])->p.isVirtual=1;
      //printf("virtual2 was created: %d rel_to %d\n", max_seen_particle, collision_queue[i].pp2); 
  
      // Create bond between the virtual particles
      bondG[0] =collision_detection_bond_marker_virtual;
      bondG[1] =max_seen_particle-1;
      local_change_bond(max_seen_particle, bondG, 0);
      //printf("virtual bond was created\n");   
    }

    on_particle_change();
  }

  // Reset the collision queue	 
  number_of_collisions = 0;
  free(collision_queue);

}





#endif
