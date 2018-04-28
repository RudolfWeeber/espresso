

import espressomd
from espressomd.magnetostatics import DipolarP3M
from espressomd.observables import MagneticDipoleMoment
from espressomd.polymer import create_polymer
from espressomd.interactions import HarmonicBond
import sys
from numpy import *
import numpy as np
from time import time
import argparse
import cPickle
import gzip

parser =argparse.ArgumentParser()
parser.add_argument("--base",type=str,required=True)
tmp=parser.parse_args()
base=tmp.base


args= cPickle.load(open(base+"_args.pcl"))
print args

if args.n_nodes<1: raise Exception()
if args.phi<=0 or args.phi>=1: raise Exception()
if args.lam<0: raise Exception()
sigma_mag=10.
sigma_pol=1.
sigma_mixed=0.5*(sigma_mag+sigma_pol)
eps=1
cut_mag=sigma_mag *2**(1./6)
cut_pol=sigma_pol *2**(1./6)
cut_mixed=sigma_mixed *2**(1./6)

kT=1.
gamma=1
gamma_mag=sigma_mag/sigma_pol
dipm=(args.lam*sigma_mag**3)**0.5




# box
l=(args.n_nodes/6. *pi *sigma_mag**3/args.phi)**(1./3.)
r=sigma_mag/2
if abs(args.n_nodes*4./3. *pi*r**3/l**3-args.phi) >1E-4: raise Exception()

s=espressomd.System(box_l=[l]*3)

s.seed=[args.seed]*s.cell_system.get_state()["n_nodes"]
np.random.seed(args.seed)
# integrator and cells
s.cell_system.skin=0.4
s.min_global_cut=cut_mixed*1.01
s.time_step=0.008

s.thermostat.set_langevin(kT=kT,gamma=gamma)

# interactions
# magnetic particles
s.non_bonded_inter[0,0].lennard_jones.set_params(epsilon=eps,sigma=sigma_mag,cutoff=cut_mag,shift="auto")
# Polymers
s.non_bonded_inter[1,1].lennard_jones.set_params(epsilon=eps,sigma=sigma_pol,cutoff=cut_pol,shift="auto")
# Mixed
s.non_bonded_inter[0,1].lennard_jones.set_params(epsilon=eps,sigma=sigma_mixed,cutoff=cut_mixed,shift="auto")
# Chain ends
s.non_bonded_inter[3,3].lennard_jones.set_params(epsilon=eps,sigma=sigma_pol,cutoff=cut_pol,shift="auto")
# Chain ends and polymer beeds
s.non_bonded_inter[1,3].lennard_jones.set_params(epsilon=eps,sigma=sigma_pol,cutoff=cut_pol,shift="auto")
# chain ends and mag particles
s.non_bonded_inter[0,3].lennard_jones.set_params(epsilon=eps,sigma=sigma_mixed,cutoff=cut_mixed,shift="auto")


bond_pol=HarmonicBond(r_0=cut_pol,k=1000)
s.bonded_inter.add(bond_pol)
bond_mag_pol=HarmonicBond(r_0=cut_mixed,k=800)
s.bonded_inter.add(bond_mag_pol)

# Load particles
parts=cPickle.load(open(base+"_equil.pcl"))
s.part.add(**parts)
s.part[0].type=2
s.part[0].type=0
print s.part[0]
print s.part[args.n_nodes+1]
print s.analysis.energy()

# Identify chain ends
ends=[]
s.part[args.n_nodes].type=3
ends.append(args.n_nodes)
for i in range(args.n_nodes+1,len(s.part)):
  if s.part[i].bonds==():
      s.part[i].type=3
      s.part[i-1].type=3
      ends.append(i-1)
      ends.append(i)
print "ends",ends
print s.analysis.energy()

# diffusion coeff for nodes
nodes =s.part[:args.n_nodes]
print "gamma"
for node in nodes:
    print node.id
    node.gamma=float(sigma_mag)/sigma_pol
    node.gamma_rot=float(sigma_mag)/sigma_pol
    node=1,1,1
print "done"

print "lambda",args.lam
if args.lam>0:
    print "Activating p3m"
    s.actors.add(DipolarP3M(prefactor=1,accuracy=1E-4,tune=True))
tick=time()
s.integrator.run(20)
print (time()-tick)

M =MagneticDipoleMoment(ids=s.part[:].id)
print M.calculate()
M_sat=args.n_nodes*dipm
print "skin",s.cell_system.tune_skin(min_skin=0.4,max_skin=3,tol=0.125,int_steps=20)
tick=time()
s.integrator.run(20)
print (time()-tick)

print "coldet"
s.collision_detection.set_params(
          mode="glue_to_surface",
          distance=cut_mixed,
          distance_glued_particle_to_vs=cut_pol, 
          bond_centers=bond_mag_pol,
          bond_vs=bond_pol,
          part_type_vs=2,
          part_type_to_attach_vs_to=0,
          part_type_to_be_glued=3,part_type_after_glueing=1)
    
tick=time()
s.integrator.run(20)
print (time()-tick)
print "skin",s.cell_system.tune_skin(min_skin=0.4,max_skin=2.5,tol=0.125,int_steps=20)
tick=time()
s.integrator.run(20)
print (time()-tick)

loops =10000
steps=100
for i in range(loops):
    tick=time()
    s.integrator.run(steps)
    E=s.analysis.energy()
    m=M.calculate()
    print time()-tick, s.analysis.energy()["total"],s.analysis.energy()["dipolar"],  np.linalg.norm(m)/M_sat,np.array(m)/M_sat
    n_vs=len(s.part)-args.n_nodes*(1+args.chains_per_node*args.chain_length)
    print "number of vs" ,n_vs
    cPickle.dump({
        "id": s.part[:].id,
        "pos": s.part[:].pos,
        "dip": s.part[:].dip, 
        "v": s.part[:].v, 
        "omega_lab": s.part[:].omega_lab, 
        "type": s.part[:].type, 
        "bonds": s.part[:].bonds,
        "virtual": s.part[:].virtual,
        "vs_relative": s.part[:].vs_relative,
        "gamma": s.part[:].gamma,
        "gamma_rot": s.part[:].gamma_rot
        },gzip.open(args.output_base+"_%0000d.pcl.gz" % (i,),"w"),-1)
    if n_vs==2*args.n_nodes*args.chains_per_node:
        print "Fully connected" 
        quit()




#





