
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

parser =argparse.ArgumentParser()
parser.add_argument("--n_nodes",type=int,required=True)
parser.add_argument("--chains_per_node",type=int,required=True)
parser.add_argument("--chain_length",type=int,required=True)
parser.add_argument("--chain_length_width",type=int,required=True)
parser.add_argument("--seed",type=int,required=True)
parser.add_argument("--lambda",dest="lam",type=float,required=True)
parser.add_argument("--phi",type=float,required=True)
parser.add_argument("--output_base",type=str,required=True)
args=parser.parse_args()


cPickle.dump(args,open(args.output_base+"_args.pcl","w"))


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
s.time_step=0.008

s.thermostat.set_langevin(kT=kT,gamma=gamma)

# interactions
s.non_bonded_inter[0,0].lennard_jones.set_params(epsilon=eps,sigma=sigma_mag,cutoff=cut_mag,shift="auto")
s.non_bonded_inter[1,1].lennard_jones.set_params(epsilon=eps,sigma=sigma_pol,cutoff=cut_pol,shift="auto")
s.non_bonded_inter[0,1].lennard_jones.set_params(epsilon=eps,sigma=sigma_mixed,cutoff=cut_mixed,shift="auto")
bond_pol=HarmonicBond(r_0=cut_pol,k=1000)
s.bonded_inter.add(bond_pol)

# place magnetic particles
for i in range(args.n_nodes):
  cos_theta=2* np.random.random() -1
  sin_theta=np.sin(np.arccos(cos_theta))
  phi =2*np.pi*np.random.random()
  dip=dipm*np.array((np.sin(phi) *sin_theta,np.cos(phi)*sin_theta,cos_theta))
  assert abs(np.linalg.norm(dip) -dipm) <1E-10
  s.part.add(pos=random.random(3)*l,dip=dip,rotation=(1,1,1))




# warmup
s.integrator.set_steepest_descent(f_max=0,gamma=0.0001,max_displacement=0.01)
while s.analysis.energy()["total"] > 2*len(s.part):
  print( s.analysis.energy()["total"] )
  s.integrator.run(10)


# add chains
for i in range(args.chains_per_node*args.n_nodes):
    monomers=args.chain_length +args.chain_length_width*np.random.normal()
    if monomers <2: monomers=2
    print "Polymer",i,monomers,len(s.part)
    create_polymer(start_pos=np.random.random(3)*l,start_id=len(s.part),N_P=1,MPC=monomers,bond_length=cut_pol,bond=bond_pol,type_poly_neutral=1,mode=2)
s.part[args.n_nodes:].type=1
cPickle.dump({
    "id": s.part[:].id,
    "pos": s.part[:].pos,
    "dip": s.part[:].dip, 
    "v": s.part[:].v, 
    "omega_lab": s.part[:].omega_lab, 
    "type": s.part[:].type, 
    "bonds": s.part[:].bonds}, 
    open(args.output_base+"_pol.pcl","w"),-1)



print s.part[:].pos
print s.non_bonded_inter[1,1].lennard_jones.get_params()
print s.analysis.energy()


s.integrator.set_steepest_descent(f_max=0,gamma=0.001,max_displacement=0.01)
while s.analysis.energy()["total"] > 0.2*len(s.part):
  print( s.analysis.energy()["total"] )
  s.integrator.run(10,recalc_forces=True)
s.integrator.set_vv()
s.time_step=0.0001
s.integrator.run(2000)
s.time_step=0.008


tick=time()
s.integrator.run(100)
print (time()-tick)
print "skin",s.cell_system.tune_skin(min_skin=0.1,max_skin=2,tol=0.1,int_steps=50)
tick=time()
s.integrator.run(100)
print (time()-tick)

print "Polymer warmup"
for i in range(1000):
    s.integrator.run(1000)
    print i


print "lambda",args.lam
if args.lam>0:
    print "Activating p3m"
    s.cell_system.skin=1.5
    s.actors.add(DipolarP3M(prefactor=1,accuracy=1E-4,tune=True))
tick=time()
s.integrator.run(10)
print (time()-tick)

M =MagneticDipoleMoment(ids=s.part[:].id)
print M.calculate()
M_sat=args.n_nodes*dipm

loops =5000
steps=100
for i in range(loops):
    tick=time()
    s.integrator.run(steps)
    E=s.analysis.energy()
    m=M.calculate()
    print time()-tick, s.analysis.energy()["total"],s.analysis.energy()["dipolar"],  np.linalg.norm(m)/M_sat,np.array(m)/M_sat
    if i%200==0:
      print "skin",s.cell_system.tune_skin(min_skin=0.1,max_skin=2,tol=0.1,int_steps=100)
        

cPickle.dump({
    "id": s.part[:].id,
    "pos": s.part[:].pos,
    "dip": s.part[:].dip, 
    "v": s.part[:].v, 
    "omega_lab": s.part[:].omega_lab, 
    "type": s.part[:].type, 
    "bonds": s.part[:].bonds}, 
    open(args.output_base+"_equil.pcl","w"),-1)










