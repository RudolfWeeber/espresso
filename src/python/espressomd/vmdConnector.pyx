cimport vmdConnector
import system
from time import sleep
from utils cimport *
cimport numpy as np
from particle_data cimport *
from globals cimport *

cdef class VmdConnector:
  """Handle for an imd connection to the vmd visualizer"""

  cdef void *sock
  cdef void *initsock


  def __init__(self, S):
    """VmdConnector(EspressoSystem)"""

  def bind(self,port=None):
    """portUsed = bind(port=None)
      Bind to the given port. automatically choose one, if unspecified
      Throughs exception if unsuccessful
      """
    

    # Destroy existing connection
    self.unbind()

    vmdsock_init()
    self.initsock = vmdsock_create()
    p=port

    if p !=None:
      # Port was specified by user. Try to bind
      if vmdsock_bind(self.initsock,p):
          vmdsock_destroy(self.initsock)
          self.initsock = NULL
          raise Exception("Could not bind to given port %d. May be already in use?".format(port))
    else:
      # Automatically find a free port
      success=0
      for p in range(10000,20000,1):
        if vmdsock_bind(self.initsock,p):
          vmdsock_destroy(self.initsock)
          self.initsock =NULL
        else:
            success=1
            break
    
      if not success:
        raise Exception("Could not find a free port")
    
    if vmdsock_listen(self.initsock):
        vmdsock_destroy(self.initsock)
        self.initsock=NULL
        raise Exception("Could not listen on port %d. In use?".format(p))

    return p


  def unbind(self):
    """Disconnect and unbind from the socket"""
    # Destroy potential previous connection
    if (self.sock):
        vmdsock_destroy(self.sock)

    if (self.initsock):
        vmdsock_destroy(self.initsock)
    self.sock=NULL
    self.initsock=NULL

  def waitForVmd(self,seconds):
    """waitForVmd(seconds)
      Wait the give3n time for vmd to connect to the socket
      """
    cdef np.int32_t length
    success=0
    i=0
    while (self.initsock != NULL and (self.sock==NULL) and i<seconds):
      i+=1
      if vmdsock_selread(self.initsock, 0) > 0:
        self.sock = vmdsock_accept(self.initsock)
        if self.sock==NULL:
          raise Exception("Socket not valid after vmdsock_accept()")
        if (imd_handshake(self.sock)):
          self.unbind()
          raise Exception("IMD handshake failed. Wrong VMD version ?")
        success=1
        break
      sleep(1)
    
    if not success:
        raise Exception("Vmd did not answer within specified waiting time")
    
    if self.sock==NULL:
          raise Exception("Socket not valid after handshake")
  
#    if (vmdsock_selread(self.sock, 0) != 1):
#      raise Exception("Read from seocket failed")
    if (imd_recv_header(self.sock, &length) != IMD_GO):
        self.unbind()
        raise Exception("No go from VMD. Wrong VMD version ?")
  
    sleep(1)
    if self.sock!=NULL:
          return
  
    raise Exception("Vmd did not connect.")


  def update(self,flag=None):
    """Send current positions to vmd"""

    cdef double shift[3] 
    shift = (0.0,0.0,0.0)
    
    cdef float *coord
    cdef int i, j
  
    # !!! TODO: Implement option for unfolded coordinates and unfolded chains (src/tcl/imd_tcl.cpp tclcommand_parse_imd_pos)
    
    # Connection checking
    if not self.initsock:
      raise Exception("No connection")
    
    if not self.sock:
        raise Exception("no connection")
    
    if not vmdsock_selwrite(self.sock, 60) > 0:
        raise Exception("could not write to IMD socket.")
  
    # Consecutive particle numbering
    if n_part != max_seen_particle + 1: 
      raise Exception("for IMD, store particles consecutively starting with 0.")
    
    # Copy particles to node 0 and pack coordinates
    updatePartCfg(1) # WIth bonds
    coord = <float*>malloc(n_part*3*sizeof(float))
    # sort partcles according to identities */
    cdef int dummy[3] 
    cdef double tmpCoord[3]
    for i in range(n_part): 
      dummy=(0,0,0)
      tmpCoord[0] = partCfg[i].r.p[0]
      tmpCoord[1] = partCfg[i].r.p[1]
      tmpCoord[2] = partCfg[i].r.p[2]
      
      # Fold
      if flag == None:    
        fold_position(tmpCoord, dummy)
      j = 3*partCfg[i].p.identity
      coord[j    ] = tmpCoord[0]
      coord[j + 1] = tmpCoord[1]
      coord[j + 2] = tmpCoord[2]
  


    # Use information from the analyse set command to fold chain molecules
    if flag == "FOLD_CHAINS":
        if analyze_fold_molecules(coord, shift ) !=0:
          raise Exception("could not fold chain.")
  
    # Send the stuff 
    if imd_send_fcoords(self.sock, n_part, coord): 
      raise Exception("could not write to IMD socket.")
    
    # Free mem
    free(coord)
    
  
  
      
  
  
  
       
