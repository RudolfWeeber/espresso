from numpy cimport int32_t

cdef extern from "vmdsock.hpp":
  int   vmdsock_init()
  void *vmdsock_create()
  int   vmdsock_bind(void *, int)
  int   vmdsock_listen(void *)
  void *vmdsock_accept(void *)
  int   vmdsock_connect(void *, const char *, int)
  int   vmdsock_write(void *, const void *, int)
  int   vmdsock_read(void *, void *, int)
  int   vmdsock_selread(void *, int)
  int   vmdsock_selwrite(void *, int)
  void  vmdsock_destroy(void *)


cdef extern from "imd.hpp":
   ctypedef enum IMDType: 
      IMD_DISCONNECT,
      IMD_ENERGIES, 
      IMD_FCOORDS,   
      IMD_GO,
      IMD_HANDSHAKE, 
      IMD_KILL,      
      IMD_MDCOMM,    
      IMD_PAUSE,
      IMD_TRATE,
      IMD_IOERROR
   
   int   imd_send_fcoords(void *, int32_t, const float *)
   int   imd_handshake(void *)
   IMDType imd_recv_header(void *, int32_t *)

cdef extern from "statistics_molecule.hpp":
    int analyze_fold_molecules(float *coord, double shift[3])

cdef extern from "grid.hpp":
  void fold_position(double pos[3],int image_box[3])


