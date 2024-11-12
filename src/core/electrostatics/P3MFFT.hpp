#pragma once\



#include "utils/Vector.hpp"
#include "heffte.h" 
#include "boost/mpi/communicator.hpp" 

template <typename FloatType>
class P3MFFT {
private:
   using backend_tag = heffte::backend::default_backend<heffte::tag::cpu>::type;
   using Box = heffte::box3d<>;
  Utils::Vector3i memory_layout;
  std::shared_ptr<Box> in_box;
  std::shared_ptr<Box> out_box;
  boost::mpi::communicator comm;
  Utils::Vector3i global_mesh;
  heffte::fft3d<backend_tag> fft3d;


public:
  P3MFFT(boost::mpi::communicator comm, 
     Utils::Vector3i global_mesh_size, 
     Utils::Vector3i rs_local_ld_index,
     Utils::Vector3i rs_local_ur_index,
     Utils::Vector3i _mem_layout) :
     comm(comm),
     memory_layout(_mem_layout),
     in_box(std::make_shared<Box>(rs_local_ld_index.as_array(), (rs_local_ur_index -Utils::Vector3i::broadcast(1)).as_array(),memory_layout.as_array())),
     out_box(std::make_shared<Box>(rs_local_ld_index.as_array(), (rs_local_ur_index-Utils::Vector3i::broadcast(1)).as_array()  ,memory_layout.as_array())),
     global_mesh(global_mesh_size),
     fft3d(*in_box,*out_box, comm) 

   {
     init_fft();
   };

   void set_preferred_kspace_decomposition() {
      int n = comm.size();
      Utils::Vector3i proc_grid{1,1,1};
      Utils::Vector3i local_size=global_mesh;
      for (int d: {2,1,0}) {
        while (n%2==0 and local_size[d]%2==0) {
          proc_grid[d]*=2;
          n/=2;
          local_size[d]/=2;
        }
      }
      if (proc_grid[1]==1 and proc_grid[0]==1 and proc_grid[2]%4==0 and local_size[1]%2==0) {
        proc_grid[2]/=2;
        proc_grid[1]*=2;
      }
      auto global_box = heffte::box3d<>({0,0,0}, {global_mesh[0]-1,global_mesh[1]-1,global_mesh[2]-1},mem_layout().as_array());
      auto all_boxes = heffte::split_world(global_box, proc_grid.as_array());
      auto new_box = all_boxes[comm.rank()];
      out_box=std::make_shared<Box>(new_box);
//      std::cout<<global_mesh<<std::endl;
//      std::cout<<Utils::Vector3i(out_box->low)<<"  | " <<Utils::Vector3i(out_box->high)<<std::endl;

//      std::cout << "reshaping to " <<proc_grid<<std::endl;
      init_fft();
   };


   void init_fft() {
    // at this stage we can manually adjust some HeFFTe options
    heffte::plan_options options = heffte::default_options<backend_tag>();

    // use strided 1-D FFT operations
    // some backends work just as well when the entries of the data are not contiguous
    // then there is no need to reorder the data in the intermediate stages which saves time
    options.use_reorder = false;

    // use point-to-point communications
    // collaborative all-to-all and individual point-to-point communications are two alternatives
    // one may be better than the other depending on
    // the version of MPI, the hardware interconnect, and the problem size
    options.algorithm = heffte::reshape_algorithm::p2p_plined;

    // in the intermediate steps, the data can be shapes as either 2-D slabs or 1-D pencils
    // for sufficiently large problem, it is expected that the pencil decomposition is better
    // but for smaller problems, the slabs may perform better (depending on hardware and backend)
   options.use_pencils = false;
   fft3d = heffte::fft3d<backend_tag>(*in_box,*out_box,comm,options);
  };

  Utils::Vector3i ks_local_ld_index() const {
     return Utils::Vector3i(out_box->low);
  };
  Utils::Vector3i ks_local_ur_index() const {
     return Utils::Vector3i(out_box->high)+Utils::Vector3i::broadcast(1);
  };
  Utils::Vector3i ks_local_size() const { 
    return ks_local_ur_index()-ks_local_ld_index();
  };
  template <typename T> 
  auto forward(T& in) { return fft3d.forward(in); };
  template <typename In,typename Out>
  void forward(In in, Out out) { fft3d.forward(in,out); };
  template <typename T> 
  auto backward(T& in) { return fft3d.backward(in); };
  template <typename T1, typename T2> 
  auto backward_batch(int n, T1 in, T2 out) { return fft3d.backward(n,in, out); };
  Utils::Vector3i mem_layout() const { return memory_layout; }

};
  
