#pragma once

namespace walberla {
template <typename FT, typename PdfStencil, lbmpy::Arch AT = lbmpy::Arch::CPU>
struct FieldTrait {
  using PdfField = field::GhostLayerField<FT, PdfStencil::Size>;
  using VectorField = field::GhostLayerField<FT, uint_t{3u}>;
  template <class Field> using PackInfo = field::communication::PackInfo<Field>;
  using PackInfoStreamingPdf = detail::KernelTrait<FT, AT>::PackInfoPdf;
  using PackInfoStreamingVec = detail::KernelTrait<FT, AT>::PackInfoVec;
  template <class Stencil>
  using RegularCommScheme =
      blockforest::communication::UniformBufferedScheme<Stencil>;
  template <class Stencil>
  using BoundaryCommScheme =
      blockforest::communication::UniformBufferedScheme<Stencil>;
};

#if defined(__CUDACC__)
template <typename FT, typename PdfStencil>
struct FieldTrait<FT, PdfStencil, lbmpy::Arch::GPU> {
private:
  static auto constexpr AT = lbmpy::Arch::GPU;
  template <class Field>
  using MemcpyPackInfo = gpu::communication::MemcpyPackInfo<Field>;

public:
  template <typename Stencil>
  class UniformGPUScheme
      : public gpu::communication::UniformGPUScheme<Stencil> {
  public:
    explicit UniformGPUScheme(auto const &bf)
        : gpu::communication::UniformGPUScheme<Stencil>(
              bf, /* sendDirectlyFromGPU */ false,
              /* useLocalCommunication */ false) {}
  };
  using PdfField = gpu::GPUField<FT>;
  using VectorField = gpu::GPUField<FT>;
  template <class Field> using PackInfo = MemcpyPackInfo<Field>;
  using PackInfoStreamingPdf = detail::KernelTrait<FT, AT>::PackInfoPdf;
  using PackInfoStreamingVec = detail::KernelTrait<FT, AT>::PackInfoVec;
  template <class Stencil> using RegularCommScheme = UniformGPUScheme<Stencil>;
  template <class Stencil>
  using BoundaryCommScheme =
      blockforest::communication::UniformBufferedScheme<Stencil>;
};
#endif
} // namespace walberla
