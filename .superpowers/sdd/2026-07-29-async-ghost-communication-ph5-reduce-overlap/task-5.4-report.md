# Task 5.4 Report: Non-zeroing CommBuf resize

## Scheme chosen

`DefaultInitAllocator<char>` — a thin `std::allocator<T>` wrapper that overrides
`construct(T*)` to use `::new (p) T` (default-initialization) instead of
`::new (p) T()` (value-initialization).  Used as the allocator for `CommBuf::buf`
(`ByteVec = std::vector<char, DefaultInitAllocator<char>>`).

**Why this over alternatives:**
- `reserve + manual m_size` tracker: breaks `std::span`, `begin()/end()`, and
  `make_span()` because `size()` would lie about capacity. Requires a separate
  `std::ptrdiff_t m_size` field and extra bookkeeping everywhere `buf.data()` is
  currently passed directly to MPI.
- `resize_and_overwrite` (C++23): not available yet (project targets C++20).
- Custom heap buffer (`std::unique_ptr<char[]>` + `m_size`/`m_cap`): much larger
  change, breaks `make_span()` and all `.data()`/`.size()` callers uniformly.
- `DefaultInitAllocator`: touches only the allocator typedef in `CommBuf`; the
  public API (`data()`, `size()`, `resize()`, `make_span()`, `bonds()`) is
  unchanged. The allocator's `allocate`/`deallocate` are inherited from
  `std::allocator<T>` unchanged.

## Written-before-read invariant

Documented in the `CommBuf` class comment (particle_packing.hpp) and in the
`resize()` docstring:

> Invariant: every byte in the non-bond buffer is written by the caller
> (pack_cells / MPI irecv) before it is read (unpack_cells / add_forces).
> Bytes added by growth are left uninitialized — the caller must write them
> before reading.

## CommBuf users audited

All callers that call `CommBuf::resize()` or directly access the underlying storage:

| Site | File | Role | Written before read? |
|------|------|------|----------------------|
| `pack_cells` | particle_packing.cpp:224 | send-side; `resize` then fills via `MemcpyOArchive` | YES — archiver writes every byte |
| `local_cell_copy` | particle_packing.cpp:326 | temp buffer per particle; resize then `MemcpyOArchive` | YES |
| `unpack_cells` | particle_packing.cpp:258 | read-only on buf already filled by MPI | N/A (not resizing here) |
| `add_forces` | particle_packing.cpp:294 | read-only on buf already filled by MPI | N/A |
| `add_rattle` | particle_packing.cpp:307 | read-only on buf already filled by MPI | N/A |
| `run_collective` (broadcast non-root) | HaloExchange.cpp:200 | `resize` then `boost::mpi::broadcast` fills buffer | YES — broadcast overwrites |
| `run_collective` (reduce root) | HaloExchange.cpp:216 | `resize` then `boost::mpi::reduce` fills recv_buf | YES — reduce overwrites |
| `halo_exchange_start` recv sizing | HaloExchange.cpp:320 | `resize` then `comm.irecv` targets `buf.data()` | YES — irecv fills on completion |
| `pack_regions` | HaloExchange.cpp:115 | delegates to `pack_cells` | covered above |

`make_span()` and `buf.data()`/`buf.size()` are used only to pass raw pointers to
`MemcpyOArchive`, `MemcpyIArchive`, or MPI — none of these read uninit bytes
without a prior write.

`bonds()` (a plain `std::vector<char>`) is NOT changed; it is serialized via
`boost::mpi` (binary archives) which requires exact byte semantics and is on the
cold resort-only path.

## Test evidence

- `make -C /ssd/weeber/comm-build -j8 check_unit_tests`: **149/149 passed**
- Parity gate (`ctest -R "^lj$|^nsquare$|^hybrid_decomposition$|^lees_edwards$"`):
  **4/4 passed**
- clang-format-19 applied; no diff after formatting.
