#ifndef VOLE_HUGE_ALLOC_H__
#define VOLE_HUGE_ALLOC_H__

// Arrays backed by 2 MiB pages where the OS offers them.
//
// The LPN encoder (lpn.h) reads d random entries of a k-entry table per
// output. With 4 KiB pages a table larger than the TLB's reach (~12 MiB on
// current x86: e.g. FP192x2 at k = 2^19 is 25 MiB) turns nearly every read
// into a page walk; on an AMD Zen 5 box that doubled the p192 receiver's LPN
// time (180M dTLB misses per 3 rounds vs 0.7M with 2 MiB pages). On Linux the
// buffer is 2 MiB aligned and marked MADV_HUGEPAGE before first touch, which
// takes effect under the default THP mode ("madvise"); elsewhere, or if the
// hint is refused, it is an ordinary aligned allocation.

#include <cstddef>
#include <cstdlib>
#include <memory>
#include <new>
#include <type_traits>
#if defined(__linux__)
#include <sys/mman.h>
#endif

namespace vole {

constexpr std::size_t kHugePage = std::size_t{2} << 20;

// n default-constructed elements; release with huge_free(p).
template <typename T>
T *huge_new(std::size_t n) {
  static_assert(std::is_trivially_destructible_v<T>, "huge_free does not run destructors");
  std::size_t bytes = n * sizeof(T);
  std::size_t align = bytes >= kHugePage ? kHugePage : alignof(std::max_align_t);
  bytes = (bytes + align - 1) / align * align;
  if (bytes == 0) bytes = align;
  void *p = std::aligned_alloc(align, bytes);
  if (p == nullptr) throw std::bad_alloc();
#if defined(__linux__) && defined(MADV_HUGEPAGE)
  if (align == kHugePage) madvise(p, bytes, MADV_HUGEPAGE);  // a hint; failure is harmless
#endif
  T *t = static_cast<T *>(p);
  std::uninitialized_default_construct_n(t, n);
  return t;
}

template <typename T>
void huge_free(T *p) {
  std::free(p);
}

}  // namespace vole

#endif  // VOLE_HUGE_ALLOC_H__
