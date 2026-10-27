#ifndef _ZLIBCOMPRESSOR_H_
#define _ZLIBCOMPRESSOR_H_

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>
#include <zlib.h>

#include "Parallel.h"

// Encode one VTK appended array with native-endian UInt64 headers.
// Each 32 KiB block is an independent zlib stream.
class ZLibCompressor {
public:
  static constexpr size_t vtkBlockSize = 32768;

  bool compress_vtk_block(const void* data, size_t sizeInBytes,
                          std::vector<uint8_t>& out) const {
    const size_t numBlocks =
        sizeInBytes / vtkBlockSize + (sizeInBytes % vtkBlockSize != 0);
    const size_t headerBytes = (3 + numBlocks) * sizeof(uint64_t);
    out.assign(headerBytes, 0);
    write_header(out, 0, numBlocks);
    write_header(out, 1, vtkBlockSize);
    // VTK uses zero when there is no partial final block.
    write_header(out, 2, sizeInBytes % vtkBlockSize);

    const auto* src = static_cast<const uint8_t*>(data);

    // The blocks are independent, so several threads compress them at the
    // same time. Thread t handles the contiguous blocks [t*per, (t+1)*per).
    // Small arrays are not worth the thread start-up cost.
    constexpr size_t minBlocksPerThread = 64;
    const size_t nThread = std::max<size_t>(
        1, std::min(converter_num_threads(), numBlocks / minBlocksPerThread));
    const size_t blocksPerThread = (numBlocks + nThread - 1) / nThread;

    std::vector<std::vector<uint8_t> > packed(nThread);
    std::vector<std::vector<uint64_t> > sizes(nThread);
    std::vector<char> hasFailed(nThread, 0);

    run_parallel(nThread, [&](size_t t) {
      const size_t bBegin = std::min(numBlocks, t * blocksPerThread);
      const size_t bEnd = std::min(numBlocks, bBegin + blocksPerThread);
      std::vector<uint8_t> compressed;
      sizes[t].reserve(bEnd - bBegin);
      for (size_t b = bBegin; b < bEnd; ++b) {
        const size_t offset = b * vtkBlockSize;
        const uLong chunkSize =
            static_cast<uLong>(std::min(vtkBlockSize, sizeInBytes - offset));
        uLongf compressedSize = ::compressBound(chunkSize);
        compressed.resize(compressedSize);
        const int status = ::compress(compressed.data(), &compressedSize,
                                      src + offset, chunkSize);
        if (status != Z_OK) {
          hasFailed[t] = 1;
          return;
        }
        sizes[t].push_back(compressedSize);
        packed[t].insert(packed[t].end(), compressed.begin(),
                         compressed.begin() + compressedSize);
      }
    });

    size_t iBlock = 0;
    for (size_t t = 0; t < nThread; ++t) {
      if (hasFailed[t]) {
        out.clear();
        return false;
      }
      for (const uint64_t size : sizes[t])
        write_header(out, 3 + iBlock++, size);
    }
    for (size_t t = 0; t < nThread; ++t)
      out.insert(out.end(), packed[t].begin(), packed[t].end());
    return true;
  }

private:
  static void write_header(std::vector<uint8_t>& out, size_t index,
                           uint64_t value) {
    // A byte vector need not provide alignment for uint64_t accesses.
    std::memcpy(out.data() + index * sizeof(value), &value, sizeof(value));
  }
};

#endif
