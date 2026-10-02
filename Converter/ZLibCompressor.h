#ifndef _ZLIBCOMPRESSOR_H_
#define _ZLIBCOMPRESSOR_H_

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>
#include <zlib.h>

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
    std::vector<uint8_t> compressed;
    for (size_t b = 0; b < numBlocks; ++b) {
      const size_t offset = b * vtkBlockSize;
      const uLong chunkSize =
          static_cast<uLong>(std::min(vtkBlockSize, sizeInBytes - offset));
      uLongf compressedSize = ::compressBound(chunkSize);
      compressed.resize(compressedSize);
      const int status = ::compress(compressed.data(), &compressedSize,
                                    src + offset, chunkSize);
      if (status != Z_OK) {
        out.clear();
        return false;
      }
      write_header(out, 3 + b, compressedSize);
      compressed.resize(compressedSize);
      out.insert(out.end(), compressed.begin(), compressed.end());
    }
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
