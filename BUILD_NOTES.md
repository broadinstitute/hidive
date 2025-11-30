# Build Notes for hidive

## System Dependencies

### sdsl-lite (vgteam fork)

The `odgi-ffi` dependency requires sdsl-lite, which has known compilation issues with newer Clang versions. The vgteam fork (https://github.com/vgteam/sdsl-lite) includes fixes for these issues.

**Installation:**

```bash
# Install the vgteam fork of sdsl-lite system-wide
# This can be done via Homebrew or by building from source:
git clone https://github.com/vgteam/sdsl-lite.git
cd sdsl-lite
mkdir build && cd build
cmake ..
make
sudo make install
```

### OpenMP (libomp)

OpenMP is required for sdsl-lite compilation.

**Installation on macOS:**

```bash
brew install libomp
```

### Build Configuration

When building, ensure the system-installed sdsl-lite and OpenMP are found by setting:

```bash
export CMAKE_PREFIX_PATH="/opt/homebrew/opt/libomp:$CMAKE_PREFIX_PATH"
export CFLAGS="-I/opt/homebrew/opt/libomp/include"
export CXXFLAGS="-I/opt/homebrew/opt/libomp/include"
export LDFLAGS="-L/opt/homebrew/opt/libomp/lib"
export CMAKE_ARGS="-Dsdsl-lite_ROOT=/usr/local/include -DOpenMP_C_INCLUDE_DIR=/opt/homebrew/opt/libomp/include -DOpenMP_CXX_INCLUDE_DIR=/opt/homebrew/opt/libomp/include -DOpenMP_C_LIB_NAMES=omp -DOpenMP_CXX_LIB_NAMES=omp -DOpenMP_omp_LIBRARY=/opt/homebrew/opt/libomp/lib/libomp.dylib -DCMAKE_C_FLAGS=-I/opt/homebrew/opt/libomp/include -DCMAKE_CXX_FLAGS=-I/opt/homebrew/opt/libomp/include"
```

Or add these to your shell profile for persistent configuration.

**Note:** The vgteam fork of sdsl-lite fixes the `louds_tree.hpp` compilation issues that occur with newer Clang versions. Make sure you have this version installed system-wide rather than relying on the bundled version in odgi-ffi.

## odgi-ffi Status

The `odgi-ffi` crate is **disabled** due to build issues with Clang 17+ and the `atomic_queue` dependency. Instead, hidive uses a modular approach:

### Design Decision

**Users prepare ODGI graphs externally** using the command-line `odgi` tool. This approach:
- Avoids build complexity and dependency issues
- Allows users to use their preferred graph construction tools
- Makes the workflow more modular and debuggable
- Enables reuse of graphs across multiple analyses

### Workflow

1. **Build pangenome graph**: `hidive build-pangenome` outputs a GFA file
2. **Convert to ODGI**: Users run `odgi build -g <gfa> -o <odgi>` (external tool)
3. **Use ODGI graph**: `hidive train-crf` and `hidive infer-haplotypes` accept ODGI files

See `PANGENOME_WORKFLOW.md` for detailed instructions.

### Dependencies

- **seqwish**: Included as Rust dependency, used by `build-pangenome`
- **odgi**: Command-line tool, must be installed separately by users

