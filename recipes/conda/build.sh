#!/bin/bash
set -euo pipefail

# Build the project using the Makefile
make INCLUDES="-I${PREFIX}/include" LDFLAGS="-L${PREFIX}/lib"

# the binaries are linked with an rpath relative to their own directory
# (@executable_path/lib on macOS, $ORIGIN/lib on Linux - see Makefile's
# R_PATH), so libhighs must land in bin/lib, not a top-level lib/.
mkdir -p "${PREFIX}/bin/lib"
cp hapscotch hapcure hictools seqtools hapcount "${PREFIX}/bin/"
cp -a lib/. "${PREFIX}/bin/lib/"
cp scripts/run_pipeline.py scripts/selfaln.py scripts/hicaln.py scripts/hicmap.py "${PREFIX}/bin/"
chmod +x "${PREFIX}"/bin/*.py