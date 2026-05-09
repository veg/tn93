# Release Notes: tn93 v1.0.16

This release introduces major performance optimizations, new screening features, and a robust automated testing infrastructure. It represents a significant step forward in scaling `tn93` to handle modern large-scale genomic datasets.

## Key Changes Since v1.0.15

### 1. Performance & Parallelism
*   **4-Way Accumulator Unrolling**: The core nucleotide counting loop in `computeTN93` has been refactored to use four independent accumulator arrays. This allows the CPU to execute multiple increments in parallel, significantly increasing throughput and achieving up to a **2x speedup** in the hot loop.
*   **Interleaved OpenMP Scheduling**: Switched the primary pairwise loop to `schedule(static, 1)`. This interleaving strategy perfectly balances the triangular workload of pairwise comparisons, ensuring that heavy early iterations and light late iterations are evenly distributed across all available CPU cores.
*   **Reduced Synchronization Contention**: Replaced global `critical` sections with `atomic` counter updates for high-frequency metrics (`pairIndex`, `foundLinks`). Independent tasks now use named critical blocks (e.g., `stats`, `fwrite`, `progress`) to prevent threads from unnecessarily blocking each other.
*   **Optimized Progress Reporting**: Minimized I/O overhead by making progress updates conditional on meaningful completion increments ($>0.1\%$).

### 2. New Features
*   **Hamming Skip (-H)**: Introduced a high-speed early-exit heuristic for rapid screening. The implementation includes tuned parameters (`WORTH_DOING=1024`, `PERIODIC_CHECK=128`) to minimize overhead while maximizing detection of highly divergent sequences.
*   **High-Performance Output Buffering**: Implemented thread-local buffering for large CSV exports. This ensures that writing millions of pairwise distances does not become a bottleneck or cause excessive lock contention.
*   **Enhanced Sequence Metadata**: Sequence descriptors now track `total_gaps` and `total_ambigs`, providing deeper insights into sequence quality and alignment characteristics.

### 3. Testing & Validation
*   **Comprehensive CTest Suite**: Added a full suite of 11 regression tests covering small and large datasets, various ambiguity resolution modes, and optimization flags.
*   **Python-based Validation Runner**: A new robust test runner that validates result identity (up to row permutation) against reference implementation baselines with configurable floating-point tolerance.
*   **Automated CI/CD**: Added a GitHub Actions workflow (`.github/workflows/tests.yml`) to automatically build and validate the codebase on every push to major branches.

### 4. Bug Fixes & Maintenance
*   **Gap-Run Skipping**: Refined the `GET_JUMP` logic to correctly handle sparse internal gaps in genomic alignments.
*   **ReadReduce Consistency**: Fixed edge cases in sequence reduction to ensure deterministic results.
*   **Thread Safety**: Hardened shared data structures and output streams for high-concurrency environments.
*   **Code Quality**: Cleaned up signed/unsigned character comparisons and improved register-level optimizations throughout the shared library.

---
**Build Requirements**: CMake 3.5+, Python 3 (for tests), OpenMP-capable compiler.
**Target Hardware**: Optimized for modern multi-core systems, with specific tuning for Apple Silicon (M-series) and high-core-count servers.
