# Release Notes: tn93 v1.0.16

This release introduces major performance optimizations, new screening features, and a robust automated testing infrastructure. It represents a significant step forward in scaling `tn93` to handle modern large-scale genomic datasets.

## Key Changes Since v1.0.15

### 1. Performance & Parallelism
*   **4-Way Accumulator Unrolling**: The core nucleotide counting loop in `computeTN93` has been refactored to use four independent accumulator arrays. This allows the CPU to execute multiple increments in parallel, significantly increasing throughput and achieving up to a **2x speedup** in the hot loop.
*   **Advanced Ambiguity Batching**: Ambiguity resolutions are stashed into thread-local buffers and resolved in bulk after the primary counting loop. This strategy prevents branch-heavy resolution logic from disrupting the CPU's instruction pipeline and code flow during the hot counting phase.
*   **Contiguous Resolved Run Optimization**: Implemented a specialized processing path for sequence pairs that share a contiguous gap- and ambiguity-free region of at least 1024bp. These regions are handled separately using the 4-way unrolled loop to maximize instruction-level parallelism.
*   **Dynamic OpenMP Scheduling**: Optimized the primary pairwise loop using `schedule(dynamic, 1)`. This strategy provides superior load balancing compared to interleaved static scheduling by allowing threads to adapt to the stochastic variability in execution times introduced by the Hamming skip early-exit heuristic.
*   **Unrolled Early-Exit Validation**: The `check_early_exit` loop has been completely unrolled to sum transition counts across all four independent accumulators directly. This minimizes branching overhead in the hot loop's periodic status checks.
*   **Performance Benchmarks**: These cumulative optimizations achieve up to a **1.6x speedup** on large-scale datasets compared to the previous stable release, while maintaining bit-perfect result parity.
*   **Refined Gap-Run Skipping**: Enhanced the `GET_JUMP` mechanism to synchronize skipping across both sequences. The loop pointer advances by the maximum jump length available, effectively bypassing contiguous gap regions without per-nucleotide branching or state disruptions.
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
