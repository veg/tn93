#include "subquadratic.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace argparse;

#define HISTOGRAM_BINS 200
#define HISTOGRAM_SLICE ((double)HISTOGRAM_BINS)
#define BUFFER_FLUSH_LENGTH 8192

struct BitPlanes
{
  vector<uint64_t> A, C, G, T, Gap;
  uint32_t first_nongap, last_nongap;
  uint32_t total_nongap;
  uint32_t ambig_count;
};

static void
build_bitplanes (const char *seq, unsigned long L, BitPlanes &bp)
{
  unsigned long num_words = (L + 63) / 64;
  bp.A.assign (num_words, 0ULL);
  bp.C.assign (num_words, 0ULL);
  bp.G.assign (num_words, 0ULL);
  bp.T.assign (num_words, 0ULL);
  bp.Gap.assign (num_words, 0ULL);
  bp.first_nongap = L;
  bp.last_nongap = 0;
  bp.total_nongap = 0;
  bp.ambig_count = 0;

  for (unsigned long i = 0; i < L; i++)
    {
      unsigned char c = (unsigned char)seq[i];
      unsigned long word = i / 64;
      uint64_t bit = 1ULL << (i % 64);

      if (IS_GAP (c, GAP))
        {
          bp.Gap[word] |= bit;
        }
      else
        {
          if (i < bp.first_nongap)
            bp.first_nongap = i;
          if (i > bp.last_nongap)
            bp.last_nongap = i;
          bp.total_nongap++;

          if (c < 4)
            {
              if (c == 0)
                bp.A[word] |= bit;
              else if (c == 1)
                bp.C[word] |= bit;
              else if (c == 2)
                bp.G[word] |= bit;
              else if (c == 3)
                bp.T[word] |= bit;
            }
          else if (c < GAP)
            {
              bp.ambig_count++;
              const long *res = resolve_char (c, false, false);
              if (res[0])
                bp.A[word] |= bit;
              if (res[1])
                bp.C[word] |= bit;
              if (res[2])
                bp.G[word] |= bit;
              if (res[3])
                bp.T[word] |= bit;
            }
        }
    }
}

static inline bool
quick_bitplane_filter (const BitPlanes &p1, const BitPlanes &p2,
                       double threshold, unsigned long min_overlap)
{
  unsigned long fn = std::max (p1.first_nongap, p2.first_nongap);
  unsigned long ln = std::min (p1.last_nongap, p2.last_nongap);
  if (fn > ln || (ln - fn + 1) < min_overlap)
    return false;

  unsigned long max_possible_overlap = ln - fn + 1;
  uint32_t max_allowed_mismatches
      = (uint32_t)std::floor (threshold * max_possible_overlap);

  unsigned long start_w = fn / 64;
  unsigned long end_w = ln / 64;

  uint32_t total_valid = 0;
  uint32_t total_mismatch = 0;

  const uint64_t *const __restrict__ a1 = p1.A.data ();
  const uint64_t *const __restrict__ c1 = p1.C.data ();
  const uint64_t *const __restrict__ g1 = p1.G.data ();
  const uint64_t *const __restrict__ t1 = p1.T.data ();
  const uint64_t *const __restrict__ gap1 = p1.Gap.data ();

  const uint64_t *const __restrict__ a2 = p2.A.data ();
  const uint64_t *const __restrict__ c2 = p2.C.data ();
  const uint64_t *const __restrict__ g2 = p2.G.data ();
  const uint64_t *const __restrict__ t2 = p2.T.data ();
  const uint64_t *const __restrict__ gap2 = p2.Gap.data ();

  for (unsigned long w = start_w; w <= end_w; w++)
    {
      uint64_t valid = ~(gap1[w] | gap2[w]);
      if (w == start_w)
        {
          valid &= ~((1ULL << (fn % 64)) - 1);
        }
      if (w == end_w && (ln % 64) != 63)
        {
          valid &= (1ULL << ((ln % 64) + 1)) - 1;
        }

      uint64_t match = (a1[w] & a2[w]) | (c1[w] & c2[w])
                       | (g1[w] & g2[w]) | (t1[w] & t2[w]);
      uint64_t mismatch = valid & ~match;

      total_mismatch += __builtin_popcountll (mismatch);
      if (total_mismatch > max_allowed_mismatches)
        return false;

      total_valid += __builtin_popcountll (valid);
    }

  if (total_valid < min_overlap)
    return false;
  return (double)total_mismatch <= threshold * total_valid;
}

static void
dump_subq_histogram (ostream *outStream, const char *tag, unsigned long *hist)
{
  if (tag)
    {
      (*outStream) << "\t\"Histogram " << tag << "\" : [";
    }
  else
    {
      (*outStream) << "\t\"Histogram\" : [";
    }
  for (unsigned long k = 0; k < HISTOGRAM_BINS; k++)
    {
      if (k)
        {
          (*outStream) << ',';
        }
      (*outStream) << '[' << (k + 1.) / HISTOGRAM_BINS << ',' << hist[k]
                   << ']';
    }
  (*outStream) << "]" << endl;
}

static inline unsigned char
canonical_base (unsigned char c, unsigned char cons_base)
{
  if (c < 4)
    return c;
  if (IS_GAP (c, GAP))
    return GAP;
  const long *res = resolve_char (c, false, false);
  if (res[cons_base])
    return cons_base;
  for (unsigned char b = 0; b < 4; b++)
    {
      if (res[b])
        return b;
    }
  return 0;
}

static vector<unsigned char>
compute_consensus (const StringBuffer &sequences, const Vector &seqLengths,
                   unsigned long L, unsigned long N)
{
  vector<vector<unsigned int> > counts (L, vector<unsigned int> (4, 0));
  for (unsigned long sid = 0; sid < N; sid++)
    {
      const char *s = stringText (sequences, seqLengths, sid);
      for (unsigned long i = 0; i < L; i++)
        {
          unsigned char c = (unsigned char)s[i];
          if (c < 4)
            {
              counts[i][c]++;
            }
        }
    }
  vector<unsigned char> consensus (L, 0);
  for (unsigned long i = 0; i < L; i++)
    {
      unsigned int max_c = 0;
      unsigned char best_base = 0;
      for (unsigned char b = 0; b < 4; b++)
        {
          if (counts[i][b] > max_c)
            {
              max_c = counts[i][b];
              best_base = b;
            }
        }
      consensus[i] = best_base;
    }
  return consensus;
}

struct LshEntry
{
  uint64_t hash;
  uint32_t seq_id;
};

int
run_subquadratic_tn93 (argparse::args_t &args, StringBuffer &sequences,
                       Vector &seqLengths, StringBuffer &names,
                       Vector &nameLengths,
                       sequence_gap_structure *sequence_descriptors,
                       long firstSequenceLength, Vector &counts,
                       int resolutionOption, unsigned long seqLengthInFile1,
                       unsigned long seqLengthInFile2)
{
  unsigned long sequenceCount = seqLengths.length () - 1;
  bool cross_comparison_only = (args.input2 != NULL);

  unsigned long pairwise
      = cross_comparison_only
            ? seqLengthInFile1 * seqLengthInFile2
            : (sequenceCount - 1) * (sequenceCount) / 2;

  long upperBound = cross_comparison_only ? seqLengthInFile1 : sequenceCount;

  double *distanceMatrix = NULL;
  if (args.format == hyphy && !args.do_count)
    {
      distanceMatrix = new double[sequenceCount * sequenceCount];
      for (unsigned long i = 0; i < sequenceCount * sequenceCount; i++)
        distanceMatrix[i] = 100.;
      for (unsigned long i = 0; i < sequenceCount; i++)
        distanceMatrix[i * sequenceCount + i] = 0.;
    }

  // Precompute bitplanes
  vector<BitPlanes> bitplanes (sequenceCount);
#pragma omp parallel for schedule(static)
  for (unsigned long sid = 0; sid < sequenceCount; sid++)
    {
      build_bitplanes (stringText (sequences, seqLengths, sid),
                       firstSequenceLength, bitplanes[sid]);
    }

  bool report_self
      = (args.input2 == NULL && args.report_self && !args.do_count);

  long foundLinks = 0;
  long actualEvals = 0;
  long pairIndex = 0;
  double percentDone = 0.0;
  double global_max_d = 0.0;
  double global_sum_d = 0.0;
  double global_weighted_links = 0.0;

  unsigned long global_hist[HISTOGRAM_BINS];
  for (unsigned long k = 0; k < HISTOGRAM_BINS; k++)
    global_hist[k] = 0;

  time_t before, after;
  time (&before);

  const unsigned long LSH_THRESHOLD = 15000;
  bool use_lsh = (!cross_comparison_only && sequenceCount > LSH_THRESHOLD);

  if (!use_lsh)
    {
      // Direct Pairwise Bitplane Screening (Subquadratic empirical time via SIMD early exit)
#pragma omp parallel shared(                                                  \
    bitplanes, sequence_descriptors, resolutionOption, foundLinks,            \
        actualEvals, pairIndex, sequences, seqLengths, sequenceCount,         \
        firstSequenceLength, args, nameLengths, names, pairwise, percentDone, \
        cerr, global_max_d, global_sum_d, global_weighted_links,              \
        distanceMatrix, upperBound, seqLengthInFile1, seqLengthInFile2,       \
        cross_comparison_only, report_self)
      {
        StringBuffer local_buffer;
        unsigned long local_hist[HISTOGRAM_BINS];
        for (unsigned long k = 0; k < HISTOGRAM_BINS; k++)
          local_hist[k] = 0;

        long local_links = 0;
        long local_evals = 0;
        double local_max_d = 0.0;
        double local_sum_d = 0.0;
        double local_weighted = 0.0;

#pragma omp for schedule(dynamic, 1)
        for (long sid1 = 0; sid1 < upperBound; sid1++)
          {
            char *n1 = stringText (names, nameLengths, sid1);
            const unsigned long n1L = stringLength (nameLengths, sid1);
            const char *s1 = stringText (sequences, seqLengths, sid1);
            long instances1 = counts.value (sid1);

            if (report_self)
              {
                if (args.format == csv)
                  {
                    char float_buf[128];
                    unsigned written = snprintf (float_buf, 128, "%g", 0.0);
                    local_buffer.appendBuffer (n1, n1L);
                    local_buffer.appendChar (args.delimiter);
                    local_buffer.appendBuffer (n1, n1L);
                    local_buffer.appendChar (args.delimiter);
                    local_buffer.appendBuffer (float_buf, written);
                    local_buffer.appendChar ('\n');
                  }
                else if (args.format == csvn)
                  {
                    char link_buf[1024];
                    int written = snprintf (link_buf, 1024, "%ld%c%ld%c%g\n",
                                            sid1, args.delimiter, sid1,
                                            args.delimiter, 0.0);
                    local_buffer.appendBuffer (link_buf, written);
                  }
              }

            long lowerBound = cross_comparison_only ? seqLengthInFile1 : sid1 + 1;

            for (unsigned long sid2 = lowerBound; sid2 < sequenceCount; sid2++)
              {
                if (!quick_bitplane_filter (bitplanes[sid1], bitplanes[sid2],
                                            args.distance, args.overlap))
                  {
                    continue;
                  }

                local_evals++;
                long weighted_count = instances1 * counts.value (sid2);
                const char *s2 = stringText (sequences, seqLengths, sid2);

                double thisD = sequence_descriptors
                                   ? computeTN93 (
                                         s1, s2, firstSequenceLength,
                                         resolutionOption, NULL, args.overlap,
                                         NULL, HISTOGRAM_SLICE, HISTOGRAM_BINS,
                                         weighted_count, 1L,
                                         &sequence_descriptors[sid1],
                                         &sequence_descriptors[sid2],
                                         args.distance)
                                   : computeTN93 (
                                         s1, s2, firstSequenceLength,
                                         resolutionOption, NULL, args.overlap,
                                         NULL, HISTOGRAM_SLICE, HISTOGRAM_BINS,
                                         weighted_count, 1L, NULL, NULL,
                                         args.distance);

                if (thisD >= args.min_distance && thisD <= args.distance)
                  {
                    local_links += weighted_count;
                    local_sum_d += thisD * weighted_count;
                    local_weighted += weighted_count;
                    if (thisD > local_max_d)
                      local_max_d = thisD;

                    unsigned long bin
                        = (unsigned long)(thisD * HISTOGRAM_SLICE);
                    if (bin >= HISTOGRAM_BINS)
                      bin = HISTOGRAM_BINS - 1;
                    local_hist[bin] += weighted_count;

                    if (!args.do_count)
                      {
                        if (args.format == csv)
                          {
                            char float_buf[128];
                            unsigned written
                                = snprintf (float_buf, 128, "%g", thisD);
                            local_buffer.appendBuffer (n1, n1L);
                            local_buffer.appendChar (args.delimiter);
                            local_buffer.appendBuffer (
                                stringText (names, nameLengths, sid2),
                                stringLength (nameLengths, sid2));
                            local_buffer.appendChar (args.delimiter);
                            local_buffer.appendBuffer (float_buf, written);
                            local_buffer.appendChar ('\n');
                          }
                        else if (args.format == csvn)
                          {
                            char link_buf[1024];
                            int written = snprintf (
                                link_buf, 1024, "%ld%c%ld%c%g\n", sid1,
                                args.delimiter, sid2, args.delimiter, thisD);
                            local_buffer.appendBuffer (link_buf, written);
                          }
                        else if (distanceMatrix)
                          {
#pragma omp critical(dist_matrix)
                            {
                              distanceMatrix[sid1 * sequenceCount + sid2]
                                  = thisD;
                              distanceMatrix[sid2 * sequenceCount + sid1]
                                  = thisD;
                            }
                          }

                        if (local_buffer.length () > BUFFER_FLUSH_LENGTH)
                          {
#pragma omp critical(fwrite)
                            {
                              fwrite (local_buffer.getString (), sizeof (char),
                                      local_buffer.length (), args.output);
                            }
                            local_buffer.resetString ();
                          }
                      }
                  }
              }

            long current_pair_index;
#pragma omp atomic capture
            current_pair_index = pairIndex += cross_comparison_only
                                                  ? seqLengthInFile2
                                                  : (sequenceCount - sid1 - 1);

            if (!args.quiet
                && (current_pair_index * 100. / pairwise - percentDone > 0.1
                    || sid1 == (long)upperBound - 1))
              {
#pragma omp critical(progress)
                {
                  if (current_pair_index * 100. / pairwise - percentDone > 0.1
                      || sid1 == (long)upperBound - 1)
                    {
                      time (&after);
                      percentDone = current_pair_index * 100. / pairwise;
                      cerr << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
                              "\b\b"
                              "\b\b\b\b\b"
                              "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
                              "\b\b"
                              "\b\b\b\b\b"
                              "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
                              "\b\b"
                              "\bProgress"
                              ":"
                           << setw (8) << percentDone << "% (" << setw (8)
                           << (foundLinks + local_links) << " links found, "
                           << setw (12) << std::setprecision (3)
                           << current_pair_index / std::max (1.0, difftime (after, before))
                           << " evals/sec)";
                      after = before;
                    }
                }
              }
          }

        if (local_buffer.length () > 0)
          {
#pragma omp critical(fwrite)
            {
              fwrite (local_buffer.getString (), sizeof (char),
                      local_buffer.length (), args.output);
            }
          }

#pragma omp critical(hist)
        {
          foundLinks += local_links;
          actualEvals += local_evals;
          global_sum_d += local_sum_d;
          global_weighted_links += local_weighted;
          if (local_max_d > global_max_d)
            global_max_d = local_max_d;
          for (unsigned long k = 0; k < HISTOGRAM_BINS; k++)
            global_hist[k] += local_hist[k];
        }
      }
    }
  else
    {
      // Windowed LSH + Outlier Fallback for massive alignments (N > 15,000)
      vector<unsigned char> consensus
          = compute_consensus (sequences, seqLengths, firstSequenceLength,
                               sequenceCount);

      // Outlier identification
      vector<uint32_t> fallback_seqs;
      vector<bool> is_fallback (sequenceCount, false);
      for (unsigned long sid = 0; sid < sequenceCount; sid++)
        {
          unsigned long span = (bitplanes[sid].last_nongap
                                >= bitplanes[sid].first_nongap)
                                   ? (bitplanes[sid].last_nongap
                                      - bitplanes[sid].first_nongap + 1)
                                   : 0;
          if (span < args.overlap || bitplanes[sid].ambig_count > span * 0.05)
            {
              fallback_seqs.push_back (sid);
              is_fallback[sid] = true;
            }
        }

      const unsigned long W = 120;
      const unsigned long S = 60;
      vector<pair<unsigned long, unsigned long> > tiles;
      for (unsigned long start = 0;
           start + W <= (unsigned long)firstSequenceLength; start += S)
        {
          tiles.push_back ({ start, start + W });
        }
      if (tiles.empty ()
          || tiles.back ().second < (unsigned long)firstSequenceLength)
        {
          unsigned long s = (firstSequenceLength > (long)W)
                                ? firstSequenceLength - W
                                : 0;
          tiles.push_back ({ s, (unsigned long)firstSequenceLength });
        }

      const int k_samples = 14;
      const int b_bands = 8;
      struct TileBand
      {
        vector<unsigned long> sample_cols;
      };

      vector<TileBand> all_bands;
      mt19937 rng (42);
      for (const auto &tile : tiles)
        {
          vector<unsigned long> cols;
          for (unsigned long p = tile.first; p < tile.second; p++)
            cols.push_back (p);
          if (cols.size () < (size_t)k_samples)
            continue;

          for (int b = 0; b < b_bands; b++)
            {
              TileBand tb;
              vector<unsigned long> pool = cols;
              shuffle (pool.begin (), pool.end (), rng);
              tb.sample_cols.assign (pool.begin (), pool.begin () + k_samples);
              all_bands.push_back (tb);
            }
        }

      vector<vector<LshEntry> > tables (all_bands.size ());
#pragma omp parallel for schedule(dynamic, 1)
      for (size_t b_idx = 0; b_idx < all_bands.size (); b_idx++)
        {
          const auto &tb = all_bands[b_idx];
          auto &table = tables[b_idx];

          for (unsigned long sid = 0; sid < sequenceCount; sid++)
            {
              if (is_fallback[sid])
                continue;
              const char *s = stringText (sequences, seqLengths, sid);

              bool has_gap = false;
              uint64_t h = 0;
              for (int k = 0; k < k_samples; k++)
                {
                  unsigned long pos = tb.sample_cols[k];
                  unsigned char c = (unsigned char)s[pos];
                  if (IS_GAP (c, GAP))
                    {
                      has_gap = true;
                      break;
                    }
                  unsigned char cb = canonical_base (c, consensus[pos]);
                  h = (h << 2) | (cb & 3);
                }
              if (!has_gap)
                table.push_back ({ h, (uint32_t)sid });
            }

          std::sort (table.begin (), table.end (),
                     [](const LshEntry &a, const LshEntry &b) {
                       return a.hash < b.hash;
                     });
        }

#pragma omp parallel shared(                                                  \
    tables, all_bands, is_fallback, fallback_seqs, bitplanes,                 \
        sequence_descriptors, resolutionOption, foundLinks, actualEvals,      \
        pairIndex, sequences, seqLengths, sequenceCount,                      \
        firstSequenceLength, args, nameLengths, names, pairwise, percentDone, \
        cerr, global_max_d, global_sum_d, global_weighted_links,              \
        distanceMatrix, upperBound, report_self)
      {
        StringBuffer local_buffer;
        vector<uint32_t> candidates;
        vector<uint32_t> seen (sequenceCount, 0);
        uint32_t epoch = 0;
        unsigned long local_hist[HISTOGRAM_BINS];
        for (unsigned long k = 0; k < HISTOGRAM_BINS; k++)
          local_hist[k] = 0;

        long local_links = 0;
        long local_evals = 0;
        double local_max_d = 0.0;
        double local_sum_d = 0.0;
        double local_weighted = 0.0;

#pragma omp for schedule(dynamic, 1)
        for (long sid1 = 0; sid1 < upperBound; sid1++)
          {
            char *n1 = stringText (names, nameLengths, sid1);
            const unsigned long n1L = stringLength (nameLengths, sid1);
            const char *s1 = stringText (sequences, seqLengths, sid1);
            long instances1 = counts.value (sid1);

            if (report_self)
              {
                if (args.format == csv)
                  {
                    char float_buf[128];
                    unsigned written = snprintf (float_buf, 128, "%g", 0.0);
                    local_buffer.appendBuffer (n1, n1L);
                    local_buffer.appendChar (args.delimiter);
                    local_buffer.appendBuffer (n1, n1L);
                    local_buffer.appendChar (args.delimiter);
                    local_buffer.appendBuffer (float_buf, written);
                    local_buffer.appendChar ('\n');
                  }
                else if (args.format == csvn)
                  {
                    char link_buf[1024];
                    int written = snprintf (link_buf, 1024, "%ld%c%ld%c%g\n",
                                            sid1, args.delimiter, sid1,
                                            args.delimiter, 0.0);
                    local_buffer.appendBuffer (link_buf, written);
                  }
              }

            epoch++;
            if (epoch == 0)
              {
                std::fill (seen.begin (), seen.end (), 0);
                epoch = 1;
              }

            candidates.clear ();
            if (!is_fallback[sid1])
              {
                for (size_t b_idx = 0; b_idx < all_bands.size (); b_idx++)
                  {
                    const auto &tb = all_bands[b_idx];
                    const auto &table = tables[b_idx];

                    bool has_gap = false;
                    uint64_t h = 0;
                    for (int k = 0; k < k_samples; k++)
                      {
                        unsigned long pos = tb.sample_cols[k];
                        unsigned char c = (unsigned char)s1[pos];
                        if (IS_GAP (c, GAP))
                          {
                            has_gap = true;
                            break;
                          }
                        unsigned char cb = canonical_base (c, consensus[pos]);
                        h = (h << 2) | (cb & 3);
                      }
                    if (has_gap)
                      continue;

                    auto it = std::lower_bound (
                        table.begin (), table.end (), LshEntry{ h, 0 },
                        [](const LshEntry &a, const LshEntry &b) {
                          return a.hash < b.hash;
                        });
                    while (it != table.end () && it->hash == h)
                      {
                        uint32_t sid2 = it->seq_id;
                        if (sid2 > (uint32_t)sid1 && seen[sid2] != epoch)
                          {
                            seen[sid2] = epoch;
                            candidates.push_back (sid2);
                          }
                        ++it;
                      }
                  }

                for (uint32_t fb_id : fallback_seqs)
                  {
                    if (fb_id > (uint32_t)sid1 && seen[fb_id] != epoch)
                      {
                        seen[fb_id] = epoch;
                        candidates.push_back (fb_id);
                      }
                  }
              }
            else
              {
                for (unsigned long sid2 = sid1 + 1; sid2 < sequenceCount; sid2++)
                  candidates.push_back (sid2);
              }

            for (uint32_t sid2 : candidates)
              {
                if (!quick_bitplane_filter (bitplanes[sid1], bitplanes[sid2],
                                            args.distance, args.overlap))
                  continue;

                local_evals++;
                long weighted_count = instances1 * counts.value (sid2);
                const char *s2 = stringText (sequences, seqLengths, sid2);

                double thisD = sequence_descriptors
                                   ? computeTN93 (
                                         s1, s2, firstSequenceLength,
                                         resolutionOption, NULL, args.overlap,
                                         NULL, HISTOGRAM_SLICE, HISTOGRAM_BINS,
                                         weighted_count, 1L,
                                         &sequence_descriptors[sid1],
                                         &sequence_descriptors[sid2],
                                         args.distance)
                                   : computeTN93 (
                                         s1, s2, firstSequenceLength,
                                         resolutionOption, NULL, args.overlap,
                                         NULL, HISTOGRAM_SLICE, HISTOGRAM_BINS,
                                         weighted_count, 1L, NULL, NULL,
                                         args.distance);

                if (thisD >= args.min_distance && thisD <= args.distance)
                  {
                    local_links += weighted_count;
                    local_sum_d += thisD * weighted_count;
                    local_weighted += weighted_count;
                    if (thisD > local_max_d)
                      local_max_d = thisD;

                    unsigned long bin
                        = (unsigned long)(thisD * HISTOGRAM_SLICE);
                    if (bin >= HISTOGRAM_BINS)
                      bin = HISTOGRAM_BINS - 1;
                    local_hist[bin] += weighted_count;

                    if (!args.do_count)
                      {
                        if (args.format == csv)
                          {
                            char float_buf[128];
                            unsigned written
                                = snprintf (float_buf, 128, "%g", thisD);
                            local_buffer.appendBuffer (n1, n1L);
                            local_buffer.appendChar (args.delimiter);
                            local_buffer.appendBuffer (
                                stringText (names, nameLengths, sid2),
                                stringLength (nameLengths, sid2));
                            local_buffer.appendChar (args.delimiter);
                            local_buffer.appendBuffer (float_buf, written);
                            local_buffer.appendChar ('\n');
                          }
                        else if (args.format == csvn)
                          {
                            char link_buf[1024];
                            int written = snprintf (
                                link_buf, 1024, "%ld%c%ld%c%g\n", sid1,
                                args.delimiter, (long)sid2, args.delimiter,
                                thisD);
                            local_buffer.appendBuffer (link_buf, written);
                          }

                        if (local_buffer.length () > BUFFER_FLUSH_LENGTH)
                          {
#pragma omp critical(fwrite)
                            {
                              fwrite (local_buffer.getString (), sizeof (char),
                                      local_buffer.length (), args.output);
                            }
                            local_buffer.resetString ();
                          }
                      }
                  }
              }

            long current_pair_index;
#pragma omp atomic capture
            current_pair_index = pairIndex += (sequenceCount - sid1 - 1);

            if (!args.quiet
                && (current_pair_index * 100. / pairwise - percentDone > 0.1
                    || sid1 == (long)upperBound - 1))
              {
#pragma omp critical(progress)
                {
                  if (current_pair_index * 100. / pairwise - percentDone > 0.1
                      || sid1 == (long)upperBound - 1)
                    {
                      time (&after);
                      percentDone = current_pair_index * 100. / pairwise;
                      cerr << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
                              "\b\b"
                              "\b\b\b\b\b"
                              "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
                              "\b\b"
                              "\b\b\b\b\b"
                              "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
                              "\b\b"
                              "\bProgress"
                              ":"
                           << setw (8) << percentDone << "% (" << setw (8)
                           << (foundLinks + local_links) << " links found, "
                           << setw (12) << std::setprecision (3)
                           << current_pair_index
                                  / std::max (1.0, difftime (after, before))
                           << " evals/sec)";
                      after = before;
                    }
                }
              }
          }

        if (local_buffer.length () > 0)
          {
#pragma omp critical(fwrite)
            {
              fwrite (local_buffer.getString (), sizeof (char),
                      local_buffer.length (), args.output);
            }
          }

#pragma omp critical(hist)
        {
          foundLinks += local_links;
          actualEvals += local_evals;
          global_sum_d += local_sum_d;
          global_weighted_links += local_weighted;
          if (local_max_d > global_max_d)
            global_max_d = local_max_d;
          for (unsigned long k = 0; k < HISTOGRAM_BINS; k++)
            global_hist[k] += local_hist[k];
        }
      }
    }

  if (distanceMatrix)
    {
      fprintf (args.output, "{");
      for (unsigned long s1 = 0; s1 < sequenceCount; s1++)
        {
          fprintf (args.output, "\n{%g", distanceMatrix[s1 * sequenceCount]);
          for (unsigned long s2 = 1; s2 < sequenceCount; s2++)
            {
              fprintf (args.output, ",%g",
                       distanceMatrix[s1 * sequenceCount + s2]);
            }
          fprintf (args.output, "}");
        }
      fprintf (args.output, "\n}\n");
      delete[] distanceMatrix;
    }

  if (!args.quiet)
    cerr << endl;

  ostream *outStream = &cout;
  if (args.output == stdout)
    outStream = &cerr;

  (*outStream) << "{" << endl;
  (*outStream) << "\t\"Note\" : \"Subquadratic mode active; comparisons and "
                  "histograms reflect candidate-filtered evaluations.\","
               << endl;
  (*outStream) << "\t\"Actual comparisons performed\" : " << actualEvals << ','
               << endl;
  (*outStream) << "\t\"Total comparisons possible\" : " << pairwise << ','
               << endl;
  (*outStream) << "\t\"Links found\" : " << foundLinks << ',' << endl;
  (*outStream) << "\t\"Maximum distance\" : " << global_max_d << ',' << endl;

  if (args.input2 == NULL)
    {
      (*outStream) << "\t\"Sequences\" : " << sequenceCount << ',' << endl;
    }
  else
    {
      (*outStream) << "\t\"Sequences in first file\" : " << seqLengthInFile1
                   << ',' << endl;
      (*outStream) << "\t\"Sequences in second file\" : " << seqLengthInFile2
                   << ',' << endl;
    }

  (*outStream) << "\t\"Mean distance\" : "
               << (global_weighted_links > 0.0
                       ? (global_sum_d / global_weighted_links)
                       : 0.0)
               << ',' << endl;
  dump_subq_histogram (outStream, NULL, global_hist);
  (*outStream) << '}' << endl;

  if (args.do_count)
    {
      fprintf (args.output, "Found %ld links among %ld pairwise comparisons\n",
               foundLinks, actualEvals);
    }

  return 0;
}
