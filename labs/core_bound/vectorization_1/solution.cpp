#include "solution.hpp"
#include <algorithm>
#include <cassert>
#include <type_traits>
#include <memory>

using transposed_t =
  std::array<std::array<uint8_t, sequence_count_v>, sequence_size_v>;

namespace {
void transpose(std::vector<sequence_t> const &sequences,
               transposed_t &transposed) {
  for (size_t i = 0; i < sequence_size_v; ++i) {
    for (size_t j = 0; j < sequence_count_v; ++j) {
      transposed[i][j] = sequences[j][i];
    }
  }
}
}

// The alignment algorithm which computes the alignment of the given sequence
// pairs.
result_t compute_alignment(std::vector<sequence_t> const &sequences1,
                           std::vector<sequence_t> const &sequences2) {
  result_t result{};
  // NB: This algorithm processes multiple sequences at once to allow
  // vectorization opportunities.
  auto transposed_sequences1 = std::make_unique<transposed_t>();
  transpose(sequences1, transposed_sequences1);
  auto transposed_sequences2 = std::make_unique<transposed_t>();
  transpose(sequences2, transposed_sequences2);

  for (size_t sequence_idx = 0; sequence_idx < sequences1.size();
       ++sequence_idx) {
    using score_t = int16_t;
    using columns_t = std::array<std::array<score_t, sequence_count_v>, sequence_size_v + 1>;

    sequence_t const &sequence1 = sequences1[sequence_idx];
    sequence_t const &sequence2 = sequences2[sequence_idx];

    /*
     * Initialise score values.
     */
    constexpr score_t gap_open{-11};
    constexpr score_t gap_extension{-1};
    constexpr score_t match{6};
    constexpr score_t mismatch{-4};

    /*
     * Setup the matrices.
     * Note we can compute each matrix with just one column of each one in
     * memory, since we are only interested in the last value of the last
     * column in each score matrix.
     */
    auto score_columns = std::make_unique<columns_t>();
    auto horizontal_gap_column = std::make_unique<columns_t>();
    std::array<score_t, sequence_count_v> last_vertical_gaps{};

    /*
     * Initialise the first column of each matrix.
     */
    std::fill_n(horizontal_gap_column.begin(), sequence_size_v + 1, gap_open);
    std::fill(last_vertical_gaps.begin(), last_vertical_gaps.end(), gap_open);

    for (size_t i = 1; i < sequence_size_v + 1; ++i) {
      for (size_t seq = 1; seq < sequence_size_v; ++i) {
        score_columns[i][seq] = last_vertical_gap[seq];
        horizontal_gap_column[i][seq] = last_vertical_gap[seq] + gap_open;
        last_vertical_gap[seq] += gap_extension;
      }
    }

    /*
     * Compute the main recursion to fill the matrices.
     */
    for (unsigned col = 1; col <= sequence2.size(); ++col) {
      // Cache last diagonal scores to compute the cells at this position.
      std::array<score_t, sequence_count_v> last_diagonal_scores{};
      std::copy(
        score_column[0].cbegin(),
        score_column[0].cend(),
        last_diagonal_scores.begin());
      std::copy(
        horizontal_gap_column[0].cbegin(),
        horizontal_gap_column[0].cend(),
        score_column.begin());
      std::transform(
        horizontal_gap_column[0].cbegin(),
        horizontal_gap_column[0].cend(),
        last_vertical_gap.begin(),
        [](auto horiz_gap) {
          return horiz_gap + gap_open;
        });
      std::for_each_n(
        horizontal_gap_column[0].begin(),
        horizontal_gap_column[0].end(),
        [](auto& horiz_gap) {
          return horiz_gap + gap_extension;
        });

      for (unsigned row = 1; row <= sequence1.size(); ++row) {
        for (size_t seq = 1; seq < sequence_size_v; ++i) {
          // Compute next score from diagonal direction with match/mismatch.
          score_t best_cell_score =
              last_diagonal_score +
              (sequence1[row - 1] == sequence2[col - 1] ? match : mismatch);
          // Determine best score from diagonal, vertical, or horizontal
          // direction.
          best_cell_score = std::max(best_cell_score, last_vertical_gap);
          best_cell_score = std::max(best_cell_score, horizontal_gap_column[row]);
          // Cache next diagonal value and store optimum in score_column.
          last_diagonal_score = score_column[row];
          score_column[row] = best_cell_score;
          // Compute the next values for vertical and horizontal gap.
          best_cell_score += gap_open;
          last_vertical_gap += gap_extension;
          horizontal_gap_column[row] += gap_extension;
          // Store optimum between gap open and gap extension.
          last_vertical_gap = std::max(last_vertical_gap, best_cell_score);
          horizontal_gap_column[row] =
              std::max(horizontal_gap_column[row], best_cell_score);
        }
      }
    }

    // Report the best score.
    result[sequence_idx] = score_column.back();
  }

  return result;
}
