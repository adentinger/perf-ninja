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
}  // namespace

// The alignment algorithm which computes the alignment of the given sequence
// pairs.
result_t compute_alignment(std::vector<sequence_t> const &sequences1,
                           std::vector<sequence_t> const &sequences2) {
  assert(sequences1.size() == sequence_count_v);
  assert(sequences2.size() == sequence_count_v);

  result_t result{};
  // NB: This algorithm processes multiple sequences at once to allow
  // vectorization opportunities.
  auto tr_sequences1 = std::make_unique<transposed_t>();
  transpose(sequences1, *tr_sequences1);
  auto tr_sequences2 = std::make_unique<transposed_t>();
  transpose(sequences2, *tr_sequences2);

  using score_t = int16_t;
  using columns_t = std::array<std::array<score_t, sequence_count_v>, sequence_size_v + 1>;

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
  auto horizontal_gap_columns = std::make_unique<columns_t>();
  std::array<score_t, sequence_count_v> last_vertical_gaps{};

  /*
    * Initialise the first column of each matrix.
    */
  std::fill(
    (*horizontal_gap_columns)[0].begin(),
    (*horizontal_gap_columns)[0].end(),
    gap_open);
  std::fill(last_vertical_gaps.begin(), last_vertical_gaps.end(), gap_open);

  for (size_t i = 1; i < sequence_size_v + 1; ++i) {
    std::copy(
      last_vertical_gaps.cbegin(),
      last_vertical_gaps.cend(),
      (*score_columns)[i].begin());
    std::transform(
      last_vertical_gaps.cbegin(),
      last_vertical_gaps.cend(),
      (*horizontal_gap_columns)[i].begin(),
      [](auto lvg) {
        return lvg + gap_open;
      }
    );
    std::for_each(
      last_vertical_gaps.begin(),
      last_vertical_gaps.end(),
      [](auto& lvg) {
        lvg += gap_extension;
      }
    );
  }

  /*
    * Compute the main recursion to fill the matrices.
    */
  for (unsigned col = 1; col <= sequence_size_v; ++col) {
    // Cache last diagonal scores to compute the cells at this position.
    std::array<score_t, sequence_count_v> last_diagonal_scores{};
    std::copy(
      (*score_columns)[0].cbegin(),
      (*score_columns)[0].cend(),
      last_diagonal_scores.begin());
    std::copy(
      (*horizontal_gap_columns)[0].cbegin(),
      (*horizontal_gap_columns)[0].cend(),
      (*score_columns)[0].begin());
    std::transform(
      (*horizontal_gap_columns)[0].cbegin(),
      (*horizontal_gap_columns)[0].cend(),
      last_vertical_gaps.begin(),
      [](auto horiz_gap) {
        return horiz_gap + gap_open;
      });
    std::for_each(
      (*horizontal_gap_columns)[0].begin(),
      (*horizontal_gap_columns)[0].end(),
      [](auto& horiz_gap) {
        horiz_gap += gap_extension;
      });

    for (unsigned row = 1; row <= sequence_size_v; ++row) {
      for (size_t seq = 0; seq < sequence_count_v; ++seq) {
        // Compute next score from diagonal direction with match/mismatch.
        score_t best_cell_score =
            last_diagonal_scores[seq] +
            ((*tr_sequences1)[row - 1][seq] == (*tr_sequences2)[col - 1][seq] ? match : mismatch);
        // Determine best score from diagonal, vertical, or horizontal
        // direction.
        best_cell_score = std::max(best_cell_score, last_vertical_gaps[seq]);
        best_cell_score = std::max(best_cell_score, (*horizontal_gap_columns)[row][seq]);
        // Cache next diagonal value and store optimum in score_column.
        last_diagonal_scores[seq] = (*score_columns)[row][seq];
        (*score_columns)[row][seq] = best_cell_score;
        // Compute the next values for vertical and horizontal gap.
        best_cell_score += gap_open;
        last_vertical_gaps[seq] += gap_extension;
        (*horizontal_gap_columns)[row][seq] += gap_extension;
        // Store optimum between gap open and gap extension.
        last_vertical_gaps[seq] = std::max(last_vertical_gaps[seq], best_cell_score);
        (*horizontal_gap_columns)[row][seq] =
            std::max((*horizontal_gap_columns)[row][seq], best_cell_score);
      }
    }
  }

  // Report the best score.
  std::copy(
    score_columns->back().cbegin(),
    score_columns->back().cend(),
    result.begin());

  return result;
}
