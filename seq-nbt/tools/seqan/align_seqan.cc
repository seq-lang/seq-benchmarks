#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <seqan3/alignment/configuration/align_config_gap_cost_affine.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <string>
#include <utility>
using namespace std;
using namespace std::chrono;

double score_only(char **argv) {
  std::ifstream queries(argv[1]);
  std::ifstream targets(argv[2]);
  auto config=
      seqan3::align_cfg::method_local{} |
      seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{
          seqan3::match_score{1}, seqan3::mismatch_score{-2}}} |
      seqan3::align_cfg::gap_cost_affine{
          seqan3::align_cfg::open_score{-4},
          seqan3::align_cfg::extension_score{-2}} |
      seqan3::align_cfg::output_score{};
  long total = 0;
  long num = 0;
  std::string query, target;
  auto t = high_resolution_clock::now();
  while (std::getline(queries, query) && std::getline(targets, target)) {
    seqan3::dna4_vector q{};
    seqan3::dna4_vector t{};
    for (char c : query)
      q.push_back(seqan3::assign_char_to(c, seqan3::dna4{}));
    for (char c : target)
      t.push_back(seqan3::assign_char_to(c, seqan3::dna4{}));
    auto results = seqan3::align_pairwise(std::tie(q, t), config);
    auto &res = *results.begin();
    total += res.score();
    ++num;
  }
  // std::cout << num << " " << total << " " << x << std::endl;
  return duration_cast<milliseconds>(high_resolution_clock::now() - t).count() /
          1000.0;
}
double backtrace(char **argv) {
  std::ifstream queries(argv[1]);
  std::ifstream targets(argv[2]);
  auto config =
    seqan3::align_cfg::method_local{} |
    seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{
        seqan3::match_score{1}, seqan3::mismatch_score{-2}}} |
    seqan3::align_cfg::gap_cost_affine{
        seqan3::align_cfg::open_score{-4},
        seqan3::align_cfg::extension_score{-2}} |
    seqan3::align_cfg::output_alignment{};
  long total = 0;
  long num = 0;
  std::string query, target;
  auto t = high_resolution_clock::now();
  while (std::getline(queries, query) && std::getline(targets, target)) {
    seqan3::dna4_vector q{};
    seqan3::dna4_vector t{};
    for (char c : query)
      q.push_back(seqan3::assign_char_to(c, seqan3::dna4{}));
    for (char c : target)
      t.push_back(seqan3::assign_char_to(c, seqan3::dna4{}));
    auto results = seqan3::align_pairwise(std::tie(q, t), config);
    auto &res = *results.begin();
    // total += res.score();
    auto a = res.alignment();
    ++num;
  }
  // std::cout << num << " " << total << " " << x << std::endl;
  return duration_cast<milliseconds>(high_resolution_clock::now() - t).count() /
          1000.0;
}

int main(int argc, char *argv[]) {
  assert(argc == 3);
  vector<double> tm;
  double tms = 0;
  for (int so = 0; so < 2; so++) {
    for (int _ = 0; _ < 3; _++) {
      double x = so == 1 ? score_only(argv) : backtrace(argv);
      tm.push_back(x);
      tms += x;
    }

    double stdev = 0, mean = tms / tm.size();
    for (auto t : tm)
      stdev += (mean - t) * (mean - t);
    cout << "[sw-time] seqan " << so << " " << mean << " " << sqrt(stdev / tm.size()) << endl;
  }
  return 0;
}