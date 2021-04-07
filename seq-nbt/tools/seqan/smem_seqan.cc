#include <algorithm>
#include <chrono>
#include <iostream>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp> // pretty printing
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>
#include <seqan3/std/filesystem>
using namespace std;
using namespace std::chrono;

constexpr int64_t min_seed = 17;
constexpr int64_t min_iwidth = 20;
constexpr int64_t min_intv = 1;

template <typename Seq> void print(const Seq &seq, int64_t from, int64_t to) {
  for (size_t i = from; i < to; i++)
    std::cout << seq[i].to_char();
}

template <typename Cursor> struct Interval {
  Cursor fwd;
  Cursor rev;
  bool valid_fwd = true;
  bool valid_rev = true;

  size_t count() const {
    return (valid_fwd ? fwd.count() : 0) + (valid_rev ? rev.count() : 0);
  }

  template <typename T> void extend_right(const T &x) {
    if (valid_fwd)
      valid_fwd = fwd.extend_right(x);
    if (valid_rev)
      valid_rev = rev.extend_left(x.complement());
  }

  template <typename T> void extend_left(const T &x) {
    if (valid_fwd)
      valid_fwd = fwd.extend_left(x);
    if (valid_rev)
      valid_rev = rev.extend_right(x.complement());
  }

  auto locate_fwd() const { return fwd.locate(); }

  auto locate_rev() const { return rev.locate(); }
};

template <typename Intv> struct SMEM {
  Intv interval;
  int64_t start;
  int64_t stop;

  int64_t len() const { return stop - start; }
  operator bool() const { return stop >= start; }
};

template <typename T> bool N(const T &c) { return c.to_char() == 'N'; }

template <typename FMIndex, typename Rec, typename Ref>
void findSMEMs(FMIndex &fmi, const Rec &rec, const Ref &ref) {
  using Cursor = typename FMIndex::cursor_type;
  using Intv = Interval<Cursor>;
  using MySMEM = SMEM<Intv>;
  auto q = seqan3::get<seqan3::field::seq>(rec);
  std::vector<MySMEM> prev, curr, mems;
  int64_t start = 0;
  int64_t l = q.size();

  while (true) {
    while (start < l && N(q[start]))
      ++start;
    if (start >= l)
      break;

    mems.clear();
    prev.clear();
    curr.clear();
    int64_t x = start;

    if (N(q[x]))
      return;

    Intv cursor = {fmi.cursor(), fmi.cursor()};
    cursor.extend_right(q[x]);
    MySMEM ik = {cursor, x, x + 1};

    // forwards search
    int64_t i = x + 1;
    for (; i < l; i++) {
      if (!N(q[i])) {
        auto ok = ik.interval;
        ok.extend_right(q[i]);
        if (ok.count() != ik.interval.count()) {
          curr.push_back(ik);
          if (ok.count() < min_intv)
            break;
        }
        ik = {ok, x, i + 1};
      } else {
        curr.push_back(ik);
        break;
      }
    }

    if (i == l)
      curr.push_back(ik);
    std::reverse(curr.begin(), curr.end());
    int64_t ret = curr[0].stop;
    prev.swap(curr);

    // backward search for MEMs
    for (i = x - 1; i >= -1; i--) {
      bool c = (i >= 0 && !N(q[i]));
      curr.clear();
      for (const auto &p : prev) {
        Intv ok;
        if (c) {
          ok = p.interval;
          seqan3::dna5 b = q[i];
          ok.extend_left(b);
        }
        if (!c || ok.count() < min_intv) {
          if (curr.empty()) {
            if (mems.empty() || i + 1 < mems.back().start) {
              MySMEM ik = {p.interval, i + 1, p.stop};
              if (ik.len() >= min_seed) {
                mems.push_back(ik);
              }
            }
          }
        } else if (curr.empty() || ok.count() != curr.back().interval.count()) {
          curr.push_back({ok, p.start, p.stop});
        }
      }
      if (curr.empty())
        break;
      prev.swap(curr);
    }

    std::reverse(mems.begin(), mems.end());
    start = ret;

    for (const auto &mem : mems) {
      auto &intv = mem.interval;
      auto offset = mem.start;
      auto match_size = mem.len();
      std::cout << seqan3::get<seqan3::field::id>(rec) << "\tEM\t" << offset
                << "\t" << (offset + match_size) << "\t" << intv.count();
      if (intv.count() <= min_iwidth) {
        if (intv.valid_fwd) {
          for (const auto &result : intv.locate_fwd()) {
            auto id = std::get<0>(result);
            auto pos = std::get<1>(result);
            std::cout << "\t" << id << ":+" << (pos + 1) << ":";
            print(q, offset, offset + match_size);
            std::cout << ":";
            print(ref[id], pos, pos + match_size);
          }
        }
        if (intv.valid_rev) {
          for (const auto &result : intv.locate_rev()) {
            auto id = std::get<0>(result);
            auto pos = std::get<1>(result);
            std::cout << "\t" << id << ":-" << (pos + 1) << ":";
            print(q, offset, offset + match_size);
            std::cout << ":";
            print(ref[id], pos, pos + match_size);
          }
        }
      } else {
        std::cout << "\t*";
      }
      std::cout << std::endl;
    }
  }
}

int main(int argc, char *argv[]) {
  std::cin.tie(nullptr);
  std::cout.sync_with_stdio(false);
  std::cerr << "Seqan\n";
  assert(argc == 3);
  seqan3::sequence_file_input fasta{argv[1]};
  std::vector<seqan3::dna5_vector> seqs;
  auto t = high_resolution_clock::now();
  for (auto &rec : fasta) {
    seqs.push_back(std::move(seqan3::get<seqan3::field::seq>(rec)));
    break; // only 1st chromosome
  }
  seqan3::bi_fm_index fmi{seqs};
  std::cerr
      << "index "
      << duration_cast<milliseconds>(high_resolution_clock::now() - t).count() /
             1000.0
      << "\n";
  t = high_resolution_clock::now();

  seqan3::sequence_file_input fastq{argv[2]};
  for (auto &rec : fastq) {
    findSMEMs(fmi, rec, seqs);
  }
  std::cerr
      << "smem "
      << duration_cast<milliseconds>(high_resolution_clock::now() - t).count() /
             1000.0
      << "\n";
}