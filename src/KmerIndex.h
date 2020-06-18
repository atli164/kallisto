#ifndef KALLISTO_KMERINDEX_H
#define KALLISTO_KMERINDEX_H


#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <stdint.h>
#include <ostream>
//#include <map>


#include "common.h"
#include "Kmer.hpp"

#include "hash.hpp"

#include <ColoredCDBG.hpp>

std::string revcomp(const std::string s);

struct TRInfo {
  int trid;
  int start;
  int stop; //exclusive [start,stop)
  bool sense; // true for sense, false for anti-sense
};

using EcMap = std::vector<std::vector<int>>;

struct SortedVectorHasher {
  size_t operator()(const std::vector<int>& v) const {
    uint64_t r = 0;
    int i=0;
    for (auto x : v) {
      uint64_t t;
      MurmurHash3_x64_64(&x,sizeof(x), 0,&t);
      t = (x>>i) | (x<<(64-i));
      r = r ^ t;
      i = (i+1)%64;
    }
    return r;
  }
};

struct ContigToTranscript {
  int trid;
  int pos; 
  bool sense; // true for sense, 
};

class KmerEntry : public Bifrost::CCDBG_Data_t<KmerEntry> {
public:
  Bifrost::Kmer kmer;
  uint32_t _pos; // 0-based forward distance to EC-junction
  int length; // number of k-mers in contig
  int ec, id;
  std::string seq; // sequence
  std::vector<ContigToTranscript> transcripts;

  KmerEntry() : _pos(0xFFFFFFF), length(0), kmer(Bifrost::Kmer()) {}
  KmerEntry(int id, int len, int pos, bool isFw, Bifrost::Kmer k) : id(id), length(len), kmer(k) {
    setPos(pos);
    setDir(isFw);
    seq = "";
    ec = -1;
    transcripts = std::vector<ContigToTranscript>();
  }

  inline int getPos() const {return (_pos & 0x0FFFFFFF);}
  inline int isFw() const  {return (_pos & 0xF0000000) == 0; }
  inline void setPos(int p) {_pos = (_pos & 0xF0000000) | (p & 0x0FFFFFFF);}
  inline void setDir(bool _isFw) {_pos = (_pos & 0x0FFFFFFF) | ((_isFw) ? 0 : 0xF0000000);}
  inline int getDist(bool fw) const {
    if (isFw() == fw) {
      return (length - 1 - getPos());
    } else {
      return getPos();
    }
  }
};

struct KmerIndex {
  KmerIndex(const ProgramOptions& opt) : k(opt.k), num_trans(0), skip(opt.skip), target_seqs_loaded(false) { }

  ~KmerIndex() {}

  void match(const char *s, int l, std::vector<std::pair<KmerEntry, int>>& v) const;
  int mapPair(const char *s1, int l1, const char *s2, int l2, int ec) const;
  std::vector<int> intersect(int ec, const std::vector<int>& v) const;

  void BuildTranscripts(const ProgramOptions& opt);
  void BuildDeBruijnGraph(const ProgramOptions& opt, const std::vector<std::string>& seqs);
  void BuildEquivalenceClasses(const ProgramOptions& opt, const std::vector<std::string>& seqs);
  void FixSplitContigs(const ProgramOptions& opt, std::vector<std::vector<TRInfo>>& trinfos);
  bool fwStep(Bifrost::Kmer km, Bifrost::Kmer& end) const;

  // output methods
  void write(const std::string& index_out, bool writeBifrostTable = true);
  void writePseudoBamHeader(std::ostream &o) const;
  
  // note opt is not const
  // load methods
  void load(ProgramOptions& opt, bool loadBifrostTable = true);
  void loadTranscriptSequences() const;
  void clear();

  // positional information
  std::pair<int,bool> findPosition(int tr, Bifrost::Kmer km, KmerEntry val, int p = 0) const;
  std::pair<int,bool> findPosition(int tr, Bifrost::Kmer km, int p) const;

  int k; // k-mer size used
  int num_trans; // number of targets
  int skip;
  int idcnt = 0; // used to assign incremental ids

  Bifrost::ColoredCDBG<KmerEntry> dbGraph;
  EcMap ecmap;
  std::unordered_map<std::vector<int>, int, SortedVectorHasher> ecmapinv;
  const size_t INDEX_VERSION = 10; // increase this every time you change the fileformat

  std::vector<int> target_lens_;

  std::vector<std::string> target_names_;
  std::vector<std::string> target_seqs_; // populated on demand
  bool target_seqs_loaded;
};

#endif // KALLISTO_KMERINDEX_H
