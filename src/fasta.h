// fasta.h - minimal FASTA reader.
#pragma once

#include <string>

struct FastaRecord {
    std::string name;      // first word of the '>' header line
    std::string sequence;  // upper-case bases, only A C G T N
};

// Reads the FIRST record of a FASTA file (a warning is printed to stderr if more records follow).
// Whitespace and line breaks are ignored and lower-case is converted to upper-case. Any character
// other than A, C, G, T, N is an error. A file with no bases yields an empty sequence.
// Throws std::runtime_error on I/O or format problems.
FastaRecord read_fasta(const std::string& path);
