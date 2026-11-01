// fasta.cpp - minimal FASTA reader.
#include "fasta.h"

#include <cctype>
#include <fstream>
#include <iostream>
#include <stdexcept>

FastaRecord read_fasta(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open '" + path + "'");

    FastaRecord rec;
    bool seen_header = false;
    std::string line;
    int line_no = 0;

    while (std::getline(in, line)) {
        ++line_no;
        if (!line.empty() && line[0] == '>') {
            if (seen_header || !rec.sequence.empty()) {
                std::cerr << "warning: '" << path << "' has more than one record; using only the first\n";
                break;
            }
            seen_header = true;
            // keep just the first word of the header as the sequence name
            size_t end = line.find_first_of(" \t", 1);
            rec.name = line.substr(1, end == std::string::npos ? std::string::npos : end - 1);
            continue;
        }
        for (char raw : line) {
            if (std::isspace(static_cast<unsigned char>(raw))) continue;
            const char c = static_cast<char>(std::toupper(static_cast<unsigned char>(raw)));
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'N') {
                throw std::runtime_error(path + ":" + std::to_string(line_no) + ": invalid base '" +
                                         std::string(1, raw) + "' (expected A, C, G, T or N)");
            }
            rec.sequence.push_back(c);
        }
    }

    if (rec.name.empty()) {
        const size_t slash = path.find_last_of("/\\");
        rec.name = (slash == std::string::npos) ? path : path.substr(slash + 1);
    }
    return rec;
}
