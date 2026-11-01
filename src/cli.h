// cli.h - tiny command-line parser shared by sw_align and sw_bench.
#pragma once

#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

class Args {
public:
    // `valued` lists the options that take a value (e.g. "--threads"). --help / -h is always accepted.
    // A value is always the next argument, so negative numbers such as "--gap -2" work.
    Args(int argc, char** argv, const std::set<std::string>& valued) {
        for (int k = 1; k < argc; ++k) {
            const std::string arg = argv[k];
            if (arg == "--help" || arg == "-h") {
                help_ = true;
            } else if (arg.rfind("--", 0) == 0) {
                if (!valued.count(arg)) throw std::runtime_error("unknown option " + arg);
                if (k + 1 >= argc) throw std::runtime_error("option " + arg + " needs a value");
                options_[arg] = argv[++k];
            } else {
                positional_.push_back(arg);
            }
        }
    }

    bool help() const { return help_; }
    const std::vector<std::string>& positional() const { return positional_; }

    std::string get(const std::string& name, const std::string& def) const {
        auto it = options_.find(name);
        return it == options_.end() ? def : it->second;
    }

    int get_int(const std::string& name, int def) const {
        auto it = options_.find(name);
        return it == options_.end() ? def : parse_int(it->second, name);
    }

    // "1,2,4,8" -> {1,2,4,8}
    std::vector<int> get_int_list(const std::string& name, const std::vector<int>& def) const {
        auto it = options_.find(name);
        if (it == options_.end()) return def;
        std::vector<int> out;
        const std::string& s = it->second;
        size_t start = 0;
        while (start <= s.size()) {
            size_t comma = s.find(',', start);
            if (comma == std::string::npos) comma = s.size();
            out.push_back(parse_int(s.substr(start, comma - start), name));
            start = comma + 1;
        }
        return out;
    }

private:
    static int parse_int(const std::string& text, const std::string& name) {
        try {
            size_t used = 0;
            const int v = std::stoi(text, &used);
            if (used != text.size()) throw std::invalid_argument("trailing characters");
            return v;
        } catch (const std::exception&) {
            throw std::runtime_error("option " + name + " expects an integer, got '" + text + "'");
        }
    }

    bool help_ = false;
    std::vector<std::string> positional_;
    std::map<std::string, std::string> options_;
};
