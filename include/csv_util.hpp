#ifndef CSV_UTIL
#define CSV_UTIL
#include <iostream>
#include <set>
#include <vector>

namespace CSV {
    inline std::vector<std::string>
    parse_line(const std::string& line, char separator, const std::set<char>& escape_chars,
               const std::set<char>& ignored_chars) {
        std::vector<std::string> result;
        std::string              word;
        bool                     in_escape = false;
        for (char c: line) {
            if (ignored_chars.find(c) != ignored_chars.end()) {
                continue;
            }
            if (in_escape) {
                if (escape_chars.find(c) != escape_chars.end()) {
                    word += c;
                    in_escape = false;
                } else {
                    word += c;
                }
            } else {
                if (escape_chars.find(c) != escape_chars.end()) {
                    in_escape = true;
                } else if (c == separator) {
                    result.push_back(word);
                    word = "";
                } else {
                    word += c;
                }
            }
        }
        result.push_back(word);
        return result;
    };

} // namespace CSV

#endif // CSV_UTIL
