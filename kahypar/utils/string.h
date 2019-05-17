#pragma once

#include <string>

namespace kahypar {
namespace string {
bool starts_with(const std::string& haystack, const std::string &needle) {
  return haystack.substr(0, needle.length()) == needle;
}

std::string trim(std::string text, const char *symbols = " \t\n\r\f\v") {
  text.erase(0, text.find_first_not_of(symbols));
  text.erase(text.find_last_not_of(symbols) + 1);
  return text;
}

std::vector<std::string> split(const std::string& text, const std::string& delimiter, bool trim_tokens = true) {
  std::vector<std::string> tokens;
  std::size_t start = 0;
  std::size_t end = text.find(delimiter);

  while (end != std::string::npos) {
    std::string token = text.substr(start, end - start);
    start = end + delimiter.length();
    end = text.find(delimiter, start);

    if (trim_tokens) token = trim(std::move(token));
    tokens.push_back(std::move(token));
  }

  // last token
  std::string token = text.substr(start, end);
  if (trim_tokens) token = trim(std::move(token));
  tokens.push_back(std::move(token));

  return tokens;
}
} // namespace string
} // namespace kahypar