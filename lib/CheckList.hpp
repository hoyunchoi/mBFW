#pragma once

#include <algorithm>
#include <array>
#include <stdexcept>
#include <utility>

namespace mBFW {

template <typename Key, typename Value, std::size_t Size>
struct CheckList {
    std::array<std::pair<Key, Value>, Size> data;

    constexpr Value at(const Key& key) const {
        const auto itr = std::find_if(begin(data), end(data), [&key](const auto& v) { return v.first == key; });
        if (itr != end(data)) {
            return itr->second;
        } else {
            throw std::range_error("Not Found");
        }
    }
};

}  // namespace mBFW
