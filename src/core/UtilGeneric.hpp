#ifndef RGM_UTILGENERIC_HPP_
#define RGM_UTILGENERIC_HPP_

#include <vector>

namespace RGM
{

namespace detail
{

template<typename T>
struct pair_less_than_ {
    template<typename T1>
    bool operator()(const std::pair<T, T1> &a, const std::pair<T, T1> &b)
    {
        return a.first < b.first;
    }
};

template<typename T>
struct pair_larger_than_ {
    template<typename T1>
    bool operator()(const std::pair<T, T1> &a, const std::pair<T, T1> &b)
    {
        return a.first > b.first;
    }
};

template<typename T>
struct equal_to_ {
    bool operator()(T i, T j)
    {
        return (i == j);
    }
};

} // namespace detail

// Fci - forward const iterator
template<typename Fci>
std::vector<std::vector<Fci> > enumerateCombinations_(Fci begin,
        Fci end, unsigned int combination_size);

template<typename T>
std::vector<T> & operator+=(std::vector<T> & a, const std::vector<T> & b);

template <typename T>
void uniqueVector_(std::vector<T> & vec);


template <typename T>
std::vector<int> sortVector_(std::vector<T> & vData, bool ascending = false, bool bSorted =false);

} // namespace RGM

#endif // RGM_UTILGENERIC_HPP_
