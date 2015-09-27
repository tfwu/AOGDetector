#include <algorithm>
#include <assert.h>

#include "UtilGeneric.hpp"
#include "Common.hpp"
#include "UtilLog.hpp"


namespace RGM
{

template<typename Fci>
std::vector<std::vector<Fci> > enumerateCombinations_(Fci begin,
        Fci end, unsigned int combination_size)
{
    DEFINE_RGM_LOGGER;

    // empty set of combinations
    std::vector<std::vector<Fci> > result;
    if (combination_size == 0u) {
        return result;    // there is exactly one combination of size 0 - empty set
    }

    std::vector<Fci> current_combination;
    current_combination.reserve(combination_size + 1u); // one additional slot

    // in my vector to store
    // the end sentinel there.
    // The code is cleaner thanks to that
    for (unsigned int i = 0u; i < combination_size && begin != end;
            ++i, ++begin) {
        current_combination.push_back(begin);           // Construction of the first combination
    }

    // Since I assume the iterators support only incrementing, I have to iterate over
    // the set to get its size, which is expensive. Here I had to iterate anyway to
    // produce the first combination, so I use the loop to also check the size.
    RGM_CHECK_GE(current_combination.size(), combination_size);

    result.push_back(current_combination);              // Store the first combination in the results set
    current_combination.push_back(end);                 // Here I add mentioned earlier sentinel to
    // simplify rest of the code. If I did it
    // earlier, previous statement would get ugly.
    while (true) {
        unsigned int i = combination_size;
        Fci tmp;                                        // Thanks to the sentinel I can find first
        do {                                            // iterator to change, simply by scaning
            // from right to left and looking for the
            tmp = current_combination[--i];             // first "bubble". The fact, that it's
            ++tmp;                                      // a forward iterator makes it ugly but I
        }                                               // can't help it.
        while (i > 0u && tmp == current_combination[i + 1u]);

        // Here is probably my most obfuscated expression.
        // Loop above looks for a "bubble". If there is no "bubble", that means, that
        // current_combination is the last combination, Expression in the if statement
        // below evaluates to true and the function exits returning result.
        // If the "bubble" is found however, the statement below has a side effect of
        // incrementing the first iterator to the left of the "bubble".
        if (++current_combination[i] == current_combination[i + 1u]) {
            return result;
        }
        // Rest of the code sets positions of the rest of the iterators
        // (if there are any), that are to the right of the incremented one,
        // to form next combination

        while (++i < combination_size) {
            current_combination[i] = current_combination[i - 1u];
            ++current_combination[i];
        }
        // Below is the ugly side of using the sentinel. Well it had to have some
        // disadvantage. Try without it.
        result.push_back(
            std::vector<Fci>(current_combination.begin(),
                             current_combination.end() - 1));
    }
}

template std::vector<std::vector<std::vector<int>::const_iterator > > enumerateCombinations_(std::vector<int>::const_iterator begin,
        std::vector<int>::const_iterator end, unsigned int combination_size);


template<typename T>
std::vector<T> & operator+=(std::vector<T> & a, const std::vector<T> & b)
{
    DEFINE_RGM_LOGGER;

    RGM_CHECK_EQ(a.size(),  b.size());
    if (a.size() == 0) {
        return a;
    }

    std::vector<T> result;
    result.reserve(a.size());

    std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result),
                   std::plus<T>());

    a.swap(result);

    return a;
}

template std::vector<int> & operator+=(std::vector<int> & a, const std::vector<int> & b);
template std::vector<float> & operator+=(std::vector<float> & a, const std::vector<float> & b);
template std::vector<double> & operator+=(std::vector<double> & a, const std::vector<double> & b);

template <typename T>
void uniqueVector_(std::vector<T> & vec)
{
    std::sort(vec.begin(), vec.end());
    typename std::vector<T>::iterator it = std::unique(vec.begin(), vec.end(), detail::equal_to_<T>());
    vec.resize(std::distance(vec.begin(), it));
}

template void uniqueVector_<int>(std::vector<int> & vec);
template void uniqueVector_<unsigned int>(std::vector<unsigned int> & vec);
template void uniqueVector_<float>(std::vector<float> & vec);
template void uniqueVector_<double>(std::vector<double> & vec);

template <typename T>
std::vector<int> sortVector_(std::vector<T> & vData, bool ascending/* = false */, bool bSorted/* =false */)
{
    int num = (int)vData.size();
    assert(num>0);

    std::vector<std::pair<T, int> > vTemp(num);
    for ( int i=0; i<num; ++i )	{
        vTemp[i] = std::pair<T, int>(vData[i], i);
    }

    if ( ascending )
        std::sort(vTemp.begin(), vTemp.end(), detail::pair_less_than_<T>());
    else
        std::sort(vTemp.begin(), vTemp.end(), detail::pair_larger_than_<T>());

    std::vector<int> vIndex(num);
    for ( int i=0; i<num; ++i ) {
        vIndex[i] = vTemp[i].second;

        if ( bSorted )
            vData[i] = vTemp[i].first;
    }

    return vIndex;
}

template std::vector<int> sortVector_(std::vector<int> & vData, bool ascending, bool bSorted);
template std::vector<int> sortVector_(std::vector<float> & vData, bool ascending, bool bSorted);
template std::vector<int> sortVector_(std::vector<double> & vData, bool ascending, bool bSorted);

} // namespace RGM
