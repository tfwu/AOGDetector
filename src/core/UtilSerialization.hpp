#ifndef RGM_UTILSERIALIZATION_HPP_
#define RGM_UTILSERIALIZATION_HPP_

#include "Common.hpp"


namespace boost
{
namespace serialization
{
/// Eigen::Array
template <class Archive, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
void save(Archive & ar, const Eigen::Array<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols> & m, const unsigned int version)
{
    int rows=m.rows(),cols=m.cols();
    ar & BOOST_SERIALIZATION_NVP(rows);
    ar & BOOST_SERIALIZATION_NVP(cols);
    ar & make_array(m.data(), rows*cols);
}

template <class Archive, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
void load(Archive & ar, Eigen::Array<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols> & m, const unsigned int version)
{
    int rows,cols;
    ar & BOOST_SERIALIZATION_NVP(rows);
    ar & BOOST_SERIALIZATION_NVP(cols);
    m.resize(rows,cols);
    ar & make_array(m.data(), rows*cols);
}

template <class Archive, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
void serialize(Archive & ar, Eigen::Array<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols> & m, const unsigned int version)
{
    split_free(ar,m,version);
}

///  Eigen::Matrix
template <class Archive, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
void save(Archive & ar, const Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols> & m, const unsigned int version)
{
    int rows=m.rows(),cols=m.cols();
    ar & BOOST_SERIALIZATION_NVP(rows);
    ar & BOOST_SERIALIZATION_NVP(cols);
    ar & make_array(m.data(), rows*cols);
}

template <class Archive, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
void load(Archive & ar, Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols> & m, const unsigned int version)
{
    int rows,cols;
    ar & BOOST_SERIALIZATION_NVP(rows);
    ar & BOOST_SERIALIZATION_NVP(cols);
    m.resize(rows,cols);
    ar & make_array(m.data(), rows*cols);
}

template <class Archive, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
void serialize(Archive & ar, Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols> & m, const unsigned int version)
{
    split_free(ar,m,version);
}

} // namespace serialziation
} // namespace boost

#endif // RGM_UTILSERIALIZATION_HPP_
