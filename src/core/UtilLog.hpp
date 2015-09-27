/// This file is modified from c++ boost log library
/*
 *          Copyright Andrey Semashev 2007 - 2013.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef RGM_UTILLOG_HPP_
#define RGM_UTILLOG_HPP_

#include <cstddef>
#include <string>
#include <string.h>
#include <ostream>
#include <fstream>
#include <iomanip>
#include <boost/format.hpp>
#include <boost/smart_ptr/shared_ptr.hpp>
#include <boost/smart_ptr/make_shared_object.hpp>
#include <boost/phoenix/bind.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/attributes.hpp>
#include <boost/log/sources/basic_logger.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/sources/severity_channel_logger.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/log/sinks/sync_frontend.hpp>
#include <boost/log/sinks/text_ostream_backend.hpp>
#include <boost/log/attributes/scoped_attribute.hpp>
#include <boost/log/utility/value_ref.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/utility/setup/console.hpp>

#include "UtilFile.hpp"

namespace RGM
{
namespace logging  = boost::log;
namespace src      = boost::log::sources;
namespace expr     = boost::log::expressions;
namespace sinks    = boost::log::sinks;
namespace attrs    = boost::log::attributes;
namespace keywords = boost::log::keywords;

/// We define our own severity levels
enum severity_level {
    normal,
    notification,
    warning,
    error,
    critical
};

/// The operator puts a human-friendly representation of the severity level to the stream
std::ostream& operator<< (std::ostream& strm, severity_level level);

/// init the logger
BOOST_LOG_ATTRIBUTE_KEYWORD(line_id, "LineID", unsigned int)
BOOST_LOG_ATTRIBUTE_KEYWORD(severity, "Severity", severity_level)
BOOST_LOG_ATTRIBUTE_KEYWORD(tag_attr, "Tag", std::string)

#ifndef RGM_NO_LOG
void log_init();
#endif

/// Get file name
#define FILENAME \
    (strrchr(__FILE__, FileUtil::FILESEP) ? \
     strrchr(__FILE__, FileUtil::FILESEP) + 1 : __FILE__)

/// Define a logger (multithread)
#ifdef RGM_NO_LOG
#define DEFINE_RGM_LOGGER
#else
#define DEFINE_RGM_LOGGER \
    src::severity_logger_mt< severity_level > slg
#endif

/// Record a log
#ifdef RGM_NO_LOG
#define RGM_LOG(level, msg)
#else
#define RGM_LOG(level, msg) \
    do { \
        BOOST_LOG_SCOPED_THREAD_TAG("Tag", FILENAME); \
        BOOST_LOG_SEV(slg, level) << msg; \
    } while (0)
#endif

/// Define CHECK macro
#ifdef RGM_NO_LOG
#define RGM_CHECK(conditions, level)
#else
#define RGM_CHECK(conditions, level) \
    do { \
        if (!(conditions)) { \
            BOOST_LOG_SCOPED_THREAD_TAG("Tag", FILENAME); \
            BOOST_LOG_SEV(slg, level) << "Check Failed: " <<  #conditions << "(" << FILENAME << ")"; \
            if ( level >= error ) {\
                exit(0);\
            }\
        } \
    } while (0)
#endif

/// Check equality
#define RGM_CHECK_EQ(lhs, rhs) \
    RGM_CHECK((lhs == rhs), error)

/// Check not equal
#define RGM_CHECK_NOTEQ(lhs, rhs) \
    RGM_CHECK((lhs != rhs), error)

/// Check greater than
#define RGM_CHECK_GT(lhs, rhs) \
    RGM_CHECK((lhs > rhs), error)

/// Check greater than or equal to
#define RGM_CHECK_GE(lhs, rhs) \
    RGM_CHECK((lhs >= rhs), error)

/// Check less than
#define RGM_CHECK_LT(lhs, rhs) \
    RGM_CHECK((lhs < rhs), error)

/// Check less than or equal to
#define RGM_CHECK_LE(lhs, rhs) \
    RGM_CHECK((lhs <= rhs), error)

/// Check not null
#define RGM_CHECK_NOTNULL(arg) \
    RGM_CHECK((arg != NULL), error)

} // namespace RGM

#endif // RGM_UTILLOG_HPP_
