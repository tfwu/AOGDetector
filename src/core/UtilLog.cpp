#include "UtilLog.hpp"

namespace  RGM
{

std::ostream& operator<< (std::ostream& strm, severity_level level)
{
    static const char* strings[] = {
        "normal",
        "notification",
        "warning",
        "error",
        "critical"
    };

    if (static_cast< std::size_t >(level) < sizeof(strings) / sizeof(*strings)) {
        strm << strings[level];
    } else {
        strm << static_cast< int >(level);
    }

    return strm;
}

#ifndef RGM_NO_LOG
void log_init()
{
    // Setup the common formatter for all sinks
    logging::formatter fmt = expr::stream
                            << expr::if_(expr::has_attr(tag_attr))
                            [
                                expr::stream << "[" << tag_attr << "] "
                            ]
                             << std::setw(6) << std::setfill('0')
                             << line_id << std::setfill(' ')
                             << ": <" << severity << ">\t"
                             << expr::smessage;

    // Initialize sinks
    typedef sinks::synchronous_sink< sinks::text_ostream_backend > text_sink;
    boost::shared_ptr< text_sink > sink = boost::make_shared< text_sink >();

    sink->locked_backend()->add_stream(
        boost::make_shared< std::ofstream >("full.log"));

    sink->set_formatter(fmt);

    logging::core::get()->add_sink(sink);

    sink = boost::make_shared< text_sink >();

    sink->locked_backend()->add_stream(
        boost::make_shared< std::ofstream >("important.log"));

    sink->set_formatter(fmt);

    sink->set_filter(severity >= warning ||
                     (expr::has_attr(tag_attr) && tag_attr == "IMPORTANT_MESSAGE"));

    logging::core::get()->add_sink(sink);

    // add console
    logging::add_console_log(std::cout,
                             keywords::auto_flush = true);

    // Add attributes
    logging::add_common_attributes();
}
#endif

} // namespace RGM
