#include "logging.h"

#include <format>
#include <iomanip>
#include <iostream>
#include <mutex>

namespace
{
    // set a mutex to enable OpenMPI writes
    std::mutex log_mutex;

    constexpr const char *info_prefix = "[Planck][INF] => ";
    constexpr const char *error_prefix = "[Planck][ERR] => ";
}

void logging(LogLevel level, const std::string &label, const std::string &message)
{
    std::lock_guard<std::mutex> lock(log_mutex);
    const char *prefix = (level == LogLevel::Info) ? info_prefix : error_prefix;
    std::ostream &out_stream = (level == LogLevel::Info) ? std::cout : std::cerr;
    out_stream << std::setw(20) << std::left << prefix << std::setw(35) << std::left << label << message << '\n';
}
