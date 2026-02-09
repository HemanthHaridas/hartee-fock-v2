#pragma once

#include <string>

enum class LogLevel
{
    Info,
    Error
};

/// Thread-safe, formatted logging
void logging(LogLevel level, const std::string &label, const std::string &message);
