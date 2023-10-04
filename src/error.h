/* @file error.h
**
** error checking macros/functions and error messages

MIT License

Copyright (c) 2018  Hasindu Gamaarachchi (hasindu@unsw.edu.au)
Copyright (c) 2018  Thomas Daniell

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


******************************************************************************/

#ifndef ERROR_H
#define ERROR_H

#include <errno.h>
#include <stdio.h>
#include <string.h>

#define log_level get_log_level()

// the level of verbosity in the log printed to the standard error
enum log_level_opt {
    LOG_OFF,      // nothing at all
    LOG_ERR,      // error messages
    LOG_WARN,     // warning and error messages
    LOG_INFO,     // information, warning and error messages
    LOG_VERB,     // verbose, information, warning and error messages
    LOG_DBUG,     // debugging, verbose, information, warning and error messages
    LOG_TRAC      // tracing, debugging, verbose, information, warning and error messages
};

enum log_level_opt get_log_level();
void set_log_level(enum log_level_opt level);

#define DEBUG_PREFIX "[DEBUG] %s: " /* TODO function before debug */
#define VERBOSE_PREFIX "[INFO] %s: "
#define INFO_PREFIX "[%s::INFO]\033[1;34m "
#define WARNING_PREFIX "[%s::WARNING]\033[1;33m "
#define ERROR_PREFIX "[%s::ERROR]\033[1;31m "
#define NO_COLOUR "\033[0m"

#define LOG_TRACE(msg, ...) { \
    if (log_level >= LOG_TRAC) { \
        fprintf(stderr, DEBUG_PREFIX msg \
                " At %s:%d\n", \
                __func__, __VA_ARGS__, __FILE__, __LINE__ - 1); \
    } \
}

#define LOG_DEBUG(msg, ...) { \
    if (log_level >= LOG_DBUG) { \
        fprintf(stderr, DEBUG_PREFIX msg \
                " At %s:%d\n", \
                __func__, __VA_ARGS__, __FILE__, __LINE__ - 1); \
    } \
}

#define VERBOSE(msg, ...) { \
    if (log_level >= LOG_VERB) { \
        fprintf(stderr, VERBOSE_PREFIX msg "\n", __func__, __VA_ARGS__); \
    } \
}

#define INFO(msg, ...) { \
    if (log_level >= LOG_INFO) { \
        fprintf(stderr, INFO_PREFIX msg NO_COLOUR "\n", __func__, __VA_ARGS__); \
    } \
}

#define WARNING(msg, ...) { \
    if (log_level >= LOG_WARN) { \
        fprintf(stderr, WARNING_PREFIX msg NO_COLOUR \
                " At %s:%d\n", \
                __func__, __VA_ARGS__, __FILE__, __LINE__ - 1); \
    } \
}

#define ERROR(msg, ...) { \
    if (log_level >= LOG_ERR) { \
        fprintf(stderr, ERROR_PREFIX msg NO_COLOUR \
                " At %s:%d\n", \
                __func__, __VA_ARGS__, __FILE__, __LINE__ - 1); \
    } \
}

#define MALLOC_CHK(ret) { \
    if ((ret) == NULL) { \
        MALLOC_ERROR() \
        exit(EXIT_FAILURE); \
    } \
}

#define MALLOC_ERROR() ERROR("Failed to allocate memory: %s", strerror(errno))

#define F_CHK(ret, file) { \
    if ((ret) == NULL) { \
        ERROR("Could not to open file %s: %s", file, strerror(errno)) \
        exit(EXIT_FAILURE); \
    } \
}

#define NULL_CHK(ret) { \
    if ((ret) == NULL) { \
        ERROR("NULL returned: %s.", strerror(errno)) \
        exit(EXIT_FAILURE); \
    } \
}

#define NEG_CHK(ret) { \
    if ((ret) < 0) { \
        ERROR("Negative value returned: %s.", strerror(errno)) \
        exit(EXIT_FAILURE); \
    } \
}

#define ASSERT(ret) { \
    if ((ret) == 0){ \
        fprintf(stderr, ERROR_PREFIX "Assertion failed." NO_COLOUR \
                " At %s:%d\nExiting.\n", \
                __func__ , __FILE__, __LINE__ - 1); \
        exit(EXIT_FAILURE); \
    } \
}

#endif
