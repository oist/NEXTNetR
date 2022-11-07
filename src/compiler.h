#pragma once

#if defined(__clang__)
    #define WARNINGS_DISABLE \
        _Pragma("clang diagnostic push") \
        _Pragma("clang diagnostic ignored \"-Weverything\"") \
        _Pragma("clang diagnostic ignored \"-Wdeprecated-declarations\"")
#elif defined(__GNUC__) || defined(__GNUG__)
    #define WARNINGS_DISABLE \
        _Pragma("GCC diagnostic push") \
        _Pragma("GCC diagnostic ignored \"-Wall\"") \
        _Pragma("GCC diagnostic ignored \"-Wsign-compare\"") \
        _Pragma("GCC diagnostic ignored \"-Wparentheses\"") \
        _Pragma("GCC diagnostic ignored \"-Wdeprecated-declarations\"")
#else
    #define WARNINGS_DISABLE
#endif


#if defined(__clang__)
    #define WARNINGS_ENABLE \
      _Pragma("clang diagnostic pop")
#elif defined(__GNUC__) || defined(__GNUG__)
    #define WARNINGS_ENABLE \
      _Pragma("GCC diagnostic pop")
#else
    #define WARNINGS_ENABLE
#endif
