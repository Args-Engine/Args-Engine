#pragma once
/**
 * @file platform.hpp
 */

#if defined(_WIN64)
 /**@def LEGION_WINDOWS
  * @brief Defined when compiling for Windows.
  */
#define LEGION_WINDOWS

#define WIN32_LEAN_AND_MEAN
#define VC_EXTRALEAN
#define NOMINMAX
#include <Windows.h>
#include <shellapi.h>
#include <shlobj.h>
#include <shlwapi.h>
#include <objbase.h>
#elif defined(__linux__)
 /**@def LEGION_LINUX
  * @brief Defined when compiling for Linux.
  */
#define LEGION_LINUX

#include <sys/prctl.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <sched.h>
#include <errno.h>
#endif

#if !defined(PROJECT_NAME)
#define PROJECT_NAME user_project
#endif

#define NARGS_(_1, _2, _3, _4, _5 , _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, N, ...) N
#define NARGS(args...) NARGS_(args..., 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1)

#define CONCAT(A, B) A ## B

#define CONCAT_DEFINE(A, B) CONCAT(A, B)

                                                               
#define  L_NAME_1(x)                                                                    #x
#define  L_NAME_2(x, x2)                                                                #x , #x2
#define  L_NAME_3(x, x2, x3)                                                            #x , #x2 , #x3
#define  L_NAME_4(x, x2, x3, x4)                                                        #x , #x2 , #x3 , #x4
#define  L_NAME_5(x, x2, x3, x4, x5)                                                    #x , #x2 , #x3 , #x4 , #x5
#define  L_NAME_6(x, x2, x3, x4, x5, x6)                                                #x , #x2 , #x3 , #x4 , #x5 , #x6
#define  L_NAME_7(x, x2, x3, x4, x5, x6, x7)                                            #x , #x2 , #x3 , #x4 , #x5 , #x6 , #x7
#define  L_NAME_8(x, x2, x3, x4, x5, x6, x7, x8)                                        #x , #x2 , #x3 , #x4 , #x5 , #x6 , #x7 , #x8
#define  L_NAME_9(x, x2, x3, x4, x5, x6, x7, x8, x9)                                    #x , #x2 , #x3 , #x4 , #x5 , #x6 , #x7 , #x8 , #x9
#define L_NAME_10(x, x2, x3, x4, x5, x6, x7, x8, x9, x10)                               #x , #x2 , #x3 , #x4 , #x5 , #x6 , #x7 , #x8 , #x9 , #x10
#define L_NAME_11(x, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11)                          #x , #x2 , #x3 , #x4 , #x5 , #x6 , #x7 , #x8 , #x9 , #x10 , #x11
#define L_NAME_12(x, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12)                     #x , #x2 , #x3 , #x4 , #x5 , #x6 , #x7 , #x8 , #x9 , #x10 , #x11 , #x12
#define L_NAME_13(x, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13)                #x , #x2 , #x3 , #x4 , #x5 , #x6 , #x7 , #x8 , #x9 , #x10 , #x11 , #x12 , #x13
#define L_NAME_14(x, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14)           #x , #x2 , #x3 , #x4 , #x5 , #x6 , #x7 , #x8 , #x9 , #x10 , #x11 , #x12 , #x13 , #x14
#define L_NAME_15(x, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15)      #x , #x2 , #x3 , #x4 , #x5 , #x6 , #x7 , #x8 , #x9 , #x10 , #x11 , #x12 , #x13 , #x14 , #x15
#define L_NAME_16(x, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16) #x , #x2 , #x3 , #x4 , #x5 , #x6 , #x7 , #x8 , #x9 , #x10 , #x11 , #x12 , #x13 , #x14 , #x15 , #x16

#define STRINGIFY_SEPERATE(args...) CONCAT_DEFINE(L_NAME_, NARGS(args))(args)

#define RULE_OF_5(type)\
type() = default;\
type(const type&) = default;\
type(type&&) = default;\
type& operator=(const type&) = default;\
type& operator=(type&&) = default;

/**@def LEGION_CPP17V
 * @brief the version number of c++17 as long
 */
#define LEGION_CPP17V 201703L

#define LEGION_DEBUG_VALUE 1
#define LEGION_RELEASE_VALUE 2

#if defined(DOXY_EXCLUDE)
#define NDOXY(...)
#define CNDOXY(...)
#else
#define NDOXY(args...) args
#define CNDOXY(args...) , args
#endif

#if defined(_DEBUG) || defined(DEBUG)
    /**@def LEGION_DEBUG
     * @brief Defined in debug mode.
     */
    #define LEGION_DEBUG
    #define LEGION_CONFIGURATION LEGION_DEBUG_VALUE
#else
    /**@def LEGION_RELEASE
     * @brief Defined in release mode.
     */
    #define LEGION_RELEASE 
    #define LEGION_CONFIGURATION LEGION_RELEASE_VALUE
#endif

#if (!defined(LEGION_LOW_POWER) && !defined(LEGION_HIGH_PERFORMANCE))
    /**@def LEGION_HIGH_PERFORMANCE
     * @brief Automatically defined if LEGION_LOW_POWER was not defined. It makes Legion ask the hardware's full attention to run as fast as possible.
     * @note Define LEGION_LOW_POWER to run Legion with minimal resources instead.
     */
    #define LEGION_HIGH_PERFORMANCE
#endif

#ifndef __FUNC__
    #define __FUNC__ __func__ 
#endif

#if defined(__clang__)
    // clang
#define LEGION_CLANG

#if defined(__GNUG__) || (defined(__GNUC__) && defined(__cplusplus))
#define LEGION_CLANG_GCC
#elif defined(_MSC_VER)
#define LEGION_CLANG_MSVC
#endif

#define L_PAUSE_INSTRUCTION __builtin_ia32_pause
#elif defined(__GNUG__) || (defined(__GNUC__) && defined(__cplusplus))
    // gcc
#define LEGION_GCC
#define L_PAUSE_INSTRUCTION __builtin_ia32_pause
#elif defined(_MSC_VER)
    // msvc
#define LEGION_MSVC
#define L_PAUSE_INSTRUCTION _mm_pause
#endif

#ifndef __FULL_FUNC__
#if defined(LEGION_CLANG) || defined(LEGION_GCC)
#define __FULL_FUNC__ __PRETTY_FUNCTION__
#elif defined(LEGION_MSVC)
#define __FULL_FUNC__ __FUNCSIG__
#else
#define __FULL_FUNC__ __func__
#endif
#endif

#if (defined(LEGION_WINDOWS) && !defined(LEGION_WINDOWS_USE_CDECL)) || defined (DOXY_INCLUDE)
    /**@def LEGION_CCONV
     * @brief the calling convention exported functions will use in the args engine
     */
    #define LEGION_CCONV __fastcall
#elif defined(LEGION_MSVC)
    #define LEGION_CCONV __cdecl
#else
    #define LEGION_CCONV
#endif

#if defined(LEGION_GCC) || defined(LEGION_CLANG)
#define L_ALWAYS_INLINE __attribute__((always_inline))
#else
#define L_ALWAYS_INLINE __forceinline
#endif

/**@def NO_MANGLING
 * @brief exports functions with C style names instead of C++ mangled names
 */
#define NO_MANGLING extern "C"

/**@def LEGION_FUNC
 * @brief export setting + calling convention used by the engine
 */
#define LEGION_FUNC LEGION_CCONV

/**@def LEGION_INTERFACE
 * @brief un-mangled function name +  export setting + calling convention used by the engine
 */
#define LEGION_INTERFACE NO_MANGLING LEGION_CCONV 

#if defined(__has_cpp_attribute)|| defined(DOXY_INCLUDE) 
/**@def L_HASCPPATTRIB
 * @brief checks if a certain attribute exists in this version of c++
 * @param x attribute you want to test for
 * @return true if attribute exists
 */
#  define L_HASCPPATTRIB(x) __has_cpp_attribute(x)
#else
#  define L_HASCPPATTRIB(x) 0
#endif

#if L_HASCPPATTRIB(fallthrough)
#define L_FALLTHROUGH [[fallthrough]]
#else
#define L_FALLTHROUGH
#endif

#if __cplusplus >= LEGION_CPP17V || L_HASCPPATTRIB(nodiscard) || defined(DOXY_INCLUDE)

/**@def L_NODISCARD
 * @brief Marks a function as "nodiscard" meaning that result must be captured and should not be discarded.
 */
#define L_NODISCARD [[nodiscard]]
#else
#define L_NODISCARD
#endif

#if __cplusplus >= LEGION_CPP17V || L_HASCPPATTRIB(maybe_unused) || defined(DOXY_INCLUDE)

/**@def L_MAYBEUNUSED
 * @brief [[maybe_unused]]
 */
#define L_MAYBEUNUSED [[maybe_unused]]
#else
#define L_MAYBEUNUSED
#endif



#if __cplusplus > LEGION_CPP17V || L_HASCPPATTRIB(noreturn) || defined(DOXY_INCLUDE)
/**@def L_NORETURN
 * @brief Marks a function as "noreturn" meaning that the function will never finish, or terminate the application
 */
#define L_NORETURN [[noreturn]]
#else
#define L_NORETURN
#endif

/**@def LEGION_PURE
 * @brief Marks a function as pure virtual.
 */
#define LEGION_PURE =0

/**@def LEGION_IMPURE
 * @brief Marks a function as overridable but default implemented.
 */
#define LEGION_IMPURE {}

/**@def LEGION_IMPURE_RETURN
 * @brief Marks a function as overridable but default implemented with certain default return value.
 * @param x value the function should return.
 */
#define LEGION_IMPURE_RETURN(x) { return (x); }

#if !defined(LEGION_MIN_THREADS)
#define LEGION_MIN_THREADS 5
#endif
