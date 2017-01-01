#if !defined(INCLUDED_PETSCCONF_H)
#define INCLUDED_PETSCCONF_H

#ifndef IS_COLORING_MAX
#define IS_COLORING_MAX 65535
#endif

#ifndef MPI_Type_create_struct
#define MPI_Type_create_struct(count,lens,displs,types,newtype) MPI_Type_struct((count),(lens),(displs),(types),(newtype))
#endif

#ifndef STDC_HEADERS
#define STDC_HEADERS 1
#endif

#ifndef MPI_Comm_create_errhandler
#define MPI_Comm_create_errhandler(p_err_fun,p_errhandler) MPI_Errhandler_create((p_err_fun),(p_errhandler))
#endif

#ifndef MPI_Comm_set_errhandler
#define MPI_Comm_set_errhandler(comm,p_errhandler) MPI_Errhandler_set((comm),(p_errhandler))
#endif

#ifndef MPIU_COLORING_VALUE
#define MPIU_COLORING_VALUE MPI_UNSIGNED_SHORT
#endif

#ifndef PETSC_UINTPTR_T
#define PETSC_UINTPTR_T uintptr_t
#endif

#ifndef PETSC_HAVE_PTHREAD
#define PETSC_HAVE_PTHREAD 1
#endif

#ifndef PETSC_DEPRECATED
#define PETSC_DEPRECATED(why) __attribute((deprecated))
#endif

#ifndef PETSC_REPLACE_DIR_SEPARATOR
#define PETSC_REPLACE_DIR_SEPARATOR '\\'
#endif

#ifndef PETSC_HAVE_SO_REUSEADDR
#define PETSC_HAVE_SO_REUSEADDR 1
#endif

#ifndef PETSC_HAVE_MPI
#define PETSC_HAVE_MPI 1
#endif

#ifndef PETSC_PREFETCH_HINT_T2
#define PETSC_PREFETCH_HINT_T2 _MM_HINT_T2
#endif

#ifndef PETSC_PREFETCH_HINT_T0
#define PETSC_PREFETCH_HINT_T0 _MM_HINT_T0
#endif

#ifndef PETSC_PREFETCH_HINT_T1
#define PETSC_PREFETCH_HINT_T1 _MM_HINT_T1
#endif

#ifndef PETSC_USE_SOCKET_VIEWER
#define PETSC_USE_SOCKET_VIEWER 1
#endif

#ifndef PETSC_DIR
#define PETSC_DIR "/root/petsc-3.7.5"
#endif

#ifndef PETSC_LIB_DIR
#define PETSC_LIB_DIR "/root/petsc-3.7.5/arch-linux2-c-opt/lib"
#endif

#ifndef PETSC_HAVE_FORTRAN
#define PETSC_HAVE_FORTRAN 1
#endif

#ifndef PETSC_USE_ISATTY
#define PETSC_USE_ISATTY 1
#endif

#ifndef PETSC_SLSUFFIX
#define PETSC_SLSUFFIX "so"
#endif

#ifndef PETSC_FUNCTION_NAME_CXX
#define PETSC_FUNCTION_NAME_CXX __func__
#endif

#ifndef PETSC_HAVE_FBLASLAPACK
#define PETSC_HAVE_FBLASLAPACK 1
#endif

#ifndef PETSC_HAVE_ATOLL
#define PETSC_HAVE_ATOLL 1
#endif

#ifndef PETSC_HAVE_DOUBLE_ALIGN_MALLOC
#define PETSC_HAVE_DOUBLE_ALIGN_MALLOC 1
#endif

#ifndef PETSC_UNUSED
#define PETSC_UNUSED __attribute((unused))
#endif

#ifndef PETSC_FUNCTION_NAME_C
#define PETSC_FUNCTION_NAME_C __func__
#endif

#ifndef PETSC_USE_SINGLE_LIBRARY
#define PETSC_USE_SINGLE_LIBRARY 1
#endif

#ifndef PETSC_HAVE_BUILTIN_EXPECT
#define PETSC_HAVE_BUILTIN_EXPECT 1
#endif

#ifndef PETSC_DIR_SEPARATOR
#define PETSC_DIR_SEPARATOR '/'
#endif

#ifndef PETSC_PATH_SEPARATOR
#define PETSC_PATH_SEPARATOR ':'
#endif

#ifndef PETSC_HAVE_XMMINTRIN_H
#define PETSC_HAVE_XMMINTRIN_H 1
#endif

#ifndef PETSC_PREFETCH_HINT_NTA
#define PETSC_PREFETCH_HINT_NTA _MM_HINT_NTA
#endif

#ifndef PETSC_Prefetch
#define PETSC_Prefetch(a,b,c) _mm_prefetch((const char*)(a),(c))
#endif

#ifndef PETSC_HAVE_BLASLAPACK
#define PETSC_HAVE_BLASLAPACK 1
#endif

#ifndef PETSC_HAVE_GZIP
#define PETSC_HAVE_GZIP 1
#endif

#ifndef PETSC_HAVE_STRING_H
#define PETSC_HAVE_STRING_H 1
#endif

#ifndef PETSC_HAVE_SYS_TYPES_H
#define PETSC_HAVE_SYS_TYPES_H 1
#endif

#ifndef PETSC_HAVE_ENDIAN_H
#define PETSC_HAVE_ENDIAN_H 1
#endif

#ifndef PETSC_HAVE_SYS_PROCFS_H
#define PETSC_HAVE_SYS_PROCFS_H 1
#endif

#ifndef PETSC_HAVE_DLFCN_H
#define PETSC_HAVE_DLFCN_H 1
#endif

#ifndef PETSC_HAVE_SCHED_H
#define PETSC_HAVE_SCHED_H 1
#endif

#ifndef PETSC_HAVE_STDINT_H
#define PETSC_HAVE_STDINT_H 1
#endif

#ifndef PETSC_HAVE_LINUX_KERNEL_H
#define PETSC_HAVE_LINUX_KERNEL_H 1
#endif

#ifndef PETSC_HAVE_TIME_H
#define PETSC_HAVE_TIME_H 1
#endif

#ifndef PETSC_HAVE_MATH_H
#define PETSC_HAVE_MATH_H 1
#endif

#ifndef PETSC_TIME_WITH_SYS_TIME
#define PETSC_TIME_WITH_SYS_TIME 1
#endif

#ifndef PETSC_HAVE_SYS_PARAM_H
#define PETSC_HAVE_SYS_PARAM_H 1
#endif

#ifndef PETSC_HAVE_PTHREAD_H
#define PETSC_HAVE_PTHREAD_H 1
#endif

#ifndef PETSC_HAVE_UNISTD_H
#define PETSC_HAVE_UNISTD_H 1
#endif

#ifndef PETSC_HAVE_STDLIB_H
#define PETSC_HAVE_STDLIB_H 1
#endif

#ifndef PETSC_HAVE_SYS_WAIT_H
#define PETSC_HAVE_SYS_WAIT_H 1
#endif

#ifndef PETSC_HAVE_SETJMP_H
#define PETSC_HAVE_SETJMP_H 1
#endif

#ifndef PETSC_HAVE_LIMITS_H
#define PETSC_HAVE_LIMITS_H 1
#endif

#ifndef PETSC_HAVE_SYS_UTSNAME_H
#define PETSC_HAVE_SYS_UTSNAME_H 1
#endif

#ifndef PETSC_HAVE_NETINET_IN_H
#define PETSC_HAVE_NETINET_IN_H 1
#endif

#ifndef PETSC_HAVE_SYS_SOCKET_H
#define PETSC_HAVE_SYS_SOCKET_H 1
#endif

#ifndef PETSC_HAVE_FLOAT_H
#define PETSC_HAVE_FLOAT_H 1
#endif

#ifndef PETSC_HAVE_SEARCH_H
#define PETSC_HAVE_SEARCH_H 1
#endif

#ifndef PETSC_HAVE_SYS_RESOURCE_H
#define PETSC_HAVE_SYS_RESOURCE_H 1
#endif

#ifndef PETSC_HAVE_SYS_TIMES_H
#define PETSC_HAVE_SYS_TIMES_H 1
#endif

#ifndef PETSC_HAVE_NETDB_H
#define PETSC_HAVE_NETDB_H 1
#endif

#ifndef PETSC_HAVE_MALLOC_H
#define PETSC_HAVE_MALLOC_H 1
#endif

#ifndef PETSC_HAVE_PWD_H
#define PETSC_HAVE_PWD_H 1
#endif

#ifndef PETSC_HAVE_FCNTL_H
#define PETSC_HAVE_FCNTL_H 1
#endif

#ifndef PETSC_HAVE_STRINGS_H
#define PETSC_HAVE_STRINGS_H 1
#endif

#ifndef PETSC_HAVE_SYS_SYSINFO_H
#define PETSC_HAVE_SYS_SYSINFO_H 1
#endif

#ifndef PETSC_HAVE_SYS_TIME_H
#define PETSC_HAVE_SYS_TIME_H 1
#endif

#ifndef PETSC_USING_F90
#define PETSC_USING_F90 1
#endif

#ifndef PETSC_USING_F2003
#define PETSC_USING_F2003 1
#endif

#ifndef PETSC_HAVE_RTLD_NOW
#define PETSC_HAVE_RTLD_NOW 1
#endif

#ifndef PETSC_HAVE_RTLD_LOCAL
#define PETSC_HAVE_RTLD_LOCAL 1
#endif

#ifndef PETSC_HAVE_RTLD_LAZY
#define PETSC_HAVE_RTLD_LAZY 1
#endif

#ifndef PETSC_C_STATIC_INLINE
#define PETSC_C_STATIC_INLINE static inline
#endif

#ifndef PETSC_HAVE_FORTRAN_UNDERSCORE
#define PETSC_HAVE_FORTRAN_UNDERSCORE 1
#endif

#ifndef PETSC_HAVE_CXX_NAMESPACE
#define PETSC_HAVE_CXX_NAMESPACE 1
#endif

#ifndef PETSC_HAVE_RTLD_GLOBAL
#define PETSC_HAVE_RTLD_GLOBAL 1
#endif

#ifndef PETSC_C_RESTRICT
#define PETSC_C_RESTRICT  __restrict__
#endif

#ifndef PETSC_CXX_RESTRICT
#define PETSC_CXX_RESTRICT  __restrict__
#endif

#ifndef PETSC_CXX_STATIC_INLINE
#define PETSC_CXX_STATIC_INLINE static inline
#endif

#ifndef PETSC_HAVE_LIBDL
#define PETSC_HAVE_LIBDL 1
#endif

#ifndef PETSC_HAVE_LIBM
#define PETSC_HAVE_LIBM 1
#endif

#ifndef PETSC_HAVE_LIBPTHREAD
#define PETSC_HAVE_LIBPTHREAD 1
#endif

#ifndef PETSC_HAVE_FENV_H
#define PETSC_HAVE_FENV_H 1
#endif

#ifndef PETSC_HAVE_LOG2
#define PETSC_HAVE_LOG2 1
#endif

#ifndef PETSC_HAVE_LIBFLAPACK
#define PETSC_HAVE_LIBFLAPACK 1
#endif

#ifndef PETSC_HAVE_ERF
#define PETSC_HAVE_ERF 1
#endif

#ifndef PETSC_HAVE_LIBFBLAS
#define PETSC_HAVE_LIBFBLAS 1
#endif

#ifndef PETSC_HAVE_LIBQUADMATH
#define PETSC_HAVE_LIBQUADMATH 1
#endif

#ifndef PETSC_HAVE_TGAMMA
#define PETSC_HAVE_TGAMMA 1
#endif

#ifndef PETSC_HAVE_LIBGFORTRAN
#define PETSC_HAVE_LIBGFORTRAN 1
#endif

#ifndef PETSC_ARCH
#define PETSC_ARCH "arch-linux2-c-opt"
#endif

#ifndef PETSC_CLANGUAGE_C
#define PETSC_CLANGUAGE_C 1
#endif

#ifndef PETSC_USE_ERRORCHECKING
#define PETSC_USE_ERRORCHECKING 1
#endif

#ifndef PETSC_BITS_PER_BYTE
#define PETSC_BITS_PER_BYTE 8
#endif

#ifndef PETSC_SIZEOF_VOID_P
#define PETSC_SIZEOF_VOID_P 8
#endif

#ifndef PETSC_RETSIGTYPE
#define PETSC_RETSIGTYPE void
#endif

#ifndef PETSC_HAVE_CXX_COMPLEX
#define PETSC_HAVE_CXX_COMPLEX 1
#endif

#ifndef PETSC_SIZEOF_LONG
#define PETSC_SIZEOF_LONG 8
#endif

#ifndef PETSC_USE_FORTRANKIND
#define PETSC_USE_FORTRANKIND 1
#endif

#ifndef PETSC_SIZEOF_SIZE_T
#define PETSC_SIZEOF_SIZE_T 8
#endif

#ifndef PETSC_HAVE_SIGINFO_T
#define PETSC_HAVE_SIGINFO_T 1
#endif

#ifndef PETSC_SIZEOF_CHAR
#define PETSC_SIZEOF_CHAR 1
#endif

#ifndef PETSC_SIZEOF_DOUBLE
#define PETSC_SIZEOF_DOUBLE 8
#endif

#ifndef PETSC_SIZEOF_FLOAT
#define PETSC_SIZEOF_FLOAT 4
#endif

#ifndef PETSC_HAVE_C99_COMPLEX
#define PETSC_HAVE_C99_COMPLEX 1
#endif

#ifndef PETSC_SIZEOF_INT
#define PETSC_SIZEOF_INT 4
#endif

#ifndef PETSC_SIZEOF_LONG_LONG
#define PETSC_SIZEOF_LONG_LONG 8
#endif

#ifndef PETSC_SIZEOF_SHORT
#define PETSC_SIZEOF_SHORT 2
#endif

#ifndef PETSC_HAVE_STRCASECMP
#define PETSC_HAVE_STRCASECMP 1
#endif

#ifndef PETSC_HAVE_GET_NPROCS
#define PETSC_HAVE_GET_NPROCS 1
#endif

#ifndef PETSC_HAVE_POPEN
#define PETSC_HAVE_POPEN 1
#endif

#ifndef PETSC_HAVE_SIGSET
#define PETSC_HAVE_SIGSET 1
#endif

#ifndef PETSC_HAVE_GETWD
#define PETSC_HAVE_GETWD 1
#endif

#ifndef PETSC_HAVE_VSNPRINTF
#define PETSC_HAVE_VSNPRINTF 1
#endif

#ifndef PETSC_HAVE_TIMES
#define PETSC_HAVE_TIMES 1
#endif

#ifndef PETSC_HAVE_DLSYM
#define PETSC_HAVE_DLSYM 1
#endif

#ifndef PETSC_HAVE_SNPRINTF
#define PETSC_HAVE_SNPRINTF 1
#endif

#ifndef PETSC_HAVE_GETHOSTBYNAME
#define PETSC_HAVE_GETHOSTBYNAME 1
#endif

#ifndef PETSC_HAVE_GETCWD
#define PETSC_HAVE_GETCWD 1
#endif

#ifndef PETSC_HAVE_DLERROR
#define PETSC_HAVE_DLERROR 1
#endif

#ifndef PETSC_HAVE_FORK
#define PETSC_HAVE_FORK 1
#endif

#ifndef PETSC_HAVE_RAND
#define PETSC_HAVE_RAND 1
#endif

#ifndef PETSC_HAVE_GETTIMEOFDAY
#define PETSC_HAVE_GETTIMEOFDAY 1
#endif

#ifndef PETSC_HAVE_DLCLOSE
#define PETSC_HAVE_DLCLOSE 1
#endif

#ifndef PETSC_HAVE_UNAME
#define PETSC_HAVE_UNAME 1
#endif

#ifndef PETSC_HAVE_GETHOSTNAME
#define PETSC_HAVE_GETHOSTNAME 1
#endif

#ifndef PETSC_HAVE_MKSTEMP
#define PETSC_HAVE_MKSTEMP 1
#endif

#ifndef PETSC_HAVE_SIGACTION
#define PETSC_HAVE_SIGACTION 1
#endif

#ifndef PETSC_HAVE_DRAND48
#define PETSC_HAVE_DRAND48 1
#endif

#ifndef PETSC_HAVE_MEMALIGN
#define PETSC_HAVE_MEMALIGN 1
#endif

#ifndef PETSC_HAVE_VA_COPY
#define PETSC_HAVE_VA_COPY 1
#endif

#ifndef PETSC_HAVE_CLOCK
#define PETSC_HAVE_CLOCK 1
#endif

#ifndef PETSC_HAVE_ACCESS
#define PETSC_HAVE_ACCESS 1
#endif

#ifndef PETSC_HAVE_SIGNAL
#define PETSC_HAVE_SIGNAL 1
#endif

#ifndef PETSC_HAVE_USLEEP
#define PETSC_HAVE_USLEEP 1
#endif

#ifndef PETSC_HAVE_GETRUSAGE
#define PETSC_HAVE_GETRUSAGE 1
#endif

#ifndef PETSC_HAVE_VFPRINTF
#define PETSC_HAVE_VFPRINTF 1
#endif

#ifndef PETSC_HAVE_NANOSLEEP
#define PETSC_HAVE_NANOSLEEP 1
#endif

#ifndef PETSC_HAVE_GETDOMAINNAME
#define PETSC_HAVE_GETDOMAINNAME 1
#endif

#ifndef PETSC_HAVE_TIME
#define PETSC_HAVE_TIME 1
#endif

#ifndef PETSC_HAVE_LSEEK
#define PETSC_HAVE_LSEEK 1
#endif

#ifndef PETSC_HAVE_SOCKET
#define PETSC_HAVE_SOCKET 1
#endif

#ifndef PETSC_HAVE_SYSINFO
#define PETSC_HAVE_SYSINFO 1
#endif

#ifndef PETSC_HAVE_READLINK
#define PETSC_HAVE_READLINK 1
#endif

#ifndef PETSC_HAVE_REALPATH
#define PETSC_HAVE_REALPATH 1
#endif

#ifndef PETSC_HAVE_DLOPEN
#define PETSC_HAVE_DLOPEN 1
#endif

#ifndef PETSC_HAVE_MEMMOVE
#define PETSC_HAVE_MEMMOVE 1
#endif

#ifndef PETSC_HAVE__GFORTRAN_IARGC
#define PETSC_HAVE__GFORTRAN_IARGC 1
#endif

#ifndef PETSC_SIGNAL_CAST
#define PETSC_SIGNAL_CAST  
#endif

#ifndef PETSC_HAVE_SLEEP
#define PETSC_HAVE_SLEEP 1
#endif

#ifndef PETSC_HAVE_VPRINTF
#define PETSC_HAVE_VPRINTF 1
#endif

#ifndef PETSC_HAVE_BZERO
#define PETSC_HAVE_BZERO 1
#endif

#ifndef PETSC_HAVE_GETPAGESIZE
#define PETSC_HAVE_GETPAGESIZE 1
#endif

#ifndef PETSC_WRITE_MEMORY_BARRIER
#define PETSC_WRITE_MEMORY_BARRIER() asm volatile("sfence":::"memory")
#endif

#ifndef PETSC_MEMORY_BARRIER
#define PETSC_MEMORY_BARRIER() asm volatile("mfence":::"memory")
#endif

#ifndef PETSC_READ_MEMORY_BARRIER
#define PETSC_READ_MEMORY_BARRIER() asm volatile("lfence":::"memory")
#endif

#ifndef PETSC_CPU_RELAX
#define PETSC_CPU_RELAX() asm volatile("rep; nop" ::: "memory")
#endif

#ifndef PETSC_BLASLAPACK_UNDERSCORE
#define PETSC_BLASLAPACK_UNDERSCORE 1
#endif

#ifndef PETSC_HAVE_MPIUNI
#define PETSC_HAVE_MPIUNI 1
#endif

#ifndef PETSC_HAVE_MPI_COMM_C2F
#define PETSC_HAVE_MPI_COMM_C2F 1
#endif

#ifndef PETSC_HAVE_MPI_IN_PLACE
#define PETSC_HAVE_MPI_IN_PLACE 1
#endif

#ifndef PETSC_HAVE_MPI_COMM_F2C
#define PETSC_HAVE_MPI_COMM_F2C 1
#endif

#ifndef PETSC_HAVE_MPI_FINT
#define PETSC_HAVE_MPI_FINT 1
#endif

#ifndef PETSC_USE_INFO
#define PETSC_USE_INFO 1
#endif

#ifndef PETSC_Alignx
#define PETSC_Alignx(a,b)   
#endif

#ifndef PETSC_USE_BACKWARD_LOOP
#define PETSC_USE_BACKWARD_LOOP 1
#endif

#ifndef PETSC_IS_COLOR_VALUE_TYPE_SIZE
#define PETSC_IS_COLOR_VALUE_TYPE_SIZE PETSC_SIZEOF_SHORT
#endif

#ifndef PETSC_USE_MALLOC_COALESCED
#define PETSC_USE_MALLOC_COALESCED 1
#endif

#ifndef PETSC_USE_LOG
#define PETSC_USE_LOG 1
#endif

#ifndef PETSC_IS_COLOR_VALUE_TYPE
#define PETSC_IS_COLOR_VALUE_TYPE short
#endif

#ifndef PETSC_USE_CTABLE
#define PETSC_USE_CTABLE 1
#endif

#ifndef PETSC_MEMALIGN
#define PETSC_MEMALIGN 16
#endif

#ifndef PETSC_USE_SCALAR_REAL
#define PETSC_USE_SCALAR_REAL 1
#endif

#ifndef PETSC_HAVE_ISINF
#define PETSC_HAVE_ISINF 1
#endif

#ifndef PETSC_HAVE_ISNAN
#define PETSC_HAVE_ISNAN 1
#endif

#ifndef PETSC_HAVE_ISNORMAL
#define PETSC_HAVE_ISNORMAL 1
#endif

#ifndef PETSC_USE_REAL_DOUBLE
#define PETSC_USE_REAL_DOUBLE 1
#endif

#ifndef PETSC_HAVE_DYNAMIC_LIBRARIES
#define PETSC_HAVE_DYNAMIC_LIBRARIES 1
#endif

#ifndef PETSC_HAVE_SHARED_LIBRARIES
#define PETSC_HAVE_SHARED_LIBRARIES 1
#endif

#ifndef PETSC_USE_SHARED_LIBRARIES
#define PETSC_USE_SHARED_LIBRARIES 1
#endif

#ifndef PETSC_LEVEL1_DCACHE_LINESIZE
#define PETSC_LEVEL1_DCACHE_LINESIZE 64
#endif

#ifndef PETSC_LEVEL1_DCACHE_SIZE
#define PETSC_LEVEL1_DCACHE_SIZE 32768
#endif

#ifndef PETSC_LEVEL1_DCACHE_ASSOC
#define PETSC_LEVEL1_DCACHE_ASSOC 8
#endif

#ifndef PETSC_HAVE_FORTRAN_GET_COMMAND_ARGUMENT
#define PETSC_HAVE_FORTRAN_GET_COMMAND_ARGUMENT 1
#endif

#ifndef PETSC_HAVE_GFORTRAN_IARGC
#define PETSC_HAVE_GFORTRAN_IARGC 1
#endif

#ifndef PETSC_USE_PROC_FOR_SIZE
#define PETSC_USE_PROC_FOR_SIZE 1
#endif

#ifndef PETSC__GNU_SOURCE
#define PETSC__GNU_SOURCE 1
#endif

#ifndef PETSC__BSD_SOURCE
#define PETSC__BSD_SOURCE 1
#endif

#ifndef PETSC__DEFAULT_SOURCE
#define PETSC__DEFAULT_SOURCE 1
#endif

#ifndef PETSC_HAVE_SCHED_CPU_SET_T
#define PETSC_HAVE_SCHED_CPU_SET_T 1
#endif

#ifndef PETSC_HAVE_PTHREAD_BARRIER_T
#define PETSC_HAVE_PTHREAD_BARRIER_T 1
#endif

#ifndef PETSC_HAVE_SYS_SYSCTL_H
#define PETSC_HAVE_SYS_SYSCTL_H 1
#endif

#endif
