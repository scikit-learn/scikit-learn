#if defined(_MSC_VER)
// MS compiler doesn't do standardized _Thread_local until very recent versions:
#define tls_int __declspec(thread) int
#else
// _Thread_local is available in C11 and later:
#define tls_int _Thread_local int
#endif
