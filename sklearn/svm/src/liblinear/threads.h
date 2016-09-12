#ifndef _LL_THREADS_H
#define _LL_THREADS_H

#ifdef __cplusplus
extern "C" {
#endif

#define _LL_THREADS

#ifdef _LL_THREADS

    // pick semaphore implementation
#ifdef __APPLE__
#define __SEM_NAMED
#else
#define __SEM_MEMORY
#endif

#ifdef __SEM_NAMED
#include "stdio.h"
#endif

#if defined(_WIN32)
#include <windows.h>
#else
#include <pthread.h>
#include <semaphore.h>
#endif

#if defined(_WIN32)
    typedef CRITICAL_SECTION   _LL_Mutex;
    typedef HANDLE             _LL_Semaphore;
    typedef DWORD              _LL_ThreadReturnValue;
    #define _LL_STDCALL        _stdcall
#else
    typedef void               *_LL_ThreadReturnValue;
    typedef pthread_mutex_t    _LL_Mutex;
    #define _LL_STDCALL
#ifdef __SEM_NAMED
    typedef sem_t              *_LL_Semaphore;
#else
    typedef sem_t              _LL_Semaphore;
#endif
#endif

    typedef void *_LL_ThreadHandle;
    typedef _LL_ThreadReturnValue (_LL_STDCALL *_LL_ThreadStartFunction) (void *startData);

    bool _LL_MutexInit(_LL_Mutex  * const mutex);
    bool _LL_MutexDestroy(_LL_Mutex  * const mutex);
    bool _LL_MutexLock(_LL_Mutex  * const mutex);
    bool _LL_MutexUnlock(_LL_Mutex  * const mutex);
    bool _LL_THREADS_SemInit(_LL_Semaphore  * const semaphore, int initialValue, int maxValue);
    bool _LL_THREADS_SemWait(_LL_Semaphore  * const semaphore);
    bool _LL_THREADS_SemPost(_LL_Semaphore  * const semaphore);
    bool _LL_THREADS_SemDestroy(_LL_Semaphore  * const semaphore);
    bool _LL_ThreadCreate(_LL_ThreadHandle          * const threadHandle,
                          _LL_ThreadStartFunction           startRoutine,
                          void                      * const startData,
                          int                               stackSize);
    bool _LL_ThreadJoin(const _LL_ThreadHandle threadHandle);


#endif // _LL_THREADS

#ifdef __cplusplus
}
#endif

#endif /* _LL_THREADS_H */

