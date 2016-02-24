/* 
 
   LibLinear Threads

   Created 2016:
   - Added threading to liblinear - Claes-Fredrik Mannby
 */

#include "threads.h"

#ifdef _LL_THREADS

bool _LL_MutexInit(_LL_Mutex  * const mutex)
{
#if defined(_WIN32)

#if (_WIN32_WINNT >= 0x0403)

    if (InitializeCriticalSectionAndSpinCount((CRITICAL_SECTION *)mutex, 0x400)) {
        return true;
    }
    else {
        return false;
    }

#else // _WIN32_WINNT

    InitializeCriticalSection((CRITICAL_SECTION *)mutex);

    return true;

#endif // _WIN32_WINNT

#else

    if (!pthread_mutex_init((pthread_mutex_t *)mutex, NULL)) {
        return true;
    }
    else {
        return false;
    }

#endif
}

bool _LL_MutexDestroy(_LL_Mutex  * const mutex)
{
#if defined(_WIN32)

    DeleteCriticalSection((CRITICAL_SECTION *)mutex);

    return true;

#else

    if (!pthread_mutex_destroy((pthread_mutex_t *)mutex)) {
        return true;
    }
    else {
        return false;
    }

#endif
}

bool _LL_MutexLock(_LL_Mutex  * const mutex)
{
#if defined(_WIN32)

    EnterCriticalSection((CRITICAL_SECTION *)mutex);

    return true;

#else

    if (!pthread_mutex_lock((pthread_mutex_t *)mutex)) {
        return true;
    }
    else {
        return false;
    }

#endif
}

bool _LL_MutexUnlock(_LL_Mutex  * const mutex)
{
#if defined(_WIN32)

    LeaveCriticalSection((CRITICAL_SECTION *)mutex);

    return true;

#else

    if (!pthread_mutex_unlock((pthread_mutex_t *)mutex)) {
        return true;
    }
    else {
        return false;
    }

#endif
}

bool _LL_THREADS_SemInit(_LL_Semaphore  * const semaphore, int initialValue, int maxValue)
{
#if defined(_WIN32)

    *semaphore = CreateSemaphore(NULL, initialValue, maxValue, NULL);
    if (*semaphore) {
        return true;
    }
    else {
        return false;
    }

#else

    (void)maxValue;
#ifdef __SEM_NAMED
    {
        const size_t bytes = sizeof("/llt-%p") + 2 * sizeof(void *);
        char semName[bytes];
        snprintf(semName, bytes, "/llt-%p", semaphore);
        (void)sem_unlink(semName);
        *semaphore = sem_open(semName, O_CREAT | O_EXCL, S_IRUSR | S_IWUSR, initialValue);
        if (*semaphore != SEM_FAILED) {
            return true;
        }
        else {
            return false;
        }
    }
#else // __SEM_NAMED
    if (!sem_init(semaphore, 0, initialValue)) {
        return true;
    }
    else {
        return false;
    }
#endif // __SEM_NAMED

#endif
}

bool _LL_THREADS_SemWait(_LL_Semaphore  * const semaphore)
{
#if defined(_WIN32)

    if (WaitForSingleObject(*(HANDLE * const)semaphore, INFINITE) == WAIT_OBJECT_0) {
        return true;
    }
    else {
        return false;
    }

#else

#ifdef __SEM_NAMED
    if (!sem_wait(*(sem_t **)semaphore)) {
        return true;
    }
    else {
        return false;
    }
#else // __SEM_NAMED
    if (!sem_wait((sem_t *)semaphore)) {
        return true;
    }
    else {
        return false;
    }
#endif /* __SEM_NAMED */

#endif
}

bool _LL_THREADS_SemPost(_LL_Semaphore  * const semaphore)
{
#if defined(_WIN32)

    if (ReleaseSemaphore(*(HANDLE * const)semaphore, 1, NULL)) {
        return true;
    }
    else {
        return false;
    }

#else

#ifdef __SEM_NAMED
    if (!sem_post(*(sem_t **)semaphore)) {
        return true;
    }
    else {
        return false;
    }
#else /* __SEM_NAMED */
    if (!sem_post((sem_t *)semaphore)) {
        return true;
    }
    else {
        return false;
    }
#endif /* __SEM_NAMED */

#endif
}

bool _LL_THREADS_SemDestroy(_LL_Semaphore  * const semaphore)
{
#if defined(_WIN32)

    if (CloseHandle(*(HANDLE * const)semaphore)) {
        return true;
    }
    else {
        return false;
    }

#else

#ifdef __SEM_NAMED
    {
        const size_t bytes = sizeof("/llt-%p") + 2 * sizeof(void *);
        char semName[bytes];
        snprintf(semName, bytes, "/llt-%p", semaphore);
        if (!sem_unlink(semName)) {
            if (!sem_close(*(sem_t **)semaphore)) {
                return true;
            }
            else {
                return false;
            }
        }
        else {
            return false;
        }
    }
#else // __SEM_NAMED
    if (!sem_destroy((sem_t *)semaphore)) {
        return true;
    }
    else {
        return false;
    }
#endif // __SEM_NAMED

#endif
}

bool _LL_ThreadCreate(_LL_ThreadHandle          * const threadHandle,
                      _LL_ThreadStartFunction           startRoutine,
                      void                      * const startData,
                      int                               stackSize)
{
#if defined(_WIN32)

    *threadHandle = (_LL_ThreadHandle)CreateThread(NULL, stackSize, startRoutine, startData, 0, NULL);

    if (*threadHandle) {
        return true;
    }
    else {
        return false;
    }

#else

    int status;
    pthread_attr_t attr;

    status = pthread_attr_init(&attr);

    if (!status) {
        status = pthread_attr_setstacksize(&attr, stackSize);
    }

    if (!pthread_create((pthread_t *)threadHandle, (status ? NULL : &attr), startRoutine, startData)) {
        return true;
    }
    else {
        return false;
    }

#endif
}

bool _LL_ThreadJoin(const _LL_ThreadHandle threadHandle)
{
    if (!threadHandle) {
        return true;
    }

#if defined(_WIN32)

    {
        DWORD status = WaitForSingleObject((HANDLE)threadHandle, INFINITE);

        if (!CloseHandle((HANDLE)threadHandle)) {
            return false;
        }

        if (status == WAIT_OBJECT_0) {
            return true;
        }
        else {
            return false;
        }
    }
    
#else
    
    if (!pthread_join((pthread_t)threadHandle, NULL)) {
        return true;
    }
    else {
        return false;
    }
    
#endif
}

#endif // _LL_THREADS
