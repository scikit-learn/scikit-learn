typedef int omp_lock_t;

void omp_init_lock(omp_lock_t *lock) {}

void omp_destroy_lock(omp_lock_t *lock) {}

void omp_set_lock(omp_lock_t *lock) {}

void omp_unset_lock(omp_lock_t *lock) {}

int omp_get_thread_num() { return 0; }

int omp_get_max_threads() { return 1; }
