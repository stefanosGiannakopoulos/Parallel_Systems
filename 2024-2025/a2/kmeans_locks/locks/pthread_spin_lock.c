#include <pthread.h>

#include "alloc.h"
#include "lock.h"
char LOCKNAME[32];

struct lock_struct {
	pthread_spinlock_t spinlock;
};

lock_t *lock_init(int nthreads)
{
	strcpy(LOCKNAME,"pthread-spinlock");
	lock_t *lock;

	XMALLOC(lock, 1);
	pthread_spin_init(&lock->spinlock, PTHREAD_PROCESS_SHARED);
	return lock;
}

void lock_free(lock_t *lock)
{
	pthread_spin_destroy(&lock->spinlock);
	XFREE(lock);
}

void lock_acquire(lock_t *lock)
{
	lock_t *l = lock;

	pthread_spin_lock(&lock->spinlock);
}

void lock_release(lock_t *lock)
{
	lock_t *l = lock;
	pthread_spin_unlock(&lock->spinlock);
}
