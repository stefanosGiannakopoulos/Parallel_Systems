#include <pthread.h>

#include "alloc.h"
#include "lock.h"
char LOCKNAME[32];

struct lock_struct {
	pthread_mutex_t mutex;
};

lock_t *lock_init(int nthreads)
{
	strcpy(LOCKNAME,"pthread-mutex");
	lock_t *lock;

	XMALLOC(lock, 1);
	pthread_mutex_init(&lock->mutex, NULL);
	return lock;
}

void lock_free(lock_t *lock)
{
	pthread_mutex_destroy(&lock->mutex);
	XFREE(lock);
}

void lock_acquire(lock_t *lock)
{
	lock_t *l = lock;

	pthread_mutex_lock(&lock->mutex);
}

void lock_release(lock_t *lock)
{
	lock_t *l = lock;
	pthread_mutex_unlock(&lock->mutex);
}
