#include "alloc.h"
#include "lock.h"
char LOCKNAME[32];

#define FALSE 0
#define TRUE  1

#define FLAG_ENTRY_SIZE 1024

struct lock_struct {
	volatile char *flag;
	char padding1[64 - sizeof(char *)];
	unsigned long long tail;
	char padding2[64 - sizeof(unsigned long long)];
	int size;
};

lock_t *lock_init(int nthreads)
{
	strcpy(LOCKNAME,"array-based");
	lock_t *lock;
	int i;

	XMALLOC(lock, 1);

	lock->size = nthreads * FLAG_ENTRY_SIZE;
	lock->tail = 0;
	XMALLOC(lock->flag, lock->size);

	for (i=0; i < lock->size; i++)
		lock->flag[i] = FALSE;
	lock->flag[0] = TRUE;

	// printf("lock->size = %d\n", lock->size);
	return lock;
}

void lock_free(lock_t *lock)
{
	lock_t *l = lock;
	XFREE((void*)l->flag);
	XFREE(l);
}

__thread int mySlot;

void lock_acquire(lock_t *lock)
{
	lock_t *l = lock;
	int slot = __sync_fetch_and_add(&l->tail, FLAG_ENTRY_SIZE) % l->size;
	
	mySlot = slot;
	// if (mySlot < 0) {
	// printf("mySlot = %d\n", mySlot);
	// printf("Tail = %d\n", l->tail);
	// }

	while (l->flag[slot] == FALSE)
		/* do nothing */ ;
}

void lock_release(lock_t *lock)
{
	lock_t *l = lock;
	int slot = mySlot;

	l->flag[slot] = FALSE;
	l->flag[(slot+FLAG_ENTRY_SIZE) % l->size] = TRUE;
}
