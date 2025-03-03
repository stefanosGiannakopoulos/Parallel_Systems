#include <stdio.h>
#include <stdlib.h> /* rand() */
#include <limits.h>
#include <pthread.h> /* for pthread_spinlock_t */

#include "../lib/alloc.h"
#include "ll.h"

typedef struct ll_node {
	int key;
	struct ll_node *next;
	pthread_spinlock_t lock;
	short int marked;
} ll_node_t;

struct linked_list {
	ll_node_t *head;
};

/**
 * Create a new linked list node.
 **/
static ll_node_t *ll_node_new(int key)
{
	ll_node_t *ret;

	XMALLOC(ret, 1);
	ret->key = key;
	ret->next = NULL;
	pthread_spin_init(&ret->lock, PTHREAD_PROCESS_SHARED);
	ret->marked = 0;

	return ret;
}

/**
 * Free a linked list node.
 **/
static void ll_node_free(ll_node_t *ll_node)
{
	XFREE(ll_node);
}

/**
 * Create a new empty linked list.
 **/
ll_t *ll_new()
{
	ll_t *ret;

	XMALLOC(ret, 1);
	ret->head = ll_node_new(-1);
	ret->head->next = ll_node_new(INT_MAX);
	ret->head->next->next = NULL;

	return ret;
}

/**
 * Free a linked list and all its contained nodes.
 **/
void ll_free(ll_t *ll)
{
	ll_node_t *next, *curr = ll->head;
	while (curr) {
		next = curr->next;
		ll_node_free(curr);
		curr = next;
	}
	XFREE(ll);
}

#define LOCK_NODE(node) pthread_spin_lock(&(node)->lock)
#define UNLOCK_NODE(node) pthread_spin_unlock(&(node)->lock)

#define TRAVERSE_LIST() \
	do { \
		curr = ll->head; \
		next = curr->next; \
		 \
		while (next->key < key) { \
			curr = next; \
			next = curr->next; \
		} \
	} while (0)

static int validate(ll_node_t *curr, ll_node_t *next)
{
	return (!curr->marked && !next->marked && curr->next == next);
}

int ll_contains(ll_t *ll, int key)
{
	ll_node_t *curr, *next;

	TRAVERSE_LIST();

	return (next->key == key && !next->marked);
}

int ll_add(ll_t *ll, int key)
{
	int ret = 0;
	ll_node_t *curr, *next;
	ll_node_t *new_node;

	do {
		ret = 0;
		curr = next = NULL;

		TRAVERSE_LIST();

		LOCK_NODE(curr);
		LOCK_NODE(next);
		if (validate(curr, next)) {
			if (key != next->key) {
				ret = 1;
				new_node = ll_node_new(key);
				new_node->next = next;
				curr->next = new_node;
				UNLOCK_NODE(curr);
				UNLOCK_NODE(next);
				break;
			} else {
				UNLOCK_NODE(curr);
				UNLOCK_NODE(next);
				break;
			}
		}
		UNLOCK_NODE(curr);
		UNLOCK_NODE(next);
	} while (1);

	return ret;
}

int ll_remove(ll_t *ll, int key)
{
	int ret = 0;
	ll_node_t *curr, *next;

	do {
		ret = 0;
		curr = next = NULL;

		TRAVERSE_LIST();

		LOCK_NODE(curr);
		LOCK_NODE(next);

		if (validate(curr, next)) {
			if (key == next->key) {
				ret = 1;
				next->marked = 1;
				curr->next = next->next;
				UNLOCK_NODE(curr);
				UNLOCK_NODE(next);
//				ll_node_free(next);
				break;
			} else {
				UNLOCK_NODE(curr);
				UNLOCK_NODE(next);
				break;
			}
		}

		UNLOCK_NODE(curr);
		UNLOCK_NODE(next);
	} while (1);

	return ret;
}

/**
 * Print a linked list.
 **/
void ll_print(ll_t *ll)
{
	ll_node_t *curr = ll->head;
	printf("LIST [");
	while (curr) {
		if (curr->key == INT_MAX)
			printf(" -> MAX");
		else
			printf(" -> %d", curr->key);
		curr = curr->next;
	}
	printf(" ]\n");
}
