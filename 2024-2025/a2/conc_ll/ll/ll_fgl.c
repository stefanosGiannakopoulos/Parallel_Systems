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
		LOCK_NODE(ll->head); \
		curr = ll->head; \
		LOCK_NODE(curr->next); \
		next = curr->next; \
		 \
		while (next->key < key) { \
			UNLOCK_NODE(curr); \
			curr = next; \
			LOCK_NODE(curr->next); \
			next = curr->next; \
		} \
	} while (0)

int ll_contains(ll_t *ll, int key)
{
	int ret = 0;
	ll_node_t *curr, *next;

	TRAVERSE_LIST();

	ret = (key == next->key);
	UNLOCK_NODE(curr);
	UNLOCK_NODE(next);
	return ret;
}

int ll_add(ll_t *ll, int key)
{
	int ret = 0;
	ll_node_t *curr, *next;
	ll_node_t *new_node;

	TRAVERSE_LIST();

	if (key != next->key) {
		ret = 1;
		new_node = ll_node_new(key);
		new_node->next = next;
		curr->next = new_node;
	}

	UNLOCK_NODE(curr);
	UNLOCK_NODE(next);
	return ret;
}

int ll_remove(ll_t *ll, int key)
{
	int ret = 0;
	ll_node_t *curr, *next;

	TRAVERSE_LIST();

	if (key == next->key) {
		ret = 1;
		curr->next = next->next;
	}

	UNLOCK_NODE(curr);
	UNLOCK_NODE(next);
	if (ret)
		free(next);
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
