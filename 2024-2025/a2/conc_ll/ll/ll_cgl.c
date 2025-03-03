#include <stdio.h>
#include <stdlib.h> /* rand() */
#include <limits.h>
#include <pthread.h> /* for pthread_spinlock_t */

#include "../lib/alloc.h"
#include "ll.h"

typedef struct ll_node {
	int key;
	struct ll_node *next;
} ll_node_t;

struct linked_list {
	ll_node_t *head;
	pthread_spinlock_t lock;
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

	return ret;
}

/**
 * Free a linked list node.
 **/
static void ll_node_free(ll_node_t *ll_node)
{
	XFREE(ll_node);
}

ll_t *ll_new()
{
	ll_t *ret;

	XMALLOC(ret, 1);
	ret->head = ll_node_new(-1);
	ret->head->next = ll_node_new(INT_MAX);
	ret->head->next->next = NULL;
	pthread_spin_init(&ret->lock, PTHREAD_PROCESS_SHARED);

	return ret;
}

void ll_free(ll_t *ll)
{
	ll_node_t *next, *curr = ll->head;
	while (curr) {
		next = curr->next;
		ll_node_free(curr);
		curr = next;
	}
	pthread_spin_destroy(&ll->lock);
	XFREE(ll);
}

int ll_contains(ll_t *ll, int key)
{
	ll_node_t *curr = ll->head;
	int ret = 0;
	
	pthread_spin_lock(&ll->lock);
	while (curr->key < key)
		curr = curr->next;

	ret = (key == curr->key);
	pthread_spin_unlock(&ll->lock);
	return ret;
}

int ll_add(ll_t *ll, int key)
{
	int ret = 0;
	ll_node_t *curr, *next;
	ll_node_t *new_node;

	pthread_spin_lock(&ll->lock);
	curr = ll->head;
	next = curr->next;

	while (next->key < key) {
		curr = next;
		next = curr->next;
	}

	if (key != next->key) {
		ret = 1;
		new_node = ll_node_new(key);
		new_node->next = next;
		curr->next = new_node;
	}
	pthread_spin_unlock(&ll->lock);

	return ret;
}

int ll_remove(ll_t *ll, int key)
{
	int ret = 0;
	ll_node_t *curr, *next;

	pthread_spin_lock(&ll->lock);
	curr = ll->head;
	next = curr->next;

	while (next->key < key) {
		curr = next;
		next = curr->next;
	}

	if (key == next->key) {
		ret = 1;
		curr->next = next->next;
		ll_node_free(next);
	}
	pthread_spin_unlock(&ll->lock);

	return ret;
}

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
