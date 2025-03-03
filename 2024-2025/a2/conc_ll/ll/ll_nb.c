#include <stdio.h>
#include <stdlib.h> /* rand() */
#include <limits.h>

#include "../lib/alloc.h"
#include "ll.h"

#define CAS_VAL(addr, old_val, new_val) \
	__sync_val_compare_and_swap((addr), (old_val), (new_val))

typedef struct ll_node {
	int key;
	struct ll_node *next;
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

static inline int is_marked_reference(void *ptr)
{
	long w = (long)ptr;
	return ((int)(w & 0x1L));
}

static inline void *get_unmarked_reference(void *ptr)
{
	long w = (long)ptr;
	return ((void *)(w & ~0x1L));
}

static inline void *get_marked_reference(void *ptr)
{
	long w = (long)ptr;
	return ((void *)(w | 0x1L));
}

static inline int physical_delete_right(ll_node_t *l, ll_node_t *r)
{
	ll_node_t *rnext, *cas_result;

	rnext = get_unmarked_reference(r->next);
	cas_result = CAS_VAL(&l->next, r, rnext);
	return (cas_result == r);
}

static inline ll_node_t *list_search(ll_t *ll, int key, ll_node_t **left)
{
	ll_node_t *l, *r; /* left, right */

retry:
	l = ll->head;
	r = ll->head->next;

	while (1) {
		if (l->next != r)
			goto retry;

		if (is_marked_reference(r->next)) {
			if (!physical_delete_right(l, r))
				goto retry;
		} else {
			if (r->key >= key)
				break;
			l = r;
		}
		r = get_unmarked_reference(r->next);
	}

	*left = l;
	return r;
}

int ll_contains(ll_t *ll, int key)
{
	int ret = 0;
	ll_node_t *l, *r;

	r = list_search(ll, key, &l);
	if (r->key == key && !is_marked_reference(r->next))
		ret = 1;

	return ret;
}

int ll_add(ll_t *ll, int key)
{
	ll_node_t *l, *r, *cas_result;
	ll_node_t *new_node;

	do {
		r = list_search(ll, key, &l);
		if (r->key == key)
			return 0;
		new_node = ll_node_new(key);
		new_node->next = r;
		cas_result = CAS_VAL(&l->next, r, new_node);
	} while (cas_result != r);

	return 1;
}

int ll_remove(ll_t *ll, int key)
{
	ll_node_t *l, *r, *cas_result;
	void *unmarked_ref, *marked_ref;

	do {
		r = list_search(ll, key, &l);
		if (r->key != key)
			return 0;

		unmarked_ref = get_unmarked_reference(r->next);
		marked_ref = get_marked_reference(unmarked_ref);
		cas_result = CAS_VAL(&r->next, unmarked_ref, marked_ref);
	} while (cas_result != unmarked_ref);

	physical_delete_right(l, r);
	return 1;
}
