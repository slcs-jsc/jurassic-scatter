#ifndef WORKQUEUE_H
#define WORKQUEUE_H

#include <stdlib.h> /* malloc */

/* Module containing all variables and functions related to a work-queue for scattering */

#define Queue_Inactive -1
#define Queue_Prepare   2 
#define Queue_Execute   8
#define Queue_Collect  32

typedef struct {
    void* input;
    void* result;
    int ir;
} queue_item_t;

typedef struct {
    queue_item_t* items;
    int capacity;
    int begin;
    int end;
} queue_t;
          
/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

int init_queue(queue_t *q, int size);
int push_queue(queue_t *q, void* inp, void* out, int ir);
int get_queue_item(queue_t *q, void** inp, void **out, int *ir, int index);
int pop_queue(queue_t *q, void** inp, void **out, int *ir);

#endif
