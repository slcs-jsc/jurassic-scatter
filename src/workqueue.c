#include <stdlib.h> /* size_t */
#include <stdio.h> /* printf */

/* #define WORK_QUEUE_DEBUG */

#include "workqueue.h"


int init_queue(queue_t *q, int size) {
    size_t mem;
#ifdef  WORK_QUEUE_DEBUG
    printf("# %s(%p, %d);\n", __func__, (void*)q, size);
#endif
    q->capacity = 0;
    if (size > 0) {
      q->capacity = size;
      mem = (size_t)q->capacity * sizeof(queue_item_t); /* memory requirement in Byte */
      q->items = malloc(mem);
      printf("# %s(%d) requires %.3f MiByte\n", __func__, size, (double)mem/1048576.);
    } else if (NULL != q->items) { 
      printf("# %s(%d) releases the queue.\n", __func__, size);
      free(q->items);
      q->items = NULL; 
      return 0;
    }
    q->begin = 0;
    q->end = 0;
    return (NULL == q->items); /* raise error if malloc failed */
}

int push_queue(queue_t *q, void* inp, void* out, int ir) {
    int index;
    index = q->end++;
    if (index >= q->capacity) return -1-q->capacity; /* raise error */
    q->items[index].input  = inp;
    q->items[index].result = out;
    q->items[index].ir = ir;
#ifdef  WORK_QUEUE_DEBUG
    printf("# %d = %s(%p, %p) ir=%d;\n", index, __func__, inp, out, ir);
#endif
    return index;
}

int get_queue_item(queue_t *q, void** inp, void **out, int *ir, int index) {
    *inp = q->items[index].input;
    *out = q->items[index].result;
    *ir  = q->items[index].ir;
#ifdef  WORK_QUEUE_DEBUG
    printf("# %d = %s(%p, %p) ir=%d;\n", index, __func__, *inp, *out, *ir);
#endif
    return index;
}

int pop_queue(queue_t *q, void** inp, void **out, int *ir) {
    int index;
    index = q->begin++;
    get_queue_item(q, inp, out, ir, index);
#ifdef  WORK_QUEUE_DEBUG
    printf("# %d = %s(%p, %p) ir=%d;\n", index, __func__, *inp, *out, *ir);
#endif
    return index;
}

