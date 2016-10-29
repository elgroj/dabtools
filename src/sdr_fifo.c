#include <stdlib.h>
#include <stdio.h>
#include "sdr_fifo.h"


/* http://en.wikipedia.org/wiki/Circular_buffer */

void error(char * msg, int ctx, int arg) {
  fprintf (stderr, "Error in sdr_fifo: %s (%d, %d)\n", msg, ctx, arg);
  exit(2);
}


void cbInit(CircularBuffer *cb, uint32_t size) {
    cb->size  = size;
    cb->start = 0;
    cb->count = 0;
    cb->elems = (uint8_t *)calloc(cb->size, sizeof(uint8_t));
}
void cbFree(CircularBuffer *cb) {
    free(cb->elems);
}

int cbIsFull(CircularBuffer *cb) {
    return cb->count == cb->size;
}

int cbIsEmpty(CircularBuffer *cb) {
    return cb->count == 0;
}

void cbWrite(CircularBuffer *cb, uint8_t *elem)
{
  if (cbIsFull(cb)) {
    error("buffer overflow on write", cb->count, 1);
  }
  uint32_t end = (cb->start + cb->count) % cb->size;
  cb->elems[end] = *elem;
  cb->count++;
}

void cbRead(CircularBuffer *cb, uint8_t *elem)
{
  if (cbIsEmpty(cb)) {
    error("buffer underflow on read", 0, 1);
  }
  *elem = cb->elems[cb->start];
  cb->start = (cb->start + 1) % cb->size;
  cb->count--;
}

// just gets some old data back into the buffer
// if you do this before you have ever added something, you will get init data
void cbUnread(CircularBuffer *cb)
{
  if (cbIsFull(cb)) {
    error("buffer overflow on unread", cb->count, 1);
  }
  cb->start = (cb->start - 1) % cb->size;
  cb->count++;
}


int32_t sdr_read_fifo(CircularBuffer * fifo, uint32_t bytes, int32_t shift, uint8_t * buffer)
{
  int32_t i=0;
  uint32_t j=0;
  uint8_t dump;

  if (shift > 0) {
    if (shift > (int)fifo->count) {
      error("fifo read underflow on shift", fifo->count, shift);
    }
    for (i=0; i<shift; i++)
      cbRead(fifo, &dump);
  }
  else {
    // shift < 0
    if ((int)fifo->count - shift > (int)fifo->size) {
      error("fifo unread overflow on shift", fifo->count, shift);
    }
    for (i=0; i<-shift; i++)
      cbUnread(fifo);
  }

  if (bytes > fifo->count) {
    error("fifo read underflow on actual read", fifo->count, bytes);
  }
  for (j=0; j<bytes; j++)
    cbRead(fifo, &buffer[j]);

  return 1;
}
