#include <stdlib.h>
#include <sys/time.h>
struct timeval tStart;

void StartTimer() {
  gettimeofday(&tStart, NULL);
}

double GetElapsedTime()
{
  struct timeval tStop, tElapsed;
  gettimeofday(&tStop, NULL);
  timersub(&tStop, &tStart, &tElapsed);
  gettimeofday(&tStart, NULL); // reset original timer
  return tElapsed.tv_sec + tElapsed.tv_usec / 1000.0 / 1000.0;
}
