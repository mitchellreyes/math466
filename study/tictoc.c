#include <sys/time.h>
#include <sys/resource.h>

static double tic_time;
void tic() {
    struct timeval tp;
    gettimeofday(&tp,0);
    tic_time=tp.tv_sec+tp.tv_usec*1.0e-6;
}
double toc() {
    struct timeval tp;
    gettimeofday(&tp,0);
    return (tp.tv_sec+tp.tv_usec*1.0e-6)-tic_time;
}
