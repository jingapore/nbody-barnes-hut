#ifndef DIM_SIZE
#define DIM_SIZE 3
#endif

#ifndef _BODY_H_N_BODY_
#define _BODY_H_N_BODY_

struct Body
{
    int idx; // to ensure grader can check output
    double pos[DIM_SIZE];
    double vel[DIM_SIZE];
    double m;
    double work; // work will be used for loadbalancing across procs
};

#endif // _BODY_H_N_BODY_
