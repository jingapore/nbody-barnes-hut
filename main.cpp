/* input/output */
#include <iostream>

/* math */
#include <cmath>

/* srand */
#include <stdlib.h>

/* STL types */
#include <vector>
#include <array>
#include <algorithm>

/* MPI library */
#include <mpi.h>

/* Custom MPI structs */
#include "misc/mpi_types.h"

/* Body class */
#include "misc/body.h"

/* Reading and writing */
#include "misc/readwrite.h"

/* Tree building */
#include "tree/orb.h"
#include "tree/tree.h"
#include "tree/build_tree.h"

/* Input parsing */
#include "misc/inputparser.h"

int main(int argc, char *argv[])
{

    int size, rank, tmax, N, nbodies, default_sampling_interval;
    std::vector<Body> bodies;
    std::vector<Body> bodies_out_of_range;
    std::vector<pair<double, int>> splits;
    vector<pair<array<double, DIM_SIZE>, array<double, DIM_SIZE>>> bounds, other_bounds;
    vector<pair<int, bool>> partners;
    vector<double> comp_time;
    double dt, min[DIM_SIZE], max[DIM_SIZE];
    double start_time, stop_time;
    InputParser ip;
    bool overwrite;

    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* Number of processes must be a power of 2 */
    // this should be fixed for future implementation so that it is more flexible
    if(fmod(log2(size), 1) != 0){
        if(rank == 0){
            std::cerr << "Error: Number of processes must be a power of 2 for now" << std::endl;
        }
        MPI_Finalize();
        return -1;
    }

    if (!ip.parse(argc, argv))
    {
        if (rank == 0)
        {
            std::cout << "Error: Invalid command line input" << std::endl;
            print_usage(argc, argv);
        }
        return -1;
    }

    // /* Initialize custom MPI structures */
    init_mpi_types();

    if (ip.read_bodies())
    {
        /* Read bodies from file */
        auto p = read_bodies(ip.in_file().c_str(), MPI_COMM_WORLD, nbodies, DIM_SIZE);
        bodies = p.first;
        N = p.second;
    }
    else
    {
        std::cout << "Error: Need to provide input file; we do not generate bodies" << std::endl;
        return -1;
    }


    overwrite = true;
    tmax = ip.n_steps(); // number of time steps
    default_sampling_interval = tmax; // for every x timesteps, we output positions
    dt = ip.time_step(); // time step

    if (ip.clock_run())
    {
        MPI_Barrier(MPI_COMM_WORLD);
        start_time = MPI_Wtime();
    }

    for (int t = 0; t < tmax; t++)
    {
        /* Reset variables */
        bounds.clear();
        other_bounds.clear();
        partners.clear();

        /* Domain composition and transfer of bodies */
        // global_minmax(bodies, min, max); // this can be softcoded in future, but for now we take lab's inputs

        for (int c=0; c<DIM_SIZE; c++) {
            min[c] = 0.0;
            max[c] = 4.0;
        }

        orb(bodies, bounds, other_bounds, partners, min, max, rank, size);


        /* Build the local tree */
        Tree tree(min, max, ip.bh_approx_constant(), ip.grav_constant(), 0.03);
        build_tree(bodies, bounds, other_bounds, partners, tree, rank);


        /* Compute forces */
        std::vector<array<double, DIM_SIZE>> forces;
        for (Body &b : bodies)
        {
            /* time the computation */
            double start_time = MPI_Wtime();
            array<double, DIM_SIZE> f = tree.compute_force(&b);
            /* update the workload for the body */
            b.work = MPI_Wtime() - start_time;
            forces.push_back(f);
        }

        /* Update positions */
        bool removal_required = false;
        for (int i = 0; i < bodies.size(); i++)
        {
            Body &b = bodies[i];
            for (int c = 0; c < DIM_SIZE; c++)
            {
                const double a = forces[i][c] / b.m;
                b.pos[c] = b.pos[c] + (b.vel[c] * dt) + (0.5 * a * dt * dt); // original code didn't dampen by 0.5, and just used 1.0
                b.vel[c] = b.vel[c] + (a * dt);
                if (b.pos[c]<0 || b.pos[c]>4) {
                    b.m = -1;
                    bodies_out_of_range.push_back(b);
                    removal_required = true;
                }
            }
        }

        if (removal_required) {
            bodies.erase(std::remove_if(bodies.begin(), bodies.end(), [](const Body& b) {return b.m==-1;}), bodies.end());
        }

        /* Output */
        /* Print time step to stdout */
        if (rank == 0 and ip.verbose())
        {
            std::cout << "\rTime step: " << t + 1 << "/" << tmax;
            if (t == tmax - 1)
            {
                std::cout << std::endl;
            }
        }

        if ((t+1) % default_sampling_interval == 0)
        {

            /* Stop the time */
            if (ip.clock_run())
            {
                MPI_Barrier(MPI_COMM_WORLD);
                stop_time = MPI_Wtime();
            }

            /* Write tree */
            if (rank == 0 and ip.write_tree())
            {
                write_tree(ip.out_tree_file().c_str(), tree, true, overwrite);
            }

            /* Write tree size*/
            if (rank == 0 and ip.write_tree_size())
            {
                write_to_file(ip.out_tree_size_file().c_str(), tree.size(), overwrite);
            }

            /* Write positions */
            if (ip.write_positions())
            {
                bodies.insert(bodies.end(), bodies_out_of_range.begin(), bodies_out_of_range.end());
                write_bodies(ip.out_file().c_str(), bodies, MPI_COMM_WORLD, false);
            }

            /* Write running time */
            if (ip.clock_run())
            {
                if (rank == 0)
                {
                    std::cout << stop_time - start_time << std::endl;
                }
            }

            if (overwrite)
            {
                overwrite = false;
            }

            if (ip.clock_run())
            {
                // start the time
                MPI_Barrier(MPI_COMM_WORLD);
                start_time = MPI_Wtime();
            }
        }
    }

    // /* Finalize */
    free_mpi_types();
    MPI_Finalize();

    if(ip.write_summary()){
        if(rank == 0){
            write_summary(ip, N, size);
        }
    }
}
