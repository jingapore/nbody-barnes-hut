#include "build_tree.h"

#include<math.h>

#include <iostream>
#include <algorithm>
#include <mpi.h>
#include "../misc/mpi_types.h"

using std::cout; using std::endl;
using std::copy; using std::fill;
using std::begin; using std::end;


void build_tree(const vector<Body > & bodies, const bound_vec & bounds, const bound_vec & other_bounds,
                const vector<pair<int, bool> > & partners, Tree & tree, int rank){

    MPI_Status status;
    int n_recv_cells, partner;
    vector<Cell *> cells_to_send;
    vector<MPICell> recv_cells;
    vector<MPICell> send_cells;
    

    /* insert cells corresponding to the ORB tree created on this process */
    // progressively, each cell from bounds (a vec) is a child of the preceding one because of how ORB was done
    // what happens is that root gets 1 subcell for the 1st bound
    // then when inserting into root the next subcell for the next bound, a child is created since that next bound 
    // will be within the boundaries of the preceding cell (tested in `bounds_in_cell`)
    for(int i = 0; i < bounds.size(); i++){
        tree.insert_emptycell(&bounds.at(i).first[0], &bounds.at(i).second[0]);
    }

    /* insert bodies residing on this process */
    for(const Body & b : bodies){
        tree.insert_body(&b);
    }

    
    /* receive and send cells to other processes */
    /* loop over each split from the ORB */
    for(int i = 0; i < bounds.size(); i++){
        /* boundaries of the ORB split */
        const auto & bound = bounds.at(i);
        const auto & other_bound = other_bounds.at(i);

        /* gets cells to send to other domain*/
        cells_to_send.clear();
        tree.cells_to_send(&other_bound.first[0], &other_bound.second[0], i, cells_to_send);

        send_cells.clear();
        recv_cells.clear();
        /* construct temporary cell objects for transmission */
        // the beauty of this is how a processor will get all of the necessary info
        // because the for-loop makes each channel (or bisection) happen in sequence instead of in parallel
        // while the processors can talk to their respective partner for a channel in parallel
        // take the example of:
        // | 0 | 2 |
        // ---------
        // | 1 | 3 |
        // 0 and 2 will exchange info first.
        // but what if 2 wants to know what is going on  in 1? 0 won't have that info.
        // however, while 0 and 2 are exchanging info, let's remember that 1 and 3 are exchanging info, 
        // and all relevant info that 2 needs from 1 would have been sent to 3 
        // (since 3 is closer to 1 than 2 is).
        // so in the next channel iteration, 3 will be able to send info from 1 that 2 needs, 
        // since 3's local tree would have been updated with 1's info.
        for(Cell * cell : cells_to_send){
            MPICell mpi_cell;
            mpi_cell.m = cell->m;
            mpi_cell.parent_idx = cell->parent_idx;
            copy(begin(cell->min_bounds), end(cell->min_bounds), begin(mpi_cell.min_bounds));         
            copy(begin(cell->max_bounds), end(cell->max_bounds), begin(mpi_cell.max_bounds));         
            copy(begin(cell->rm), end(cell->rm), begin(mpi_cell.rm));         
            send_cells.push_back(mpi_cell);
        }
        
        /* get communication partner (same as during ORB) */
        partner = partners.at(i).first;
        
        /* send and receive cells from other domain */
        if(partners.at(i).second){ //`.second` is bool that is set by above
            MPI_Probe(partner, 0, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, mpi_cell_type, &n_recv_cells);
            recv_cells.resize(n_recv_cells);
            MPI_Recv(&recv_cells[0], n_recv_cells, mpi_cell_type, partner, 0, MPI_COMM_WORLD, &status);
            MPI_Send(&send_cells[0], send_cells.size(), mpi_cell_type, partner, 0, MPI_COMM_WORLD);
        }
        else{
            MPI_Send(&send_cells[0], send_cells.size(), mpi_cell_type, partner, 0, MPI_COMM_WORLD);
            MPI_Probe(partner, 0, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, mpi_cell_type, &n_recv_cells);
            recv_cells.resize(n_recv_cells);
            MPI_Recv(&recv_cells[0], n_recv_cells, mpi_cell_type, partner, 0, MPI_COMM_WORLD, &status);
        }
        

        /* construct subrees of received cells */
        vector<Cell *> root_cells = construct_received_trees(recv_cells);

        /* insert received cells into tree */
        for(Cell * cell : root_cells){
            tree.insert_cell(cell);
        } 

        /* prune the tree so as to remove subtrees not needed for evaluating the force
           on each body residing on this process */
        // Must prune tree after inserting since we always send at least one cell
        tree.prune_tree(&bound.first[0], &bound.second[0]);
    }
}


vector<Cell *> construct_received_trees(const vector<MPICell> & recv_cells){

    vector<Cell *> root_cells;
    vector<Cell *> all_cells;

    for (const MPICell & mpicell : recv_cells){
        Cell * cell = new Cell;
        cell->m = mpicell.m;
        copy(begin(mpicell.min_bounds), end(mpicell.min_bounds), begin(cell->min_bounds));         
        copy(begin(mpicell.max_bounds), end(mpicell.max_bounds), begin(cell->max_bounds));         
        copy(begin(mpicell.rm), end(mpicell.rm), begin(cell->rm));         
        fill(cell->subcells, &cell->subcells[8], nullptr);

        if(mpicell.parent_idx == -1){
            root_cells.push_back(cell); //unclear why we need vec, since i assume there should only be 1 cell.
        }
        else{
            Cell * parent = all_cells.at(mpicell.parent_idx);
            for(int i = 0; i < 8; i++){
                if(parent->subcells[i] == nullptr){
                    parent->subcells[i] = cell; 
                    break;
                }
            }
        }

        all_cells.push_back(cell);
        
    }
    return root_cells;
}
