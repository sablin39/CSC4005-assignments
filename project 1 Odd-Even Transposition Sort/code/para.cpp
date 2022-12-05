#include <mpi.h>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <chrono>
#include "math.h"

int odd_even_sort(int my_array[],int num_my_element,int rank,int world_size,MPI_Comm comm){

    int tmp_rx;
    int tmp_tx;
    int buf_tx;
    int buf_rx;
    int local_indx;
    int worst=num_my_element*world_size;
    for (int i=0;i<worst;i++){

        if (i%2==0){
            for (local_indx=num_my_element-1;local_indx>0;local_indx-=2){
                if (my_array[local_indx] < my_array[local_indx - 1]) {
                    std::swap(my_array[local_indx - 1], my_array[local_indx]);
                }
            }
        }else{
            for (local_indx=num_my_element-2;local_indx>0;local_indx-=2){
                if (my_array[local_indx] < my_array[local_indx - 1]) {
                    std::swap(my_array[local_indx - 1], my_array[local_indx]);
                }
            }
            if(rank!=0){
                tmp_tx=my_array[0];
                MPI_Send(&tmp_tx,1,MPI_INT,rank-1,0,comm);
                MPI_Recv(&tmp_rx,1,MPI_INT,rank-1,0,comm,MPI_STATUS_IGNORE);
                if (tmp_rx>my_array[0])
                    my_array[0]=tmp_rx;
            }
            if (rank!=world_size-1){
                buf_tx=my_array[num_my_element-1];
                MPI_Recv(&tmp_rx,1,MPI_INT,rank+1,0,comm,MPI_STATUS_IGNORE);
                MPI_Send(&buf_tx,1,MPI_INT,rank+1,0,comm);
                if (tmp_rx<my_array[num_my_element-1])
                    my_array[num_my_element-1]=tmp_rx;
            }
        }
    }
    return 0;
}

int main (int argc, char **argv){

    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int num_elements; // number of elements to be sorted

    num_elements = atoi(argv[1]); // convert command line argument to num_elements

    int elements[num_elements]; // store elements
    int sorted_elements[num_elements]; // store sorted elements

    if (rank == 0) { // read inputs from file (master process)
        std::ifstream input(argv[2]);
        int element;
        int i = 0;
        while (input >> element) {
            elements[i] = element;
            i++;
        }
        std::cout << "actual number of elements:" << i << std::endl;
    }

    std::chrono::high_resolution_clock::time_point t1;
    std::chrono::high_resolution_clock::time_point t2;
    std::chrono::duration<double> time_span;
    if (rank == 0){
        t1 = std::chrono::high_resolution_clock::now(); // record time
    }

    /* TODO BEGIN
        Implement parallel odd even transposition sort
        Code in this block is not a necessary.
        Replace it with your own code.
        Useful MPI documentation: https://rookiehpc.github.io/mpi/docs

        MPI_Scatter(elements, num_my_element, MPI_INT, my_element, num_my_element, MPI_INT, 0, MPI_COMM_WORLD); // distribute elements to each process
        MPI_Gather(my_element, num_my_element, MPI_INT, sorted_elements, num_my_element, MPI_INT, 0, MPI_COMM_WORLD); // collect result from each process
    */
    int num_my_element = num_elements / world_size; // number of elements allocated to each process
    //int spare=num_elements-num_my_element*world_size;
    int my_element[num_my_element]; // store elements of each process
    MPI_Scatter(elements, num_my_element, MPI_INT, my_element, num_my_element, MPI_INT, 0, MPI_COMM_WORLD); // distribute elements to each process
    odd_even_sort(my_element,num_my_element,rank,world_size,MPI_COMM_WORLD);
    MPI_Gather(my_element, num_my_element, MPI_INT, sorted_elements, num_my_element, MPI_INT, 0, MPI_COMM_WORLD); // collect result from each process
    if (rank==0){
        //std::cout<<elements[num_elements-1]<<"|"<<elements[num_elements-2];
        int spare=num_elements-num_my_element*world_size;
        if (spare!=0){
            for (int k=num_my_element*world_size;k<num_elements;k++){
                sorted_elements[k]=elements[k];
            }

            for (int j=0;j<num_elements;j++)
                sorted_elements[j]=elements[j];
            bool is_sorted= false;
            while (!is_sorted){
                is_sorted= true;
                for (int i=0;i<num_elements-1;i+=2){

                    if (sorted_elements[i]>sorted_elements[i+1]){
                        std::swap(sorted_elements[i],sorted_elements[i+1]);
                        is_sorted= false;
                    }
                }
                for (int i=1;i<num_elements-2;i+=2){

                    if (sorted_elements[i]>sorted_elements[i+1]){
                        std::swap(sorted_elements[i],sorted_elements[i+1]);
                        is_sorted= false;
                    }
                }
            }
        }
    }

    /* TODO END */

    if (rank == 0){ // record time (only executed in master process)
        t2 = std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        std::cout << "Student ID: " << "121090327" << std::endl; // replace it with your student id
        std::cout << "Name: " << "Zijuan Lin" << std::endl; // replace it with your name
        std::cout << "Assignment 1" << std::endl;
        std::cout << "Run Time: " << time_span.count() << " seconds" << std::endl;
        //std::cout <<time_span.count()<< std::endl;
        std::cout << "Input Size: " << num_elements << std::endl;
        std::cout << "Process Number: " << world_size << std::endl;
    }

    if (rank == 0){ // write result to file (only executed in master process)
        std::ofstream output(argv[2]+std::string(".parallel.out"), std::ios_base::out);
        for (int i = 0; i < num_elements; i++) {
            output << sorted_elements[i] << std::endl;
        }
    }

    MPI_Finalize();

    return 0;
}


