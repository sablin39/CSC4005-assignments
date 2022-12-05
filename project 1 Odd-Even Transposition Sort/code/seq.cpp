#include <cstdlib>
#include <fstream>
#include <iostream>
#include <chrono>
#include "vector"
#include "math.h"

int odd_even_sort(int arr[],int *res,int size){
    bool is_sorted= false;
    while (!is_sorted){
        is_sorted= true;
        std::cout<<"enter loop\n";
        for (int i=0;i<size-1;i+=2){

            if (arr[i]>arr[i+1]){
                std::swap(arr[i],arr[i+1]);
                is_sorted= false;
            }
        }
        for (int i=1;i<size-2;i+=2){

            if (arr[i]>arr[i+1]){
                std::swap(arr[i],arr[i+1]);
                is_sorted= false;
            }
        }
    }
    *res=*arr;
    for (int i=0;i<size;i++)
        std::cout<<res<<"\n";
    return 0;
}

int main (int argc, char **argv){

    int num_elements; // number of elements to be sorted
    num_elements = atoi(argv[1]); // convert command line argument to num_elements
    int elements[num_elements]; // store elements
    int sorted_elements[num_elements]; // store sorted elements
    std::string dir=argv[2];
    std::ifstream input(dir);
    int element;
    int i = 0;
    while (input >> element) {
        elements[i] = element;
        i++;
    }
    std::cout << "actual number of elements:" << i << std::endl;

    std::chrono::high_resolution_clock::time_point t1;
    std::chrono::high_resolution_clock::time_point t2;
    std::chrono::duration<double> time_span;
    t1 = std::chrono::high_resolution_clock::now(); // record time

    /* TODO BEGIN
        Implement sequential odd even transposition sort
        Code in this block is not a necessary.
        Replace it with your own code.
    */

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

    /* TODO END */

    t2 = std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    //std::cout << "Student ID: " << "121090327" << std::endl; // replace it with your student id
    //std::cout << "Name: " << "Zijuan Lin" << std::endl; // replace it with your name
    //std::cout << "Assignment 1" << std::endl;
    //std::cout << "Run Time: " << time_span.count() << " seconds" << std::endl;
    std::cout  << time_span.count()<< std::endl;
    std::cout << "Input Size: " << num_elements << std::endl;
    std::cout << "Process Number: " << 1 << std::endl;

    //std::ofstream output(argv[2]+std::string(".seq.out"), std::ios_base::out);
    std::ofstream output(argv[2]+std::string(".seq.out"), std::ios_base::out);
    for (int i = 0; i < num_elements; i++) {
        output << sorted_elements[i] << std::endl;
    }

    return 0;
}


