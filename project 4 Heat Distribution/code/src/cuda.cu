#include <cuda.h>
#include <cuda_runtime.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>

#ifdef GUI
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "./headers/physics.h"
/*
__device__ __managed__ float Threshold=0.00001f;
__device__ __managed__ float Fire_temp=90.0f;
__device__ __managed__ float Wall_temp=0.0f;
__device__ __managed__ float Fire_size=100;
__device__ __managed__ int Resolution=800;
*/

int block_size = 512; // cuda thread block size
int size; // problem size
int iteration=1000;

__global__ void initialize(float *data,int size) {
    // TODO: intialize the temperature distribution (in parallelized way)
    float Wall_temp=0.0f;
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= size*size) {
        return;
    }
    data[i] = Wall_temp;
}


__global__ void generate_fire_area(bool *fire_area,int size){
    float Fire_size=100;
    // TODO: generate the fire area (in parallelized way)
    int idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx >= size*size) {
        return; 
    }
    fire_area[idx]=0;
    float fire1_r2 = Fire_size * Fire_size;
    float fire2_r2 = (Fire_size / 2) * (Fire_size / 2);
    int i=idx/size;
    int j=idx%size;
    int a = i - size / 2;
    int b = j - size / 2;
    int r2 = 0.5 * a * a + 0.8 * b * b - 0.5 * a * b;
    if (r2 < fire1_r2) fire_area[idx] = 1;
    
    int c = i - 1 * size / 3;
    int d = j - 1 * size / 3;
    int r3 = c * c + d * d;
    if (r3 < fire2_r2) fire_area[i * size + j] = 1;
}


__global__ void update(float *data, float *new_data,int size) {
    // TODO: update temperature for each point  (in parallelized way)
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >=size*size) {
        return;
    }
    float total=0;
        int count=0;
        if (i-size>=0){
            total+=data[i-size];
            count++;
        }
        if (i+size<size*size){
            total+=data[i+size];
            count++;
        }
        if (i>=1){
            total+=data[i-1];
            count++;
        }
        if (i+1<size*size){
            total+=data[i+1];
            count++;
        }
        new_data[i]=total/count;
        
    
}


__global__ void maintain_wall(float *data,int size) {
    float Wall_temp=0.0f;
    // TODO: maintain the temperature of the wall (sequential is enough)
    for (int i=0;i<size-1;i++){
        data[i]=Wall_temp;
        data[(i+1)*(size)]=Wall_temp;
        data[(i+1)*(size)+1]=Wall_temp;
        data[(size)*(size-1)+i]=Wall_temp;
    }
}


__global__ void maintain_fire(float *data, bool *fire_area,int size) {
    float Fire_temp=90.0f;
    // TODO: maintain the temperature of the fire (in parallelized way)
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >=size*size) { 
        return; 
    }
    if (fire_area[i])
        data[i]=Fire_temp;
}


#ifdef GUI
__global__ void data2pixels(float *data, GLubyte* pixels,int size){
    int Resolution=800;
    float Fire_temp=90.0f;
    
    // TODO: convert rawdata (large, size^2) to pixels (small, Resolution^2) for faster rendering speed (in parallelized way)
    int idx = blockDim.x * blockIdx.x + threadIdx.x;
    
    if (idx >= size*size) {
        return; 
    }
    
    int x=idx/Resolution;
    int y=idx%Resolution;
    float factor_data_pixel = (float) size / Resolution;
    float factor_temp_color = (float) 255 / Fire_temp;
    
    
    int idx_pixel = idx * 3;
    int x_raw = x * factor_data_pixel;
    int y_raw = y * factor_data_pixel;
    int idx_raw = y_raw * size + x_raw;
    float temp = data[idx_raw];
    int color =  ((int) temp / 5 * 5) * factor_temp_color;
    pixels[idx_pixel] = color;
    pixels[idx_pixel + 1] = 255 - color;
    pixels[idx_pixel + 2] = 255 - color;
    
}


void plot(GLubyte* pixels){
    // visualize temprature distribution
    #ifdef GUI
    int Resolution=800;
    glClear(GL_COLOR_BUFFER_BIT);
    glDrawPixels(Resolution, Resolution, GL_RGB, GL_UNSIGNED_BYTE, pixels);
    glutSwapBuffers();
    #endif
}
#endif


void master() {
    int Resolution=800;
    float *data_odd;
    float *data_even;
    bool *fire_area;

    cudaMalloc(&data_odd, size * size * sizeof(float));
    cudaMalloc(&data_even, size * size * sizeof(float));
    cudaMalloc(&fire_area, size * size * sizeof(bool));
    //cudaMemcpyToSymbol(size,&size,sizeof(float));
    #ifdef GUI
    GLubyte *pixels;
    GLubyte *host_pixels;
    host_pixels = new GLubyte[Resolution * Resolution * 3];
    cudaMalloc(&pixels, resolution * resolution * 3 * sizeof(GLubyte));
    #endif

    int n_block_size = size * size / block_size + 1;
    int n_block_resolution = Resolution * Resolution / block_size + 1;

    initialize<<<n_block_size, block_size>>>(data_odd,size);
    generate_fire_area<<<n_block_size, block_size>>>(fire_area,size);
    
    int count = 1;
    double total_time = 0;

    while (true){
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

        // TODO: modify the following lines to fit your need.
        if (count % 2 == 1) {
            update<<<n_block_size, block_size>>>(data_odd, data_even,size);
            maintain_fire<<<n_block_size, block_size>>>(data_even, fire_area,size);
            maintain_wall<<<1, 1>>>(data_even,size);
        } else {
            update<<<n_block_size, block_size>>>(data_even, data_odd,size);
            maintain_fire<<<n_block_size, block_size>>>(data_odd, fire_area,size);
            maintain_wall<<<1, 1>>>(data_odd,size);
        }

        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        double this_time = std::chrono::duration<double>(t2 - t1).count();
        total_time += this_time;
        //printf("Iteration %d, elapsed time: %.6f\n", count, this_time);
        count++;
        
        #ifdef GUI
        if (count % 2 == 1) {
            data2pixels<<<n_block_resolution, block_size>>>(data_even, pixels,size);
        } else {
            data2pixels<<<n_block_resolution, block_size>>>(data_odd, pixels,size);
        }
        if (count==iteration)
            break;
        cudaMemcpy(host_pixels, pixels, resolution * resolution * 3 * sizeof(GLubyte), cudaMemcpyDeviceToHost);
        //printf("%d",host_pixels[100]);
        plot(host_pixels);
        #endif

    }

    //printf("Converge after %d iterations, elapsed time: %.6f, average computation time: %.6f\n", count-1, total_time, (double) total_time / (count-1));
    printf("%d %d %.6f %.6f\n",size, count-1, total_time, (double) total_time / (count-1));

    cudaFree(data_odd);
    cudaFree(data_even);
    cudaFree(fire_area);

    #ifdef GUI
    cudaFree(pixels);
    delete[] host_pixels;
    #endif
    
}


int main(int argc, char *argv[]){
    
    size = atoi(argv[1]);
    
    #ifdef GUI
    int Resolution=800;
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(Resolution, Resolution);
    glutCreateWindow("Heat Distribution Simulation Sequential Implementation");
    gluOrtho2D(0, Resolution, 0, Resolution);
    #endif

    master();

    //printf("Student ID: 121090327\n"); // replace it with your student id
    //printf("Name: Zijuan Lin\n"); // replace it with your name
    //printf("Assignment 4: Heat Distribution CUDA Implementation\n");

    return 0;

}


