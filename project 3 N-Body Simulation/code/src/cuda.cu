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
#include "./headers/checkpoint.h"


int block_size = 1024;


int n_body;
int n_iteration;


__global__ void update_position(double *x, double *y, double *vx, double *vy, int n) {
    //TODO: update position 
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    double dt = 0.0001f;
    
    if (i < n) {
        x[i]+=vx[i]*dt;
        y[i]+=vy[i]*dt;
    }
}

__global__ void update_velocity(double *m, double *x, double *y, double *vx, double *vy, int n) {
    //TODO: calculate force and acceleration, update velocity
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    // if (i < n) {  
    // }
    int bound_x = 4000;
    int bound_y = 4000;
    double error = 1e-9f;
    double dt = 0.0001f;
    double gravity_const = 100000.0f;
    double radius2 = 2.0f;

    double tmpax=0,tmpay=0;
    if (i<n){
        for (int j=0;j<n;j++){
        if (i==j)
            continue;
        double dist=pow((x[i]-x[j]),2)+pow((y[i]-y[j]),2);
        if (dist<=radius2){//collision
            continue;
            vx[i]=(vx[i]*(m[i]-m[j])+2*m[j]*vx[j])/(m[i]+m[j]);
            vy[i]=(vy[i]*(m[i]-m[j])+2*m[j]*vy[j])/(m[i]+m[j]);
            vx[j]=(vx[j]*(m[j]-m[i])+2*m[i]*vx[i])/(m[i]+m[j]);
            vy[j]=(vy[j]*(m[j]-m[i])+2*m[i]*vy[i])/(m[i]+m[j]);
            //continue;
        }
        dist+=error;
        double ax=((gravity_const*m[j])/(dist))*((x[j]-x[i])/pow(dist,0.5));
        double ay=((gravity_const*m[j])/(dist))*((y[j]-y[i])/pow(dist,0.5));
        tmpax+=ax;
        tmpay+=ay;
        if (isnan(vx[i]))
            printf("(%d,%lf)",j,dist);
        }
        if (x[i]<=0.1)
            vx[i]=-vx[i];
        if (y[i]<=0.1)
            vy[i]=-vy[i];
        if (x[i]>=0.99*bound_x)
            vx[i]=-vx[i];
        if (y[i]>=0.99*bound_y)
            vy[i]=-vy[i];
        vx[i]+=tmpax*dt;
        vy[i]+=tmpay*dt;
    }

}


void generate_data(double *m, double *x,double *y,double *vx,double *vy, int n) {
    // TODO: Generate proper initial position and mass for better visualization
    srand((unsigned)time(NULL));
    for (int i = 0; i < n; i++) {
        m[i] = rand() % max_mass + 1.0f;
        x[i] = 2000.0f + rand() % (bound_x / 4);
        y[i] = 2000.0f + rand() % (bound_y / 4);
        vx[i] = 0.0f;
        vy[i] = 0.0f;
    }
}


double master() {
    double* m = new double[n_body];
    double* x = new double[n_body];
    double* y = new double[n_body];
    double* vx = new double[n_body];
    double* vy = new double[n_body];
    double totaltime=0;
    generate_data(m, x, y, vx, vy, n_body);

    Logger l = Logger("cuda", n_body, bound_x, bound_y);

    double *device_m;
    double *device_x;
    double *device_y;
    double *device_vx;
    double *device_vy;

    cudaMalloc(&device_m, n_body*sizeof(double));
    cudaMalloc(&device_x, n_body*sizeof(double));
    cudaMalloc(&device_y, n_body*sizeof(double));
    cudaMalloc(&device_vx, n_body*sizeof(double));
    cudaMalloc(&device_vy, n_body*sizeof(double));

    cudaMemcpy(device_m, m, n_body*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(device_x, x, n_body*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(device_y, y, n_body*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(device_vx, vx, n_body*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(device_vy, vy, n_body*sizeof(double), cudaMemcpyHostToDevice);

    int n_block = n_body / block_size + 1;

    for (int i = 0; i < n_iteration; i++){
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

        update_velocity<<<n_block, block_size>>>(device_m, device_x, device_y, device_vx, device_vy, n_body);
        update_position<<<n_block, block_size>>>(device_x, device_y, device_vx, device_vy, n_body);

        cudaMemcpy(x, device_x, n_body*sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(y, device_y, n_body*sizeof(double), cudaMemcpyDeviceToHost);

        l.save_frame(x, y);

        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_span = t2 - t1;
        totaltime+=time_span.count();
        printf("Iteration %d, elapsed time: %.3f\n", i, time_span.count());

        #ifdef GUI
        glClear(GL_COLOR_BUFFER_BIT);
        glColor3f(1.0f, 0.0f, 0.0f);
        glPointSize(2.0f);
        glBegin(GL_POINTS);
        double xi;
        double yi;
        for (int i = 0; i < n_body; i++){
            xi = x[i];
            yi = y[i];
            glVertex2f(xi, yi);
        }
        glEnd();
        glFlush();
        glutSwapBuffers();
        #else

        #endif

    }

    cudaFree(device_m);
    cudaFree(device_x);
    cudaFree(device_y);
    cudaFree(device_vx);
    cudaFree(device_vy);

    cudaFree(device_m);
    cudaFree(device_x);
    cudaFree(device_y);
    cudaFree(device_vx);
    cudaFree(device_vy);

    delete m;
    delete x;
    delete y;
    delete vx;
    delete vy;
    return totaltime;
}


int main(int argc, char *argv[]){
    
    n_body = atoi(argv[1]);
    n_iteration = atoi(argv[2]);

    #ifdef GUI
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(500, 500);
    glutCreateWindow("N Body Simulation CUDA Implementation");
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    gluOrtho2D(0, bound_x, 0, bound_y);
    #endif

    double time=master();

    printf("Student ID: 119010001\n"); // replace it with your student id
    printf("Name: Your Name\n"); // replace it with your name
    printf("Assignment 2: N Body Simulation CUDA Implementation\n");
    printf("totaltime=%lf",time);
    return 0;

}


