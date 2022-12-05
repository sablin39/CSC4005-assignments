#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <mpi.h>
#include <omp.h>
#ifdef GUI
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "./headers/physics.h"


int size; // problem size


int my_rank;
int world_size;
int n_omp_threads;
int iteration=1000;

void initialize(float *data) {
    // intialize the temperature distribution
    int len = size * size;
    for (int i = 0; i < len; i++) {
        data[i] = wall_temp;
    }
}


void generate_fire_area(bool *fire_area){
    // generate the fire area
    int len = size * size;
    for (int i = 0; i < len; i++) {
        fire_area[i] = 0;
    }

    float fire1_r2 = fire_size * fire_size;
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            int a = i - size / 2;
            int b = j - size / 2;
            int r2 = 0.5 * a * a + 0.8 * b * b - 0.5 * a * b;
            if (r2 < fire1_r2) fire_area[i * size + j] = 1;
        }
    }

    float fire2_r2 = (fire_size / 2) * (fire_size / 2);
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            int a = i - 1 * size / 3;
            int b = j - 1 * size / 3;
            int r2 = a * a + b * b;
            if (r2 < fire2_r2) fire_area[i * size + j] = 1;
        }
    }
}


void update(float *data, float *new_data, int begin, int end) {
    // TODO: update the temperature of each point, and store the result in `new_data` to avoid data racin
    omp_set_num_threads(n_omp_threads);
    #pragma omp parallel for
    for (int i = begin; i < end; i++){
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
}


void maintain_fire(float *data, bool* fire_area, int begin, int end) {
    // TODO: maintain the temperature of fire
    omp_set_num_threads(n_omp_threads);
    #pragma omp parallel for
    for (int i = begin; i < end; i++){
        if (fire_area[i]) data[i] = fire_temp;
    }
}


void maintain_wall(float *data, int begin, int end) {
    // TODO: maintain the temperature of the wall
    omp_set_num_threads(16);
    #pragma omp parallel for
    for (int i=0;i<size-1;i++){
        data[i]=wall_temp;
        data[(i+1)*(size)]=wall_temp;
        data[(i+1)*(size)+1]=wall_temp;
        data[(size)*(size-1)+i]=wall_temp;
    }
}


#ifdef GUI
void data2pixels(float *data, GLubyte* pixels){
    // convert rawdata (large, size^2) to pixels (small, resolution^2) for faster rendering speed
    float factor_data_pixel = (float) size / resolution;
    float factor_temp_color = (float) 255 / fire_temp;
    for (int x = 0; x < resolution; x++){
        for (int y = 0; y < resolution; y++){
            int idx = x * resolution + y;
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
    }
}

void plot(GLubyte* pixels){
    // visualize temprature distribution
    #ifdef GUI
    glClear(GL_COLOR_BUFFER_BIT);
    glDrawPixels(resolution, resolution, GL_RGB, GL_UNSIGNED_BYTE, pixels);
    glutSwapBuffers();
    #endif
}
#endif



void work() {
    // TODO: MPI routine (one possible solution, you can use another partition method)
    float* data_odd = new float[size * size];
    float* data_even = new float[size * size];
    bool* fire_area = new bool[size * size];
    
    int my_size=size*size/world_size+1;
    int start=my_size*my_rank;
    int end=(my_rank+1)*my_size>size*size ? size*size : (my_rank+1)*my_size;

    initialize(data_odd);
    generate_fire_area(fire_area);

    #ifdef GUI
    GLubyte* pixels;
    pixels = new GLubyte[resolution * resolution * 3];
    #endif

    int count = 1;
    double total_time = 0;

    // TODO: Send initial distribution to each slave process

    while (true) {
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

        // TODO: Computation of my part
        float *odd, *even;
        if (count%2==1){
            odd=data_odd;
            even=data_even;
        }else{
            odd=data_even;
            even=data_odd;
        }
        update(odd,even,start,end);
        maintain_fire(even,fire_area,start,end);
        maintain_wall(even,start,end);
        if (my_rank==0){
            MPI_Gather(MPI_IN_PLACE,my_size,MPI_FLOAT,even,my_size,MPI_FLOAT,0,MPI_COMM_WORLD);
        }else{
            MPI_Gather(even+start,my_size,MPI_FLOAT,even,my_size,MPI_FLOAT,0,MPI_COMM_WORLD);
        }
        
        if (my_rank!=0){
            //MPI_Request send_req;
            MPI_Send(even+start,size,MPI_FLOAT,my_rank-1,0,MPI_COMM_WORLD);
        }
        if (my_rank!=world_size-1){
            MPI_Send(even+end-size,size,MPI_FLOAT,my_rank+1,0,MPI_COMM_WORLD);
        } 
        if (my_rank!=0){
            //MPI_Request send_req;
            MPI_Recv(even+start-size,size,MPI_FLOAT,my_rank-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        if (my_rank!=world_size-1){
            MPI_Recv(even+end,size,MPI_FLOAT,my_rank+1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        } 
        

        // TODO: Send border row to neighbours
        if (my_rank==0){
            
            std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
            double this_time = std::chrono::duration<double>(t2 - t1).count();
            total_time += this_time;
            //printf("Iteration %d, elapsed time: %.6f\n", count, this_time);
        }
        count++;
        if (my_rank==0){
        #ifdef GUI
        if (count % 2 == 1) {
            // TODO: Gather pixels of slave processes
            data2pixels(data_even, pixels);
        } else {
            // TODO: Gather pixels of slave processes
            data2pixels(data_odd, pixels);
        }
        plot(pixels);
        #endif
        }
        if (count==iteration){
            MPI_Barrier(MPI_COMM_WORLD);
            if (my_rank==0){
                //printf("Converge after %d iterations, elapsed time: %.6f, average computation time: %.6f\n", count-1, total_time, (double) total_time / (count-1));
                printf("%d %d %.6f %.6f\n",world_size, count-1, total_time, (double) total_time / (count-1));
            }
            break;
        }
        
    }

    delete[] data_odd;
    delete[] data_even;
    delete[] fire_area;

    #ifdef GUI
    delete[] pixels;
    #endif
}




int main(int argc, char *argv[]) {
    size = atoi(argv[1]);
    n_omp_threads = atoi(argv[2]);
	MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);


	if (my_rank == 0) {
        #ifdef GUI
        int window_size=800;
        glutInit(&argc, argv);
        glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
        glutInitWindowPosition(0, 0);
        glutInitWindowSize(window_size, window_size);
        glutCreateWindow("Heat Distribution Simulation Sequential Implementation");
        gluOrtho2D(0, resolution, 0, resolution);
        #endif
	} 
    work();
    

	if (my_rank == 0){
		//printf("Student ID: 121090327\n"); // replace it with your student id
		//printf("Name: Zijuan Lin\n"); // replace it with your name
		//printf("Assignment 4: Heat Distribution Simulation MPI Implementation\n");
	}

	MPI_Finalize();

	return 0;
}

