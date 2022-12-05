#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include<pthread.h>

#ifdef GUI
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "./headers/physics.h"


int size; // problem size
int n_thd;

int iteration=1000;

typedef struct {
    int my_thd;
    int n_thd;
    int size;
    int count;
    float *odd;
    float *even;
    bool *fire_area;
    
} Args;
Args *arg;
pthread_barrier_t bar;

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
    for (int i = begin; i < end; i++){
        if (fire_area[i]) data[i] = fire_temp;
    }
}


void maintain_wall(float *data, int begin, int end) {
    // TODO: maintain the temperature of the wall
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



void* worker(void* args){
    Args* my_arg = (Args*) args;
    int thd=my_arg->my_thd;
    int n_thd=my_arg->n_thd;
    int size=my_arg->size;
    int my_size=size*size/n_thd+1;
    int start=my_size*thd;
    int end=(thd+1)*my_size>size*size ? size*size : (thd+1)*my_size;
    float *odd=my_arg->odd;
    float *even=my_arg->even;
    bool *fire_area=my_arg->fire_area;
    
    update(odd,even,start,end);
    maintain_fire(even,fire_area,start,end);
    pthread_barrier_wait(&bar);
    if (thd==0)
        maintain_wall(even,start,end);
    pthread_exit(NULL);

}

void work() {
    float* data_odd = new float[size * size];
    float* data_even = new float[size * size];
    bool* fire_area = new bool[size * size];
    // TODO: MPI routine (one possible solution, you can use another partition method)
    printf("start");
    initialize(data_odd);
    generate_fire_area(fire_area);
    pthread_barrier_init(&bar,NULL,n_thd);
    printf("Initialized");
    #ifdef GUI
    GLubyte* pixels;
    pixels = new GLubyte[resolution * resolution * 3];
    #endif

    int count = 1;
    double total_time = 0;

    
    // TODO: Send initial distribution to each slave process

    while (true) {
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        arg= new Args[n_thd];
        
        // TODO: Computation of my part
        printf("creating threads");
        pthread_t thds[n_thd];
        for (long i=0;i<n_thd;i++){
            arg[i].count=count;
            arg[i].my_thd=i;
            arg[i].n_thd=n_thd;
            arg[i].size=size;
            if (count%2==1){
            arg[i].odd=data_odd;
            arg[i].even=data_even;
            }else{
                arg[i].odd=data_even;
                arg[i].even=data_odd;
            }
            arg[i].fire_area=fire_area;
            pthread_create(&thds[i],NULL,worker,&arg[i]);
        }
            
        for (int i=0;i<n_thd;i++)
            pthread_join(thds[i],NULL);

        delete[] arg;
            
        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        double this_time = std::chrono::duration<double>(t2 - t1).count();
        total_time += this_time;
        //printf("Iteration %d, elapsed time: %.6f\n", count, this_time);
        
        count++;
        
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
        
        if (count==iteration){
            //printf("Converge after %d iterations, elapsed time: %.6f, average computation time: %.6f\n", count-1, total_time, (double) total_time / (count-1));
            printf("%d %d %.6f %.6f\n",n_thd, count-1, total_time, (double) total_time / (count-1));
            
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
    n_thd= atoi(argv[2]);
    
    #ifdef GUI
    int window_size=800;
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(window_size, window_size);
    glutCreateWindow("Heat Distribution Simulation Sequential Implementation");
    gluOrtho2D(0, resolution, 0, resolution);
    #endif
	
    work();
    

	
		//printf("Student ID: 121090327\n"); // replace it with your student id
		//printf("Name: Zijuan Lin\n"); // replace it with your name
		//printf("Assignment 4: Heat Distribution Simulation MPI Implementation\n");
	

	

	return 0;
}

