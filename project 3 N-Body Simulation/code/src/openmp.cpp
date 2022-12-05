#include <omp.h>
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


int n_body;
int n_iteration;

int n_omp_threads;


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



void update_position(double *x, double *y, double *vx, double *vy, int i) {
    //TODO: update position
    //int slice=n_body/n_omp_threads;
    x[i]+=vx[i]*dt;
    y[i]+=vy[i]*dt;


}

void update_velocity(double *m, double *x, double *y, double *vx, double *vy, int i) {
    //TODO: calculate force and acceleration, update velocity
    double tmpax=0,tmpay=0;
    for (int j=0;j<n_body;j++){
            if (i==j)
                continue;
            double dist=pow((x[i]-x[j]),2)+pow((y[i]-y[j]),2);
            if (dist<=radius2){//collision
                continue;
                vx[i]=(vx[i]*(m[i]-m[j])+2*m[j]*vx[j])/(m[i]+m[j]);
                vy[i]=(vy[i]*(m[i]-m[j])+2*m[j]*vy[j])/(m[i]+m[j]);
                vx[j]=(vx[j]*(m[j]-m[i])+2*m[i]*vx[i])/(m[i]+m[j]);
                vy[j]=(vy[j]*(m[j]-m[i])+2*m[i]*vy[i])/(m[i]+m[j]);
                
            }
            dist+=error;
            double ax=((gravity_const*m[j])/(dist))*((x[j]-x[i])/pow(dist,0.5));
            double ay=((gravity_const*m[j])/(dist))*((y[j]-y[i])/pow(dist,0.5));
            tmpax+=ax;
            tmpay+=ay;
            if (isnan(vx[i]))
                printf("%lf,%lf",ax,ay);
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


double master() {
    double* m = new double[n_body];
    double* x = new double[n_body];
    double* y = new double[n_body];
    double* vx = new double[n_body];
    double* vy = new double[n_body];
    double totaltime=0.0;
    generate_data(m, x, y, vx, vy, n_body);

    Logger l = Logger("openmp", n_body, bound_x, bound_y);

    for (int i = 0; i < n_iteration; i++){
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        
        //TODO: choose better threads configuration
        omp_set_num_threads(n_omp_threads);
        #pragma omp parallel for
        for (int i = 0; i < n_body; i++) {
            update_velocity(m, x, y, vx, vy, i);
        }

        //omp_set_num_threads(8);
        omp_set_num_threads(8);
        #pragma omp parallel for
        for (int i = 0; i < n_body; i++) {
            update_position(x, y, vx, vy, i);
        }

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
    n_omp_threads = atoi(argv[3]);
    int windowsize=500;
    double totaltime;
    #ifdef GUI
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(windowsize, windowsize);
    glutCreateWindow("N Body Simulation Sequential Implementation");
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    gluOrtho2D(0, bound_x, 0, bound_y);
    #endif
    totaltime=master();

    printf("Student ID: 121090327\n"); // replace it with your student id
    printf("Name: Zijuan Lin\n"); // replace it with your name
    printf("Assignment 2: N Body Simulation OpenMP Implementation\n");
    printf("total=%lfs",totaltime);
    return 0;

}


