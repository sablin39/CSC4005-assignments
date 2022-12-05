#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <pthread.h>

#ifdef GUI
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "./headers/physics.h"
#include "./headers/checkpoint.h"

int n_thd; // number of threads

int n_body;
int n_iteration;
double timecost;
int *slices;
int *pile;

double *x;
double *y;
double *vx;
double *vy;
double *m;

pthread_barrier_t bar;

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



void update_position(double *x, double *y, double *vx, double *vy, int n) {
    //TODO: update position 
    for (int i=pile[n];i<pile[n]+slices[n];i++){
        x[i]+=vx[i]*dt;
        y[i]+=vy[i]*dt;
    }
}

void update_velocity(double *m, double *x, double *y, double *vx, double *vy, int n) {
    //TODO: calculate force and acceleration, update velocity
    
    for (int i=pile[n];i<pile[n]+slices[n];i++){
        for (int j=0;j<n_body;j++){
            if (i==j)
                continue;
            double dist=pow((x[i]-x[j]),2.0)+pow((y[i]-y[j]),2.0);
            if (dist<=10){//collision
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
            vx[i]+=ax*dt;
            vy[i]+=ay*dt;
            if (isnan(vx[i]))
                printf("(%lf,%lf)",ax,ay);
        }
        if (x[i]<=0.1)
            vx[i]=-vx[i];
        if (y[i]<=0.1)
            vy[i]=-vy[i];
        if (x[i]>=0.99*bound_x)
            vx[i]=-vx[i];
        if (y[i]>=0.99*bound_y)
            vy[i]=-vy[i];
        
    }
}


typedef struct {
    //TODO: specify your arguments for threads
    //int a;
    //int b;
    //TODO END
} Args;


void* worker(void* args) {
    long thd=(long)args;
    update_velocity(m,x,y,vx,vy,thd);
    pthread_barrier_wait(&bar);
    update_position(x,y,vx,vy,thd);
    pthread_exit(NULL);
    //TODO: procedure in each threads

    // Args* my_arg = (Args*) args;
    // int a = my_arg->a;
    // int b = my_arg->b;

    // TODO END
}


double master(){
    m = new double[n_body];
    x = new double[n_body];
    y = new double[n_body];
    vx = new double[n_body];
    vy = new double[n_body];
    double totaltime=0;
    generate_data(m, x, y, vx, vy, n_body);
    slices=new int[n_thd];
    pile=new int[n_thd];
    for (int i=0;i<n_thd;++i){
        if (i<n_body%n_thd){
            slices[i]=n_body/n_thd+1;
        }else
            slices[i]=n_body/n_thd;
    }
    pile[0]=0;
    for (int i=0;i<n_thd-1;++i){
        pile[i+1]=pile[i]+slices[i];
    }
    pthread_barrier_init(&bar,NULL,n_thd);



    Logger l = Logger("pthread", n_body, bound_x, bound_y);

    for (int i = 0; i < n_iteration; i++){
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        //TODO: assign jobs
        pthread_t thds[n_thd];
        for (long i=0;i<n_thd;i++)
            pthread_create(&thds[i],NULL,worker,(void *)i);
        for (int i=0;i<n_thd;i++)
            pthread_join(thds[i],NULL);
        //TODO End

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


int main(int argc, char *argv[]) {
    n_body = atoi(argv[1]);
    n_iteration = atoi(argv[2]);
    n_thd = atoi(argv[3]);

    #ifdef GUI
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(800, 800);
	glutInitWindowPosition(0, 0);
	glutCreateWindow("Pthread");
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glMatrixMode(GL_PROJECTION);
	gluOrtho2D(0, bound_x, 0, bound_y);
    #endif
    timecost=master();
    printf("Student ID: 121090327\n"); // replace it with your student id
    printf("Name: Zijuan Lin\n"); // replace it with your name
    printf("Assignment 2: N Body Simulation Pthread Implementation\n");
    printf("totaltime=%lf",timecost);
	return 0;
}

