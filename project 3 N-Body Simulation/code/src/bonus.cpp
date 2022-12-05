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
#include "./headers/checkpoint.h"


int n_body;
int n_iteration;
int n_omp_threads;
int *slices;
int *pile;
double totaltime;
int my_rank;
int world_size;


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
    omp_set_num_threads(n_omp_threads);
    #pragma omp parallel for
    for (int i=n;i<n+slices[my_rank];i++){
        x[i]+=vx[i]*dt;
        y[i]+=vy[i]*dt;
    }
}


void update_velocity(double *m, double *x, double *y, double *vx, double *vy, int n) {
    //TODO: calculate force and acceleration, update velocity
    omp_set_num_threads(n_omp_threads);
    #pragma omp parallel for
    
    for (int i=n;i<n+slices[my_rank];i++){
        for (int j=0;j<n;j++){
            if (i==j)
                continue;
            double dist=pow((x[i]-x[j]),2)+pow((y[i]-y[j]),2);
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


void slave(){
    // TODO: MPI routine
    double* total_m = new double[n_body];
    double* total_x = new double[n_body];
    double* total_y = new double[n_body];
    double* total_vx = new double[n_body];
    double* total_vy = new double[n_body];
    MPI_Bcast(total_m,n_body,MPI_DOUBLE,0,MPI_COMM_WORLD);

    double* local_m = new double[slices[my_rank]];
    double* local_x = new double[slices[my_rank]];
    double* local_y = new double[slices[my_rank]];
    double* local_vx = new double[slices[my_rank]];
    double* local_vy = new double[slices[my_rank]];
    printf("slice %d.",slices[my_rank]);
    printf("pile %d", pile[my_rank]);


    for (int i = 0; i < n_iteration; i++){
        
        MPI_Bcast(total_x,n_body,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(total_y,n_body,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(total_vx,n_body,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(total_vy,n_body,MPI_DOUBLE,0,MPI_COMM_WORLD);
        update_velocity(total_m,total_x,total_y,total_vx,total_vy,pile[my_rank]);
        update_position(total_x,total_y,total_vx,total_vy,pile[my_rank]);
        for (int i=0;i<slices[my_rank];i++){

            local_x[i]=total_x[i+pile[my_rank]];
            local_y[i]=total_y[i+pile[my_rank]];
            local_vx[i]=total_vx[i+pile[my_rank]];
            local_vy[i]=total_vy[i+pile[my_rank]];
        }
        MPI_Gatherv(local_x,slices[my_rank],MPI_DOUBLE,total_x,slices,pile,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Gatherv(local_y,slices[my_rank],MPI_DOUBLE,total_y,slices,pile,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Gatherv(local_vx,slices[my_rank],MPI_DOUBLE,total_vx,slices,pile,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Gatherv(local_vy,slices[my_rank],MPI_DOUBLE,total_vy,slices,pile,MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
    delete local_x;
    delete local_y;
    delete local_vx;
    delete local_vy;
    delete total_m;
    delete total_x;
    delete total_y;
    delete total_vx;
    delete total_vy;
    // TODO End
    
}

double master() {
    double* total_m = new double[n_body];
    double* total_x = new double[n_body];
    double* total_y = new double[n_body];
    double* total_vx = new double[n_body];
    double* total_vy = new double[n_body];
    double totaltime=0;
    generate_data(total_m, total_x, total_y, total_vx, total_vy, n_body);

    Logger l = Logger("mpi_openmp", n_body, bound_x, bound_y);
    MPI_Bcast(total_m,n_body,MPI_DOUBLE,0,MPI_COMM_WORLD);
    
    double* local_m = new double[slices[my_rank]];
    double* local_x = new double[slices[my_rank]];
    double* local_y = new double[slices[my_rank]];
    double* local_vx = new double[slices[my_rank]];
    double* local_vy = new double[slices[my_rank]];



    for (int i = 0; i < n_iteration; i++){
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        MPI_Bcast(total_x,n_body,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(total_y,n_body,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(total_vx,n_body,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(total_vy,n_body,MPI_DOUBLE,0,MPI_COMM_WORLD);
        update_velocity(total_m,total_x,total_y,total_vx,total_vy,pile[my_rank]);
        update_position(total_x,total_y,total_vx,total_vy,pile[my_rank]);
        for (int i=0;i<slices[my_rank];i++){
            local_x[i]=total_x[i+pile[my_rank]];
            local_y[i]=total_y[i+pile[my_rank]];
            local_vx[i]=total_vx[i+pile[my_rank]];
            local_vy[i]=total_vy[i+pile[my_rank]];
        }
        MPI_Gatherv(local_x,slices[my_rank],MPI_DOUBLE,total_x,slices,pile,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Gatherv(local_y,slices[my_rank],MPI_DOUBLE,total_y,slices,pile,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Gatherv(local_vx,slices[my_rank],MPI_DOUBLE,total_vx,slices,pile,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Gatherv(local_vy,slices[my_rank],MPI_DOUBLE,total_vy,slices,pile,MPI_DOUBLE,0,MPI_COMM_WORLD);
        // TODO: MPI routine
        

        // TODO End
        l.save_frame(total_x, total_y);

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
            xi = total_x[i];
            yi = total_y[i];
            glVertex2f(xi, yi);
        }
        glEnd();
        glFlush();
        glutSwapBuffers();
        #else

        #endif
    }

    delete local_x;
    delete local_y;
    delete local_vx;
    delete local_vy;
    delete total_m;
    delete total_x;
    delete total_y;
    delete total_vx;
    delete total_vy;
    return totaltime;
}




int main(int argc, char *argv[]) {
    n_body = atoi(argv[1]);
    n_iteration = atoi(argv[2]);
    n_omp_threads=n_omp_threads = atoi(argv[3]);
	MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    slices=new int[world_size];
    pile=new int[world_size];
    for (int i=0;i<world_size;++i){
        if (i<n_body%world_size){
            slices[i]=n_body/world_size+1;
        }else
            slices[i]=n_body/world_size;
    }
    pile[0]=0;
    for (int i=0;i<world_size-1;++i){
        pile[i+1]=pile[i]+slices[i];
    }
    

	if (my_rank == 0) {
		#ifdef GUI
		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
		glutInitWindowSize(500, 500); 
		glutInitWindowPosition(0, 0);
		glutCreateWindow("N Body Simulation MPI Implementation");
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glMatrixMode(GL_PROJECTION);
		gluOrtho2D(0, bound_x, 0, bound_y);
		#endif
        totaltime=master();
	} else {
        printf("enter slave");
        slave();
        
    }

	if (my_rank == 0){
		printf("Student ID: 1210909327\n"); // replace it with your student id
		printf("Name: Zijuan Lin\n"); // replace it with your name
		printf("Assignment 2: N Body Simulation MPI Implementation\n");
        printf("totaltime=%lf",totaltime);
	}

	MPI_Finalize();

	return 0;
}

