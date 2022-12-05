#include "asg2.h"
#include <stdio.h>
#include <mpi.h>


int rank;
int world_size;
void slave();

void master() {
	//TODO: procedure run in master process
	slave();
	//TODO END
}


void slave() {
	//TODO: procedure run in slave process
	int slice;
	if (rank==0){
		
		(total_size%world_size==0)?(slice=total_size/world_size):(slice=(total_size/world_size)+1);
	}
	MPI_Bcast(&slice,1,MPI_INT,0,MPI_COMM_WORLD);
	float *res=new float[slice];
	Point *p=data;
	p+=rank*slice;

	for (int i=0;i<slice;i++){
		compute(p);
		res[i]=p->color;
		p++;
	}
	
	//printf("%d;",rank);	
	float *result=new float[slice*world_size];

	MPI_Gather(res,slice,MPI_FLOAT,result,slice,MPI_FLOAT,0,MPI_COMM_WORLD);
		
	if (rank==0){
		for (int j=0;j<total_size;j++){
			
			data[j].color=result[j];
		}
	}
	//TODO END
}


int main(int argc, char *argv[]) {
	if ( argc == 4 ) {
		X_RESN = atoi(argv[1]);

		Y_RESN = atoi(argv[2]);
		max_iteration = atoi(argv[3]);
	} else {
		X_RESN = 1000;
		Y_RESN = 1000;
		max_iteration = 100;
	}

	if (rank == 0) {
		#ifdef GUI
		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
		glutInitWindowSize(500, 500); 
		glutInitWindowPosition(0, 0);
		glutCreateWindow("MPI");
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glMatrixMode(GL_PROJECTION);
		gluOrtho2D(0, X_RESN, 0, Y_RESN);
		glutDisplayFunc(plot);
		#endif
	}

	/* computation part begin */
	MPI_Init(&argc, &argv);
    	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	initData();
	if (rank == 0){
		
		t1 = std::chrono::high_resolution_clock::now();
	}
	/*
	if (rank == 0) {
		// you may change this part
		//master();
	} else {
		// you may change this part
		slave();
	}
	*/
	slave();
	if (rank == 0){
		t2 = std::chrono::high_resolution_clock::now();  
		time_span = t2 - t1;
	}

	if (rank == 0){
		printf("Student ID: 121090327\n"); // replace it with your student id
		printf("Name: Zijuan Lin\n"); // replace it with your name
		printf("Assignment 2 MPI\n");
		printf("Run Time: %f seconds\n", time_span.count());
		printf("Problem Size: %d * %d, %d\n", X_RESN, Y_RESN, max_iteration);
		printf("Process Number: %d\n", world_size);
	}

	MPI_Finalize();
	/* computation part end */

	if (rank == 0){
		#ifdef GUI
		glutMainLoop();
		#endif
	}

	return 0;
}
