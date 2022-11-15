#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <mpi.h>

#ifdef GUI
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "./headers/physics.h"
#include "./headers/logger.h"


int n_body;
int n_iteration;

int my_rank;
int world_size;

MPI_Group group_world;
MPI_Group slave_world;
MPI_Comm slaves;


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


void update_position(double *x, double *y, double *vx, double *vy, int i, int start_point, double* my_x, double* my_y, double* my_vx, double* my_vy) {
    //TODO: update position 
    double x_new = x[i] + vx[i] * dt;
    double y_new = y[i] + vy[i] * dt;
    if (x_new >= bound_x || x_new <= 0) {
        vx[i] = - vx[i] / 2;
        my_vx[i - start_point] = vx[i];
    }
    if (y_new >= bound_y || y_new <= 0) {
        vy[i] = - vy[i] / 2;
        my_vy[i - start_point] = vy[i];
    }
    x[i] += vx[i] * dt;
    y[i] += vy[i] * dt;
    my_x[i - start_point] = x[i];
    my_y[i - start_point] = y[i];
}


void update_velocity(double *m, double *x, double *y, double *vx, double *vy, int n, int ith, int start_point, double *my_vx, double *my_vy) {
    //TODO: calculate force and acceleration, update velocity
    double fx = 0.0;
    double fy = 0.0;

    for (int j = 0; j < n; j++) {
        if (ith == j) continue;
        double distance = sqrt(pow(x[j] - x[ith], 2) + pow(y[j] - y[ith], 2));
        if (distance < 2 * sqrt(radius2)) {
            vx[ith] = - vx[ith] / 2;
            vy[ith] = - vy[ith] / 2;
        }
        fx += ((gravity_const * m[ith] * m[j] * (x[j] - x[ith])) / pow(distance + err, 3));
        fy += ((gravity_const * m[ith] * m[j] * (y[j] - y[ith])) / pow(distance + err, 3));
    }

    double vx_new = vx[ith] + (fx * dt / m[ith]);
    double vy_new = vy[ith] + (fy * dt / m[ith]);
    vx[ith] = vx_new;
    vy[ith] = vy_new;
    my_vx[ith - start_point] = vx[ith];
    my_vy[ith - start_point] = vy[ith];
}


void slave(){
    // TODO: MPI routine
    double* local_m = new double[n_body];
    double* local_x = new double[n_body];
    double* local_y = new double[n_body];
    double* local_vx = new double[n_body];
    double* local_vy = new double[n_body];
    int local_rank;
    int local_size;
    MPI_Comm_rank(slaves, &local_rank);
    MPI_Comm_size(slaves, &local_size);
    MPI_Recv(local_m, n_body, MPI_DOUBLE, 0, my_rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // int my_num = n_body / world_size;
    // if (my_rank <= n_body % world_size) my_num += 1;
    // int start_point = (my_rank - 1) * my_num;
    // if (my_rank - 1 <= n_body % world_size) start_point = (my_rank - 1) * (int(n_body / world_size) + 1);
    // int end_point = start_point + my_num;
    int my_num = n_body / local_size;
    if (n_body % local_size != 0) my_num = (n_body + (n_body % local_size)) / local_size;
    int start_point = local_rank * my_num;
    int end_point = start_point + my_num;

    double* my_x = new double[my_num];
    double* my_y = new double[my_num];
    double* my_vx = new double[my_num];
    double* my_vy = new double[my_num];

    int term = 0;
    while (term == 0) {
        MPI_Barrier(slaves);
        MPI_Recv(&term, 1, MPI_INT, 0, my_rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(local_x, n_body, MPI_DOUBLE, 0, my_rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(local_y, n_body, MPI_DOUBLE, 0, my_rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(local_vx, n_body, MPI_DOUBLE, 0, my_rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(local_vy, n_body, MPI_DOUBLE, 0, my_rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (int i = start_point; i < end_point; i++) {
            update_velocity(local_m, local_x, local_y, local_vx, local_vy, n_body, i, start_point, my_vx, my_vy);
        }
        MPI_Barrier(slaves);
        for (int i = start_point; i < end_point; i++) {
            update_position(local_x, local_y, local_vx, local_vy, i, start_point, my_x, my_y, my_vx, my_vy);
        }

        MPI_Barrier(slaves);
        // if (local_rank == 1) {
        //     printf("%f, %f, %d \n", my_x[0], my_y[0], start_point);
        // }
        MPI_Gather(my_x, my_num, MPI_DOUBLE, local_x, my_num, MPI_DOUBLE, 0, slaves);
        MPI_Gather(my_y, my_num, MPI_DOUBLE, local_y, my_num, MPI_DOUBLE, 0, slaves);
        MPI_Gather(my_vx, my_num, MPI_DOUBLE, local_vx, my_num, MPI_DOUBLE, 0, slaves);
        MPI_Gather(my_vy, my_num, MPI_DOUBLE, local_vy, my_num, MPI_DOUBLE, 0, slaves);

        if (local_rank == 0) {
            MPI_Send(local_x, n_body, MPI_DOUBLE, 0, my_rank, MPI_COMM_WORLD);
            MPI_Send(local_y, n_body, MPI_DOUBLE, 0, my_rank, MPI_COMM_WORLD);
            MPI_Send(local_vx, n_body, MPI_DOUBLE, 0, my_rank, MPI_COMM_WORLD);
            MPI_Send(local_vy, n_body, MPI_DOUBLE, 0, my_rank, MPI_COMM_WORLD);
        }
        
    }
    
    // TODO End
}

void master() {
    double* total_m = new double[n_body];
    double* total_x = new double[n_body];
    double* total_y = new double[n_body];
    double* total_vx = new double[n_body];
    double* total_vy = new double[n_body];

    generate_data(total_m, total_x, total_y, total_vx, total_vy, n_body);

    Logger l = Logger("mpi", n_body, bound_x, bound_y);
    for (int i = 1; i < world_size; i++) MPI_Send(total_m, n_body, MPI_DOUBLE, i, i, MPI_COMM_WORLD);
    
    for (int i = 0; i < n_iteration; i++){
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

        // TODO: MPI routine
        int term = 0;
        if (i == n_iteration - 1) term = 1;
        for (int i = 1; i < world_size; i++) MPI_Send(&term, 1, MPI_INT, i, i, MPI_COMM_WORLD);
        for (int i = 1; i < world_size; i++) MPI_Send(total_x, n_body, MPI_DOUBLE, i, i, MPI_COMM_WORLD);
        for (int i = 1; i < world_size; i++) MPI_Send(total_y, n_body, MPI_DOUBLE, i, i, MPI_COMM_WORLD);
        for (int i = 1; i < world_size; i++) MPI_Send(total_vx, n_body, MPI_DOUBLE, i, i, MPI_COMM_WORLD);
        for (int i = 1; i < world_size; i++) MPI_Send(total_vy, n_body, MPI_DOUBLE, i, i, MPI_COMM_WORLD);
        
        MPI_Recv(total_x, n_body, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(total_y, n_body, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(total_vx, n_body, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(total_vy, n_body, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // TODO End

        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_span = t2 - t1;

        printf("Iteration %d, elapsed time: %.3f\n", i, time_span);

        l.save_frame(total_x, total_y);

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

    delete total_m;
    delete total_x;
    delete total_y;
    delete total_vx;
    delete total_vy;

}




int main(int argc, char *argv[]) {
    n_body = atoi(argv[1]);
    n_iteration = atoi(argv[2]);

	MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    MPI_Comm_group(MPI_COMM_WORLD, &group_world);
    int* slave_member = new int[world_size - 1];
    for (int i = 0; i < world_size - 1; i++) {
        slave_member[i] = i + 1;
    }
    MPI_Group_incl(group_world, world_size - 1, slave_member, &slave_world);
    MPI_Comm_create(MPI_COMM_WORLD, slave_world, &slaves);

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
        master();
	} else {
        slave();
    }

	if (my_rank == 0){
		printf("Student ID: 119010369\n"); // replace it with your student id
		printf("Name: Your Name\n"); // replace it with your name
		printf("Assignment 2: N Body Simulation MPI Implementation\n");
	}

	MPI_Finalize();

	return 0;
}

