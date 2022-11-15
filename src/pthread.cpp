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
#include "./headers/logger.h"

int n_thd; // number of threads

int n_body;
int n_iteration;

pthread_barrier_t my_barrier;
pthread_mutex_t my_lock;
pthread_mutex_t my_lock1;
int current_v;
int current_p;
std::chrono::duration<double> total_time;

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
    double x_new = x[i] + vx[i] * dt;
    double y_new = y[i] + vy[i] * dt;
    if (x_new >= bound_x || x_new <= 0) vx[i] = - vx[i] / 2;
    if (y_new >= bound_y || y_new <= 0) vy[i] = - vy[i] / 2;
    x[i] += vx[i] * dt;
    y[i] += vy[i] * dt;
}

void update_velocity(double *m, double *x, double *y, double *vx, double *vy, int n, int ith) {
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
}

typedef struct {
    //TODO: specify your arguments for threads
    double* m;
    double* x;
    double* y;
    double* vx;
    double* vy;
    int id;
    //TODO END
} Args;

void* worker(void* args) {
    //TODO: procedure in each threads
    
    Args* my_arg = (Args*) args;
    // printf("%d start v\n", my_arg->id);
    // int a = my_arg->a;
    // int b = my_arg->b;
    // update velocity
    int ith = 0;
    while (current_v < n_body) {
        pthread_mutex_lock(&my_lock);
        if (current_v >= n_body) {
            pthread_mutex_unlock(&my_lock);
            break;
        };
        ith = current_v;
        current_v = current_v + 1;
        pthread_mutex_unlock(&my_lock);

        update_velocity(my_arg->m, my_arg->x, my_arg->y, my_arg->vx, my_arg->vy, n_body, ith);
    }
    ith = 0;
    // printf("%d complete v\n", my_arg->id);

    pthread_barrier_wait(&my_barrier);

    // printf("%d start p\n", my_arg->id);

    while (current_p < n_body) {
        pthread_mutex_lock(&my_lock1);
        if (current_p >= n_body) {
            pthread_mutex_unlock(&my_lock1);
            break;
        }
        ith = current_p;
        current_p = current_p + 1;
        pthread_mutex_unlock(&my_lock1);

        update_position(my_arg->x, my_arg->y, my_arg->vx, my_arg->vy, ith);
    }
    // printf("%d complete p\n", my_arg->id);
    // TODO END
}

void master(){
    double* m = new double[n_body];
    double* x = new double[n_body];
    double* y = new double[n_body];
    double* vx = new double[n_body];
    double* vy = new double[n_body];

    generate_data(m, x, y, vx, vy, n_body);

    Logger l = Logger("pthread", n_body, bound_x, bound_y);

    for (int i = 0; i < n_iteration; i++){
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        //TODO: assign jobs
        pthread_t thds[n_thd]; // thread pool
        Args args[n_thd]; // arguments for all threads
        for (int thd = 0; thd < n_thd; thd++) {
            args[thd].m = m;
            args[thd].x = x;
            args[thd].y = y;
            args[thd].vx = vx;
            args[thd].vy = vy;
            args[thd].id = thd;
        }
        current_p = 0;
        current_v = 0;
        for (int thd = 0; thd < n_thd; thd++) pthread_create(&thds[thd], NULL, worker, &args[thd]);
        // printf("creat_complete\n");
        for (int thd = 0; thd < n_thd; thd++) pthread_join(thds[thd], NULL);
        
        //TODO End

        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_span = t2 - t1;
        total_time += time_span;

        printf("Iteration %d, elapsed time: %.3f\n", i, time_span);

        l.save_frame(x, y);

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


}

int main(int argc, char *argv[]) {
    n_body = atoi(argv[1]);
    n_iteration = atoi(argv[2]);
    n_thd = atoi(argv[3]);
    pthread_mutex_init(&my_lock, NULL);
    pthread_mutex_init(&my_lock1, NULL);
    pthread_barrier_init(&my_barrier, NULL, n_thd);

    #ifdef GUI
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutCreateWindow("Pthread");
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glMatrixMode(GL_PROJECTION);
	gluOrtho2D(0, bound_x, 0, bound_y);
    #endif
    master();

    printf("Student ID: 119010369\n");
    printf("Name: Bodong Yan\n");
    printf("Assignment 3: N Body Simulation pthread Implementation\n");
    printf("total computation time: %.3f\n", total_time);

	return 0;
}

