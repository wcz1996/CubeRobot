////#include <stdio.h>
////#include <stdlib.h>
////#include <math.h>
////#include "GLTools.h"
////#include <GLUT/GLUT.h>
////#include <time.h>
////#include <iostream>
////#include <cmath>
////#include <fstream>
////#include <string>
////#include <sstream>
////#include <random>
////
////using namespace std;
////
////const double g = 9.8;
////const double m = 0.1;
//////const int spring_constant = 10000;
////const int ground_constant = 10000;
////const double dt = 0.001;
////const double damping = 0.999;
////const double initial_height = 0.0;
////const double pi = 3.1415926;
////const double w = 4*pi;
////const double mutation_rate = 1.0;
////const int recordNumber = 1;
////const int gen = 2;
////const int T = 5000;
////const int N = 5;
////const int length = N*N*N;
////
////float final_result[2][30000][512][3] = {0};//!!!!!!!!!!!!!!!!!!!!!!!!!!
////static int day = 0;
////int recordFlag = 1;
////
////int final_l0_type[512][512] = {0};
////int final_ball_num = 0;
////
////static int spin = 0;
////static GLdouble ex = 0.0, ey = 0.0, ez = 5.0;
////static GLfloat position[] = { 0.0, 0.0, 0.0, 1.0 };
////
////
//////////////////////////////////////////////////////////////////////////////////////////////////
////double norm(double array[]);
//////////////////////////////////////////////////////////////////////////////////////////////////
////class cube
////{
////public:
////    float initial_position[512][3] = { 0 };
////    float center[512][3] = { 0 };
////    float initial_l0[512][512] = { 0 };//!!!!!!!!!!!!!!!!!!!!!!!!!!
////    int l_0_type[512][512] = { 0 };//!!!!!!!!!!!!!!!!!!!!!!!!!!
////    double position[512][3] = { 0 };
////    double velocity[512][3] = { 0 };
////    //double acceration[1000][3] = { 0 };
////    double force[512][3] = { 0 };
////    int ball_num = 0;
////    cube()
////    {
////        // Set obstacle
////        for(int i=0; i<512; i++)
////        {
////            for(int j=0; j<3; j++)
////            {
////                initial_position[i][j] = 999;
////                center[i][j] = 999;
////            }
////        }
////        // Initial cube center
////        center[0][0] = 0.05;
////        center[0][1] = 0.05;
////        center[0][2] = 0.05;
////        // Initial position and initial_position
////        int i = 0;
////        for (double x = center[0][0]-0.05; x <= center[0][0]+0.05; x+=0.1)
////        {
////            for (double y = center[0][1]-0.05; y <= center[0][1]+0.05; y+=0.1)
////            {
////                for (double z = center[0][2]-0.05; z <= center[0][2]+0.05; z+=0.1)
////                {
////                    position[i][0] = x;
////                    position[i][1] = y;
////                    position[i][2] = z + initial_height;
////                    initial_position[i][0] = x;
////                    initial_position[i][1] = y;
////                    initial_position[i][2] = z;
////                    i++;
////                }
////            }
////        }
////    }
////};
//////////////////////////////////////////////////////////////////////////////////////////////////
////
////double norm(double array[])
////{
////    double sum = 0;
////    for (int i = 0; i < 3; i++)
////    {
////        sum += pow(array[i], 2);
////    }
////
////    sum = sqrt(sum);
////    return sum;
////}
////
////bool is_nan(double dVal)
////{
////    if (dVal==dVal)
////        return false;
////    return true;
////}
////
//////////////////////////////////////////////////////////////////////////////////////////////////
////void varify(int *p, int i, int j, int k)
////{
////    if ((i<0)||(i>=N)||(j<0)||(j>=N)||(k<0)||(k>=N)) return;
////    if ((p[i+N*j+N*N*k] > 10)||(p[i+N*j+N*N*k] == 0)) return;
////    p[i+N*j+N*N*k] += 10;    //Flag the cube has been approached.
////    varify(p,i+1,j,k);
////    varify(p,i-1,j,k);
////    varify(p,i,j+1,k);
////    varify(p,i,j-1,k);
////    varify(p,i,j,k+1);
////    varify(p,i,j,k-1);
////    return;
////}
//////////////////////////////////////////////////////////////////////////////////////////////////
////double get_fitness(int parent[], int loop)
////{
////    //    cout<<center_type[0]<<' '<<"abc"<<endl;
////    cube  c;
////
////    int center_type[512];
////    for (int i = 0; i<N*N*N; i++)
////    {
////        center_type[i] = parent[i];
////    }
////
////    varify(center_type, 0, 0, 0);
////    for (int i = 0; i<N*N*N; i++){
////        if (center_type[i] >10)
////            center_type[i]-=10;
////        else center_type[i]=0;
////    }
////
////    double links[5][3] = {0.0};
////    double fitness = 0.0;
////    // Initialize a cub
////    // Get ball num
////    c.ball_num = 0;
////    for(int i=0; i<512; i++)
////    {
////        if(c.initial_position[i][0] != 999)
////            c.ball_num += 1;
////        else
////            break;
////    }
////    // Initial l0 of first cube
////    for(int i=0; i<8; i++)
////    {
////        for(int j=0; j<8; j++)
////        {
////            double delta_l[3] = {0.0};
////            delta_l[0] = c.initial_position[i][0] - c.initial_position[j][0];
////            delta_l[1] = c.initial_position[i][1] - c.initial_position[j][1];
////            delta_l[2] = c.initial_position[i][2] - c.initial_position[j][2];
////            c.initial_l0[i][j] = norm(delta_l);
////            c.l_0_type[i][j] = center_type[0];
////        }
////    }
////
////    // Set 4 kinds of type
////    links[1][0] = 10000;
////    links[1][1] = 0;
////    links[1][2] = 0;
////    links[2][0] = 10000;
////    links[2][1] = 0.1;
////    links[2][2] = 0;
////    links[3][0] = 10000;
////    links[3][1] = 0.1;
////    links[3][2] = pi;
////    links[4][0] = 1000;
////    links[4][1] = 0;
////    links[4][2] = 0;
////
////    // Reconstruct center list
////    for (int x = 0; x < N; x++)
////    {
////        for (int y = 0; y < N; y++)
////        {
////            for (int z = 0; z < N; z++)
////            {
////                if(center_type[x+N*y+N*N*z] != 0)
////                {
////                    c.center[x+N*y+N*N*z][0] = x*0.1+0.05;
////                    c.center[x+N*y+N*N*z][1] = y*0.1+0.05;
////                    c.center[x+N*y+N*N*z][2] = z*0.1+0.05;
////
////                }
////            }
////        }
////    }
////
////    // Reconstruct Robot according to center
////    double temp_initial_pos[8][3] = {0};
////    for(int i=1; i<N*N*N; i++)
////    {
////        if(c.center[i][0] == 999)
////            continue;
////        else
////        {
////            int temp_i = 0;
////            for (double x = c.center[i][0]-0.05; x <= c.center[i][0]+0.06; x+=0.1)
////            {
////                for (double y = c.center[i][1]-0.05; y <= c.center[i][1]+0.06; y+=0.1)
////                {
////                    for (double z = c.center[i][2]-0.05; z <= c.center[i][2]+0.06; z+=0.1)
////                    {
////                        temp_initial_pos[temp_i][0] = x;
////                        temp_initial_pos[temp_i][1] = y;
////                        temp_initial_pos[temp_i][2] = z;
////                        temp_i += 1;
////                    }
////                }
////            }
////        }
////        int last_ball_num = c.ball_num;
////        int ball_list[8] = {0};
////        for(int ii=0; ii<8; ii++)
////        {
////            int repeat = 0; // 0 for new point
////            for(int jj=0; jj<last_ball_num; jj++)
////            {
////                if(c.initial_position[jj][0]-temp_initial_pos[ii][0] < 0.01 && c.initial_position[jj][0]-temp_initial_pos[ii][0] > -0.01)
////                {
////                    if(c.initial_position[jj][1]-temp_initial_pos[ii][1] < 0.01 && c.initial_position[jj][1]-temp_initial_pos[ii][1] > -0.01)
////                    {
////                        if(c.initial_position[jj][2]-temp_initial_pos[ii][2] < 0.01 && c.initial_position[jj][2]-temp_initial_pos[ii][2] > -0.01)
////                        {
////                            ball_list[ii] = jj;
////                            repeat = 1;
////                            break;
////                        }
////                    }
////                }
////            }
////            if(repeat == 0)
////            {
////                //cout << i << ' '<< temp_initial_pos[i][0] << ' ' << temp_initial_pos[i][1] << ' ' << temp_initial_pos[i][2] << endl;
////                ball_list[ii] = c.ball_num;
////                c.initial_position[c.ball_num][0] = temp_initial_pos[ii][0];
////                c.initial_position[c.ball_num][1] = temp_initial_pos[ii][1];
////                c.initial_position[c.ball_num][2] = temp_initial_pos[ii][2];
////                c.position[c.ball_num][0]=temp_initial_pos[ii][0];
////                c.position[c.ball_num][1]=temp_initial_pos[ii][1];
////                c.position[c.ball_num][2]=temp_initial_pos[ii][2] + initial_height;
////                c.ball_num += 1;
////            }
////        }
////        for(int ii=0; ii<8; ii++)
////        {
////            for(int jj=0; jj<8; jj++)
////            {
////                int iii = ball_list[ii];
////                int jjj = ball_list[jj];
////                double delta_l[3] = {0.0};
////                delta_l[0] = c.initial_position[iii][0] - c.initial_position[jjj][0];
////                delta_l[1] = c.initial_position[iii][1] - c.initial_position[jjj][1];
////                delta_l[2] = c.initial_position[iii][2] - c.initial_position[jjj][2];
////                c.initial_l0[iii][jjj] = norm(delta_l);
////                c.l_0_type[iii][jjj] = center_type[i];
////            }
////        }
////    }
////
////    for(int i=0; i<c.ball_num; i++)
////    {
////        for(int j=0; j<c.ball_num; j++)
////        {
////            final_l0_type[i][j] = c.l_0_type[i][j];
////        }
////    }
////    final_ball_num = c.ball_num;
////    //    cout<<c.ball_num<<endl;
////    //    cout << c.position[0][0] << ' ' << c.position[0][1] << ' ' << c.position[0][2] << endl;
////    //    for(int i=0; i<16; i++)
////    //        cout << i << ' '<< c.position[i][0] << ' ' << c.position[i][1] << ' ' << c.position[i][2] << endl;
////
////    // Simulator
////    int t = 0;
////    double fx = 0.0;
////    double fy = 0.0;
////    double fz = 0.0;
////    double fxy = 0.0;
////    double vxy = 0.0;
////    double delta_p[3] = { 0 };
////    double delta_initial_p[3] = { 0 };
////    double l_0 = 0.1;
////    double delta_l = 0;
////    double direction[3] = { 0 };
////    double fl[3] = { 0 };
////    double final_path[3] = {0};
////    double temp_path[3] = {0};
////
////    while (t <= loop)
////    {
////        //        cout << c.force[0][0] << ' ' << c.force[0][1] << ' ' << c.force[0][2] << endl;
////        //cout << c.position[0][0] << ' ' << c.position[0][1] << ' ' << c.position[0][2] << endl;
////        cout<<t<<endl;
////        for (int i = 0; i < c.ball_num; i++)
////        {
////            //cout<<c.force[0][2]<<endl;
////            // Add gravity
////            for (int j = 0; j < 3; j++)
////            {
////                if (j != 2)
////                    c.force[i][j] = 0.0;
////                else
////                    c.force[i][j] = -1 * m * g;
////            }
////            //Sprint force
////            for (int j = 0; j < c.ball_num; j++)
////            {
////                if ((c.initial_l0[i][j] - 0.1 < 0.01 && c.initial_l0[i][j] - 0.1 > -0.01)||
////                    (c.initial_l0[i][j] - 0.1 * sqrt(2) < 0.01 && c.initial_l0[i][j] - 0.1 * sqrt(2) > -0.01)||
////                    (c.initial_l0[i][j] - 0.1 * sqrt(3) < 0.01 && c.initial_l0[i][j] - 0.1 * sqrt(3) > -0.01))
////                {
////                    for (int k = 0; k < 3; k++)
////                    {
////                        delta_p[k] = c.position[i][k] - c.position[j][k];
////                    }
////
////                    //cout<<l_0<<endl;
////                    //l_0 *= (1 + parameter_b[i][j] * sin(w*t*dt + parameter_c[i][j]));
////                    //                    cout<<c.l_0_type[i][j]<<endl;
////                    //                    cout<<c.l_0_type[i][j]<<' '<<c.l_0_type[j][i]<<endl;
////                    delta_l = norm(delta_p) - c.initial_l0[i][j]*(1 + links[c.l_0_type[i][j]][1] * sin(w*t*dt + links[c.l_0_type[i][j]][2]));
////                    //                    delta_l = norm(delta_p) - c.initial_l0[i][j]*(1 + 0.1 * sin(w*t*dt + pi));
////
////                    for (int k = 0; k < 3; k++)
////                    {
////                        //cout<<norm(delta_p)<<endl;
////                        direction[k] = delta_p[k] / norm(delta_p);
////                        //fl[k] = -1 * parameter_a[i][j] * delta_l * direction[k];
////                        fl[k] = -1 * links[c.l_0_type[i][j]][0] * delta_l * direction[k];
////                        c.force[i][k] += fl[k];
////                    }
////                }
////            }
////
////            // Ground force
////            if (c.position[i][2] < 0)
////            {
////                // Vertical
////                fz = -1* ground_constant * c.position[i][2];
////
////                //Horizontal
////                fxy = sqrt(pow(c.force[i][0], 2) + pow(c.force[i][1], 2));
////                vxy = sqrt(pow(c.velocity[i][0], 2) + pow(c.velocity[i][1], 2));
////                //cout << fxy << endl;
////                // Static friction
////                if (vxy < 0.00001 && fxy <= fz * 1)
////                {
////                    c.force[i][0] += -1 * c.force[i][0];
////                    c.force[i][1] += -1 * c.force[i][1];
////                    c.force[i][2] += fz;
////                }
////                else if (vxy < 0.00001)
////                {
////                    c.force[i][0] += -1 * fz * 1;
////                    c.force[i][1] += -1 * fz * 1;
////                    c.force[i][2] += fz;
////                }
////                // Dynamic friction
////                else
////                {
////                    fx = -0.8 * fz * c.velocity[i][0] / vxy;
////                    fy = -0.8 * fz * c.velocity[i][1] / vxy;
////                    c.force[i][0] += fx;
////                    c.force[i][1] += fy;
////                    c.force[i][2] += fz;
////                }
////            }
////            //            if(isnan(c.force[i][0]) ||isnan(c.force[i][1]) ||isnan(c.force[i][2]) )
////            //            {
////            //
////            //            }
////
////        }
////        // Update parameters
////        for (int i = 0; i < c.ball_num; i++)
////        {
////            for (int j = 0; j < 3; j++)
////            {
////                //                c.acceration[i][j] = c.force[i][j] / m;
////                c.velocity[i][j] += c.force[i][j] / m * dt;
////                c.velocity[i][j] *= damping;
////                c.position[i][j] += c.velocity[i][j] * dt;
////                if (recordFlag > 0){
////                    final_result[recordFlag][t][i][j] = c.position[i][j];
////                }
////            }
////        }
////        if(t == 2000)
////        {
////            for (int i = 0; i < c.ball_num; i++)
////            {
////                temp_path[0] = (c.position[i][0] - c.initial_position[i][0]);
////                temp_path[1] = (c.position[i][1] - c.initial_position[i][1]);
////            }
////        }
////        t++;
////    }
////
////    for (int i = 0; i < c.ball_num; i++)
////    {
////        final_path[0] = (c.position[i][0] - c.initial_position[i][0]);
////        final_path[1] = (c.position[i][1] - c.initial_position[i][1]);
////        delta_p[0] = final_path[0] - temp_path[0];
////        delta_p[1] = final_path[1] - temp_path[1];
////        fitness += norm(delta_p)/c.ball_num;
////    }
////
////    return fitness;
////}
//////////////////////////////////////////////////////////////////////////////////////////////////
////void init(void)
////{
////    GLfloat mat_specular[] = { 1.0, 0.0, 1.0, 1.0 };
////    GLfloat mat_shininess[] = { 128.0 };
////    GLfloat sun_mat_ambient[] = {0.0f, 0.0f, 0.0f, 1.0f};
////    GLfloat sun_mat_diffuse[] = {0.0f, 0.0f, 0.0f, 1.0f};
////    GLfloat sun_mat_specular[] = {0.0f, 0.0f, 0.0f, 1.0f};
////    GLfloat sun_mat_emission[] = {0.5f, 0.0f, 0.0f, 1.0f};
////    GLfloat sun_mat_shininess = 0.0f;
////    glMaterialfv(GL_FRONT, GL_AMBIENT, sun_mat_ambient);
////    glMaterialfv(GL_FRONT, GL_DIFFUSE, sun_mat_diffuse);
////    glMaterialfv(GL_FRONT, GL_SPECULAR, sun_mat_specular);
////    glMaterialfv(GL_FRONT, GL_EMISSION, sun_mat_emission);
////    glMaterialf (GL_FRONT, GL_SHININESS, sun_mat_shininess);
////
////
////
////    //光的颜色照到物体上,即物体呈现出的颜色
////    GLfloat light_position[] = { 1.0, 1.0, 1.0, 0.0 }; //最后一个参数为0表示该光源是directional的
////    GLfloat light_ambient[] = { 0.0, 0.0, 0.0, 1.0 };
////    GLfloat light_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
////    GLfloat light_specular[] = { 1.0, 1.0, 1.0, 1.0 };
////
////    //环境光与光源无关，当关闭LIGHT0时，显现环境光颜色
////    //GLfloat Light_Model_Ambient[] = { 0.0 , 0.0 , 0.0 , 1.0 }; //全局环境光设为黑色
////    //GLfloat Light_Model_Ambient[] = { 1.0 ,1.0 , 1.0 , 1.0 };//这时物体呈现灰色，也即是环境光的颜色（当LIGHT0关闭时）
////    GLfloat Light_Model_Ambient[] = { 0.2 ,0.2 , 0.2 , 1.0 }; //默认的全局环境光
////    glLightModelfv( GL_LIGHT_MODEL_AMBIENT , Light_Model_Ambient ); //设置全局环境光的方法
////
////    glClearColor (0.0, 0.0, 0.0, 0.0);
////    glShadeModel (GL_SMOOTH);
////
//////    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
//////    glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
////    glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
////    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
////    glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
////    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
////
////    glEnable(GL_LIGHTING);
////    glEnable(GL_LIGHT0); //把全局环境光设置为黑色,然后注释掉本句可以发现"黑屏"
////    glDepthFunc(GL_LESS); //glDepthFunc()来设置深度测试函数
////    glEnable(GL_DEPTH_TEST);
////}
//////////////////////////////////////////////////////////////////////////////////////////////////
////void myDisplay(void)
////{
////
////    glEnable(GL_DEPTH_TEST);
////    // 清除屏幕及深度缓存
////    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
////    // 深度测试开启，实现遮挡关系
////    glDepthFunc(GL_LESS);
////    glHint(GL_PERSPECTIVE_CORRECTION_HINT,GL_NICEST);
////
////    //Background color
////    glClearColor(0.0f,0.0f,0.0f,0.0f);
////
////    //Sight Seting
////    glEnable(GL_DEPTH_TEST);
////    glMatrixMode(GL_PROJECTION);
////    glLoadIdentity();
////    gluPerspective(40, 1, 1, 40000000);
////    glMatrixMode(GL_MODELVIEW);
////    glLoadIdentity();
////    gluLookAt(-1, -2, 2, 0, 0, 0, 0, 0, 1);
////
////    // Ground
////    int k=0;
////    for (float j=2.0; j>0; j=j-0.1){
////        if (k==0) {
////            glBegin(GL_POLYGON);
////            glColor3f(0.1f, 0.1f, 0.1f);
////            for(int i=0; i<30; ++i) {
////                glVertex3f(j*cos(2*pi/30*i), j*sin(2*pi/30*i), 0.001*(2-j));
////            }
////            k=1;
////            glEnd();
////        }
////        else {
////            glBegin(GL_POLYGON);
////            glColor3f(0.5f, 0.5f, 0.5f);
////            for(int i=0; i<30; ++i) {
////                glVertex3f(j*cos(2*pi/30*i), j*sin(2*pi/30*i), 0.001*(2-j));
////            }
////            k=0;
////            glEnd();
////        }
////    }
////
////    //XYZ axis
////    glLineWidth(3.0f);
////    glBegin(GL_LINES);
////    glColor3f(1.0f, 0.0f, 0.0f);
////    glVertex3f(0.0f, 0.0f, 0.0f);
////    glVertex3f(0.0f, 0.0f,2.0f);
////    glVertex3f(0.0f, 0.0f, 0.0f);
////    glVertex3f(0.0f, 2.0f,0.0f);
////    glVertex3f(0.0f, 0.0f, 0.0f);
////    glVertex3f(2.0f, 0.0f, 0.0f);
////    glEnd();
////
////    glLineWidth(3.0f);
////    glBegin(GL_LINES);
////    glColor3f(1.0f, 1.0f, 0.0f);
////    for (int k = 1; k<=recordNumber; k++){
////        for ( int i = 0; i < final_ball_num; i++){
////            for (int j = i+1; j < final_ball_num; j ++){
////                if (final_l0_type[i][j] == 1)
////                {
////                    glColor3f(0.2f*k/recordNumber, 1.0f*k/recordNumber, 1.0f*k/recordNumber);
////                    glVertex3d(final_result[k][day][i][0]-0.05,final_result[k][day][i][1]-0.05,final_result[k][day][i][2]);
////                    glVertex3d(final_result[k][day][j][0]-0.05,final_result[k][day][j][1]-0.05,final_result[k][day][j][2]);
////                }
////                else if (final_l0_type[i][j] == 2)
////                {
////                    glColor3f(1.0f*k/recordNumber, 0.392f*k/recordNumber, 1.0f*k/recordNumber);
////                    glVertex3d(final_result[k][day][i][0]-0.05,final_result[k][day][i][1]-0.05,final_result[k][day][i][2]);
////                    glVertex3d(final_result[k][day][j][0]-0.05,final_result[k][day][j][1]-0.05,final_result[k][day][j][2]);
////                }
////                else if (final_l0_type[i][j] == 3)
////                {
////                    glColor3f(0.4f*k/recordNumber, 1.0f*k/recordNumber, 0.0f*k/recordNumber);
////                    glVertex3d(final_result[k][day][i][0]-0.05,final_result[k][day][i][1]-0.05,final_result[k][day][i][2]);
////                    glVertex3d(final_result[k][day][j][0]-0.05,final_result[k][day][j][1]-0.05,final_result[k][day][j][2]);
////                }
////                else if (final_l0_type[i][j] == 4)
////                {
////                    glColor3f(1.0f*k/recordNumber, 1.0f*k/recordNumber, 0.0f*k/recordNumber);
////                    glVertex3d(final_result[k][day][i][0]-0.05,final_result[k][day][i][1]-0.05,final_result[k][day][i][2]);
////                    glVertex3d(final_result[k][day][j][0]-0.05,final_result[k][day][j][1]-0.05,final_result[k][day][j][2]);
////                }
////            }
////        }
////    }
////    glEnd();
////
//////    //遮挡测试
//////    glBegin(GL_QUADS);
//////    glColor3f(1.0f, 0.5f, 0.5f);
//////    glVertex3i(1, 1,-0.1);
//////    glVertex3i(1+1, 1, -0.1);
//////    glVertex3i(1+1, 1+1, -0.1);
//////    glVertex3i(1, 1+1, -0.1);
//////    glEnd();
////
////    //glutSolidSphere( 1.0, 20.0, 20.0);
////
////    glFlush();
////    glutSwapBuffers();
////}
////
//////////////////////////////////////////////////////////////////////////////////////////////////
////void myIdle(void)
////{
////    day++;
////    while (day % 1 > 0)
////        day++;
////    if (day>10*T)
////        day = 0;
////    myDisplay();
////}
/////////////////////////////////////////////////////////////////////////////////
//////void reshape (int w, int h)
//////{
//////    glViewport (0, 0, (GLsizei) w, (GLsizei) h);
//////    glMatrixMode (GL_PROJECTION);
//////    glLoadIdentity();
//////    gluPerspective(40.0, (GLfloat) w/(GLfloat) h, 1.0, 20.0);
//////    glMatrixMode(GL_MODELVIEW);
//////    glLoadIdentity();
//////
//////    //glLightfv (GL_LIGHT0, GL_POSITION, position);
//////}
//////
//////void mouse(int button, int state, int x, int y)
//////{
//////    switch (button) {
//////        case GLUT_LEFT_BUTTON:
//////            if (state == GLUT_DOWN) {
//////                spin = (spin + 30) % 360;
//////                glutPostRedisplay();
//////            }
//////            break;
//////        default:
//////            break;
//////    }
//////}
//////void keyboard(unsigned char key, int x, int y)
//////{
//////    switch (key) {
//////        case 'w': ey += 0.2;  glutPostRedisplay(); break;
//////        case 's': ey -= 0.2;  glutPostRedisplay(); break;
//////        case 'd': ex += 0.2;  glutPostRedisplay(); break;
//////        case 'a': ex -= 0.2;  glutPostRedisplay(); break;
//////        case 27:
//////            exit(0);
//////            break;
//////    }
//////}
////
//////////////////////////////////////////////////////////////////////////////////////////////////
////int main(int argc, char *argv[])
////{
////
////    int spring[512] = {0};
////
////
////    char data[100];
////
////    ifstream infile;
////    infile.open("100.txt");
////
////    for (int i=0; i<512; i++){
////        infile >> data;
////        spring[i] = atoi(data);
////        cout<< spring[i] <<endl;
////    }
////    // 关闭打开的文件
////    infile.close();
////    double fitnessAtTheEnd = get_fitness(spring, T);
////
////    glutInit(&argc, argv);
////    glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE | GLUT_DEPTH);
////    glutInitWindowPosition(100, 100);
////    glutInitWindowSize(1080, 1080);
////    glutCreateWindow("Cube Robot");
////    //init ();
////    glutDisplayFunc(myDisplay);
////    glutIdleFunc(myIdle);
//////    glutReshapeFunc(reshape);
//////    glutMouseFunc(mouse);
//////    glutKeyboardFunc(keyboard);
////    glutMainLoop();
////    return 0;
////}
////
//// Licenced by Chengzhang Wang & Xipeng Xie
//// Used for Evolutionary Computation Coure in Columbia University
//// Created in Dec. 2018
//
//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>
//#include "GLTools.h"
//#include <GLUT/GLUT.h>
//#include <time.h>
//#include <iostream>
//#include <cmath>
//#include <fstream>
//#include <string>
//#include <sstream>
//
//using namespace std;
//
//const double g = 9.8;
//const double m = 0.1;
////const int spring_constant = 10000;
//const int ground_constant = 10000;
//const double dt = 0.001;
//const double damping = 0.999;
//const double initial_height = 0.0;
//const double pi = 3.1415926;
//const double w = 4*pi;
//const double mutation_rate = 1.0;
//const int recordNumber = 1;
//const int gen = 100;
//const int T = 30000;
//const int N = 5;
//const int length = N*N*N;
//
//double final_result[10][30001][130][3] = {0};
//static int day = 0;
//int recordFlag = 1;
//
//double final_l0_type[1000][1000] = {0};
//int final_ball_num = 0;
//
//
////////////////////////////////////////////////////////////////////////////////////////////////
//double norm(double array[]);
////////////////////////////////////////////////////////////////////////////////////////////////
//class cube
//{
//public:
//    double initial_position[1000][3] = { 0 };
//    double center[1000][3] = { 0 };
//    double initial_l0[130][130] = { 0 };
//    int l_0_type[1000][1000] = { 0 };
//    double position[1000][3] = { 0 };
//    double velocity[1000][3] = { 0 };
//    double acceration[1000][3] = { 0 };
//    double force[1000][3] = { 0 };
//    int ball_num = 0;
//    cube()
//    {
//        // Set obstacle
//        for(int i=0; i<1000; i++)
//        {
//            for(int j=0; j<3; j++)
//            {
//                initial_position[i][j] = 999;
//                center[i][j] = 999;
//            }
//        }
//        // Initial cube center
//        center[0][0] = 0.05;
//        center[0][1] = 0.05;
//        center[0][2] = 0.05;
//
//        int i = 0;
//        for (double x = center[0][0]-0.05; x <= center[0][0]+0.05; x+=0.1)
//        {
//            for (double y = center[0][1]-0.05; y <= center[0][1]+0.05; y+=0.1)
//            {
//                for (double z = center[0][2]-0.05; z <= center[0][2]+0.05; z+=0.1)
//                {
//                    position[i][0] = x;
//                    position[i][1] = y;
//                    position[i][2] = z + initial_height;
//                    initial_position[i][0] = x;
//                    initial_position[i][1] = y;
//                    initial_position[i][2] = z;
//                    i++;
//                }
//            }
//        }
//    }
//};
////////////////////////////////////////////////////////////////////////////////////////////////
//
//double norm(double array[])
//{
//    double sum = 0;
//    for (int i = 0; i < 3; i++)
//    {
//        sum += pow(array[i], 2);
//    }
//
//    sum = sqrt(sum);
//    return sum;
//}
//
//bool is_nan(double dVal)
//{
//    if (dVal==dVal)
//        return false;
//    return true;
//}
//
////////////////////////////////////////////////////////////////////////////////////////////////
//void varify(int *p, int i, int j, int k)
//{
//    if ((i<0)||(i>=N)||(j<0)||(j>=N)||(k<0)||(k>=N)) return;
//    if ((p[i+N*j+N*N*k] > 10)||(p[i+N*j+N*N*k] == 0)) return;
//    p[i+N*j+N*N*k] += 10;    //Flag the cube has been approached.
//    varify(p,i+1,j,k);
//    varify(p,i-1,j,k);
//    varify(p,i,j+1,k);
//    varify(p,i,j-1,k);
//    varify(p,i,j,k+1);
//    varify(p,i,j,k-1);
//    return;
//}
////////////////////////////////////////////////////////////////////////////////////////////////
//double get_fitness(int parent[], int loop)
//{
//    //    cout<<center_type[0]<<' '<<"abc"<<endl;
//    cube  c;
//
//    int center_type[1000];
//    for (int i = 0; i<N*N*N; i++)
//    {
//        center_type[i] = parent[i];
//    }
//
//    varify(center_type, 0, 0, 0);
//    for (int i = 0; i<N*N*N; i++){
//        if (center_type[i] >10)
//            center_type[i]-=10;
//        else center_type[i]=0;
//    }
//
//    double links[5][3] = {0.0};
//    double temp_fitness = 0.0;
//    double fitness = 0.0;
//    // Initialize a cub
//    // Get ball num
//    c.ball_num = 0;
//    for(int i=0; i<1000; i++)
//    {
//        if(c.initial_position[i][0] != 999)
//            c.ball_num += 1;
//        else
//            break;
//    }
//    // Initial l0 of first cube
//    for(int i=0; i<8; i++)
//    {
//        for(int j=0; j<8; j++)
//        {
//            double delta_l[3] = {0.0};
//            delta_l[0] = c.initial_position[i][0] - c.initial_position[j][0];
//            delta_l[1] = c.initial_position[i][1] - c.initial_position[j][1];
//            delta_l[2] = c.initial_position[i][2] - c.initial_position[j][2];
//            c.initial_l0[i][j] = norm(delta_l);
//            c.l_0_type[i][j] = center_type[0];
//        }
//    }
//
//    // Set 4 kinds of type
//    links[1][0] = 10000;
//    links[1][1] = 0;
//    links[1][2] = 0;
//    links[2][0] = 10000;
//    links[2][1] = 0.1;
//    links[2][2] = 0;
//    links[3][0] = 10000;
//    links[3][1] = 0.1;
//    links[3][2] = pi;
//    links[4][0] = 1000;
//    links[4][1] = 0;
//    links[4][2] = 0;
//
//    // Reconstruct center list
//    for (int x = 0; x < N; x++)
//    {
//        for (int y = 0; y < N; y++)
//        {
//            for (int z = 0; z < N; z++)
//            {
//                if(center_type[x+N*y+N*N*z] != 0)
//                {
//                    c.center[x+N*y+N*N*z][0] = x*0.1+0.05;
//                    c.center[x+N*y+N*N*z][1] = y*0.1+0.05;
//                    c.center[x+N*y+N*N*z][2] = z*0.1+0.05;
//
//                }
//            }
//        }
//    }
//
//    // Reconstruct Robot according to center
//    double temp_initial_pos[8][3] = {0};
//    for(int i=1; i<N*N*N; i++)
//    {
//        if(c.center[i][0] == 999)
//            continue;
//        else
//        {
//            int temp_i = 0;
//            for (double x = c.center[i][0]-0.05; x <= c.center[i][0]+0.06; x+=0.1)
//            {
//                for (double y = c.center[i][1]-0.05; y <= c.center[i][1]+0.06; y+=0.1)
//                {
//                    for (double z = c.center[i][2]-0.05; z <= c.center[i][2]+0.06; z+=0.1)
//                    {
//                        temp_initial_pos[temp_i][0] = x;
//                        temp_initial_pos[temp_i][1] = y;
//                        temp_initial_pos[temp_i][2] = z;
//                        temp_i += 1;
//                    }
//                }
//            }
//        }
//        int last_ball_num = c.ball_num;
//        int ball_list[8] = {0};
//        for(int ii=0; ii<8; ii++)
//        {
//            int repeat = 0; // 0 for new point
//            for(int jj=0; jj<last_ball_num; jj++)
//            {
//                if(c.initial_position[jj][0]-temp_initial_pos[ii][0] < 0.01 && c.initial_position[jj][0]-temp_initial_pos[ii][0] > -0.01)
//                {
//                    if(c.initial_position[jj][1]-temp_initial_pos[ii][1] < 0.01 && c.initial_position[jj][1]-temp_initial_pos[ii][1] > -0.01)
//                    {
//                        if(c.initial_position[jj][2]-temp_initial_pos[ii][2] < 0.01 && c.initial_position[jj][2]-temp_initial_pos[ii][2] > -0.01)
//                        {
//                            ball_list[ii] = jj;
//                            repeat = 1;
//                            break;
//                        }
//                    }
//                }
//            }
//            if(repeat == 0)
//            {
//                //cout << i << ' '<< temp_initial_pos[i][0] << ' ' << temp_initial_pos[i][1] << ' ' << temp_initial_pos[i][2] << endl;
//                ball_list[ii] = c.ball_num;
//                c.initial_position[c.ball_num][0] = temp_initial_pos[ii][0];
//                c.initial_position[c.ball_num][1] = temp_initial_pos[ii][1];
//                c.initial_position[c.ball_num][2] = temp_initial_pos[ii][2];
//                c.position[c.ball_num][0]=temp_initial_pos[ii][0];
//                c.position[c.ball_num][1]=temp_initial_pos[ii][1];
//                c.position[c.ball_num][2]=temp_initial_pos[ii][2] + initial_height;
//                c.ball_num += 1;
//            }
//        }
//        for(int ii=0; ii<8; ii++)
//        {
//            for(int jj=0; jj<8; jj++)
//            {
//                int iii = ball_list[ii];
//                int jjj = ball_list[jj];
//                double delta_l[3] = {0.0};
//                delta_l[0] = c.initial_position[iii][0] - c.initial_position[jjj][0];
//                delta_l[1] = c.initial_position[iii][1] - c.initial_position[jjj][1];
//                delta_l[2] = c.initial_position[iii][2] - c.initial_position[jjj][2];
//                c.initial_l0[iii][jjj] = norm(delta_l);
//                c.l_0_type[iii][jjj] = center_type[i];
//            }
//        }
//    }
//
//    for(int i=0; i<c.ball_num; i++)
//    {
//        for(int j=0; j<c.ball_num; j++)
//        {
//            final_l0_type[i][j] = c.l_0_type[i][j];
//        }
//    }
//    final_ball_num = c.ball_num;
//    //    cout<<c.ball_num<<endl;
//    //    cout << c.position[0][0] << ' ' << c.position[0][1] << ' ' << c.position[0][2] << endl;
//    //    for(int i=0; i<16; i++)
//    //        cout << i << ' '<< c.position[i][0] << ' ' << c.position[i][1] << ' ' << c.position[i][2] << endl;
//
//    // Simulator
//    int t = 0;
//    double fx = 0.0;
//    double fy = 0.0;
//    double fz = 0.0;
//    double fxy = 0.0;
//    double vxy = 0.0;
//    double delta_p[3] = { 0 };
//    double delta_initial_p[3] = { 0 };
//    double l_0 = 0.1;
//    double delta_l = 0;
//    double direction[3] = { 0 };
//    double fl[3] = { 0 };
//    double final_path[3] = {0};
//    double temp_path[3] = {0};
//
//    while (t <= loop)
//    {
//        //        cout << c.force[0][0] << ' ' << c.force[0][1] << ' ' << c.force[0][2] << endl;
//        //cout << c.position[0][0] << ' ' << c.position[0][1] << ' ' << c.position[0][2] << endl;
//        for (int i = 0; i < c.ball_num; i++)
//        {
//            //cout<<c.force[0][2]<<endl;
//            // Add gravity
//            for (int j = 0; j < 3; j++)
//            {
//                if (j != 2)
//                    c.force[i][j] = 0.0;
//                else
//                    c.force[i][j] = -1 * m * g;
//            }
//            //Sprint force
//            for (int j = 0; j < c.ball_num; j++)
//            {
//                if ((c.initial_l0[i][j] - 0.1 < 0.01 && c.initial_l0[i][j] - 0.1 > -0.01)||
//                    (c.initial_l0[i][j] - 0.1 * sqrt(2) < 0.01 && c.initial_l0[i][j] - 0.1 * sqrt(2) > -0.01)||
//                    (c.initial_l0[i][j] - 0.1 * sqrt(3) < 0.01 && c.initial_l0[i][j] - 0.1 * sqrt(3) > -0.01))
//                {
//                    for (int k = 0; k < 3; k++)
//                    {
//                        delta_p[k] = c.position[i][k] - c.position[j][k];
//                    }
//
//                    //cout<<l_0<<endl;
//                    //l_0 *= (1 + parameter_b[i][j] * sin(w*t*dt + parameter_c[i][j]));
//                    //                    cout<<c.l_0_type[i][j]<<endl;
//                    //                    cout<<c.l_0_type[i][j]<<' '<<c.l_0_type[j][i]<<endl;
//                    delta_l = norm(delta_p) - c.initial_l0[i][j]*(1 + links[c.l_0_type[i][j]][1] * sin(w*t*dt + links[c.l_0_type[i][j]][2]));
//                    //                    delta_l = norm(delta_p) - c.initial_l0[i][j]*(1 + 0.1 * sin(w*t*dt + pi));
//
//
//                    //cout << t << ' ' << delta_l << endl;
//
//                    for (int k = 0; k < 3; k++)
//                    {
//                        //cout<<norm(delta_p)<<endl;
//                        direction[k] = delta_p[k] / norm(delta_p);
//                        //fl[k] = -1 * parameter_a[i][j] * delta_l * direction[k];
//                        fl[k] = -1 * links[c.l_0_type[i][j]][0] * delta_l * direction[k];
//                        c.force[i][k] += fl[k];
//                    }
//                }
//            }
//
//            // Ground force
//            if (c.position[i][2] < 0)
//            {
//                // Vertical
//                fz = -1* ground_constant * c.position[i][2];
//
//                //Horizontal
//                fxy = sqrt(pow(c.force[i][0], 2) + pow(c.force[i][1], 2));
//                vxy = sqrt(pow(c.velocity[i][0], 2) + pow(c.velocity[i][1], 2));
//                //cout << fxy << endl;
//                // Static friction
//                if (vxy < 0.00001 && fxy <= fz * 1)
//                {
//                    c.force[i][0] += -1 * c.force[i][0];
//                    c.force[i][1] += -1 * c.force[i][1];
//                    c.force[i][2] += fz;
//                }
//                else if (vxy < 0.00001)
//                {
//                    c.force[i][0] += -1 * fz * 1;
//                    c.force[i][1] += -1 * fz * 1;
//                    c.force[i][2] += fz;
//                }
//                // Dynamic friction
//                else
//                {
//                    fx = -0.8 * fz * c.velocity[i][0] / vxy;
//                    fy = -0.8 * fz * c.velocity[i][1] / vxy;
//                    c.force[i][0] += fx;
//                    c.force[i][1] += fy;
//                    c.force[i][2] += fz;
//                }
//            }
//            //            if(isnan(c.force[i][0]) ||isnan(c.force[i][1]) ||isnan(c.force[i][2]) )
//            //            {
//            //
//            //            }
//
//        }
//        // Update parameters
//        for (int i = 0; i < c.ball_num; i++)
//        {
//            for (int j = 0; j < 3; j++)
//            {
//                c.acceration[i][j] = c.force[i][j] / m;
//                c.velocity[i][j] += c.acceration[i][j] * dt;
//                c.velocity[i][j] *= damping;
//                c.position[i][j] += c.velocity[i][j] * dt;
//                if (recordFlag > 0){
//                    final_result[recordFlag][t][i][j] = c.position[i][j];
//                }
//            }
//        }
//        if(t == 10)
//        {
//            for (int i = 0; i < c.ball_num; i++)
//            {
//                temp_path[0] += c.position[i][0];
//                temp_path[1] += c.position[i][1];
//            }
//        }
//        t++;
//    }
//
//    for (int i = 0; i < c.ball_num; i++)
//    {
//        final_path[0] += c.position[i][0];
//        final_path[1] += c.position[i][1];
//    }
//
//    delta_p[0] = final_path[0] - temp_path[0];
//    delta_p[1] = final_path[1] - temp_path[1];
//    fitness += norm(delta_p)/c.ball_num;
//    cout<< final_path[0]<< ' '<< final_path[1] <<endl;
//    cout<< temp_path[0]<< ' '<< temp_path[1] <<endl;
//    cout<< c.ball_num<<endl;
//    return fitness;
//}
//
////////////////////////////////////////////////////////////////////////////////////////////////
//void myDisplay(void)
//{
//
//    glEnable(GL_DEPTH_TEST);
//    // 清除屏幕及深度缓存
//    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//    // 深度测试开启，实现遮挡关系
//    glDepthFunc(GL_LEQUAL);
//    glHint(GL_PERSPECTIVE_CORRECTION_HINT,GL_NICEST);
//
//    //Background color
//    glClearColor(0.0f,0.0f,0.0f,0.0f);
//
//    //Anti-aliasing
//    glEnable(GL_BLEND);
//    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
//    glEnable(GL_POINT_SMOOTH);
//    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
//    glEnable(GL_LINE_SMOOTH);
//    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
//
//    //Sight Seting
//    glEnable(GL_DEPTH_TEST);
//    glMatrixMode(GL_PROJECTION);
//    glLoadIdentity();
//    gluPerspective(40, 2.5, 1, 40000000);
//    glMatrixMode(GL_MODELVIEW);
//    glLoadIdentity();
//    gluLookAt(1, 2, 1, 0, 0, 0, 0, 0, 1);
//
//    int k=0;
//    for (float j=2.0; j>2.9; j=j-0.1){
//        if (k==0) {
//            glBegin(GL_POLYGON);
//            glColor3f(0.1f, 0.1f, 0.1f);
//            for(int i=0; i<30; ++i) {
//                glVertex3f(j*cos(2*pi/30*i), j*sin(2*pi/30*i), 0.0f);
//            }
//            k=1;
//            glEnd();
//        }
//        else {
//            glBegin(GL_POLYGON);
//            glColor3f(0.5f, 0.5f, 0.5f);
//            for(int i=0; i<30; ++i) {
//                glVertex3f(j*cos(2*pi/30*i), j*sin(2*pi/30*i), 0.0f);
//            }
//            k=0;
//            glEnd();
//        }
//    }
//
//    //XYZ axis
//    glLineWidth(1.0f);
//    glBegin(GL_LINES);
//    glColor3f(0.0f, 1.0f, 1.0f);
//    glVertex3f(0.0f, 0.0f, 0.0f);
//    glVertex3f(0.0f, 0.0f,2.0f);
//    glVertex3f(0.0f, 0.0f, 0.0f);
//    glVertex3f(0.0f, 1000.0f,0.0f);
//    glVertex3f(0.0f, 0.0f, 0.0f);
//    glVertex3f(100.0f, 0.0f, 0.0f);
//    glEnd();
//
//    double xx = 0;
//    double yy = 0;
//    glLineWidth(3.0f);
//    glBegin(GL_LINES);
//    glColor3f(1.0f, 1.0f, 0.0f);
//    for (int k = 1; k<=recordNumber; k++){
//        for ( int i = 0; i < final_ball_num; i++){
//            for (int j = i+1; j < final_ball_num; j ++){
//                if (final_l0_type[i][j] == 1)
//                {
//                    glColor3f(0.2f*k/recordNumber, 1.0f*k/recordNumber, 1.0f*k/recordNumber);
//                    glVertex3d(final_result[k][day][i][0]-0.05-xx,final_result[k][day][i][1]-0.05-yy,final_result[k][day][i][2]);
//                    glVertex3d(final_result[k][day][j][0]-0.05-xx,final_result[k][day][j][1]-0.05-yy,final_result[k][day][j][2]);
//                }
//                else if (final_l0_type[i][j] == 2)
//                {
//                    glColor3f(1.0f*k/recordNumber, 0.392f*k/recordNumber, 1.0f*k/recordNumber);
//                    glVertex3d(final_result[k][day][i][0]-0.05-xx,final_result[k][day][i][1]-0.05-yy,final_result[k][day][i][2]);
//                    glVertex3d(final_result[k][day][j][0]-0.05-xx,final_result[k][day][j][1]-0.05-yy,final_result[k][day][j][2]);
//                }
//                else if (final_l0_type[i][j] == 3)
//                {
//                    glColor3f(0.4f*k/recordNumber, 1.0f*k/recordNumber, 0.0f*k/recordNumber);
//                    glVertex3d(final_result[k][day][i][0]-0.05-xx,final_result[k][day][i][1]-0.05-yy,final_result[k][day][i][2]);
//                    glVertex3d(final_result[k][day][j][0]-0.05-xx,final_result[k][day][j][1]-0.05-yy,final_result[k][day][j][2]);
//                }
//                else if (final_l0_type[i][j] == 4)
//                {
//                    glColor3f(1.0f*k/recordNumber, 1.0f*k/recordNumber, 0.0f*k/recordNumber);
//                    glVertex3d(final_result[k][day][i][0]-0.05-xx,final_result[k][day][i][1]-0.05-yy,final_result[k][day][i][2]);
//                    glVertex3d(final_result[k][day][j][0]-0.05-xx,final_result[k][day][j][1]-0.05-yy,final_result[k][day][j][2]);
//                }
//            }
//        }
//    }
//    glEnd();
//    glFlush();
//    glutSwapBuffers();
//}
//
////////////////////////////////////////////////////////////////////////////////////////////////
//void myIdle(void)
//{
//    day++;
//    if (day < 10000) day =10000;
//    while (day % 1 > 0)
//        day++;
//    if (day>= T)
//        day = 10000;
//    myDisplay();
//}
//
////////////////////////////////////////////////////////////////////////////////////////////////
//int main(int argc, char *argv[])
//{
//
//    int spring[512] = {0};
//
//
//        char data[100];
//
//        ifstream infile;
//        infile.open("5.txt");
//
//        for (int i=0; i<512; i++){
//            infile >> data;
//            spring[i] = atoi(data);
//            cout<< spring[i] <<endl;
//        }
//        // 关闭打开的文件
//        infile.close();
//    double fitnessAtTheEnd = get_fitness(spring, T);
//
//    glutInit(&argc, argv);
//    glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE | GLUT_DEPTH);
//    glutInitWindowPosition(100, 100);
//    glutInitWindowSize(1250, 500);
//    glutCreateWindow("第一个 OpenGL 程序");
//    glutDisplayFunc(&myDisplay);
//    glutIdleFunc(&myIdle);
//    glutMainLoop();
//    return 0;
//}
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "GLTools.h"
#include <GLUT/GLUT.h>
#include <time.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <random>

using namespace std;

const double g = 9.8;
const double m = 0.1;
//const int spring_constant = 10000;
const int ground_constant = 10000;
const double dt = 0.001;
const double damping = 0.999;
const double initial_height = 1.0;
const double pi = 3.1415926;
const double w = 4*pi;
const double mutation_rate = 1.0;
const int recordNumber = 1;
const int gen = 2;
const int T = 30000;
const int N = 3;
const int length = N*N*N;

float final_result[2][30000][512][3] = {0};//!!!!!!!!!!!!!!!!!!!!!!!!!!
static int day = 0;
int recordFlag = 1;

int final_l0_type[512][512] = {0};
int final_ball_num = 0;

//////////////////////////////////////////////////////////////////////////////////////////////
double norm(double array[]);
//////////////////////////////////////////////////////////////////////////////////////////////
class cube
{
public:
    float initial_position[512][3] = { 0 };
    float center[512][3] = { 0 };
    float initial_l0[512][512] = { 0 };//!!!!!!!!!!!!!!!!!!!!!!!!!!
    int l_0_type[512][512] = { 0 };//!!!!!!!!!!!!!!!!!!!!!!!!!!
    double position[512][3] = { 0 };
    double velocity[512][3] = { 0 };
    //double acceration[1000][3] = { 0 };
    double force[512][3] = { 0 };
    int ball_num = 0;
    cube()
    {
        // Set obstacle
        for(int i=0; i<512; i++)
        {
            for(int j=0; j<3; j++)
            {
                initial_position[i][j] = 999;
                center[i][j] = 999;
            }
        }
        // Initial cube center
        center[0][0] = 0.05;
        center[0][1] = 0.05;
        center[0][2] = 0.05;
        // Initial position and initial_position
        int i = 0;
        for (double x = center[0][0]-0.05; x <= center[0][0]+0.05; x+=0.1)
        {
            for (double y = center[0][1]-0.05; y <= center[0][1]+0.05; y+=0.1)
            {
                for (double z = center[0][2]-0.05; z <= center[0][2]+0.05; z+=0.1)
                {
                    position[i][0] = x;
                    position[i][1] = y;
                    position[i][2] = z + initial_height;
                    initial_position[i][0] = x;
                    initial_position[i][1] = y;
                    initial_position[i][2] = z;
                    i++;
                }
            }
        }
    }
};
//////////////////////////////////////////////////////////////////////////////////////////////

double norm(double array[])
{
    double sum = 0;
    for (int i = 0; i < 3; i++)
    {
        sum += pow(array[i], 2);
    }
    
    sum = sqrt(sum);
    return sum;
}

bool is_nan(double dVal)
{
    if (dVal==dVal)
        return false;
    return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////
void varify(int *p, int i, int j, int k)
{
    if ((i<0)||(i>=N)||(j<0)||(j>=N)||(k<0)||(k>=N)) return;
    if ((p[i+N*j+N*N*k] > 10)||(p[i+N*j+N*N*k] == 0)) return;
    p[i+N*j+N*N*k] += 10;    //Flag the cube has been approached.
    varify(p,i+1,j,k);
    varify(p,i-1,j,k);
    varify(p,i,j+1,k);
    varify(p,i,j-1,k);
    varify(p,i,j,k+1);
    varify(p,i,j,k-1);
    return;
}
//////////////////////////////////////////////////////////////////////////////////////////////
double get_fitness(int parent[], int loop)
{
    //    cout<<center_type[0]<<' '<<"abc"<<endl;
    cube  c;
    
    int center_type[512];
    for (int i = 0; i<N*N*N; i++)
    {
        center_type[i] = parent[i];
    }
    
    varify(center_type, 0, 0, 0);
    for (int i = 0; i<N*N*N; i++){
        if (center_type[i] >10)
            center_type[i]-=10;
        else center_type[i]=0;
    }
    
    double links[5][3] = {0.0};
    double fitness = 0.0;
    // Initialize a cub
    // Get ball num
    c.ball_num = 0;
    for(int i=0; i<512; i++)
    {
        if(c.initial_position[i][0] != 999)
            c.ball_num += 1;
        else
            break;
    }
    // Initial l0 of first cube
    for(int i=0; i<8; i++)
    {
        for(int j=0; j<8; j++)
        {
            double delta_l[3] = {0.0};
            delta_l[0] = c.initial_position[i][0] - c.initial_position[j][0];
            delta_l[1] = c.initial_position[i][1] - c.initial_position[j][1];
            delta_l[2] = c.initial_position[i][2] - c.initial_position[j][2];
            c.initial_l0[i][j] = norm(delta_l);
            c.l_0_type[i][j] = center_type[0];
        }
    }
    
    // Set 4 kinds of type
    links[1][0] = 10000;
    links[1][1] = 0;
    links[1][2] = 0;
    links[2][0] = 10000;
    links[2][1] = 0.1;
    links[2][2] = 0;
    links[3][0] = 10000;
    links[3][1] = 0.1;
    links[3][2] = pi;
    links[4][0] = 1000;
    links[4][1] = 0;
    links[4][2] = 0;
    
    // Reconstruct center list
    for (int x = 0; x < N; x++)
    {
        for (int y = 0; y < N; y++)
        {
            for (int z = 0; z < N; z++)
            {
                if(center_type[x+N*y+N*N*z] != 0)
                {
                    c.center[x+N*y+N*N*z][0] = x*0.1+0.05;
                    c.center[x+N*y+N*N*z][1] = y*0.1+0.05;
                    c.center[x+N*y+N*N*z][2] = z*0.1+0.05;
                    
                }
            }
        }
    }
    
    // Reconstruct Robot according to center
    double temp_initial_pos[8][3] = {0};
    for(int i=1; i<N*N*N; i++)
    {
        if(c.center[i][0] == 999)
            continue;
        else
        {
            int temp_i = 0;
            for (double x = c.center[i][0]-0.05; x <= c.center[i][0]+0.06; x+=0.1)
            {
                for (double y = c.center[i][1]-0.05; y <= c.center[i][1]+0.06; y+=0.1)
                {
                    for (double z = c.center[i][2]-0.05; z <= c.center[i][2]+0.06; z+=0.1)
                    {
                        temp_initial_pos[temp_i][0] = x;
                        temp_initial_pos[temp_i][1] = y;
                        temp_initial_pos[temp_i][2] = z;
                        temp_i += 1;
                    }
                }
            }
        }
        int last_ball_num = c.ball_num;
        int ball_list[8] = {0};
        for(int ii=0; ii<8; ii++)
        {
            int repeat = 0; // 0 for new point
            for(int jj=0; jj<last_ball_num; jj++)
            {
                if(c.initial_position[jj][0]-temp_initial_pos[ii][0] < 0.01 && c.initial_position[jj][0]-temp_initial_pos[ii][0] > -0.01)
                {
                    if(c.initial_position[jj][1]-temp_initial_pos[ii][1] < 0.01 && c.initial_position[jj][1]-temp_initial_pos[ii][1] > -0.01)
                    {
                        if(c.initial_position[jj][2]-temp_initial_pos[ii][2] < 0.01 && c.initial_position[jj][2]-temp_initial_pos[ii][2] > -0.01)
                        {
                            ball_list[ii] = jj;
                            repeat = 1;
                            break;
                        }
                    }
                }
            }
            if(repeat == 0)
            {
                //cout << i << ' '<< temp_initial_pos[i][0] << ' ' << temp_initial_pos[i][1] << ' ' << temp_initial_pos[i][2] << endl;
                ball_list[ii] = c.ball_num;
                c.initial_position[c.ball_num][0] = temp_initial_pos[ii][0];
                c.initial_position[c.ball_num][1] = temp_initial_pos[ii][1];
                c.initial_position[c.ball_num][2] = temp_initial_pos[ii][2];
                c.position[c.ball_num][0]=temp_initial_pos[ii][0];
                c.position[c.ball_num][1]=temp_initial_pos[ii][1];
                c.position[c.ball_num][2]=temp_initial_pos[ii][2] + initial_height;
                c.ball_num += 1;
            }
        }
        for(int ii=0; ii<8; ii++)
        {
            for(int jj=0; jj<8; jj++)
            {
                int iii = ball_list[ii];
                int jjj = ball_list[jj];
                double delta_l[3] = {0.0};
                delta_l[0] = c.initial_position[iii][0] - c.initial_position[jjj][0];
                delta_l[1] = c.initial_position[iii][1] - c.initial_position[jjj][1];
                delta_l[2] = c.initial_position[iii][2] - c.initial_position[jjj][2];
                c.initial_l0[iii][jjj] = norm(delta_l);
                c.l_0_type[iii][jjj] = center_type[i];
            }
        }
    }
    
    for(int i=0; i<c.ball_num; i++)
    {
        for(int j=0; j<c.ball_num; j++)
        {
            final_l0_type[i][j] = c.l_0_type[i][j];
        }
    }
    final_ball_num = c.ball_num;
    //    cout<<c.ball_num<<endl;
    //    cout << c.position[0][0] << ' ' << c.position[0][1] << ' ' << c.position[0][2] << endl;
    //    for(int i=0; i<16; i++)
    //        cout << i << ' '<< c.position[i][0] << ' ' << c.position[i][1] << ' ' << c.position[i][2] << endl;
    
    // Simulator
    int t = 0;
    double fx = 0.0;
    double fy = 0.0;
    double fz = 0.0;
    double fxy = 0.0;
    double vxy = 0.0;
    double delta_p[3] = { 0 };
    double delta_initial_p[3] = { 0 };
    double l_0 = 0.1;
    double delta_l = 0;
    double direction[3] = { 0 };
    double fl[3] = { 0 };
    double final_path[3] = {0};
    double temp_path[3] = {0};
    
    while (t <= loop)
    {
        //        cout << c.force[0][0] << ' ' << c.force[0][1] << ' ' << c.force[0][2] << endl;
        //cout << c.position[0][0] << ' ' << c.position[0][1] << ' ' << c.position[0][2] << endl;
        cout<<t<<endl;
        for (int i = 0; i < c.ball_num; i++)
        {
            //cout<<c.force[0][2]<<endl;
            // Add gravity
            for (int j = 0; j < 3; j++)
            {
                if (j != 2)
                    c.force[i][j] = 0.0;
                else
                    c.force[i][j] = -1 * m * g;
            }
            //Sprint force
            for (int j = 0; j < c.ball_num; j++)
            {
                if ((c.initial_l0[i][j] - 0.1 < 0.01 && c.initial_l0[i][j] - 0.1 > -0.01)||
                    (c.initial_l0[i][j] - 0.1 * sqrt(2) < 0.01 && c.initial_l0[i][j] - 0.1 * sqrt(2) > -0.01)||
                    (c.initial_l0[i][j] - 0.1 * sqrt(3) < 0.01 && c.initial_l0[i][j] - 0.1 * sqrt(3) > -0.01))
                {
                    for (int k = 0; k < 3; k++)
                    {
                        delta_p[k] = c.position[i][k] - c.position[j][k];
                    }
                    
                    //cout<<l_0<<endl;
                    //l_0 *= (1 + parameter_b[i][j] * sin(w*t*dt + parameter_c[i][j]));
                    //                    cout<<c.l_0_type[i][j]<<endl;
                    //                    cout<<c.l_0_type[i][j]<<' '<<c.l_0_type[j][i]<<endl;
                    delta_l = norm(delta_p) - c.initial_l0[i][j]*(1 + links[c.l_0_type[i][j]][1] * sin(w*t*dt + links[c.l_0_type[i][j]][2]));
                    //                    delta_l = norm(delta_p) - c.initial_l0[i][j]*(1 + 0.1 * sin(w*t*dt + pi));
                    
                    for (int k = 0; k < 3; k++)
                    {
                        //cout<<norm(delta_p)<<endl;
                        direction[k] = delta_p[k] / norm(delta_p);
                        //fl[k] = -1 * parameter_a[i][j] * delta_l * direction[k];
                        fl[k] = -1 * links[c.l_0_type[i][j]][0] * delta_l * direction[k];
                        c.force[i][k] += fl[k];
                    }
                }
            }
            
            // Ground force
            if (c.position[i][2] < 0)
            {
                // Vertical
                fz = -1* ground_constant * c.position[i][2];
                
                //Horizontal
                fxy = sqrt(pow(c.force[i][0], 2) + pow(c.force[i][1], 2));
                vxy = sqrt(pow(c.velocity[i][0], 2) + pow(c.velocity[i][1], 2));
                //cout << fxy << endl;
                // Static friction
                if (vxy < 0.00001 && fxy <= fz * 1)
                {
                    c.force[i][0] += -1 * c.force[i][0];
                    c.force[i][1] += -1 * c.force[i][1];
                    c.force[i][2] += fz;
                }
                else if (vxy < 0.00001)
                {
                    c.force[i][0] += -1 * fz * 1;
                    c.force[i][1] += -1 * fz * 1;
                    c.force[i][2] += fz;
                }
                // Dynamic friction
                else
                {
                    fx = -0.8 * fz * c.velocity[i][0] / vxy;
                    fy = -0.8 * fz * c.velocity[i][1] / vxy;
                    c.force[i][0] += fx;
                    c.force[i][1] += fy;
                    c.force[i][2] += fz;
                }
            }

        }
        // Update parameters
        for (int i = 0; i < c.ball_num; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                //                c.acceration[i][j] = c.force[i][j] / m;
                c.velocity[i][j] += c.force[i][j] / m * dt;
                c.velocity[i][j] *= damping;
                c.position[i][j] += c.velocity[i][j] * dt;
                if (recordFlag > 0){
                    final_result[recordFlag][t][i][j] = c.position[i][j];
                }
            }
        }
        if(t == 10)
        {
            for (int i = 0; i < c.ball_num; i++)
            {
                temp_path[0] += c.position[i][0];
                temp_path[1] += c.position[i][1];
            }
        }
        t++;
    }
    
    for (int i = 0; i < c.ball_num; i++)
    {
        final_path[0] += c.position[i][0];
        final_path[1] += c.position[i][1];
    }
    
    delta_p[0] = final_path[0] - temp_path[0];
    delta_p[1] = final_path[1] - temp_path[1];
    fitness += norm(delta_p)/c.ball_num;
    cout<< final_path[0]<< ' '<< final_path[1] <<endl;
    cout<< temp_path[0]<< ' '<< temp_path[1] <<endl;
    cout<< c.ball_num<<endl;

    return fitness;
}

//////////////////////////////////////////////////////////////////////////////////////////////
void myDisplay(void)
{
    
    glEnable(GL_DEPTH_TEST);
    // 清除屏幕及深度缓存
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    // 深度测试开启，实现遮挡关系
    glDepthFunc(GL_LEQUAL);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT,GL_NICEST);
    
    //Background color
    glClearColor(0.0f,0.0f,0.0f,0.0f);
    
    //Anti-aliasing
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_POINT_SMOOTH);
    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    
    //Sight Seting
    glEnable(GL_DEPTH_TEST);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(35, 2.5, 1.4, 40000000);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(1.0, 2.5, 1, 0, 0, 0, 0, 0, 1);
    
    int k=0;
    for (float j=2.0; j>0; j=j-0.1){
        if (k==0) {
            glBegin(GL_POLYGON);
            glColor3f(0.1f, 0.1f, 0.1f);
            for(int i=0; i<30; ++i) {
                glVertex3f(j*cos(2*pi/30*i), j*sin(2*pi/30*i), 0.0f);
            }
            k=1;
            glEnd();
        }
        else {
            glBegin(GL_POLYGON);
            glColor3f(0.5f, 0.5f, 0.5f);
            for(int i=0; i<30; ++i) {
                glVertex3f(j*cos(2*pi/30*i), j*sin(2*pi/30*i), 0.0f);
            }
            k=0;
            glEnd();
        }
    }
    
    //XYZ axis
    glLineWidth(1.0f);
    glBegin(GL_LINES);
    glColor3f(0.0f, 1.0f, 1.0f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 0.0f,2.0f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 2.0f,0.0f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(2.0f, 0.0f, 0.0f);
    glEnd();
    
    double xx = 0;
    double yy = 0;
    glLineWidth(3.0f);
    glBegin(GL_LINES);
    glColor3f(1.0f, 1.0f, 0.0f);
    for (int k = 1; k<=recordNumber; k++){
        for ( int i = 0; i < final_ball_num; i++){
            for (int j = i+1; j < final_ball_num; j ++){
                if (final_l0_type[i][j] == 1)
                {
                    glColor3f(0.2f*k/recordNumber, 1.0f*k/recordNumber, 1.0f*k/recordNumber);
                    glVertex3d(final_result[k][day][i][0]-0.05-xx,final_result[k][day][i][1]-0.05-yy,final_result[k][day][i][2]);
                    glVertex3d(final_result[k][day][j][0]-0.05-xx,final_result[k][day][j][1]-0.05-yy,final_result[k][day][j][2]);
                }
                else if (final_l0_type[i][j] == 2)
                {
                    glColor3f(1.0f*k/recordNumber, 0.392f*k/recordNumber, 1.0f*k/recordNumber);
                    glVertex3d(final_result[k][day][i][0]-0.05-xx,final_result[k][day][i][1]-0.05-yy,final_result[k][day][i][2]);
                    glVertex3d(final_result[k][day][j][0]-0.05-xx,final_result[k][day][j][1]-0.05-yy,final_result[k][day][j][2]);
                }
                else if (final_l0_type[i][j] == 3)
                {
                    glColor3f(0.4f*k/recordNumber, 1.0f*k/recordNumber, 0.0f*k/recordNumber);
                    glVertex3d(final_result[k][day][i][0]-0.05-xx,final_result[k][day][i][1]-0.05-yy,final_result[k][day][i][2]);
                    glVertex3d(final_result[k][day][j][0]-0.05-xx,final_result[k][day][j][1]-0.05-yy,final_result[k][day][j][2]);
                }
                else if (final_l0_type[i][j] == 4)
                {
                    glColor3f(1.0f*k/recordNumber, 1.0f*k/recordNumber, 0.0f*k/recordNumber);
                    glVertex3d(final_result[k][day][i][0]-0.05-xx,final_result[k][day][i][1]-0.05-yy,final_result[k][day][i][2]);
                    glVertex3d(final_result[k][day][j][0]-0.05-xx,final_result[k][day][j][1]-0.05-yy,final_result[k][day][j][2]);
                }
            }
        }
    }
    glEnd();
    glFlush();
    glutSwapBuffers();
}

//////////////////////////////////////////////////////////////////////////////////////////////
void myIdle(void)
{
    day++;
    while (day % 4 > 0)
        day++;
    if (day>T)
        day = 0;
    myDisplay();
}

//////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
    int spring[512] = {0};
    
    
    char data[100];
    
    ifstream infile;
    infile.open("100.txt");
    
    for (int i=0; i<512; i++){
        infile >> data;
        spring[i] = atoi(data);
        cout<< spring[i] <<endl;
    }
    // 关闭打开的文件
    infile.close();
    double fitnessAtTheEnd = get_fitness(spring, T);

    //double fitness = get_fitness(1000);
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(1250, 700);
    glutCreateWindow("第一个 OpenGL 程序");
    glutDisplayFunc(&myDisplay);
    glutIdleFunc(&myIdle);
    glutMainLoop();
    return 0;
}




