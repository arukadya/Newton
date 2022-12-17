//
//  main.cpp
//  Newton
//
//  Created by 須之内俊樹 on 2022/11/18.
//
#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <sstream>
#include <iomanip>
#include "Eigen/Core"

Eigen::Vector2d f1(Eigen::Vector2d &x){
    return Eigen::Vector2d(
        x.x() + x.y(),
        x.x()*x.x() + x.y()*x.y() - 1
    );
}
Eigen::Vector2d f2(Eigen::Vector2d &x){
    return Eigen::Vector2d(
        x.y()*x.y()*x.y() - 3*x.x()*x.x()*x.y(),//y^3 - 3(x^2)y
        x.x()*x.x()*x.x() - 3*x.x()*x.y()*x.y() - 1 // x^3 - 3x(y^2)
    );
}
Eigen::Matrix2d f1prime(Eigen::Vector2d &x){
    Eigen::Matrix2d A;
        A <<
    1 , 1,
    2*x.x(), 2*x.y();
    return A;
}
Eigen::Matrix2d f2prime(Eigen::Vector2d &x){
    Eigen::Matrix2d A;
        A <<
    -6*x.x()*x.y() , 3*x.y()*x.y() - 3*x.x()*x.x(),
    3*x.x()*x.x() - 3*x.y()*x.y(), -6*x.x()*x.y();
    return A;
}
Eigen::Matrix2d invercef1prime(Eigen::Vector2d &x){
    Eigen::Matrix2d A = f1prime(x);
    Eigen::Matrix2d B;
    B << A(1,1),-A(0,1),
    -A(1,0),A(0,0);
    return B/(A(0,0)*A(1,1) - A(0,1)*A(1,0));
}
Eigen::Matrix2d invercef2prime(Eigen::Vector2d &x){
    Eigen::Matrix2d A = f2prime(x);
    Eigen::Matrix2d B;
    B << A(1,1),-A(0,1),
    -A(1,0),A(0,0);
    return B/(A(0,0)*A(1,1) - A(0,1)*A(1,0));
    
}
Eigen::Vector2d Newton1(Eigen::Vector2d &x){
    return x - invercef1prime(x)*f1(x);
}
Eigen::Vector2d Newton2(Eigen::Vector2d &x){
    return x - invercef2prime(x)*f2(x);
}

void mapping2(int Nx,int Ny,double eps,int maxRep){
    for(int i=-Nx/2;i<=Nx/2;i++)for(int j=-Ny/2;j<=Ny/2;j++){
        int cnt = 0;
        Eigen::Vector2d now_x((double)i/Nx*2,(double)j/Ny*2);
        std::cout << now_x.x() << " " << now_x.y() << " ";
        Eigen::Vector2d startVector = invercef2prime(now_x)*f2(now_x);
        startVector.normalize();
        unsigned int startVectorColor = startVector.norm();
        std::cout << startVector.x() << " " << startVector.y() << " " << startVectorColor << std::endl;
        Eigen::Vector2d next_x = Newton2(now_x);
        double res = (next_x - now_x).transpose()*(next_x - now_x);
        cnt++;
        while(res > eps && cnt < maxRep){
            now_x = next_x;
            next_x = Newton2(now_x);
            res = (next_x - now_x).transpose()*(next_x - now_x);
            cnt++;
        }
        std::cout << cnt << " ";
        unsigned int color = 0;
        color += (maxRep/cnt)*(1<<24);
        if((next_x.y() < sqrt(3)/2 + eps) && (next_x.y() > sqrt(3)/2 - eps) &&(next_x.x() < -0.5 + eps) && (next_x.x() > -0.5 - eps))color = (((1<<24) - (1<<16))/cnt)&((1<<24) - (1<<16));
        else if((next_x.y() < -sqrt(3)/2 + eps) && (next_x.y() > -sqrt(3)/2 - eps) && (next_x.x() < -0.5 + eps) && (next_x.x() > -0.5 - eps))color = (((1<<16) - (1<<8))/cnt)&((1<<16) - (1<<8));
        else if((next_x.y() < eps) && (next_x.y() > -eps) && (next_x.x() < 1 + eps) && (next_x.x() > 1 - eps))color = ((1<<8)-1)/cnt;
        else {
            color = (1<<24)-1;
        }
        //printf("%x ",color);
        //std::cout << color << std::endl;
    }
}

void mapping1(int Nx,int Ny,double eps,int maxRep){
    for(int i=-Nx/2;i<=Nx/2;i++)for(int j=-Ny/2;j<=Ny/2;j++){
        int cnt = 0;
        Eigen::Vector2d now_x((double)i/Nx*2,(double)j/Ny*2);
        std::cout << now_x.x() << " " << now_x.y() << " ";
        Eigen::Vector2d next_x = Newton1(now_x);
        double res = (next_x - now_x).transpose()*(next_x - now_x);
        cnt++;
        while(res > eps && cnt < maxRep){
            now_x = next_x;
            next_x = Newton1(now_x);
            res = (next_x - now_x).transpose()*(next_x - now_x);
            cnt++;
        }
        std::cout << cnt << " ";
        unsigned int color = 0;
        color += (maxRep/cnt)*(1<<24);
        if((next_x.x() < 1/sqrt(2) + eps) && (next_x.x() > 1/sqrt(2) - eps)) color =(((1<<24) - (1<<16))/cnt)&((1<<24) - (1<<16));
        else if((next_x.x() < -1/sqrt(2) + eps) && (next_x.x() > -1/sqrt(2) - eps)) color = (((1<<16) - (1<<8))/cnt)&((1<<16) - (1<<8));
        else color = (1<<24)-1;
        printf("%x ",color);
        std::cout << color << std::endl;
    }
}

int main(int argc, const char * argv[]) {
    mapping2(1000,1000,1.0e-4,20);
    return 0;
}
