#include "roots.hpp"
#include <cmath>

bool bisection(std::function<double(double)> f,
               double a, double b,
               double *root){

    //check that the function evaluated at both domain bounds are opposite signs            
    if( f(a) * f(b) >= 0){
        return false;
    }       
    
    double c;
    double tolerance = 1e-6;
    
    while( std::abs(b-a) >= tolerance ){

        c = (a + b)/2;

        if( f(c) == 0.0){
            break;
        }

        //checks if f(c) and f(a) are same sign
        if( f(c) * f(a) >=0 ){
            a = c;
        } else {
            b = c;
        }

    }//while

    *root = c;
    return true;

}// bisection

bool regula_falsi(std::function<double(double)> f,
                  double a, double b,
                  double *root){
    
    //check that the function evaluated at both domain bounds are opposite signs 
    if( f(a) * f(b) >= 0){
        return false;
    } 

    double c;
    double tolerance = 1e-6;

    while( std::abs(b-a) >= tolerance ){
        
        //c is the slope
        c = ( a*f(b) - b*f(a) ) / ( f(b) - f(a) );
        
        if( f(c) == 0.0){
            break;
        }

        //checks if f(c) and f(a) are same sign
        if( f(c) * f(a) >=0 ){
            a = c;
        } else {
            b = c;
        }

    }//while

    *root = c;
    return true;

}//regular_falsi

bool newton_raphson(std::function<double(double)> f,
                    std::function<double(double)> g,
                    double a, double b, double c,
                    double *root){
    
    double tolerance = 1e-6;
    int max_iterations = 1e6;
    double x0 = c;

    for(int i=0; i <= max_iterations; i++){

        //checks for if derivative is zero
        if( g(x0) == 0 ){
            return false;
        }

        double x1 = x0 - (f(x0)/g(x0));

        if (x1 < a || x1 > b) {
            return false; 
        }

        if( std::abs(x1- x0) < tolerance ){
            *root = x1;
            return true;
        }

        x0 = x1;

    }//for

    return false;
        
}//newton_raphson

bool secant(std::function<double(double)> f,
            double a, double b, double c,
            double *root){

    double tolerance = 1e-6;
    int max_iterations = 1e6;
    double x2;
    double x1 = c;
    double x0 = b;

    for(int i=0; i <= max_iterations; i++){

        //checks for if derivative is zero
        if( ( f(x0) - f(x1) ) == 0 ){
            return false;
        }

        x2 = x1 - f(x1) * ( ( x0 - x1 ) / ( f(x0) - f(x1) ) );

        if (x2 < a || x2 > b) {
            return false; 
        }

        if( std::abs(x2 - x1) < tolerance ){
            *root = x2;
            return true;
        }

        x0 = x1;
        x1 = x2;

    }//for

    return false;   

}//secant


