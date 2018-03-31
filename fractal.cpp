//written by Blake Lockett,
//the lack of functions is by design as each call on a large scale
//makes a real impact on performance. Some sectioning of the code has happened compared
//to original for ease of reading. 
#include <iostream>
#include <cmath>
#include <complex>
#include <omp.h>
#include <mpi.h>

// defines the use of the cimg png output, turn off live displaying
#define cimg_display 0
#define cimg_use_png 1
#include "CImg.h"

//break if previous if is less than this in |x1-x| this is a linear approximation
#define ACCURACY 0.00000001

using namespace std;
using namespace cimg_library;

/*
 * prototypes
 */
void getPortion(FILE *fptr, char *buffer);
double convertToDouble(char *num);
void readValues(FILE *fptr, double values[8]);
void generateFrames(double values[8], int nodes, int rank);
complex<double> newtonsMethod(complex<double> x);
complex<double> f(complex<double> x);
complex<double> fd(complex<double> x);

/*
 * implementation the function f(z)=z^5+3z^2+4z  f'(z)=5z^4+6z^2+4
 */
int main(int argc, char* argv[]) {
    //mpi setup
    MPI_Init(NULL, NULL);
    //number of nodes and which rank the current context is running in. 
    int nodes,rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //no input validation! it's slow.
    //0 rezw, 1 rezh, 2 real num(a), 3 imagine num(a), 4 real nums(b), 5 imagine nums(b), 6 frames, 7 zoom
    double values[] = {0,0,0,0,0,0,0,0};
    //passed in file
    FILE *fptr = fopen(argv[1],"r");
    //saves values from file into vars
    readValues(fptr, values);
    fclose(fptr);

    generateFrames(values, nodes,rank);
    MPI_Finalize();
    return 0;
}

//sets out the X by Y image and changes (0,0) from top left to the correct position.
void generateFrames(double values[8], int nodes,int rank) {
    double incrx,incry,rel,imga,x1,x2,y1,y2;
    int loopcheckerx = (int) values[0], loopcheckery = (int)values[1], maxxy = loopcheckerx*loopcheckery;
    int maxFrame = (int) values[6];
    for (int zoomFrame = rank; zoomFrame < maxFrame; zoomFrame += nodes) {
        double factor = (pow(values[7],zoomFrame));
        x1 = values[2]*factor;
        x2 = values[4]*factor;
        y1 = values[3]*factor;
        y2 = values[5]*factor;
        CImg<unsigned char> img((int) values[0], (int) values[1], 1, 3);
        incrx = (fabs(x1) + fabs(x2)) / values[0];
        incry = (fabs(y1) + fabs(y2)) / values[1];
		int x,y;	
#pragma omp parallel private(rel,imga) shared(img)
{
    #pragma omp for 
    for (x = 0; x < (loopcheckerx); x++) {
        rel = (x1 + (incrx * x));
        //needs to decrease not increase as we start at top
        for(y=0; y<loopcheckery; y++){
            imga = (y1 - (incry * y));
            complex<double> start(rel, imga);
            start = newtonsMethod(start);
            //this makes each root a unique colour
            if (start.real() != NAN && start.imag() != NAN) {
                img(x, y, 0) = ((int) (19 * start.imag() * 911) + 911) % 256;
                img(x, y, 1) = ((int) (19 * start.real() * 953) + 953) % 256;
                img(x, y, 2) = ((int) (19 * 953 + start.real() * 953 + start.imag())) % 256;
            }
        }
    }
}
        string result;
        string png = ".png";
        // 0000-9999 + room to change up to 999999999999999999999 on larger problems
        char numstr[21]; 
        sprintf(numstr, "%04d", zoomFrame);
        result = numstr + png;

        img.save_png(result.c_str());
    }
}

// reading values from an input file
void readValues(FILE *fptr, double values[8]) {
    int size = 28;
    //assumed size of file plus position in the line
    char* buffer = (char *)malloc(size * sizeof(char));

    if(fptr) {

        //first line
        getPortion(fptr, buffer);
        //rezW
        values[0] = convertToDouble(buffer);

        getPortion(fptr, buffer);
        //rezH
        values[1] = convertToDouble(buffer);

        //second line
        getPortion(fptr, buffer);
        //real a
        values[2] = convertToDouble(buffer);

        getPortion(fptr, buffer);
        //imag a
        values[3] = convertToDouble(buffer);

        //third line
        getPortion(fptr, buffer);
        //real b
        values[4]= convertToDouble(buffer);

        getPortion(fptr, buffer);
        //imag b
        values[5] = convertToDouble(buffer);

        //fourth line
        getPortion(fptr, buffer);
        values[7] = convertToDouble(buffer);

        //fithline
        getPortion(fptr, buffer);
        //frames
        values[6] = convertToDouble(buffer);

    }
    free(buffer);
}

//custom get int
void getPortion(FILE *fptr, char *buffer) {
    int  pos = 0;
    int c;//will be converted into char
    c = getc(fptr);
    while(c != ' ' && c != '\n' && c != EOF) {
        buffer[pos] = (char) c;
        c = getc(fptr);
        pos++;
    }
    //tells my atoi end of string
    buffer[pos] = 0;
    pos = 0;
}


//Planned to be faster than alternatives
double convertToDouble(char *num) {
    double val = 0;
    bool flag= false;
    int places = 0;
    int neg = 1;
    if(num[0] == '-'){
        neg = -1;
        *num++;
    }

    while(*num) {
        if (*num == '.'){
            *num++;
            flag = true;
        }
        //moves numbers over by *10 and calcs the number using the value of the string 0
        val = val*10 + (*num++ - '0');
        if (flag){
            places++;
        }

    }
    return (val/(pow(10, places)))*neg;
}


//used for calculating the root more accurately of the guessed x.
complex<double> newtonsMethod(complex<double> x) {
    //f0 is the inital value, x1 will be calculated x1 initialized as .1 to avoid trigering first if
    complex<double> x1(0.1,0.1);
    complex<double> x0 = x;
    //fx and the f'x of the function
    complex<double> fx = f(x0);
    complex<double> fdx = fd(x0);

    for (int i = 0; i < REPS; ++i) {

        //we are at the turning point of the function and could never find a tangent.
        if (fdx.real() == 0){
            x1 = (NAN,NAN);
            i = REPS;
        } else{
            x1 = ((x0) - (fx/fdx));

            if(fabs(x1.real() - x0.real())<= ACCURACY && fabs(x1.imag() - x0.imag())<= ACCURACY){
                i = REPS;
            } else {
                //set up the next iteration of newtons method
                x0 = x1;
                fx = f(x0);
                fdx = fd(x0);

                if (x1.real() == 0 && x1.imag()==0 && fdx.real() ==0 && fdx.imag() ==0){
                    x1 = (0,0);
                    i =REPS;
                }
            }
        }
    }
    return x1;
}


//f(z)=z^5+3z^2+4z  
complex<double> f(complex<double> x) {
    complex<double> a(pow(x,5));
    complex<double> b(3.0*pow(x,2));
    complex<double> c(4.0*(x));
    return (a+b+c);
}

//f'(z)=5z^4+6z^2+4
complex<double> fd(complex<double> x) {
    complex<double> a = 5.0*pow(x,4),b = 6.0*(x),c = (4.0,0.0);
    return (a+b+c);
}
