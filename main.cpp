#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;
//Write a program to explore this system by computing the effective gravity (a vector)
// and potential (a scalar) at specific grid points in the xy-plane. Use units such that
// the gravitational constant G ≡ 1. The program should take as input the mass of both
// bodies and their separation. You may “hardwire” the grid dimensions into your code
// if you wish. The output should be the potential and x and y components of the
// acceleration at each grid point.
double m1,m2,d;
const float G=1;
const float FEPS = 1e-6;
const double XMIN = -2.0;
const double XMAX = 2.0;
const double DX = 0.08;

const double YMIN = -2.0;
const double YMAX = 2.0;
const double DY = 0.08;

double m = m1 + m2;
double w2 = m / (d ** 3);
double x1 = -m2 * d / m;
double x2 = m1 * d / m;
double xc = 0.5 * (x1 + x2);


double gravPotential(double m1, double m2, double d, double x, double y){
    double r1, r2;
    r1 = sqrt((x+d/2)*(x+d/2)+y*y);
    r2 = sqrt((x-d/2)*(x-d/2)+y*y);
    float potential;
    potential = (-(G*m1/r1 + G*m2/r2));
    return potential;
}

//
double rootFinder(double right, double left, double d){
    while (right - left > FEPS){
        double middle, midmid, lPoint;
        midmid = ((d/2), 0);
        middle = (left+right)/2;
        if (midmid < 0) { left = middle; }
        else { right = middle; }
        lPoint = (right+left)/2;
        return lPoint;
    }
}
double gravx(double ir13, double ir23, double x, double y){
    return -m1 * (x - x1) * ir13 - m2 * (x - x2) * ir23 + w2 * x;
}
double  gravy(double ir13, double ir23, double x, double y) {
    return -m1 * y * ir13 - m2 * y * ir23 + w2 * y;
}

double Fx(double x, double y){
    double r1, r2, a, b;
    r1 = sqrt(pow(y,2) + pow((x+d/2),2)); //Square root ( y^2 + (x+d/2)^2) )
    r2 = sqrt(pow(y,2) + pow((x-d/2),2)); //Square root ( y^2 + (x-d/2)^2) )
    a = (G * m1 * (x-d/2) ) / (pow(r1,3)); //(G*m1*(x-d)/2)/r1^3
    b = (G * m2 * (x+d/2) ) / (pow(r2,3)); //(G*m1*(x-d)/2)/r1^3
    return a - b;
}

double Fy(double x, double y){
    double r1, r2, j, k;
    r1 = sqrt(pow((x+d/2),2)) + pow(y,2); // Square root ((x+d/2)^2) + y^2
    r2 = sqrt(pow((x-d/2),2)) + pow(y,2); // Square root ((x-d/2)^2) + y^2
    j = (G * m1 / pow(r1,3)); // G*m1/r1^3
    k = (G * m2 / pow(r2,3)); // G*m2/r2^3
    return y*(j+k);
}
pair<double,double> gravAcc(double x, double y){
    return {Fx(x, y),Fy(x, y)};
}

double rtbis( double(*f)(double), double xlo, double xhi, double eps ) {
    double flo = f(xlo);
    int iter = 0;
    double xm;
    while ((xhi - xlo) > (0.5 * eps * fabs(xlo + xhi))){
        iter++;
        xm = 0.5 * ( xlo + xhi );
        double fm = f(xm);
        if ( fabs(fm) < 1e-20 ) {
            break;
        }
        if ( flo < 0 ) {
            if ( fm < 0 ){
                xlo = xm;
            }else{
                xhi = xm;
            }
        }else{
            if ( fm < 0 ){
                xhi = xm;
            }else{
                xlo = xm;
            }
        }
        if ( iter > 2000 ){
            printf( "maxiter reached in rtbis: %g %g %g %g %g %g", xlo, xhi, xm, xhi - xlo, 0.5 * eps * fabs(xhi + xlo), fm );
            break;
        }
    }
    return xm;
}

int main() {
    //init Variables
    double x, y, L1, L2, L3, L4, L5;

    //User Prompts
    cout << "This is a Program computes the effective gravity (a vector)\n"
            "and potential (a scalar) at specific grid points in the xy-plane." << endl << endl;
    cout << "Below, you are going to the prompted to add the mass of two bodies, " << endl
        << "m1 and m2 as well as the distance." << endl << endl;
    cout << "What is the mass of the first body, m1: ";
    cin >> m1;
    cout << "What is the mass of the second body, m2: ";
    cin >> m2;
    cout << "What is the distance between the two masses, d: ";
    cin >> d;

    float m = m1+m2/2;
    float x2 = m1*d/m ;
    float x1 = -m2*d/m;


    //file output
    ofstream outFile;
    ofstream Lpoints;
    ofstream mass;
    outFile.open("lagrange_grid.dat");
    Lpoints.open("lagrange_root.dat");
    mass.open("lagrange_mass.dat");

    //Lagrangian Points


    for (x = XMIN; x <= XMAX; x += 0.1){
        for (y = YMIN; y <= YMAX; y += 0.1){
            double gravp, middle;
            cout << "Writing grad data to Lagrange_Grid.dat." << endl;

            gravp = gravPotential(m1, m2,x,y,d);
            pair<double,double> gravAcc(x, y);

            L1 = rtbis(-d/2,0,d);
            L2 = rtbis(0,d/2,d);
            L3 = rtbis(d/2,2*d,d);
            middle = (L1+L2)/2;
            L4 = (middle -L1) * sqrt(3); //y coordinate
            L5 = -1*L4;  // x coordinate

            cout << x << "," << y << "\n";
            cout << "Gravity Potential: " << gravp << endl;
            cout << "Gravity Acceleration: " << gravAcc.first*10 << ", " << gravAcc.second*10 << endl;


            double y2 = pow(y, 2.0);
            double r2 = (x*x + y*y);
            double r12 = pow((x - x1), 2.0) + y2;
            double r22 = pow((x - x2), 2.0) + y2;
            r12 = fmax(r12, FEPS);
            r22 = fmax(r22, FEPS);
            double isqrtr12 = 1 / sqrt(r12);
            double isqrtr22 = 1 / sqrt(r22);
            double ir13 = isqrtr12 / r12;
            double ir23 = isqrtr22 / r22;
            double gx = gravx(ir13, ir23, x, y);
            double gy = gravy(ir13, ir23, x, y);
            double p = -m1 * isqrtr12 - m2 * isqrtr22 - 0.5 * w2 * r2;
            printf("%15.8E %15.8E %15.8E %15.8E %15.8E\n", x, y, p, gx, gy);
        //    cout << x << ", " << y << " "  << gravp << "  "
           //         <<  gravAcc.first*10 << ", " << gravAcc.second*10 << endl;
           // outFile << x << ", " << y << " "  << gravp << "  "
            //        << gravAcc.first*10 << ", " << gravAcc.second*10 << endl;
        }

    }
    */
    cout << "writing mass data to mass file\n";
    mass << m1 << endl << m2 << x1 << x2 << endl;

    Lpoints << "Your calculated Lagrangian points are: \n"
            << "L1 = (" << L1 <<", 0)\n"
            << "L2 = (" << L3 <<", 0)\n"
            << "L3 = (" << L3 <<", 0)\n"
             //double check if the value of middle works
            << "L4 = (" << middle << ", "<< L4 << ")\n" // essentially in the y component
            << "L5 = ("<< middle << ", "<< L5 << ")\n"; // essentially in the y component

    return 0;
}





/*


#include <iostream>
#include <cmath>

const int num_points = 10; // grids of points
const double G = 1; // gravitational constant

// function for potential and effective gravity (vector)
 void find_acc_poten(double m1, double m2, double sep, double &poten, double &xacc, double &yacc) {
    double x1 = -sep/2;  // x-position of M1
    double x2 = sep/2; // x-position of M2
    double y1=0,y2=0; // both at the same y position

    double x,y;
    double rad;

    poten= 0;
    xacc=0;
    yacc=0;

    // sum each of the points' contribution
    for (int i=-num_points; i<=num_points; i++) {
        for (int j=-num_points; j<=num_points; j++) {
            x= (double)i/(double)num_points;
            y= (double)j/(double)num_points;
            rad = sqrt(pow(x-x1,2)+pow(y-y1,2));
            xacc+=-m1*(x-x1)/pow(rad,3);
            yacc+=-m1*(y-y1)/pow(rad,3);
            poten+= -m1/rad;

            rad = sqrt(pow(x-x2,2)+pow(y-y2,2));
            xacc+=-m2*(x-x2)/pow(rad,3);
            yacc+=-m2*(y-y2)/pow(rad,3);
            poten+= -m2/rad;
        }
    }

}

int main () {

    double m1, m2, sep;
    std::cout << "Mass1: ";
    std::cin >> m1;
    std::cout << "Mass2: ";
    std::cin >> m2;
    std::cout << "Separation: ";
    std::cin >> sep;

    double potential, xacc, yacc;

    find_acc_poten(m1,m2,sep,potential,xacc,yacc);

    std::cout << "Potential: " << potential << std::endl;
    std::cout << "x-acceleration: " << xacc << std::endl;
    std::cout << "y-acceleration: " << yacc << std::endl;

    return 0;

}
////////////////////////////////////////////////

#file
GRID_FILE = 'lagrange_grid.dat'
ROOT_FILE = 'lagrange_root.dat'

XMIN = -2.0
XMAX = 2.0
DX = 0.1

YMIN = -2.0
YMAX = 2.0
DY = 0.1

FEPS = 1e-6

m1 = float(sys.argv[1])
m2 = float(sys.argv[2])
d = float(sys.argv[3])

m = m1+m2
w2 = m/(d**3)
x1 = -m2*d/m
x2 = m1*d/m

xc = 0.5*(x1+x2)

def rtbis(f, xlo, xhi, eps):
    flo = f(xlo)
    iter = 0
    while(xhi-xlo > 0.5*eps*abs(xlo+xhi)):
        iter = iter+1
        xm = 0.5*(xlo+xhi)
        fm = f(xm)
        if abs(fm) < 1e-20:
            break
        if (flo < 0):
            if fm < 0:
                xlo = xm
            else:
                xhi = xm
        else:
            if fm < 0:
                xhi = xm
            else:
                xlo = xm
        if iter > 2000:
            print('maxiter reached in rtbis: %g %g %g %g %g %g'%(xlo,xhi,xm,xhi-xlo,0.5*eps*abs(xhi+xlo),fm))
            break
    return xm

def gravx(ir13, ir23, x, y):
    return -m1*(x-x1)*ir13 - m2*(x-x2)*ir23 + w2*x

def gravy(ir13, ir23, x, y):
    return -m1*y*ir13 - m2*y*ir23 + w2*y

def fx(x):
    r12 = (x-x1)**2
    r22 = (x-x2)**2
    return gravx(1/r12**1.5, 1/r22**1.5, x, 0.0)

def fy(y):
    r12 = (xc-x1)**2 + y**2
    r22 = (xc-x2)**2 + y**2
    return gravy(1/r12**1.5, 1/r22**1.5, xc, y)

print('writing grid data to %s'%GRID_FILE)
xg = arange(XMIN, XMAX+0.1*DX, DX)
yg = arange(YMIN, YMAX+0.1*DY, DY)
with open(GRID_FILE,'w') as f:
    for x in xg:
        for y in yg:
            y2 = y**2
            r2 = x*x + y*y
            r12 = (x-x1)**2 + y2
            r22 = (x-x2)**2 + y2
            r12 = max(r12, FEPS)
            r22 = max(r22, FEPS)
            isqrtr12 = 1/sqrt(r12)
            isqrtr22 = 1/sqrt(r22)
            ir13 = isqrtr12/r12
            ir23 = isqrtr22/r22
            gx = gravx(ir13, ir23, x, y)
            gy = gravy(ir13, ir23, x, y)
            p = - m1*isqrtr12 - m2*isqrtr22 - 0.5*w2*r2
            f.write('%15.8E %15.8E %15.8E %15.8E %15.8E\n'%(x,y,p,gx,gy))

print('writting lagrange points to %s'%ROOT_FILE)
with open(ROOT_FILE, 'w') as f:
    #L1: on x axis, between x1 and x2
    print('L1')
    lp = rtbis(fx, x1+FEPS, x2-FEPS, FEPS)
    f.write('%15.8E %15.8E\n'%(lp,0.))
    #L2: on x axis, to the right of x2
    print('L2')
    lp = rtbis(fx, x2+FEPS, XMAX, FEPS)
    f.write('%15.8E %15.8E\n'%(lp,0.))
    #L3: on x axis, to the left of x1
    print('L3')
    lp = rtbis(fx, XMIN, x1-FEPS, FEPS)
    f.write('%15.8E %15.8E\n'%(lp,0.))
    #L4: on positive x=xc line
    print('L4')
    lp = rtbis(fy, FEPS, YMAX, FEPS)
    f.write('%15.8E %15.8E\n'%(xc, lp))
    #L5: on negative x=xc line
    print('L5')
    lp = rtbis(fy, YMIN, -FEPS, FEPS)
    f.write('%15.8E %15.8E\n'%(xc, lp))







 */