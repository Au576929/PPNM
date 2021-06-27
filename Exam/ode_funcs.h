void rkstep12 (
        void f(int n, double complex x, double complex* y, double complex* dydx),
        int n,
        double complexx,
        double complex* yx,
        double complex h,
        double complex* yh,
        double complex* dy
);

void driver ( 
        void f(int n, double complex x, double complex* y, double complex* dydx),
        int n,
        double complex a,
        double complex* ya,
        double complex b,
        double t,
        double acc,
        double eps,
        char* name
);







