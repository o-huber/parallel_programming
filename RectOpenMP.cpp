
#include <omp.h>
#include <iostream>
#include <cmath>
#include <time.h>
#include <chrono>

using namespace std;
class SimpleTimer
{
public:
    SimpleTimer()
    {
        start = std::chrono::high_resolution_clock::now();
    }
    ~SimpleTimer() {
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<float> duration = end - start;
        cout << "DURATION " << duration.count() << "s" << endl;
    }
private:
    std::chrono::time_point<std::chrono::steady_clock> start, end;
};


double f(double x) {
    return 2 * cos(x) / (2 + x * x);
}

double Rect(double a, double b, int n) {

    double h = (b - a) / n;
    double sum = f(a) + f(b);
  
   #pragma parallel for num_threads(4)
    for (int i = 1; i <= n - 1; i++) {
        sum += 2 * f(a + i * h);
    }
    sum *= h / 2;

    return sum;
    }


int main()
{
    SimpleTimer st;
    double a, b;
    int n;
    a = 0;
    b = 100;
    n = 400;

    cout << "Function integral  F (x): " << Rect(a, b, n) << endl;

    return 0;
}

