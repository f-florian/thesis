#include <gsl/gsl_integration.h>

class Integral
{
public:
    Integral(const double start, const double end, const int minsize);
    double get(const double start, const double end) const;
private:
    double *data;
    double max;
    size_t size;
};
