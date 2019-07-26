#ifndef LINEAR_HPP
#define LINEAR_HPP

#include <map>
#include <vector>
#include "data/sample.hpp"

namespace Anaquin
{
    struct LinearModel
    {
        // Constant coefficient
        double c = NAN;
        
        // Least-squared slope coefficient
        double m = NAN;
        
        // Adjusted R2
        double R2 = NAN;
        
        // Pearson correlation
        double r = NAN;
        
        // Adjusted R2
        double aR2 = NAN;
        
        double F = NAN;
        double p = NAN;
    };

    struct Point
    {
        Point(double x = 0.0, double y = 0.0) : x(x), y(y) {}
        
        double x, y;
    };

    struct SequinStats : public std::map<SequinID, Point>
    {
        struct Data
        {
            std::vector<SequinID> ids;
            std::vector<double> x, y;

            std::map<SequinID, double> id2x;
            std::map<SequinID, double> id2y;
        };
        
        inline void add(const SequinID &id, double x, double y)
        {
            (*this)[id] = Point(x, y);
        }
        
        inline void sum(const SequinID &id, double x, double y)
        {
            if (!count(id)) { add(id, x, y); } else { (*this)[id].y += y; }
        }

        Limit limitQuant() const;
        
        // Return the values after filtering
        Data data(bool shouldLog, bool ignoreZero) const;
        
        LinearModel linear(bool shouldLog = true, bool ignoreZero = false) const;
    };
}

#endif
