#ifndef POLYSHAPE_H
#define POLYSHAPE_H

//enclosedPoint function makes use of:
//https://wrf.ecse.rpi.edu/Research/Short_Notes/pnpoly.html

//getOrientation function makes use of:
//http://www.science.smith.edu/~jorourke/Code/polyorient.C

#include <clipper.hpp>
#include <helper.h>
#include <array>

using namespace ClipperLib;
namespace Poly
{
    class Polyshape
    {
    public:

        Polyshape();
        Polyshape polyBuffer(double d);
        void addHole(Path pgon);
        void addRegion(Path pgon);
        Path getHole(uint8_t nr);
        Path getRegion(uint8_t nr);
        Path getPgon(uint8_t nr);
        bool enclosedPoint(IntPoint p, uint8_t nr);
        
        bool lineIntersect(IntPoint &p1, IntPoint dir1, uint8_t nr, const bool &length, Path &points);
        
        uint8_t numRegions = 0;
        uint8_t numHoles = 0;
        std::vector< Path > regions;
        std::vector< Path > holes;

    private:
            
    };

    bool isHole( Path pgon );
    bool isRegion( Path pgon );  
    bool isCollinear(IntPoint p1, IntPoint p2, IntPoint p3, IntPoint p4);
    bool lineIntersect(IntPoint &p1, IntPoint dir1, Path &path, const bool &length, Path &points);
    bool enclosedPoint(IntPoint p, Path &pgon);
    int findDistance(IntPoint p1, IntPoint dir1, IntPoint p3);
    bool getOrientation(Path &pgon);
    Path splitOutline(IntPoint p1, IntPoint p2, Path A);
    std::array<float, 2> intersect(IntPoint &p1, IntPoint &p2, IntPoint &p3, IntPoint &p4);
    IntPoint sum_IntPoints(const std::vector<IntPoint> &v);
    IntPoint mean_IntPoint(const std::vector<IntPoint> &v);
    
}
#endif
