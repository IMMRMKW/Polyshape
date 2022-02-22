#include <polyshape.h>
using namespace ClipperLib;
using namespace helper;

namespace Poly
{

Polyshape::Polyshape()
{
    regions.reserve(5);
    holes.reserve(5);
}

/*
Returns a polyshape object that is offset d mm from the polyshape object on which
the function is executed. A negative value d means that regions will be offset with 
negative d, and holes with a positive d.
inputs: buffer distance d in mm.
returns: polyshape object containing the buffered elements
*/
Polyshape Polyshape::polyBuffer(double d)
{
    Polyshape pgon;
    Paths solution;
    ClipperOffset co;
    co.AddPaths(regions, jtMiter, etClosedPolygon);
    co.AddPaths(holes, jtMiter, etClosedPolygon);
    co.Execute(solution, d);
    
    for (size_t i = 0; i < solution.size(); i++)
    {

        bool region = true;
        for (size_t j = 0; j < solution.size(); j++)
        {
            
            if (j != i )
            {
                if (Poly::enclosedPoint(solution[i][0], solution[j]))
                {
                    region = false;
                    break;
                }
            }
        }
        if (region)
        {
            pgon.addRegion(solution[i]);
        }
        else
        {
            pgon.addHole(solution[i]);
        }
        //Somehow, the resulting polygons are wrongly oriented?!
        /*if ( isRegion( solution[i] ) )
        {
            helper::printPath(solution[i]);
            pgon.addRegion(solution[i]);
        }
        else
        {
            pgon.addHole(solution[i]);
        }*/
    }

    return pgon;
}

bool Polyshape::enclosedPoint(IntPoint p, uint8_t nr)
{
    if( nr < numRegions)
    {
        return Poly::enclosedPoint(p, regions[nr]);
    }
    else if ( nr < numRegions + numHoles )
    {
        return Poly::enclosedPoint(p, holes[nr-numRegions]);
    }
    else //nr is too high, not so many polygons in polyshape
    {
        return false;
    }
}

bool Polyshape::lineIntersect(IntPoint &p1, IntPoint dir1, uint8_t nr, const bool &length, Path &points)
{
    bool partOf = false;
    if( nr < numRegions)
    {
        partOf = Poly::lineIntersect(p1, dir1, regions[nr], length, points);
    }
    else
    {
        partOf = Poly::lineIntersect(p1, dir1, holes[nr-numRegions], length, points);
    }

    return partOf;

}

Path Polyshape::getHole(uint8_t nr)
{
    return holes[nr];
}

Path Polyshape::getRegion(uint8_t nr)
{
    return regions[nr];
}

Path Polyshape::getPgon(uint8_t nr)
{
    if(nr < numRegions)
    {
        return regions[nr];
    }
    else{
        Serial.println(nr-numRegions);
        return holes[nr-numRegions];
    }
}

/*
Adds Path pgon as a hole, regardless of its orientation. Will 
reorient pgon to fit the convention.
Input: Path pgon with IntPoints with the coordinates of the points
*/
void Polyshape::addHole(Path pgon)
{
    if(!getOrientation(pgon))
    {
        std::reverse(pgon.begin(), pgon.end());
        getOrientation(pgon);
    }
    holes.push_back(pgon);
    this->numHoles++;
}

/*
Adds Path pgon as a region, regardless of its orientation. Will 
reorient pgon to fit the convention.
Input: Path pgon with IntPoints with the coordinates of the points
*/
void Polyshape::addRegion( Path pgon )
{
    if(getOrientation(pgon))
    {
        std::reverse(pgon.begin(), pgon.end());
        getOrientation(pgon);
    }
    regions.push_back(pgon);
    this->numRegions++;
}

/*
Checks whether Path pgon is a hole.
Input: Path pgon with IntPoints
Output: bool. True if pgon is a hole.
*/
bool isHole(Path pgon)
{
    if(!isRegion(pgon))
    {
        return true;
    }
    else
    {
        return false;
    }
}

/*
Checks whether Path pgon is a region.
Input: Path pgon with IntPoints
Output: bool. True if pgon is a region.
*/
bool isRegion( Path pgon )
{
    if( getOrientation(pgon) )
    {
        return false;
    }
    else
    {
        return true;
    }
}

/*
Determines if two line segments are collinear.
inputs: line segment 1, defined by IntPoint p1 and IntPoint p2, line segment 2, defined by 
IntPoint p3 and IntPoint p4.
returns: boolean, true if collinear, false if not collinear
*/
bool isCollinear(IntPoint p1, IntPoint p2, IntPoint p3, IntPoint p4)
{
    IntPoint dir1, dir2;
    
    dir1 = 100*(p1 - p2) / (p1 - p2).norm();
    dir2 = 100*(p4 - p3) / (p4 - p3).norm();
    
    if ( (dir1 - dir2).norm() == 0 || (dir1 + dir2).norm() == 0 )
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool lineIntersect(IntPoint &p1, IntPoint dir1, Path &pgon, const bool &length, Path &points)
{
    bool partOf = false;
    IntPoint dir2;
    dir2.X = dir1.Y;
    dir2.Y = -dir1.X;

    float t, min_t, max_t;
    
    for (size_t i = 0; i < pgon.size(); i++)
    {
        t = (float) (dir1.X*(p1.Y-pgon.at(i).Y)-dir1.Y*(p1.X-pgon.at(i).X))/
        (dir2.Y*dir1.X-dir2.X*dir1.Y);
        if (i == 0)
        {
            min_t = t;
            max_t = t;
        }
        else if (t < min_t)
        {
            min_t = t;
        }
        else if (t > max_t)
        {
            max_t = t;
        }
    }
    
    if (sgn(min_t) != sgn(max_t))
    {
        int j = pgon.size()-1;
        IntPoint p2 = p1 + dir1;
        for (size_t i = 0; i < pgon.size(); i++)
        {
            
            IntPoint p3 = pgon.at(i);
            IntPoint p4 = pgon.at(j);

            if (!isCollinear(p1, p2, p3, p4))
            {
                
               std::array<float, 2> t = Poly::intersect(p1, p2, p3, p4);

               if (t[1] >= 0 && t[1] <= 1)
               {
                   if (length)
                   {
                       if (t[0] >= 0 && t[0] <= 1)
                       {
                           points.push_back(p1 + (p2-p1) * t[0]);
                       }
                   }
                   else
                   {
                       points.push_back(p1 + (p2-p1) * t[0]);
                   }
               }
            }
            else
            {
                if (length)
                {
                    if ((p3-p1).norm() == 0 && (p4-p2).norm() == 0)
                    {
                        partOf = true;
                    }
                    else if ((p3-p2).norm() == 0 && (p4-p1).norm() == 0)
                    {
                        partOf = true;
                    }
                }
            }
            j = i;
        }
    }
    
    return partOf;

}

bool enclosedPoint(IntPoint p, Path &pgon)
{

    uint32_t i, j = 0;
    uint8_t c = 0;
  
    for (i = 0, j = pgon.size() - 1; i < pgon.size(); j = i++) 
    {
        if ( ((pgon.at(i).Y > p.Y) != (pgon.at(j).Y > p.Y)) &&
	        (p.X < (pgon.at(j).X - pgon.at(i).X) * (p.Y - pgon.at(i).Y) / (pgon.at(j).Y - pgon.at(i).Y) + pgon.at(i).X) )
        {
            c = !c;
        }   
    }

    if (c)
    {
        return true;
    }
    else
    {
        return false;
    }
}

/*
Finds the distance of a point, p3, with respect to a line defined by 
origin p1 and direction dir1. 
Inputs:
IntPoint p1, which is the origin of the line
IntPoint dir1 representing the direction of the line
IntPoint p3, the coordinate of which we want to know the distance
Output: integer containing the distance in centimeters.
*/
int findDistance(IntPoint p1, IntPoint dir1, IntPoint p3)
{
    IntPoint p2 = p1 + dir1;
    IntPoint dir2(dir1.Y, -dir1.X);
    IntPoint p4 = p3 + dir2;
    std::array<float, 2> t = Poly::intersect(p1, p2, p3, p4);
    return helper::ftoi(dir1.norm()*t[1]);
}

/*
Determines whether the polygon defined by pgon is oriented (anti)clockwise
inputs: Path pgon, a polygon defined by intPoints defining the x- and y-coordinates
returns: boolean, true if anticlockwise, false if clockwise
*/
bool getOrientation(Path &pgon)
{
    size_t m,n;
    IntPoint min;

    n = pgon.size();

    min = pgon[0];
    m = 0;

    for (size_t i = 0; i < n; i++ ) 
    {
        if ( (pgon[i].Y < min.Y) || ( (pgon[i].Y == min.Y) && (pgon[i].X > min.X) ) )
        {
            m = i;
            min = pgon[i];
        }
    }
  
    int area;
    
    IntPoint a, b, c; //Just renaming

    size_t m1 = (m + (n-1)) % n; // = m - 1
  
    //asign the coordinates:
    a = pgon[m1];
    b = pgon[m];
    c = pgon[(m + 1) % n];
     
    area =
                a.X * b.Y - a.Y * b.X +
                a.Y * c.X - a.X * c.Y +
                b.X * c.Y - c.X * b.Y;

    //Serial.print("Area: ");
    //Serial.println(area);
    if (area > 0)
    {
        return true;
    }
    else
    {
        return false;
    }
}

/*
Splits an outline with two intersection points laying somewhere on the loop. It locates between which two points
each intersection point lays, and adds them to the loop. Simultaneously, if it comes across p1, it starts removing
points from the loop, and adding it to another, until it comes across p2. The result are two polygons, which share 
the line p1-p2.
Inputs:
IntPoint p1, the first intersection point
Intpoint p2, the second intersection point
Path A, a clockwise oriented polygon indicating an outline of the area
Output: the anticlockwise polygon, which is a result from splitting.
*/
Path splitOutline(IntPoint p1, IntPoint p2, Path A)
{
    float epsilon = 1;
    size_t i = 0;
    size_t j = 0;
    Path pgon;
    pgon << p1;
    while (j < 2)
    {
        IntPoint p3 = A[i % A.size()];
        IntPoint p4 = A[(i+1) % A.size()];
        if (j == 0 )
        {
            int t = findDistance(p3, p4-p3, p1);
            if(abs(t) <= epsilon)
            {
                A.insert(A.begin()+(i+1) % A.size(), p1);
                //helper::printPath(A);
                j++;
                i++;
            }
            i++;
        }
        else if ( j == 1)
        {
            pgon << p3;
            i = i % A.size();
            A.erase(A.begin() + i % A.size());
            //Serial.println(i%A.size());
            float t = findDistance(p3, p4-p3, p2);
            if(abs(t) <= epsilon)
            {
                j++;
            }
            
        }
    }
    A << p2;
    pgon << p2;
    
    //An outline is oriented clockwise, the hole anticlockwise
    //the larger section is automatically oriented clockwise
    //whereas the smaller is anticlockwise
    if (getOrientation(A))
    {
        return A;
    }
    else
    {
        return pgon;
    }
}

/*
Determines where two line segments intersect by solving the following equation p1 + t[0] * (p2-p1) = p3 + t[1] * (p4-p3)
input: line segment 1, defined by IntPoint p1 and IntPoint p2, line segment 2, defined by 
IntPoint p3 and IntPoint p4.
returns: an array with values for t
*/
std::array<float, 2> intersect(IntPoint &p1, IntPoint &p2, IntPoint &p3, IntPoint &p4)
{
    std::array<float, 2> t;
    IntPoint dir1 = p2-p1;
    IntPoint dir2 = p4-p3;

    t[0] = (float) (dir2.X*(p1.Y-p3.Y)-dir2.Y*(p1.X-p3.X))/
        (dir2.Y*dir1.X-dir2.X*dir1.Y);

    t[1] = (float) (dir1.X*(p1.Y-p3.Y)-dir1.Y*(p1.X-p3.X))/
        (dir2.Y*dir1.X-dir2.X*dir1.Y);

    return t;
}

/*
Sums all IntPoints in a vector.
Input: vector of IntPoints representing coordinates.
Output: IntPoint representing the sum of the coordinate
*/
IntPoint sum_IntPoints(const std::vector<IntPoint> &v)
{
    IntPoint result;
    result.X = std::accumulate(v.begin(), v.end(), 0,
    [&](int sum, const IntPoint& curr) { return sum + curr.X; });
    result.Y = std::accumulate(v.begin(), v.end(), 0,
    [&](int sum, const IntPoint& curr) { return sum + curr.Y; });
    return result;
}

/*
Returns the mean coordinate of a set of coordinates
Input: vector of IntPoints representing coordinates.
Output: IntPoint representing the mean coordinate
*/

IntPoint mean_IntPoint(const std::vector<IntPoint> &v)
{
    IntPoint p_mean = sum_IntPoints(v);
    p_mean = p_mean / v.size(); //operator / correctly rounds the outcome
    return p_mean;
}

} //namespace Poly

