//
// Created by 可乐 on 2017/10/28.
//

#ifndef POLYGONINTERSECTION_POLYGONINTERSECTION_H
#define POLYGONINTERSECTION_POLYGONINTERSECTION_H

#endif //POLYGONINTERSECTION_POLYGONINTERSECTION_H

#include <vector>
using namespace std;

enum pointType{vetex, cross};
typedef struct Point
{
    double x;
    double y;
    Point(double xx, double yy){x = xx; y = yy;}
    Point(const Point& pt){x = pt.x; y = pt.y;}
    Point(){x = 0; y = 0;}
}Point;
class Node
{
public:
    Node(Point* pt, pointType ptT, bool isTraverse = false, Node* next = NULL)
    {
        this->pt = pt;
        this->ptType = ptT;
        this->isTraverse = isTraverse;
        this->next = next;
    }

    Point* pt;
    pointType ptType;
    bool isTraverse;
    Node* next;
};

typedef vector<Point> Polygon;

vector<vector<vector<Polygon>>> PolygonSplit(Polygon bigPolygon, double side);
int ClockWise(const Polygon& p);
vector<Polygon> PolygonIntersection(Polygon poly1, Polygon poly2);
bool isPointOnPolygonEdge(const Point& pt, const Polygon& poly);
bool IsPointOnSegment(const Point& p, const Point& p1, const Point& p2);
bool detectIntersect(const Point& p1, const Point& p2, const Point& p3, const Point& p4, double& line_x, double& line_y);
bool PointCmp(Point* poly1, Point* poly2);
