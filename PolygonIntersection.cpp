//
// Created by 可乐 on 2017/10/28.
//

#include <iostream>
#include <cmath>
#include "PolygonIntersection.h"


vector<vector<vector<Polygon>>> PolygonSplit(Polygon bigPolygon, double side)
{
    vector<vector<vector<Polygon>>> polySplit(4);
    double left = (numeric_limits<double>::max)(),
            down = (numeric_limits<double>::max)(),
            right = -(numeric_limits<double>::max)(),
            up = -(numeric_limits<double>::max)();
    for (int i = 0; i < bigPolygon.size(); ++i)
    {
        if (bigPolygon[i].x < left)
        {
            left = bigPolygon[i].x;
        }
        if (bigPolygon[i].x > right)
        {
            right = bigPolygon[i].x;
        }
        if (bigPolygon[i].y < down)
        {
            down = bigPolygon[i].y;
        }
        if (bigPolygon[i].y > up)
        {
            up = bigPolygon[i].y;
        }
    }
    //外接矩形外扩1米，尽可能避免小正方形和多边形出现顶点重合的情况
    left--;right++;down--;up++;
    int m = ceil((up-down)/side);
    int n = ceil((right-left)/side);
    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            // 生成小正方形
            Polygon square;
            square.push_back(Point(left+j*side, up-i*side));
            square.push_back(Point(left+(j+1)*side, up-i*side));
            square.push_back(Point(left+(j+1)*side, up-(i+1)*side));
            square.push_back(Point(left+j*side, up-(i+1)*side));
            vector<Polygon> intersect = PolygonIntersection(square, bigPolygon);
            if (!intersect.empty())
            {
                int race = i%2*2+j%2;
                polySplit[race].push_back(intersect);
            }
        }
    }
    return polySplit;
}

// a是否在区间[x0, x1)内
bool between(double a, double x0, double x1)
{
    double temp1= a-x0;
    double temp2= a-x1;
    if ( ( temp1<1e-6 && temp2>1e-6 ) || ( temp2<-1e-6 && temp1>-1e-6 ) )
    {
        return true;
    }
    else
    {
        return false;
    }
}

// 判断两条直线段是否有交点，有则计算交点的坐标
// 如果两条直线平行或者重合，认为不相交，交点在一条线段上，也认为不相交
// p1,p2是直线一的端点坐标
// p3,p4是直线二的端点坐标
bool detectIntersect(const Point& p1, const Point& p2, const Point& p3, const Point& p4, double& line_x, double& line_y)
{
    //double line_x,line_y; //交点
    if ( (fabs(p1.x-p2.x)<1e-6) && (fabs(p3.x-p4.x)<1e-6) )
    {
        return false;
    }
    else if ( (fabs(p1.x-p2.x)<1e-6) ) //如果直线段p1p2垂直与y轴
    {
        if (between(p1.x,p3.x,p4.x))
        {
            double k = (p4.y-p3.y)/(p4.x-p3.x);
            line_x = p1.x;
            line_y = k*(line_x-p3.x)+p3.y;

            if (between(line_y,p1.y,p2.y))
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        else
        {
            return false;
        }
    }
    else if ( (fabs(p3.x-p4.x)<1e-6) ) //如果直线段p3p4垂直与y轴
    {
        if (between(p3.x,p1.x,p2.x))
        {
            double k = (p2.y-p1.y)/(p2.x-p1.x);
            line_x = p3.x;
            line_y = k*(line_x-p2.x)+p2.y;

            if (between(line_y,p3.y,p4.y))
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        else
        {
            return false;
        }
    }
    else
    {
        double k1 = (p2.y-p1.y)/(p2.x-p1.x);
        double k2 = (p4.y-p3.y)/(p4.x-p3.x);

        if (fabs(k1-k2)<1e-6)
        {
            return false;
        }
        else
        {
            line_x = ((p3.y - p1.y) - (k2*p3.x - k1*p1.x)) / (k1-k2);
            line_y = k1*(line_x-p1.x)+p1.y;
        }

        if (between(line_x,p1.x,p2.x)&&between(line_x,p3.x,p4.x))
        {
            return true;
        }
        else
        {
            return false;
        }
    }
}

//  The function will return YES if the point x,y is inside the polygon, or
//  NO if it is not.  If the point is exactly on the edge of the polygon,
//  then the function return YES.
bool IsPointInPolygon(const Point& pt, Polygon& poly)
{
    if (isPointOnPolygonEdge(pt, poly))
    {
        return true;
    }
    int i,j;
    bool c = false;
    for (i = 0,j = poly.size() - 1;i < poly.size();j = i++)
    {
        if ((((poly[i].y <= pt.y) && (pt.y < poly[j].y)) ||
             ((poly[j].y <= pt.y) && (pt.y < poly[i].y)))
            && (pt.x < (poly[j].x - poly[i].x) * (pt.y - poly[i].y)/(poly[j].y - poly[i].y) + poly[i].x))
        {
            c = !c;
        }
    }
    return c;
}


bool PointCmp(Point* poly1, Point* poly2)
{
    if (poly1->x > poly2->x)
    {
        return true;
    }
    if (poly1->x == poly2->x)
    {
        return poly1->y > poly2->y;
    }
    return false;
}

//将交点和顶点组合成一个环
Node* crossWithPoly(const vector<vector<Point*>>& crossPoints, Polygon& poly1, int poly2Size, bool transpose)
{
    Node* crossWithPoly1 = NULL;
    Node* pNode = crossWithPoly1;
    for (int i = 0; i < poly1.size(); ++i)
    {
        vector<Point*> ptInLine;
        for (int j = 0; j < poly2Size; ++j)
        {
            Point* curCross= transpose ? crossPoints[j][i] : crossPoints[i][j];
            if (curCross != NULL)
            {
                ptInLine.push_back(curCross);
            }
        }
        sort(ptInLine.begin(), ptInLine.end(), PointCmp);
        if (!PointCmp(&poly1[i], &poly1[(i+1)%poly1.size()]))
        {
            reverse(ptInLine.begin(), ptInLine.end	());
        }
        if (crossWithPoly1 == NULL)
        {
            crossWithPoly1 = new Node(&poly1[i], vetex, false, NULL);
            pNode = crossWithPoly1;
        }
        else
        {
            pNode->next = new Node(&poly1[i], vetex, false, NULL);
            pNode = pNode->next;
        }
        for (int j = 0; j < ptInLine.size(); ++j)
        {
            pNode->next = new Node(ptInLine[j], cross, false, NULL);
            pNode = pNode->next;
        }
    }
    pNode->next = crossWithPoly1;
    return crossWithPoly1;
}

// 求两个简单多边形的交集多边形，支持交集是多个多边形的情况，不支持多边形的交点为顶点或交点落在多边形边上的情况
// poly用形参数传递，不改变原来的poly
vector<Polygon> PolygonIntersection(Polygon poly1, Polygon poly2)
{
    vector<Polygon> polyInter;
    if (ClockWise(poly1) < 0)
    {
        reverse(poly1.begin(), poly1.end());
    }
    if (ClockWise(poly2) < 0)
    {
        reverse(poly2.begin(), poly2.end());
    }
    vector<vector<Point*>> crossPoints(poly1.size(), vector<Point*>(poly2.size(), NULL));
    int crossNum = 0;
    for (int i = 0; i < poly1.size(); ++i)
    {
        for (int j = 0; j < poly2.size(); ++j)
        {
            double x, y;
            if (detectIntersect(poly1[i], poly1[(i+1)%poly1.size()], poly2[j], poly2[(j+1)%poly2.size()], x, y))
            {
                crossPoints[i][j] = new Point(x, y);
                crossNum++;
            }
        }
    }
//考虑一个多边形包含另一个多边形的情况
    if (crossNum == 0)
    {
        int inPolyNum = 0;
        for (int i = 0; i < poly1.size(); ++i)
        {
            if (IsPointInPolygon(poly1[i], poly2))
            {
                inPolyNum++;
            }
        }
        if (inPolyNum > 0)
        {
            if (inPolyNum != poly1.size())
            {
                cout << "impossible!" << endl;
                return polyInter;
            }
            polyInter.push_back(poly1);
            return polyInter;
        }
        inPolyNum = 0;
        for (int i = 0; i < poly2.size(); ++i)
        {
            if (IsPointInPolygon(poly2[i], poly1))
            {
                inPolyNum++;
            }
        }
        if (inPolyNum > 0)
        {
            if (inPolyNum != poly2.size())
            {
                cout << "impossible!" << endl;
                return polyInter;
            }
            polyInter.push_back(poly2);
            return polyInter;
        }
        return polyInter;
    }
// 考虑有交点的情况
    Node* crossWithPoly1 = crossWithPoly(crossPoints, poly1, poly2.size(), false);
    Node* crossWithPoly2 = crossWithPoly(crossPoints, poly2, poly1.size(), true);

    int traver1 = 0;
//对于简单多边形，poly1遍历完了，所有交集部分就都出来了
    while (traver1 < poly1.size()+crossNum)
    {
// 找出交集多边形的第一个点
        Node* poly1First = crossWithPoly1;
        while(true)
        {
            if (poly1First->isTraverse == false)
            {
                poly1First->isTraverse = true;
                traver1++;
                if (poly1First->ptType == cross ||
                    (poly1First->ptType == vetex && IsPointInPolygon(*poly1First->pt, poly2)))
                {
                    break;
                }
            }
            poly1First = poly1First->next;
            if (traver1 >= poly1.size() + crossNum)
            {
                return polyInter;
            }
        }
// 交替遍历crossWithPoly1和crossWithPoly2
        Node* cur = poly1First;
        Polygon intersection;
        bool traverseSwitch = true;
        do
        {
            intersection.push_back(*(cur->pt));
            cur->isTraverse = true;
            if (traverseSwitch)
            {
                traver1++;
            }
            if (cur->ptType == cross && cur != poly1First)
            {
                Node* otherLinkStart = (traverseSwitch?crossWithPoly2:crossWithPoly1);
                while(otherLinkStart->pt != cur->pt)
                {
                    otherLinkStart = otherLinkStart->next;
                }
                cur = otherLinkStart;
                traverseSwitch = !traverseSwitch;
                cur->isTraverse = true;
                if (traverseSwitch)
                {
                    traver1++;
                }
            }
            cur = cur->next;
        }while (cur->pt != poly1First->pt);

        polyInter.push_back(intersection);
    }

    return polyInter;
}

// 判断一个多边形点集是否是顺时针
/*
   Return the clockwise status of a curve, clockwise or counterclockwise
   n vertices making up curve p
   return 0 for incomputables eg: colinear points
          CLOCKWISE == 1
          COUNTERCLOCKWISE == -1
   It is assumed that
   - the polygon is closed
   - the last point is not repeated.
   - the polygon is simple (does not intersect itself or have holes)
*/
int ClockWise(const Polygon& p)
{
    int n = p.size();
    const int colckWise = 1;
    const int counterClockWise = -1;
    int i,j,k;
    int count = 0;
    double z;

    if (n < 3)
        return(0);

    for (i=0;i<n;i++) {
        j = (i + 1) % n;
        k = (i + 2) % n;
        z  = (p[j].x - p[i].x) * (p[k].y - p[j].y);
        z -= (p[j].y - p[i].y) * (p[k].x - p[j].x);
        if (z < 0)
            count--;
        else if (z > 0)
            count++;
    }
    if (count > 0)
        return(counterClockWise);
    else if (count < 0)
        return(colckWise);
    else
        return(0);
}

// 点是否落在线段内，若点落在线段的端点上也返回true
bool IsPointOnSegment(const Point& p, const Point& p1, const Point& p2)
{
    if ((p1.x-p.x)*(p2.y-p.y)-(p2.x-p.x)*(p1.y-p.y)!=0)
    {
        return false;
    }
    if ((p.x>p1.x&&p.x>p2.x) || (p.x<p1.x &p.x<p2.x))
    {
        return false;
    }
    if ((p.y>p1.y&&p.y>p2.y) || (p.y<p1.y &p.y<p2.y))
    {
        return false;
    }
    return true;
}
bool isPointOnPolygonEdge(const Point& pt, const Polygon& poly)
{
    for (int i = 0; i < poly.size(); ++i)
    {
        if (IsPointOnSegment(pt, poly[i], poly[(i+1)%poly.size()]))
        {
            return true;
        }
    }
    return false;
}