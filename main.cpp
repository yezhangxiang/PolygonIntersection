#include <iostream>
#include <fstream>
#include "PolygonIntersection.h"

int main()
{
    ofstream out("out.txt");
    ifstream in("in.txt");
    if (!in.is_open() || !out.is_open())
    {
        cout << "Error opening file";
        exit(1);
    }
    Polygon p1, p2, p_inter;
    string buffer;
    int p1Size, p2Size;
    in >> p1Size;
    for (int i = 0; i < p1Size; ++i)
    {
        double tmpX, tmpY;
        in >> tmpX >> tmpY;
        p1.push_back(Point(tmpX, tmpY));
    }
    in >> p2Size;
    for (int i = 0; i < p2Size; ++i)
    {
        double tmpX, tmpY;
        in >> tmpX >> tmpY;
        p2.push_back(Point(tmpX, tmpY));
    }
    cout << ClockWise(p1) << endl;
    cout << ClockWise(p2) << endl;
    vector<Polygon> polyInter = PolygonIntersection(p1, p2);
    for (int i = 0; i < polyInter.size(); ++i)
    {
        for (int j = 0; j < polyInter[i].size(); ++j)
        {
            out << polyInter[i][j].x << "\t" << polyInter[i][j].y << endl;
        }
        out << endl;
    }
    vector<vector<vector<Polygon>>> polySplit = PolygonSplit(p1, 3000);
    int polyIndex = 1;
    for (auto it1 = polySplit.begin(); it1 != polySplit.end(); ++it1)
    {
        for (auto it2 = it1->begin(); it2 != it1->end(); ++it2)
        {
            for (auto it3 = it2->begin(); it3 != it2->end(); ++it3)
            {
                out << "poly(" << polyIndex << ").points = [";
                for (auto it4 = it3->begin(); it4 != it3->end(); ++it4)
                {
                    out << it4->x << "\t" << it4->y << endl;
                }
                polyIndex++;
                out << "];" << endl;
            }
            out << endl;
        }
        out << endl;
    }
    out << endl;
}
