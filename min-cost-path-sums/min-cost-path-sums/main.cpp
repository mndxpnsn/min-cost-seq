//
//  main.cpp
//  min-cost-path-sums
//
//  Created by Derek Harrison on 11/10/2025.
//

#include <iostream>
#include <math.h>
#include <vector>

using namespace std;

typedef struct point {
    double x;
    double y;
} pnt;

typedef vector<point> vecPoint;
typedef vector<vector<double>> vec2D;
typedef vector<double> vec1D;

bool comparePoint(point pa, point pb);
double len(point pa, point pb);
double lenPath(vecPoint & points, int i, int j);
double minCostBottomUp(vecPoint p);
double minCost(vecPoint & points, int i, int j, vec2D & dp);
double minCostTopDown(vecPoint & points);

int main(int argc, const char * argv[]) {
    
    int n = 16;
    
    vector<point> points;
    
    srand((unsigned int) time(NULL));
    
    for(int i = 0; i < n; i++) {
        point p;
        p.x = rand() % 100;
        p.y = rand() % 100;
        points.push_back(p);
    }
   
    sort(points.begin(), points.end(), comparePoint);
    
    double minCost1 = minCostBottomUp(points);
    double minCost2 = minCostTopDown(points);
    
    cout << "min cost bottom up : " << minCost1 << endl;
    cout << "min cost top down: " << minCost2 << endl;
    
    return EXIT_SUCCESS;
}

double minCostTopDown(vecPoint & points) {
    int n = (int) points.size();
    
    vec2D dp(n, vec1D(n, -1));
    
    double res = INT_MAX;
    
    for(int k = 0; k < n - 1; k++) {
        res = min(res, minCost(points, k, n - 1, dp));
    }
    
    return res;
}

double minCost(vecPoint & points, int i, int j, vec2D & dp) {

    if(j <= 1) {
        return 0;
    }
    
    if(dp[i][j] > -1) {
        return dp[i][j];
    }
    
    if(i < j - 1) {
        return minCost(points, i, i + 1, dp) + lenPath(points, i + 1, j);
    }
    
    double res = INT_MAX;
    
    for(int q = 0; q < i; q++) {
        res = min(res, minCost(points, q, i, dp) + len(points[q], points[j]));
    }
  
    return dp[i][j] = min(res, lenPath(points, 0, i));
}

double minCostBottomUp(vecPoint p) {
    int n = (int) p.size();
    
    double res = INT_MAX;
    
    // Allocate space for table
    vec2D C = vec2D(n + 1, vec1D(n + 1, 0));
    
    // Base case
    C[0][1] = 0.0;
    
    // Partition P1 ends at 0
    for(int i = 2; i <= n; ++i) {
        double sum = 0.0;
        for(int j = 1; j < i; ++j) {
            sum += len(p[j - 1], p[j]);
        }
    
        C[0][i] = sum;
    }
    
    for(int j = 2; j <= n; ++j) {
        
        // Partition P1 ends at j - 1 and P2 ends at j
        double cost = INT_MAX;
        for(int i = 1; i < j - 1; ++i) {
            cost = min(cost, C[i][j - 1] + len(p[i - 1], p[j - 1]));
        }
        
        C[j - 1][j] = min(cost, C[0][j - 1]);
        
        // Point was added to P2
        for(int i = j + 1; i <= n; ++i) {
            C[j - 1][i] = C[j - 1][i - 1] + len(p[i - 2], p[i - 1]);
        }
    }
    
    // Compute minimum partition cost
    for(int i = 0; i < n; ++i) {
        res = min(res, C[i][n]);
    }

    return res;
}

bool comparePoint(point pa, point pb) {
    return pa.x < pb.x;
}

double len(point pa, point pb) {
    return sqrt((pa.x - pb.x) * (pa.x - pb.x) + (pa.y - pb.y) * (pa.y - pb.y));
}

double lenPath(vecPoint & points, int i, int j) {
    double res = 0;
    
    for(int k = i; k < j; k++) {
        res += len(points[k], points[k + 1]);
    }
    
    return res;
}
