//
//  main.cpp
//  min-cost-seq
//
//  Created by mndx on 10/06/2022.
//

#include <iostream>
#include <math.h>
#include <random>
#include <time.h>

enum set {
    P1,
    P2
};

typedef struct point {
    double x;
    double y;
} point_t;

typedef struct optimal_solution {
    std::vector<set> s_vec;
    double min_cost;
} opt_data;

opt_data min(opt_data o1, opt_data o2) {
    opt_data res = {};
    res.min_cost = 0.0;
    
    if(o1.min_cost < o2.min_cost) res = o1;
    else res = o2;
    
    return res;
}

double len(point_t p1, point_t p2) {
    
    double dx = (p1.x - p2.x);
    double dy = (p1.y - p2.y);
    
    return sqrt(dx * dx + dy * dy);
}

double len_p(std::vector<set> v, point_t * p, int n) {

    std::vector<point_t> p1;
    std::vector<point_t> p2;
    
    for(int i = 0; i < n; ++i) {
        if(v[i] == P1) {
            p1.push_back(p[i]);
        }
        else {
            p2.push_back(p[i]);
        }
    }
    
    int n1 = (int) p1.size();
    double sum1 = 0.0;
    
    for(int i = 1; i < n1; ++i) {
        sum1 += len(p1[i], p1[i - 1]);
    }
    
    int n2 = (int) p2.size();
    double sum2 = 0.0;
    
    for(int i = 1; i < n2; ++i) {
        sum2 += len(p2[i], p2[i - 1]);
    }
    
    return sum1 + sum2;
}

void init_p(int n, point_t * p) {
    
    double min_x = 0;
    double max_x = 0.5;
    double min_b = 0.0;
    double min_y = 0.0;
    double max_y = 1.0;
    
    srand((unsigned) time(NULL));
    
    for(int i = 0; i < n; ++i) {
        double f1 = (double) rand() / RAND_MAX;
        double x = min_b + f1 * (max_x - min_x);
        double f2 = (double) rand() / RAND_MAX;
        double y = min_y + f2 * (max_y - min_y);
        
        p[i].x = x;
        p[i].y = y;
        
        min_b += x;
    }
}

opt_data ** init_dp(int n) {
    
    opt_data ** dp = new opt_data * [n];
    
    for(int i = 0; i < n; ++i) {
        dp[i] = new opt_data[2];
        dp[i][0].min_cost = 0.0;
        dp[i][1].min_cost = 0.0;
    }
    
    return dp;
}

void free_dp(opt_data ** dp, int n) {
    
    for(int i = 0; i < n; ++i) {
        delete [] dp[i];
    }
    
    delete [] dp;
}

opt_data min_cost_rec(point_t * p, int n, bool in_p1, int d, std::vector<set> v, opt_data ** dp) {
    opt_data res = {};
    res.min_cost = 0.0;
    
    // Get size current data vector
    int s = (int) v.size();
    
    // Get data from memo table if available
    if(dp[d][in_p1].min_cost != 0.0) {
        return dp[d][in_p1];
    }
    
    // Compute minimum length cost
    if(s < n) {
        // Pick element d and store it in P1
        std::vector<set> v1 = v;
        v1.push_back(P1);
        opt_data val1 = min_cost_rec(p, n, true, d + 1, v1, dp);
        
        // Pick element d and store it in P2
        std::vector<set> v2 = v;
        v2.push_back(P2);
        opt_data val2 = min_cost_rec(p, n, false, d + 1, v2, dp);
        res = min(val1, val2);
    }
    
    // Compute length cost
    if(s == n) {
        res.min_cost = len_p(v, p, n);
        res.s_vec = v;
        return res;
    }
    
    // Store data in memo table
    dp[d][in_p1] = res;
    
    return res;
}

opt_data min_cost_seq(point_t * p, int n) {
    
    std::vector<set> v;
    
    // Allocate space for memo table
    opt_data ** dp = init_dp(n + 1);
    
    // Compute minimum cost when first element is stored in P1
    opt_data val1 = min_cost_rec(p, n, true, 0, v, dp);
    
    // Compute minimum cost when first element is stored in P2
    opt_data val2 = min_cost_rec(p, n, false, 0, v, dp);
    
    // Free memo table space
    free_dp(dp, n + 1);
    
    return min(val1, val2);
}

int main(int argc, const char * argv[]) {

    // Number of points
    int n = 10;
    
    // Declare and initialize points and solution array
    point_t * p = new point_t[n];
    
    init_p(n, p);
    
    // Compute minimum cost
    opt_data odata = min_cost_seq(p, n);

    // Verify cost
    double len_ver = len_p(odata.s_vec, p, n);
    
    // Print results
    std::cout << "minimum length cost: " << odata.min_cost << std::endl;
    std::cout << "length verification: " << len_ver << std::endl;
    
    // Free allocated space
    delete [] p;
    
    return 0;
}