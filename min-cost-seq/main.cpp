//
//  main.cpp
//  min-cost-seq
//
//  Created by mndx on 10/06/2022.
//

#include <iostream>
#include <map>
#include <math.h>
#include <random>
#include <time.h>

typedef struct point {
    double x;
    double y;
    int index;
} point_t;

typedef struct optimal_solution {
    std::vector<bool> s_vec;
    double min_cost;
    bool is_set;
} opt_data;

void merge(std::vector<point_t> & arr, int i, int k, int j) {

    int n1 = k - i + 1;
    int n2 = j - k;
    point_t * arr1 = new point_t[n1 + 1];
    point_t * arr2 = new point_t[n2 + 1];
    int arr_index = i;
    for(int arr1_index = 0; arr1_index < n1; ++arr1_index) {
        arr1[arr1_index].x = arr[arr_index].x;
        arr1[arr1_index].y = arr[arr_index].y;
        arr1[arr1_index].index = arr[arr_index].index;
        arr_index++;
    }

    for(int arr2_index = 0; arr2_index < n2; ++arr2_index) {
        arr2[arr2_index].x = arr[arr_index].x;
        arr2[arr2_index].y = arr[arr_index].y;
        arr2[arr2_index].index = arr[arr_index].index;
        arr_index++;
    }

    arr1[n1].x = 3e8;
    arr2[n2].x = 3e8;

    int arr1_index = 0;
    int arr2_index = 0;
    for(int arr_index = i; arr_index <= j; ++arr_index) {
        if(arr1[arr1_index].x >= arr2[arr2_index].x) {
            arr[arr_index].x = arr2[arr2_index].x;
            arr[arr_index].y = arr2[arr2_index].y;
            arr[arr_index].index = arr2[arr2_index].index;
            arr2_index++;
        }
        else {
            arr[arr_index].x = arr1[arr1_index].x;
            arr[arr_index].y = arr1[arr1_index].y;
            arr[arr_index].index = arr1[arr1_index].index;
            arr1_index++;
        }
    }

    //Free data
    delete [] arr1;
    delete [] arr2;
}

void merge_sort(std::vector<point_t> & arr, int i, int j) {

    if(j - i > 0) {
        int k = (i + j) / 2;
        merge_sort(arr, i, k);
        merge_sort(arr, k + 1, j);
        merge(arr, i, k, j);
    }
}

void merge_sort_wrap(std::vector<point_t> & input_arr) {

    int n = (int) input_arr.size();

    //Perform sort
    merge_sort(input_arr, 0, n - 1);
}

double ** mat2D(int n) {
    
    double ** mat = new double * [n];
    
    for(int i = 0; i < n; ++i) {
        mat[i] = new double[n];
        
        for(int j = 0; j < n; ++j) {
            mat[i][j] = 0.0;
        }
    }
    
    return mat;
}

void free_mat2D(double ** mat, int n) {
    
    for(int i = 0; i < n; ++i) {
        delete [] mat[i];
    }
    
    delete [] mat;
}

opt_data min(opt_data o1, opt_data o2) {
    opt_data res = {};
    res.min_cost = 0.0;
    
    if(o1.min_cost < o2.min_cost) res = o1;
    else res = o2;
    
    return res;
}

double min_d(double x, double y) {
    double res = 0.0;
    
    if(x == 0.0 && y > 0) {
        return y;
    }
    
    if(x > 0 && y == 0.0) {
        return x;
    }
    
    if(x < y) res = x;
    else res = y;
    
    return res;
}

double len(point_t p1, point_t p2) {
    
    double dx = (p1.x - p2.x);
    double dy = (p1.y - p2.y);
    
    return sqrt(dx * dx + dy * dy);
}

double len_p(std::vector<bool> v, point_t * p, int n) {

    std::vector<point_t> p1;
    std::vector<point_t> p2;
    
    for(int i = 0; i < n; ++i) {
        if(v[i] == true) {
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

double len_v(std::vector<point_t> p) {

    int n = (int) p.size();
    double sum = 0.0;

    for(int i = 1; i < n; ++i) {
        sum += len(p[i], p[i - 1]);
    }

    return sum;
}

void init_p(int n, point_t * p) {
    
    double min_x = 0;
    double max_x = 10.5;
    double min_b = 0.0;
    double min_y = -50.0;
    double max_y = 50.0;
    
    srand((unsigned) time(NULL));

    for(int i = 0; i < n; ++i) {
        double f1 = (double) rand() / RAND_MAX;
        double x = min_b + f1 * (max_x - min_x);
        double f2 = (double) rand() / RAND_MAX;
        double y = min_y + f2 * (max_y - min_y);

        p[i].x = x;
        p[i].y = y;
        p[i].index = i + 1;

        min_b = x;
    }
}

void print_partition(std::vector<bool> v) {

    int n = (int) v.size();
    
    if(v[n - 1] == false) {
        for(int i = 0; i < n; ++i) {
            v[i] = !v[i];
        }
    }

    for(int i = 0; i < n; ++i) {
        if(v[i] == false) {
            std::cout << "element " << (i + 1) << " is in P" << (1) << std::endl;
        }
        if(v[i] == true) {
            std::cout << "element " << (i + 1) << " is in P" << (2) << std::endl;
        }
    }
}

void print_partition_ref(std::vector<bool> v, point_t * p, int n) {

    if(v[n - 1] == false) {
        for(int i = 0; i < n; ++i) {
            v[i] = !v[i];
        }
    }

    for(int i = 0; i < n; ++i) {
        if(v[i] == false) {
            std::cout << "P1, x: " << p[i].x << ", y: " << p[i].y << std::endl;
        }
    }
    for(int i = 0; i < n; ++i) {
        if(v[i] == true) {
            std::cout << "P2, x: " << p[i].x << ", y: " << p[i].y << std::endl;
        }
    }
}

opt_data min_cost_rec(point_t * p, int n, std::vector<bool> v) {
    opt_data res = {};
    res.min_cost = 0.0;
    
    // Get size current data vector
    int s = (int) v.size();
    
    // Compute minimum length cost
    if(s < n) {
        // Pick element and store it in P1
        std::vector<bool> v1 = v;
        v1.push_back(true);
        opt_data val1 = min_cost_rec(p, n, v1);
        
        // Pick element and store it in P2
        std::vector<bool> v2 = v;
        v2.push_back(false);
        opt_data val2 = min_cost_rec(p, n, v2);
        res = min(val1, val2);
    }
    
    // Compute length cost
    if(s == n) {
        res.min_cost = len_p(v, p, n);
        res.s_vec = v;
        
        return res;
    }
    
    return res;
}

opt_data min_cost_seq(point_t * p, int n) {
    
    std::vector<bool> v;
    
    return min_cost_rec(p, n, v);
}

int ops = 0;

void make_sol(double ** C, int n, point_t * p, int lf, int j, bool flipped, std::vector<point_t> & P1, std::vector<point_t> & P2) {

    int cut = j;

    ops++;

    for(int i = j; i > lf; --i) {
        if(!flipped) {
            P2.push_back(p[i]);
        }
        if(flipped) {
            P1.push_back(p[i]);
        }
    }

    double bounds = 3e8;

    int jj = lf + 1;
    for(int i = 1; i < jj - 1; ++i) {
        double val = C[i][jj - 1] + len(p[i], p[jj]);
        if(val < bounds) {
            bounds = val;
            cut = i;
        }
    }

    if(bounds > C[0][jj - 1]) {
        bounds = C[0][jj -1];
        cut = 0;
    }

    if(cut == 0 && !flipped) {

        for(int i = 1; i <= jj - 1; ++i) {
            P1.push_back(p[i]);
        }

        return;
    }

    if(cut == 0 && flipped) {

        for(int i = 1; i <= jj - 1; ++i) {
            P2.push_back(p[i]);
        }

        return;
    }

    bool is_flipped = !flipped;
    make_sol(C, n, p, cut, jj - 1, is_flipped, P1, P2);
}

void make_sol_wrap(double ** C, int n, point_t * p, std::vector<point_t> & P1, std::vector<point_t> & P2) {

    double bounds = 3e8;
    int cut = n - 1;

    for(int i = 0; i <= n - 1; ++i) {
        if(bounds > C[i][n]) {
            bounds = C[i][n];
            cut = i;
        }
    }

    make_sol(C, n, p, cut, n, false, P1, P2);
}

void get_partition(double ** C, int n, point_t * p, std::vector<point_t> & P1, std::vector<point_t> & P2) {

    point_t * p_ref = new point_t[n + 1];

    for(int i = 0; i < n + 1; ++i) {
        if(i == 0) {
            p_ref[i].x = p[0].x;
            p_ref[i].y = p[0].y;
            p_ref[i].index = p[0].index;
        }
        if(i > 0) {
            p_ref[i].x = p[i - 1].x;
            p_ref[i].y = p[i - 1].y;
            p_ref[i].index = p[i - 1].index;
        }
    }

    make_sol_wrap(C, n, p_ref, P1, P2);

    merge_sort_wrap(P1);
    merge_sort_wrap(P2);

    delete [] p_ref;
}

double min_cost_bottom_up(point_t * p, int n, std::vector<point_t> & P1, std::vector<point_t> & P2) {
    double res = 0.0;
    
    // Allocate space for table
    double ** C = mat2D(n + 1);
    
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
        double cost = 0.0;
        for(int i = 1; i < j - 1; ++i) {
            double val = C[i][j - 1] + len(p[i - 1], p[j - 1]);
            cost = min_d(cost, val);
        }
        
        cost = min_d(cost, C[0][j - 1]);
        
        C[j - 1][j] = cost;
        
        // Point was added to P2
        for(int i = j + 1; i <= n; ++i) {
            double sum = C[j - 1][i - 1] + len(p[i - 2], p[i - 1]);
            C[j - 1][i] = sum;
        }
    }
    
    // Compute minimum partition cost
    for(int i = 0; i < n; ++i) {
        res = min_d(res, C[i][n]);
    }

    get_partition(C, n, p, P1, P2);

    // Free space used for table
    free_mat2D(C, n + 1);

    return res;
}

void print_points(point_t * p, int n) {
    
    for(int i = 0; i < n; ++i) {
        std::cout << i << ", x: " << p[i].x << ", y: " << p[i].y << std::endl;
    }
}

void print_partition_elems(std::vector<point_t> P, int id) {

    std::cout << "Printing P" << id << std::endl;

    for(auto e : P) {
        std::cout << "element " << e.index << " is in P" << id << std::endl;
    }
}

int main(int argc, const char * argv[]) {

    // Number of points
    int n = 31;
    
    // Declare and initialize points and solution array
    point_t * p = new point_t[n];
    
    init_p(n, p);
    
    std::vector<point_t> P1;
    std::vector<point_t> P2;

    // Compute minimum cost and partition bottom-up
    double len_bottom_up = min_cost_bottom_up(p, n, P1, P2);
    
    // Verify minimum cost
    double len_verification = len_v(P1) + len_v(P2);

    // Print results
    print_partition_elems(P1, 1);
    print_partition_elems(P2, 2);

    std::cout << "length computed bottom-up: " << len_bottom_up << std::endl;
    std::cout << "length verification: " << len_verification << std::endl;
    
    // Free allocated space
    delete [] p;
    
    return 0;
}
