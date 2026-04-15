#include <vector>
#include <algorithm>
#include <math.h>
#include <functional>
#include <iostream>

#include "linearize.hpp"

bool check_condition(ToleranceCondition cond,
                     double abs_diff, double abs_tol,
                     double rel_diff, double rel_tol,
                     double y_1, double y_2, double y_point)
{
    switch (cond) {
        case ToleranceCondition::Or:
            return abs_diff > abs_tol || rel_diff > rel_tol;
        case ToleranceCondition::And:
            return abs_diff > abs_tol && rel_diff > rel_tol;
        case ToleranceCondition::AbsOnly:
            return abs_diff > abs_tol;
        case ToleranceCondition::RelOnly:
            return rel_diff > rel_tol;
        case ToleranceCondition::Adaptive: {
            double local_scale = std::max({std::abs(y_1), std::abs(y_2), std::abs(y_point)});
            double scaled_tol = rel_tol * local_scale;
            return abs_diff > std::max(abs_tol, scaled_tol);
        }
        case ToleranceCondition::AdaptiveOr: {
            double local_scale = std::max({std::abs(y_1), std::abs(y_2), std::abs(y_point)});
            double scaled_tol = rel_tol * local_scale;
            return abs_diff > abs_tol || abs_diff > scaled_tol;
        }
        case ToleranceCondition::AdaptiveAnd: {
            double local_scale = std::max({std::abs(y_1), std::abs(y_2), std::abs(y_point)});
            double scaled_tol = rel_tol * local_scale;
            return abs_diff > abs_tol && abs_diff > scaled_tol;
        }
        default:
            return false; // fallback (shouldn't happen)
    }
}

void linearize(std::vector<double>& x_points, std::vector<double>& y_points, std::function<double(double)> get_new_y, const double absolute_tolerance, const double relative_tolerance, const ToleranceCondition condition){
    std::vector<std::pair<int, int>> stack;
    stack.reserve(x_points.size());
    for(int i=0; i<x_points.size()-1; i++){
        stack.push_back(std::make_pair(i, i+1));
    }
    while (!stack.empty())
    {
        auto [left, right] = stack.back();
        stack.pop_back();

        double& x_1 = x_points[left];
        double& x_2 = x_points[right];
        double& y_1 = y_points[left];
        double& y_2 = y_points[right];

        double interval = (x_2 - x_1)/2;
        if (interval < 1e-15){
            break;
        }
        
        double x_point = x_1 + interval;
        double y_point = get_new_y(x_point);
        double interp_y = y_1 + (y_2 - y_1)*(x_point - x_1)/(x_2 - x_1);

        double abs_diff = abs(y_point - interp_y);
        double y_magnitude = std::max({std::abs(y_1), std::abs(y_2), std::abs(y_point)});
        double rel_diff = (y_magnitude > 1e-15) ? abs_diff / y_magnitude : 0;


        // if (abs_diff>absolute_tolerance && rel_diff>relative_tolerance){
        if (check_condition(condition, abs_diff, absolute_tolerance, 
                   rel_diff, relative_tolerance, y_1, y_2, y_point)) {
            x_points.push_back(x_point);
            y_points.push_back(y_point);
            stack.push_back(std::make_pair(left, x_points.size()-1));
            stack.push_back(std::make_pair(x_points.size()-1, right));
        }
    }
    std::vector<std::pair<double, double>> data;
    data.resize(x_points.size());
    for (size_t i = 0; i < x_points.size(); i++){
        data[i] = std::make_pair(x_points[i], y_points[i]);
    }
    std::sort(data.begin(), data.end());
    for (size_t i = 0; i<x_points.size(); i++){
        x_points[i] = data[i].first;
        y_points[i] = data[i].second;
    }
}