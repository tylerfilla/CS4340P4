/*
 * Tyler Filla
 * CS 4340 - Project 4
 * November 29, 2018
 */

#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

// https://github.com/lava/matplotlib-cpp
#define WITHOUT_NUMPY
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

/**
 * A single data point.
 *
 * Horizontal and vertical components are labeled h and v, respectively, to
 * avoid confusion when reading.
 */
struct Point
{
    double h;
    double v;
};

/**
 * A line in slope-intercept form.
 */
struct Line
{
    double slope;
    double yint;
};

/**
 * A single data set:
 *   D = { x1, x2 }
 */
struct D_t
{
    Point x1;
    Point x2;
};

/**
 * Compute the line between two points.
 *
 * @param a The first point
 * @param b The second point
 * @return The computed line
 */
static Line compute_line(Point a, Point b)
{
    // Find slope
    auto slope = (b.v - a.v) / (b.h - a.h);

    // Find y-intercept
    // y - y1 = m(x - x1)
    // y - y1 = mx - m(x1)
    // y = mx - m(x1) - y1
    // y = mx + (-m(x1) - y1)
    // y = mx + b where b = -m(x1) + y1
    auto yint = -slope * a.h + a.v;

    return {slope, yint};
}

int main(int argc, char* argv[])
{
    // Create a simple uniform [-1, 1] RNG with a fixed seed and double precision
    std::uniform_real_distribution rng {-1.0, 1.0 + std::numeric_limits<double>::epsilon()};
    std::minstd_rand0 random {123};

    // An algorithm to generate a point
    // The point is of the form (h, h^2) where h an instance of random variable x
    auto gen_point = [&]() -> Point
    {
        auto h = rng(random);
        auto v = std::pow(h, 2);

        return {h, v};
    };

    // An algorithm to generate a data set D
    // The set is of the form {x_1, x_2} where x_i is a point generated above
    auto gen_D = [&]() -> D_t
    {
        auto x1 = gen_point();
        auto x2 = gen_point();

        return {x1, x2};
    };

    // A collection of all generated hypotheses
    // This is not the hypothesis set H (it is a subset, though)
    std::vector<Line> h;

    // The g-bar
    // This is the "average hypothesis" that may not actually be in H
    Line gbar {};

    // Generate five hundred data sets
    // For each data set, make a hypothesis
    for (int i = 0; i < 500; ++i)
    {
        // Generate data set
        auto D = gen_D();

        std::cout << "D = {(" << D.x1.h << ", " << D.x1.v << "), (" << D.x2.h << ", " << D.x2.v << ")}\n";

        // Compute the line between both points in D
        // This *IS* the hypothesis g(x) for D
        auto g = compute_line(D.x1, D.x2);
        h.push_back(g);

        // Plot the generated hypothesis in grey
        auto gx1 = -1.05;
        auto gy1 = g.slope * gx1 + g.yint;
        auto gx2 = 1.05;
        auto gy2 = g.slope * gx2 + g.yint;
        plt::plot(std::vector {gx1, gx2}, std::vector {gy1, gy2}, "gray");

        std::cout << "g(x) = " << g.slope << "x + " << g.yint << "\n";

        // Accumulate slopes and y-intercepts
        // We'll divide by the number of g's later to get g-bar
        gbar.slope += g.slope;
        gbar.yint += g.yint;
    }

    // Finish computing g-bar
    // h.size() is the number of hypotheses generated
    gbar.slope /= h.size();
    gbar.yint /= h.size();

    std::cout << "gbar(x) = " << gbar.slope << "x + " << gbar.yint << "\n";

    // Set up results plot bounds
    plt::xlim(-1.05, 1.05);
    plt::ylim(-1.05, 1.05);

    // Plot y = x^2 in black
    std::vector<double> yx2x {};
    std::vector<double> yx2y {};
    double lx = -10;
    double ly = 10;
    for (double x = -1.1; x <= 1.1; x += 0.01)
    {
        auto y = std::pow(x, 2);

        yx2x.push_back(lx);
        yx2y.push_back(ly);
        yx2x.push_back(x);
        yx2y.push_back(y);

        lx = x;
        ly = y;
    }
    plt::named_plot("f(x)", yx2x, yx2y, "black");

    // Plot g-bar over [-1.05, 1.05] in blue
    auto gbarx1 = -1.05;
    auto gbary1 = gbar.slope * gbarx1 + gbar.yint;
    auto gbarx2 = 1.05;
    auto gbary2 = gbar.slope * gbarx2 + gbar.yint;
    plt::named_plot("gbar(x)", std::vector {gbarx1, gbarx2}, std::vector {gbary1, gbary2}, "blue");

    // Show the plot
    plt::legend();
    plt::show();

    return 0;
}
