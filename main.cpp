/*
 * Tyler Filla
 * CS 4340 - Project 4
 * November 29, 2018
 */

#include <cmath>
#include <ctime>
#include <iostream>
#include <iterator>
#include <limits>
#include <random>
#include <vector>

// https://github.com/lava/matplotlib-cpp
#define WITHOUT_NUMPY
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

/**
 * The number of data sets, and also hypotheses, to generate.
 */
constexpr static int NUM_DATA_SETS = 100;

/**
 * The number of test points to use after generating hypotheses.
 */
constexpr static int NUM_TEST_POINTS = 100;

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

template<class InputIter>
static double compute_mean(InputIter begin, InputIter end)
{
    double mean = 0;
    for (auto i = begin; i != end; ++i)
    {
        mean += *i;
    }
    mean /= std::distance(begin, end);

    return mean;
}

template<class InputIter>
static double compute_variance(InputIter begin, InputIter end)
{
    double mean = compute_mean(begin, end);

    double variance = 0;
    for (auto i = begin; i != end; ++i)
    {
        variance += (*i - mean) * (*i - mean);
    }
    variance /= std::distance(begin, end);

    return variance;
}

int main(int argc, char* argv[])
{
    // Create a simple uniform [-1, 1] RNG with a varying seed and double precision
    std::uniform_real_distribution rng {-1.0, 1.0 + std::numeric_limits<double>::epsilon()};
    std::minstd_rand0 random {static_cast<unsigned long>(time(nullptr))};

    // An algorithm to generate a point
    // The point is of the form (h, h^2) where h an instance of random variable x
    auto gen_point = [&]() -> Point
    {
        auto h = rng(random);
        auto v = std::pow(h, 2);

        return {h, v};
    };

    //
    // TRAINING SECTION
    //

    std::cout << "Begin training\n";

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

    // Generate data sets
    // For each data set, make a hypothesis
    for (int i = 0; i < NUM_DATA_SETS; ++i)
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

    std::cout << "Close the plot to continue the program\n";

    // Show the plot of f(x) and gbar(x)
    plt::legend();
    plt::show();

    //
    // TESTING
    //

    std::cout << "Begin testing\n";

    // Generate test points
    std::vector<Point> pts {};
    for (int i = 0; i < NUM_TEST_POINTS; ++i)
    {
        auto pt = gen_point();

        /*
        // Plug point's x-coordinate into all hypotheses
        std::vector<double> c {};
        for (auto&& g : h)
        {
            c.push_back(g.slope * pt.h + g.yint);
        }

        // Find mean of computed results
        auto mean = compute_mean(c.begin(), c.end());

        // Find variance of computed results
        auto var = compute_variance(c.begin(), c.end());
        */

//      std::cout << "(" << pt.h << ", " << pt.v << ") -> bias = " << bias << ", var = " << var << "\n";

        pts.push_back(pt);
    }

    std::cout << "Generated " << pts.size() << " test point(s)\n";

    // MSE terms for calculation of bias
    std::vector<double> bias_mse_terms {};

    // For each test point
    for (auto&& pt : pts)
    {
        // The empirical (against g-bar)
        auto y = gbar.slope * pt.h + gbar.yint;

        // The expected
        auto yhat = std::pow(pt.h, 2);

        // Compute SSE term
        bias_mse_terms.push_back(std::pow(y - yhat, 2));
    }

    // Compute and print bias
    // This is the expected value over all test points
    std::cout << "bias = " << compute_mean(bias_mse_terms.begin(), bias_mse_terms.end()) << "\n";

    // MSE terms for eventual calculation of E[E_out]
    std::vector<double> E_E_out_mse_terms {};

    // For each generated hypothesis
    // This is also, by extension, a loop over the generated data sets
    for (auto&& g : h)
    {
        // SSE terms
        // This is just an intermediate container for processing
        // This doesn't correspond to anything in my printed explanation
        std::vector<double> sse_terms {};

        // For each test point
        for (auto&& pt : pts)
        {
            // The empirical (against the current hypothesis g)
            auto y = g.slope * pt.h + g.yint;

            // The expected
            auto yhat = std::pow(pt.h, 2);

            // The SSE term
            sse_terms.push_back(std::pow(y - yhat, 2));
        }

        // Compute MSE for test point against all generated hypothesis
        // This is E_out for this particular data set D
        E_E_out_mse_terms.push_back(compute_mean(sse_terms.begin(), sse_terms.end()));
    }

    // Compute and print E[E_out]
    // This is the expected value over all hypotheses
    std::cout << "E[E_out] = " << compute_mean(E_E_out_mse_terms.begin(), E_E_out_mse_terms.end()) << "\n";

    return 0;
}
