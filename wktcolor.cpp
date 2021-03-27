/*
   wktcolor.cpp

   Copyright (c) 2021 by Daniel Kelley

*/

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>

#include <boost/type_index.hpp>
#include <boost/optional.hpp>
#include <boost/program_options.hpp>

#define USE_BOOST

#ifdef USE_BOOST
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/polygon.hpp>

typedef boost::geometry::model::d2::point_xy<double> point_type;
#endif

#ifdef USE_CGAL
#include <CGAL/Simple_cartesian.h> // analyze OK (3:14)
#include <CGAL/IO/WKT.h> // analyze OK (3:14)
#include <CGAL/Exact_predicates_exact_constructions_kernel.h> // analyze ??
#endif

using boost::typeindex::type_id_with_cvr;

#ifdef USE_CGAL
typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

typedef CGAL::Point_2<Kernel> Point;
typedef std::vector<Point> MultiPoint;

typedef std::vector<Point> LineString;
typedef std::vector<LineString> MultiLineString;

typedef CGAL::Polygon_with_holes_2<Kernel> Polygon;
typedef std::vector<Polygon> MultiPolygon;
#endif

class WktColor {
public:
    int colors = 4;
    void run(std::string file);
private:
#ifdef USE_CGAL
    MultiPoint points;
    MultiLineString polylines;
    MultiPolygon polygons;
#endif
    void color();
    void read(std::string file);
};

void WktColor::run(std::string file)
{
    read(file);
    color();
}

void WktColor::color()
{
}

void WktColor::read(std::string file)
{
    std::ifstream input(file);
    std::cout << "colors: " << colors << std::endl;
    std::cout << "file: " << file << std::endl;
#if 0
    // Note: doesn't handle WKT GEOMETRYCOLLECTION - gets stuck in an
    // infinite loop.
    CGAL::read_WKT(input,points,polylines,polygons);

    if (points.size() > 0) {
        std::cout << "points" << std::endl;
    }
    for(auto p : points) {
        std::cout<<p<<std::endl;
    }

    if (polylines.size() > 0) {
        std::cout << "lines" << std::endl;
    }
    for(auto ls : polylines) {
        for(auto p : ls) {
            std::cout<<p<<std::endl;
        }
    }

    if (polygons.size() > 0) {
        std::cout << "polygons" << std::endl;
    }
    for(auto p : polygons) {
        std::cout<<p<<std::endl;
    }
#endif
}

static void usage()
{
    std::cout
        << "usage: wktcolors [-cN] [-h] file.wkt"
        << std::endl;
}

int main(int argc, char *argv[])
{
    int rc = EXIT_FAILURE;
    WktColor wktcolor;

    try {
        boost::optional<int> colors;
        boost::program_options::options_description
            options("Options");
        options.add_options()
            ("colors,c", boost::program_options::value(&colors), "colors")
            ("help,h", "help");

        boost::program_options::options_description desc;
        desc.add(options);

        boost::program_options::variables_map cli;
        auto parser = boost::program_options::command_line_parser(argc, argv);
        auto args = parser.options(desc).run();
        boost::program_options::store(args, cli);
        boost::program_options::notify(cli);

        if (cli.count("help")) {
            usage();
            std::cout << desc << std::endl;
        }

        if (colors) {
            wktcolor.colors = *colors;
        }

        auto v = args.options;
#if 0
        std::cout
            << type_id_with_cvr<decltype(v)>().pretty_name()
            << std::endl;
#endif
#if 0
        for_each( v.begin(), v.end(), [] (auto val) {
            // boost::program_options::basic_option<char>
            std::cout
                << val.position_key
                << ":"
                << val.value[0]
                << std::endl;
        });
#endif
        if (v.size() == 1) {
            wktcolor.run(v[0].value[0]);
            rc = EXIT_SUCCESS;
        } else {
            usage();
            rc = EXIT_FAILURE;
        }

    }

    catch(std::exception& e)
    {
        std::cout << e.what() << std::endl;
        rc = EXIT_FAILURE;
    }

    return rc;
}
