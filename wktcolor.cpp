/*
   wktcolor.cpp

   Copyright (c) 2021 by Daniel Kelley

   Color faces of WKT Polygons

   Output is a GraphML file with the following attributes:
     color: unique color index 0-n with perhaps n < 4
     x: polygonal face centriod X
     y: polygonal face centriod Y

*/

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <utility>

#include <boost/type_index.hpp>
#include <boost/optional.hpp>
#include <boost/program_options.hpp>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

#define GEOS_USE_ONLY_R_API
#include <geos_c.h>
#include <wkt.h>

#include <igraph/igraph.h>

#include <ColPack/ColPackHeaders.h>

#include <gsl/gsl-lite.hpp>

using boost::typeindex::type_id_with_cvr;

using MyMesh = OpenMesh::PolyMesh_ArrayKernelT<>;

using vertex = std::pair<double,double>;

struct polyinfo {
    unsigned int points;        // WKT points
    unsigned int faces;         // WKT polygons
};

//
// Polygon colorizer
//
class WktColor {
public:
    void run(std::string file);
    inline void set_verbose() { verbose=1; }
    inline void set_geometry(std::string s) { geometry=s; }
    inline void set_algorithm(std::string s) { algorithm=s; }
    inline void set_ordering(std::string s) { ordering=s; }
    inline void set_adj_file(std::string s) { adj_file=s; }
    inline void set_gml_file(std::string s) { gml_file=s; }
    inline void request_clique_number() { clique=1; }
private:
    MyMesh mesh;
    struct wkt wkt = {};
    struct polyinfo info = {};
    igraph_t contact = {};
    std::map<vertex,size_t> uvertex;
    std::vector<vertex> svertex;
    std::vector<GEOSGeometry *> centroid;
    int verbose = 0;
    int clique = 0;
    std::string geometry = {};
    std::string algorithm = "igraph"; // igraph or various ColPack
    std::string ordering = "NATURAL"; // ColPack ordering
    std::string adj_file = {};
    std::string gml_file = {};
    bool adj_tmp = false;
    void init();
    void open();
    void read(std::string file);
    void count();
    void create_mesh();
    void save_mesh_geometry();
    void create_contact_graph();
    void igraph_color();
    void colpack_create_adj_file();
    void colpack_color();
    void color();
    void show_clique_number();
    void close();
};

//
// Initialize colorizing sequence
//
void WktColor::init()
{
    // must be called before other igraph functions
    // according to comment in igraph cattributes.c
    igraph_set_attribute_table(&igraph_cattribute_table);
}

//
// Open given WKT file
//
void WktColor::open()
{
    memset(&wkt, 0, sizeof(wkt));
    wkt.reader = WKT_IO_ASCII;
    auto err = wkt_open(&wkt);
    assert(!err);
}

//
// Read WKT file
//
void WktColor::read(std::string file)
{
    auto err = wkt_read(&wkt, file.c_str());
    assert(!err);
}

//
// Run algorithm over file
//
void WktColor::run(std::string file)
{
    init();
    open();
    read(file);
    count();
    create_mesh();
    if (geometry.size() > 0) {
        save_mesh_geometry();
    }
    create_contact_graph();
    if (clique) {
        show_clique_number();
    } else {
        color();
    }
    close();
}

//
// Count points and polygons (faces)
//
//   Subsequent operations need an idea of how may points and polygons
//   to expect. Make sure geometry is what is expected.  No handling
//   anything weird like non-polygons or polygons with holes.
//
void WktColor::count()
{
    auto n = GEOSGetNumGeometries_r(wkt.handle, wkt.geom);
    for (int i=0; i<n; i++) {
        auto geom = GEOSGetGeometryN_r(wkt.handle, wkt.geom, i);
        assert(geom);

        unsigned int dim = 0;
        unsigned int size = 0;

        // Only Polygons with no holes
        auto g = GEOSGetExteriorRing_r(wkt.handle, geom);
        assert(g != nullptr);
        auto h = GEOSGetNumInteriorRings_r(wkt.handle, geom);
        assert(h == 0);

        auto seq = GEOSGeom_getCoordSeq_r(wkt.handle, g);
        assert(seq != nullptr);

        // Only two dimensions
        auto ok = GEOSCoordSeq_getDimensions_r(wkt.handle, seq, &dim);
        assert(ok);
        assert(dim == 2);

        ok = GEOSCoordSeq_getSize_r(wkt.handle, seq, &size);
        assert(ok);
        assert(size != 0);

        info.faces += 1;
        info.points += size;

    }
}

//
// Create mesh data structure from polygons.
//
//   This leverages the nice half-edge data structure that OpenMesh
//   provides.
//
void WktColor::create_mesh()
{
    std::vector<MyMesh::VertexHandle> vhandle;
    std::vector<MyMesh::VertexHandle> face_vhandles;

    vhandle.resize(info.points);
    int n = GEOSGetNumGeometries_r(wkt.handle, wkt.geom);
    for (int i=0; i<n; ++i) {
        unsigned int size = 0;

        face_vhandles.clear();
        auto poly = GEOSGetGeometryN_r(wkt.handle, wkt.geom, i);
        assert(poly != nullptr);
        auto ring = GEOSGetExteriorRing_r(wkt.handle, poly);
        assert(ring != nullptr);
        auto seq = GEOSGeom_getCoordSeq_r(wkt.handle, ring);
        assert(seq != nullptr);
        auto ok = GEOSCoordSeq_getSize_r(wkt.handle, seq, &size);
        assert(ok);
        centroid.push_back(GEOSGetCentroid_r(wkt.handle, ring));

        if (verbose) {
            std::cout << "poly\n";
        }

        // GOES polygons are closed (last point == first point be definition)
        // so don't bother looking at the last point in the polygon.
        for (unsigned int j = 0; j < size-1; ++j) {
            double x = 0.0;
            double y = 0.0;
            ok = GEOSCoordSeq_getXY_r(wkt.handle, seq, j, &x, &y);
            assert(ok);
            if (verbose) {
                std::cout << "  (" << x << "," << y << ")";
            }
            vertex v(x, y);
            // idx 0 means "not found" so uvertex values are idx+1
            // other vertex vectors are zero based.
            auto idx = uvertex[v];
            if (!idx) {

                assert(svertex.size() < info.points);
                idx = svertex.size() + 1;
                uvertex[v] = idx;
                svertex.push_back(v);
                if (verbose) {
                    std::cout << "+";
                }
                vhandle[idx-1] = mesh.add_vertex(MyMesh::Point(x, y, 0));
            } else if (verbose) {
                std::cout << "@";
            }
            face_vhandles.push_back(vhandle[idx-1]);
            if (verbose) {
                std::cout << idx << "\n";
            }
        }
        mesh.add_face(face_vhandles);
    }

}

//
// Create a contact graph from the mesh. Each face is a vertex in the
// graph, and there is an edge between faces that share a common edge.
//
void WktColor::create_contact_graph()
{
    igraph_vector_int_t edge;
    auto edge_limit = info.points * 2;
    auto edge_id =  info.points * 0; // * 0 to derive type
    // edge_limit is an over-estimate
    igraph_vector_int_init(&edge, edge_limit);

    for (auto face = mesh.faces_begin(); face != mesh.faces_end(); ++face) {
        auto fh = *face;
        for (auto dart = mesh.cfh_iter(fh); dart.is_valid(); ++dart) {
            auto heh = *dart;
            auto ofh = mesh.opposite_face_handle(heh);
            if (verbose) {
                std::cout
                << "dart("
                << type_id_with_cvr<decltype(fh)>().pretty_name()
                << "("
                << fh
                << "),"
                << type_id_with_cvr<decltype(heh)>().pretty_name()
                << "("
                << heh
                << ")->"
                << type_id_with_cvr<decltype(ofh)>().pretty_name()
                << "("
                << ofh
                << ")"
                << ")\n";
            }
            if (ofh.is_valid()) {
                // If the opposite half edge face is valid,
                // i.e. opposite half edge is not a boundary, then
                // create and edge between the two faces in the
                // contact graph.
                assert(edge_id < edge_limit);
                VECTOR(edge)[edge_id] = fh.idx();
                ++edge_id;
                assert(edge_id < edge_limit);
                VECTOR(edge)[edge_id] = ofh.idx();
                ++edge_id;
            }
        }
    }

    igraph_create(&contact, &edge, 0, IGRAPH_UNDIRECTED);
    // get rid of self loops and duplicate edges
    igraph_simplify(&contact, true, true, nullptr);

    // Add centroid attributes
    auto faces = igraph_vcount(&contact);
    for (int i=0; i<faces; i++) {
        auto g = centroid[i];
        double x = 0.0;
        double y = 0.0;

        auto ok = GEOSGeomGetX_r(wkt.handle, g, &x);
        assert(ok);
        ok = GEOSGeomGetY_r(wkt.handle, g, &y);
        assert(ok);
        auto err = igraph_cattribute_VAN_set(&contact, "x", i, x);
        assert(!err);
        err = igraph_cattribute_VAN_set(&contact, "y", i, y);
        assert(!err);
    }

    igraph_vector_int_destroy(&edge);
}

//
// Colorize the contact graph with igraph.
//
void WktColor::igraph_color()
{
    igraph_vector_int_t cv;
    auto faces = igraph_vcount(&contact);
    igraph_vector_int_init(&cv, 0);
    auto err = igraph_vertex_coloring_greedy(
        &contact,
        &cv,
        IGRAPH_COLORING_GREEDY_COLORED_NEIGHBORS);
    assert(!err);

    for (int i=0; i<faces; i++) {
        auto color = VECTOR(cv)[i];
        auto ncolor = gsl::narrow_cast<double>(color);
        if (verbose) {
            std::cout
                << "color( "
                << i
                << ","
                << color
                << ")"
                << std::endl;
        }

        err = igraph_cattribute_VAN_set(&contact, "color", i, ncolor);
        assert(!err);

    }

    igraph_vector_int_destroy(&cv);
}

//
// The ColPack interface uses an adjacency matrix file format
// as input. Only the row and column of each edge should be present.
//

void WktColor::colpack_create_adj_file()
{
    if (adj_file.size() == 0) {
        // Use a temporary file if no specific output was requested.
        // We just care about a unique name, so close the returned
        // file handle.
        std::string file("wktcolorXXXXXX");
        ::close(mkstemp(data(file)));
        adj_file = file;
        adj_tmp = true;
    }

    std::ofstream adj(adj_file);
    assert(adj);

    auto faces = igraph_vcount(&contact);
    igraph_eit_t eit;
    auto err = igraph_eit_create(
        &contact,
        igraph_ess_all(IGRAPH_EDGEORDER_ID),
        &eit);
    assert(!err);

    // The MatrixMarket matrix must be square with both forward and
    // reverse edges represented as matrix coordinates. The associated
    // value is ignored.
    int edges = IGRAPH_EIT_SIZE(eit) * 2;

    adj << "%%MatrixMarket matrix coordinate real general\n"
        << faces << " "
        << faces << " "
        << edges << "\n";

    while (!IGRAPH_EIT_END(eit)) {
        igraph_integer_t v1 = 0;
        igraph_integer_t v2 = 0;
        auto eid = IGRAPH_EIT_GET(eit);
        err = igraph_edge(&contact, eid, &v1, &v2);
        // MM is One based - just like FORTRAN
        ++v1;
        ++v2;
        assert(!err);
        adj << v1 << " " << v2 << " 0.0\n"; // forward edge
        adj << v2 << " " << v1 << " 0.0\n"; // reverse edge
        IGRAPH_EIT_NEXT(eit);
    }

    igraph_eit_destroy(&eit);
    adj.close();
}

//
// Colorize the contact graph with colpack.
//
// Algorithm and Ordering is from GraphColoringInterface.h
//

void WktColor::colpack_color()
{
    auto g = ColPack::GraphColoringInterface(
        SRC_FILE,
        adj_file.c_str(),
        "MM");
    g.Coloring(ordering, algorithm);
    if (verbose) {
        g.PrintVertexColoringMetrics();
    }
    vector<int> color;
    g.GetVertexColors(color);
    for (unsigned int i=0; i<color.size(); ++i) {
        auto err = igraph_cattribute_VAN_set(&contact, "color", i, color[i]);
        assert(!err);
    }

    if (adj_tmp) {
        ::unlink(adj_file.c_str());
    }
}

void WktColor::color()
{
    if (algorithm == "igraph") {
        igraph_color();
    } else {
        colpack_create_adj_file();
        colpack_color();
    }

    gsl::owner<FILE *> output = nullptr;

    output = gml_file.size() ?
        gsl::owner<FILE *>(::fopen(gml_file.c_str(), "w")) :
        gsl::owner<FILE *>(stdout);
    assert(output);

    auto err = igraph_write_graph_graphml(&contact, output, true);
    assert(!err);

    ::fclose(output);
}

//
// Show the Clique Number for the graph
// int igraph_clique_number(const igraph_t *graph, igraph_integer_t *no);
void WktColor::show_clique_number()
{
    igraph_integer_t clique_number = 0;
    auto err = igraph_clique_number(&contact, &clique_number);

    assert(!err);
    std::cout << clique_number << std::endl;
}

//
// Save the mesh geometry to an 'OFF' file.
//
void WktColor::save_mesh_geometry()
{
    try {
        auto ok = OpenMesh::IO::write_mesh(mesh, geometry);
        assert(ok);
    }
    catch( std::exception& x ) {
        std::cerr << x.what() << std::endl;
    }
}

//
// Close files, release memory.
//
void WktColor::close()
{
    for (auto & g : centroid) {
        GEOSGeom_destroy_r(wkt.handle, g);
    }

    igraph_destroy(&contact);
    wkt_close(&wkt);
}

//
// Print a usage message
//
namespace {
void usage()
{
    std::cout
        << "usage: wktcolors [-gOFF] [-hv] file.wkt"
        << std::endl;
}
}

//
// main
//
int main(int argc, char *argv[])
{
    int rc = EXIT_FAILURE;
    WktColor wktcolor;

    try {
        boost::optional<std::string> geometry;
        boost::optional<std::string> algorithm;
        boost::optional<std::string> ordering;
        boost::optional<std::string> adj_file;
        boost::optional<std::string> gml_file;
        boost::program_options::options_description
            options("Options");
        options.add_options()
            ("verbose,v", "verbose")
            ("geometry,g", boost::program_options::value(&geometry), "geometry")
            ("clique-number,N", "Clique number")
            ("algorithm,a",
             boost::program_options::value(&algorithm),
             "Coloring algorithm")
            ("ordering,O",
             boost::program_options::value(&ordering),
             "Color ordering")
            ("adj-file,m",
             boost::program_options::value(&adj_file),
             "ADJ file")
            ("gml-file,o",
             boost::program_options::value(&gml_file),
             "GML output file")
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
            exit(0);
        }

        if (cli.count("verbose")) {
            wktcolor.set_verbose();
        }

        if (cli.count("geometry")) {
            wktcolor.set_geometry(*geometry);
        }

        if (cli.count("clique-number")) {
            wktcolor.request_clique_number();
        }

        if (cli.count("algorithm")) {
            wktcolor.set_algorithm(*algorithm);
        }

        if (cli.count("ordering")) {
            wktcolor.set_ordering(*ordering);
        }

        if (cli.count("adj-file")) {
            wktcolor.set_adj_file(*adj_file);
        }

        if (cli.count("gml-file")) {
            wktcolor.set_gml_file(*gml_file);
        }

        auto v = args.options;
        auto fidx = v.size()-1;
        if (v[fidx].position_key == 0) {
            wktcolor.run(v[fidx].value[0]);
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
