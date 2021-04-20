/*
   wktcolor.cpp

   Copyright (c) 2021 by Daniel Kelley

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

using boost::typeindex::type_id_with_cvr;

using MyMesh = OpenMesh::PolyMesh_ArrayKernelT<>;

using vertex = std::pair<double,double>;

struct polyinfo {
    unsigned int points;
    unsigned int faces;
};

class WktColor {
public:
    void run(std::string file);
    inline void set_verbose() { verbose=1; }
    inline void set_geometry(std::string s) { geometry=s; }
private:
    MyMesh mesh;
    struct wkt wkt = {};
    struct polyinfo info = {};
    igraph_t contact = {};
    std::map<vertex,int> uvertex;
    std::vector<vertex> svertex;
    std::vector<GEOSGeometry *> centroid;
    int verbose = 0;
    std::string geometry = {};

    void open();
    void read(std::string file);
    void count();
    void create_mesh();
    void save_mesh_geometry();
    void create_contact_graph();
    void color();
    void close();
};

void WktColor::open()
{
    memset(&wkt, 0, sizeof(wkt));
    wkt.reader = WKT_IO_ASCII;
    auto err = wkt_open(&wkt);
    assert(!err);
    // must be called before other igraph functions
    // according to comment in igraph cattributes.c
    igraph_i_set_attribute_table(&igraph_cattribute_table);
}

void WktColor::read(std::string file)
{
    auto err = wkt_read(&wkt, file.c_str());
    assert(!err);
}

void WktColor::run(std::string file)
{
    open();
    read(file);
    count();
    create_mesh();
    if (geometry.size() > 0) {
        save_mesh_geometry();
    }
    create_contact_graph();
    color();
    close();
}

static int geom_counter(struct wkt *wkt,
                        const GEOSGeometry *geom,
                        const char *gtype,
                        void *user_data)
{
    auto info = (struct polyinfo *)user_data;
    unsigned int dim = 0;
    unsigned int size = 0;

    // Only Polygons with no holes
    assert(!strcmp(gtype, "Polygon"));
    auto g = GEOSGetExteriorRing_r(wkt->handle, geom);
    assert(g != nullptr);
    auto n = GEOSGetNumInteriorRings_r(wkt->handle, geom);
    assert(n == 0);

    auto seq = GEOSGeom_getCoordSeq_r(wkt->handle, g);
    assert(seq != nullptr);

    // Only two dimensions
    auto ok = GEOSCoordSeq_getDimensions_r(wkt->handle, seq, &dim);
    assert(ok);
    assert(dim == 2);

    ok = GEOSCoordSeq_getSize_r(wkt->handle, seq, &size);
    assert(ok);
    assert(size != 0);

    info->faces += 1;
    info->points += size;

    return 0;
}

// Count points and faces
void WktColor::count()
{
    auto err = wkt_iterate(&wkt, geom_counter, &info);
    assert(!err);
}

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

                assert(svertex.size() < (unsigned int)info.points);
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

/*
 * for each face
 *   create graph node for face
 *   for each dart
 *     create edge to face on partner dart if not boundary
 *   end
 * end
 *  FaceHandle opposite_face_handle(HalfedgeHandle _heh)
 */
void WktColor::create_contact_graph()
{
    igraph_vector_t edge;
    auto edge_limit = info.points * 2;
    auto edge_id =  info.points * 0; // * 0 to derive type
    // edge_limit is an over-estimate
    igraph_vector_init(&edge, edge_limit);

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
    igraph_vector_destroy(&edge);
}

void WktColor::color()
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
        auto g = centroid[i];
        double x = 0.0;
        double y = 0.0;

        auto ok = GEOSGeomGetX_r(wkt.handle, g, &x);
        assert(ok);
        ok = GEOSGeomGetY_r(wkt.handle, g, &y);
        assert(ok);
        auto color = VECTOR(cv)[i];
        if (verbose) {
            std::cout
                << "color( "
                << i
                << ","
                << x
                << ","
                << y
                << ","
                << color
                << ")"
                << std::endl;
        }

        err = igraph_cattribute_VAN_set(&contact, "color", i, color);
        assert(!err);
        err = igraph_cattribute_VAN_set(&contact, "x", i, x);
        assert(!err);
        err = igraph_cattribute_VAN_set(&contact, "y", i, y);
        assert(!err);

    }

    err = igraph_write_graph_graphml(&contact, stdout, 1);
    assert(!err);

    igraph_vector_int_destroy(&cv);
}

void WktColor::save_mesh_geometry()
{
    try {
        auto err = OpenMesh::IO::write_mesh(mesh, geometry);
        assert(!err);
    }
    catch( std::exception& x ) {
        std::cerr << x.what() << std::endl;
    }
}

void WktColor::close()
{
    for (auto & g : centroid) {
        GEOSGeom_destroy_r(wkt.handle, g);
    }

    igraph_destroy(&contact);
    wkt_close(&wkt);
}

static void usage()
{
    std::cout
        << "usage: wktcolors [-cN] [-hv] file.wkt"
        << std::endl;
}

int main(int argc, char *argv[])
{
    int rc = EXIT_FAILURE;
    WktColor wktcolor;

    try {
        boost::optional<std::string> geometry;
        boost::program_options::options_description
            options("Options");
        options.add_options()
            ("verbose,v", "verbose")
            ("geometry,g", boost::program_options::value(&geometry), "geometry")
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
