/*
   wktcolor.cpp

   Copyright (c) 2021 by Daniel Kelley

*/

#include <assert.h>
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

typedef OpenMesh::PolyMesh_ArrayKernelT<>  MyMesh;

typedef std::pair<double,double> vertex;

struct polyinfo {
    int points;
    int faces;
};

class WktColor {
public:
    int colors = 4;
    int verbose = 0;
    std::string geometry = {};
    void run(std::string file);
private:
    MyMesh mesh;
    struct wkt wkt = {};
    struct polyinfo info = {};
    igraph_t contact = {};
    std::map<vertex,int> uvertex;
    std::vector<vertex>  svertex;

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
    int err;
    memset(&wkt, 0, sizeof(wkt));
    wkt.reader = WKT_IO_ASCII;
    err = wkt_open(&wkt);
    assert(!err);
}

void WktColor::read(std::string file)
{
    int err;

    err = wkt_read(&wkt, file.c_str());
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
    struct polyinfo *info = (struct polyinfo *)user_data;
    const GEOSGeometry *g;
    const GEOSCoordSequence *seq;
    int ok;
    int n;
    unsigned int dim;
    unsigned int size;

    // Only Polygons with no holes
    assert(!strcmp(gtype, "Polygon"));
    g = GEOSGetExteriorRing_r(wkt->handle, geom);
    assert(g != NULL);
    n = GEOSGetNumInteriorRings_r(wkt->handle, geom);
    assert(n == 0);

    seq = GEOSGeom_getCoordSeq_r(wkt->handle, g);
    assert(seq != NULL);

    // Only two dimensions
    ok = GEOSCoordSeq_getDimensions_r(wkt->handle, seq, &dim);
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
    int err;
    err = wkt_iterate(&wkt, geom_counter, &info);
    assert(!err);
}

void WktColor::create_mesh()
{
    MyMesh::VertexHandle vhandle[info.points];
    std::vector<MyMesh::VertexHandle> face_vhandles;
    int n;
    int ok;

    n = GEOSGetNumGeometries_r(wkt.handle, wkt.geom);
    for (int i=0; i<n; ++i) {
        const GEOSGeometry *poly;
        const GEOSGeometry *ring;
        const GEOSCoordSequence *seq;
        unsigned int size;

        face_vhandles.clear();
        poly = GEOSGetGeometryN_r(wkt.handle, wkt.geom, i);
        assert(poly != NULL);
        ring = GEOSGetExteriorRing_r(wkt.handle, poly);
        assert(ring != NULL);
        seq = GEOSGeom_getCoordSeq_r(wkt.handle, ring);
        assert(seq != NULL);
        ok = GEOSCoordSeq_getSize_r(wkt.handle, seq, &size);
        assert(ok);
        if (verbose) {
            std::cout << "poly\n";
        }
        for (unsigned int j = 0; j < size-1; ++j) {
            double x;
            double y;
            ok = GEOSCoordSeq_getXY_r(wkt.handle, seq, j, &x, &y);
            assert(ok);
            if (verbose) {
                std::cout << "  (" << x << "," << y << ")";
            }
            vertex v(x, y);
            // idx 0 means "not found" so uvertex values are idx+1
            // other vertex vectors are zero based.
            int idx = uvertex[v];
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
 *     create edge to face on partner dart
 *   end
 * end
 */
void WktColor::create_contact_graph()
{
    igraph_vector_t edge;

    // points is an over-estimate
    igraph_vector_init(&edge, info.points);
    for (MyMesh::ConstFaceIter face = mesh.faces_begin();
         face != mesh.faces_end();
         ++face) {
    }

    igraph_create(&contact, &edge, 0, IGRAPH_UNDIRECTED);
}

void WktColor::color()
{
}

void WktColor::save_mesh_geometry()
{
    int err;
    try {
        (void)OpenMesh::IO::write_mesh(mesh, geometry);
        assert(!err);
    }
    catch( std::exception& x ) {
        std::cerr << x.what() << std::endl;
    }
}

void WktColor::close()
{
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
        boost::optional<int> colors;
        boost::optional<std::string> geometry;
        boost::program_options::options_description
            options("Options");
        options.add_options()
            ("verbose,v", "verbose")
            ("colors,c", boost::program_options::value(&colors), "colors")
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

        if (colors) {
            wktcolor.colors = *colors;
        }

        if (cli.count("verbose")) {
            wktcolor.verbose = 1;
        }

        if (cli.count("geometry")) {
            wktcolor.geometry = *geometry;
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
        int fidx = v.size()-1;
        if (fidx >= 0 && v[fidx].position_key == 0) {
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
