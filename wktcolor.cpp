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

#include <wkt.h>

using boost::typeindex::type_id_with_cvr;

typedef OpenMesh::PolyMesh_ArrayKernelT<>  MyMesh;

/*

  Strategy:

  wktlib for file reading
  geos_c for basic geometry
  openmesh for polygonalization
  igraph for graph analysis

  #include <utility>
  typedef std::pair<double,double> vertex;
  std::map<vertex,int> uvertex;
  std::vector<vertex>  svertex;

  #include <OpenMesh/Core/IO/MeshIO.hh>
  #include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
  typedef OpenMesh::PolyMesh_ArrayKernelT<>  MyMesh;

  wkt_open();
  // count points and faces
  wkt_iterate() {
    error if ! polygon
    for each point P {
      ++points;
    }
    ++faces;
  }

  MyMesh mesh;
  MyMesh::VertexHandle vhandle[points];
  std::vector<MyMesh::VertexHandle>  face_vhandles;
  wkt_iterate() {
    face_vhandles.clear();
    for each point P {
      vertex v(P.x, P.y);
      int idx = uvertex[vertex];
      if (!idx) {
        idx = svertex.size()+1
        uvertex[v] = idx;
        svertex.push_back(v);
        vhandle[idx-1] = mesh.add_vertex(MyMesh::Point(v[0], v[1], 0));
        face_vhandles.push_back(vhandle[idx-1]);
      }
    }
    mesh.add_face(face_vhandles);
  }

  wkt_close();

 */

typedef std::pair<double,double> vertex;

struct polyinfo {
    int points;
    int faces;
};

class WktColor {
public:
    int colors = 4;
    void run(std::string file);
private:
    MyMesh mesh;
    struct wkt wkt;
    struct polyinfo info;
    std::map<vertex,int> uvertex;
    std::vector<vertex>  svertex;

    void open();
    void read(std::string file);
    void count();
    void create_mesh();
    void color();
    void close();
};

void WktColor::open()
{
    int err = wkt_open(&wkt);
    assert(!err);
}

void WktColor::read(std::string file)
{
    int err = wkt_read(&wkt, file.c_str());
    assert(!err);
}

void WktColor::run(std::string file)
{
    open();
    read(file);
    create_mesh();
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
    int err;
    int n;
    unsigned int dim;
    unsigned int size;

    // Only Polygons with no holes
    assert(!strcmp(gtype, "Polygon"));
    g = GEOSGetExteriorRing_r(wkt->handle, geom);
    assert(g != NULL);
    n = GEOSGetNumInteriorRings_r(wkt->handle, geom);
    assert(n == 0);

    seq = GEOSGeom_getCoordSeq(g);
    assert(seq != NULL);

    // Only two dimensions
    err = GEOSCoordSeq_getDimensions(seq, &dim);
    assert(!err);
    assert(dim == 2);

    err = GEOSCoordSeq_getSize(seq, &size);
    assert(!err);
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
  std::vector<MyMesh::VertexHandle>  face_vhandles;
  int n;
  int err;

  n = GEOSGetNumGeometries_r(wkt.handle, wkt.geom);
  for (int i=0; i<n; ++i) {
      const GEOSGeometry *g;
      const GEOSCoordSequence *seq;
      unsigned int size;

      face_vhandles.clear();
      g = GEOSGetGeometryN_r(wkt.handle, wkt.geom, i);
      assert(g != NULL);
      seq = GEOSGeom_getCoordSeq(g);
      assert(seq != NULL);
      err = GEOSCoordSeq_getSize(seq, &size);
      assert(!err);
      for (unsigned int j = 0; j < size; ++j) {
          double x;
          double y;
          GEOSCoordSeq_getXY(seq, j, &x, &y);
          vertex v(x, y);
          // idx 0 means "not found" so uvertex values are idx+1
          // other vertex vectors are zero based.
          int idx = uvertex[v];
          if (!idx) {

              assert(svertex.size() < (unsigned int)info.points);
              idx = svertex.size() + 1;
              uvertex[v] = idx;
              svertex.push_back(v);
              vhandle[idx-1] = mesh.add_vertex(MyMesh::Point(x, y, 0));
              face_vhandles.push_back(vhandle[idx-1]);
          }
      }
      mesh.add_face(face_vhandles);
  }

}

void WktColor::color()
{
    int err;
    try {
        err = OpenMesh::IO::write_mesh(mesh, "output.off");
        assert(!err);
    }
    catch( std::exception& x ) {
        std::cerr << x.what() << std::endl;
    }
}

void WktColor::close()
{
    int err = wkt_close(&wkt);
    assert(!err);
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
