#include "defines.hpp"
#include "fem.hpp"
#include "file_io.hpp"
#include <array>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

void parse_nodes(std::vector<char> *buffer, int number_of_lines,
                 std::vector<fem::Node> *nodes_list) {
  std::string _;
  std::stringstream buffer_stream;
  // use the vector buffer as a buffer for the string stream
  buffer_stream.rdbuf()->pubsetbuf(reinterpret_cast<char *>(&(*buffer)[0]),
                                   buffer->size());
  int node_number;
  double x_coord;
  double y_coord;
  for (int i = 0; i < number_of_lines; i++) {
    buffer_stream >> node_number >> _ >> x_coord >> y_coord;
    (*nodes_list)[node_number] = fem::Node{x_coord, y_coord};
  }
};

void parse_edges(std::vector<char> *buffer, int number_of_lines,
                 std::vector<fem::Element<fem::FEM_EDGE>> *edges_list) {
  std::string _;
  std::stringstream buffer_stream;
  // use the vector buffer as a buffer for the string stream
  buffer_stream.rdbuf()->pubsetbuf(reinterpret_cast<char *>(&(*buffer)[0]),
                                   buffer->size());
  int edge_number;
  int first_node;
  int second_node;
  for (int i = 0; i < number_of_lines; i++) {
    buffer_stream >> edge_number >> _ >> first_node >> second_node;
    (*edges_list)[edge_number] = fem::Element<fem::FEM_EDGE>{
        std::array<int, fem::FEM_EDGE>({first_node, second_node})};
  }
};

void fem::read_mesh_file(const char *filename) {
  std::string _;                     // use to discard string
  std::ifstream mesh_file(filename); // open file stream

  /*
  Here we will be reading and creating the nodes, the problem will be created
  when we will have the type of mesh used
  */
  int number_of_nodes;
  mesh_file >> _ >> _ >> _ >> number_of_nodes; // parse number of nodes

  std::vector<fem::Node> nodes_list(number_of_nodes);

  int lines_per_thread = number_of_nodes / (MAX_THREADS - 1);
  int remains = number_of_nodes % (MAX_THREADS - 1);
  std::array<std::vector<char>, MAX_THREADS - 1> buffers_list;
  std::array<std::thread, MAX_THREADS - 1> thread_list;
  // first thread is special since it will pickup the remaining lines
  // create first buffer
  buffers_list[0] =
      std::vector<char>(NODES_LINE_LENGTH * (lines_per_thread + remains));
  // fill it
  mesh_file.read(buffers_list[0].data(), sizeof(char) * NODES_LINE_LENGTH *
                                             (lines_per_thread + remains));
  // thread it
  thread_list[0] = std::thread(parse_nodes, &(buffers_list[0]),
                               (lines_per_thread + remains), &nodes_list);

  // all other threads
  for (int i = 1; i < MAX_THREADS - 1; i++) {
    buffers_list[i] = std::vector<char>(NODES_LINE_LENGTH * lines_per_thread);
    mesh_file.read(buffers_list[i].data(),
                   sizeof(char) * NODES_LINE_LENGTH * lines_per_thread);
    thread_list[i] = std::thread(parse_nodes, &(buffers_list[i]),
                                 lines_per_thread, &nodes_list);
  }

  /*Now the edges*/
  fem::Mesh<FEM_EDGE> edges_list;
  int number_of_edges;
  mesh_file >> _ >> _ >> _ >> number_of_edges;
  edges_list.elements_lists = std::vector<Element<FEM_EDGE>>(number_of_edges);

  lines_per_thread = number_of_edges / (MAX_THREADS - 1);
  remains = number_of_edges % (MAX_THREADS - 1);
  // first thread is special since it will pickup the remaining lines
  for (int i = 0; i < MAX_THREADS - 1; i++) {
    thread_list[i].join();
  }
  // refill it
  mesh_file.read(buffers_list[0].data(), sizeof(char) * EDGES_LINE_LENGTH *
                                             (lines_per_thread + remains));
  // thread it
  thread_list[0] = std::thread(parse_nodes, &(buffers_list[0]),
                               (lines_per_thread + remains), &edges_list);

  // all other threads
  for (int i = 1; i < MAX_THREADS - 1; i++) {
    buffers_list[i] = std::vector<char>(EDGES_LINE_LENGTH * lines_per_thread);
    mesh_file.read(buffers_list[i].data(),
                   sizeof(char) * EDGES_LINE_LENGTH * lines_per_thread);
    thread_list[i] = std::thread(parse_nodes, &(buffers_list[i]),
                                 lines_per_thread, &nodes_list);
  }
};
