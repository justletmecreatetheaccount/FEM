#include "defines.hpp"
#include "fem.hpp"
#include "file_io.hpp"
#include <algorithm>
#include <array>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <istream>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

/* This function parses the buffer and inits the nodes
 * ======
 * IN : buffer -> the buffer that contains the data to be parsed
 *    : number_of_lines
 *    : nodes_list -> the vector that will hold the parsed nodes (Output)
 */
void parse_nodes(char *buffer, int number_of_lines,
                 std::vector<fem::Node> &nodes_list) {

  int node_number;
  double x_coord;
  double y_coord;
  std::istringstream buffer_stream;
  buffer_stream.rdbuf()->pubsetbuf(buffer, number_of_lines * NODES_LINE_LENGTH);
  for (int i = 0; i < number_of_lines; i++) {
    buffer_stream >> node_number;
    buffer_stream.ignore(4);
    buffer_stream >> x_coord >> y_coord;
    nodes_list[node_number].x = x_coord;
    nodes_list[node_number].y = y_coord;
  }
};

/* This function parses the buffer and inits the elements
 * ======
 * IN : buffer -> the buffer that contains the data to be parsed
 *    : number_of_lines
 *    : element_list -> the vector that will hold the parsed elements (Output)
 */
template <fem::ElementType T>
void parse_elements(char *buffer, int number_of_lines,
                    std::vector<fem::Element<T>> &element_list,
                    fem::System system,
                    std::vector<std::vector<int>> &connections_list,
                    std::vector<std::mutex> &mutexes_list);

/* Specialisation for FEM_EDGE
 */
template <>
void parse_elements<fem::FEM_EDGE>(
    char *buffer, int number_of_lines,
    std::vector<fem::Element<fem::FEM_EDGE>> &element_list, fem::System system,
    std::vector<std::vector<int>> &connections_list,
    std::vector<std::mutex> &mutexes_list) {
  std::stringstream buffer_stream;
  // use the vector buffer as a buffer for the string stream
  buffer_stream.rdbuf()->pubsetbuf(buffer, number_of_lines * EDGES_LINE_LENGTH);
  int edge_number;
  int nodes[2];
  for (int i = 0; i < number_of_lines; i++) {
    buffer_stream >> edge_number;
    buffer_stream.ignore(6);
    buffer_stream >> nodes[0] >> nodes[1];
    std::sort(nodes, nodes + 2);
    element_list[edge_number].element_nodes[0] = nodes[0];
    element_list[edge_number].element_nodes[1] = nodes[1];

    // INIT SYSTEM
    bool already_linked = false;
    mutexes_list[nodes[0]].lock();
    for (size_t j = 0; j < connections_list[nodes[0]].size(); j++) {
      if (connections_list[nodes[0]][j] == nodes[1]) {
        already_linked = true;
      }
    }
    if (!already_linked) {
      connections_list[nodes[0]].push_back(nodes[1]);
    }
    mutexes_list[nodes[0]].unlock();
  }
};

/* Specialisation for FEM_TRIANGLE
 */
template <>
void parse_elements<fem::FEM_TRIANGLE>(
    char *buffer, int number_of_lines,
    std::vector<fem::Element<fem::FEM_TRIANGLE>> &element_list,
    fem::System system, std::vector<std::vector<int>> &connections_list,
    std::vector<std::mutex> &mutexes_list) {
  std::stringstream buffer_stream;
  // use the vector buffer as a buffer for the string stream
  buffer_stream.rdbuf()->pubsetbuf(buffer,
                                   number_of_lines * TRIANGLES_LINE_LENGTH);
  int triangle_number;
  int nodes[3];
  for (int i = 0; i < number_of_lines; i++) {
    buffer_stream >> triangle_number;
    buffer_stream.ignore(6);
    buffer_stream >> nodes[0] >> nodes[1] >> nodes[2];
    std::sort(nodes, nodes + 3);

    // very ugly but element_nodes already exists
    element_list[triangle_number].element_nodes[0] = nodes[0];
    element_list[triangle_number].element_nodes[1] = nodes[1];
    element_list[triangle_number].element_nodes[2] = nodes[2];

    // INIT SYSTEM
    // first node
    bool already_linked_1 = false;
    bool already_linked_2 = false;
    mutexes_list[nodes[0]].lock();
    for (size_t j = 0; j < connections_list[nodes[0]].size(); j++) {
      if (connections_list[nodes[0]][j] == nodes[1]) {
        already_linked_1 = true;
      }
      if (connections_list[nodes[0]][j] == nodes[2]) {
        already_linked_2 = true;
      }
    }
    if (!already_linked_1) {
      connections_list[nodes[0]].push_back(nodes[1]);
    }
    if (!already_linked_2) {
      connections_list[nodes[0]].push_back(nodes[2]);
    }
    mutexes_list[nodes[0]].unlock();
    // second node
    already_linked_2 = false;
    mutexes_list[nodes[1]].lock();

    for (size_t j = 0; j < connections_list[nodes[1]].size(); j++) {
      if (connections_list[nodes[1]][j] == nodes[2]) {
        already_linked_1 = true;
      }
    }
    if (!already_linked_2) {
      connections_list[nodes[1]].push_back(nodes[2]);
    }
    mutexes_list[nodes[1]].unlock();
  }
};

void init_system(fem::System &system,
                 std::vector<std::vector<int>> &connections_list,
                 int number_of_nodes, int number_of_edges,
                 int number_of_triangles) {
  // System Init
  int number_internal_links = (number_of_triangles * 3 - number_of_edges) / 2;
  int links_per_element = (3 - 1) * 3 / 2;
  int data_size =
      (3 * number_of_nodes +
       4 * (links_per_element * number_of_triangles - number_internal_links));
  system.row_ptr.resize(number_of_nodes * 2 + 1);
  system.data.resize(data_size);
  system.column.resize(data_size);
  system.B.resize(number_of_nodes * 2);
  system.size = 2 * number_of_nodes * 2;

  for (int i = 1; i < number_of_nodes * 2; i += 2) {
    // the node itself
    system.row_ptr[i + 1] +=
        system.row_ptr[i] + connections_list[i].size() * 2 + 2;
    system.row_ptr[i + 2] +=
        system.row_ptr[i + 1] + connections_list[i].size() * 2 + 1;
    // first row
    system.column[system.row_ptr[i]] = i * 2;
    system.column[system.row_ptr[i] + 1] = i * 2 + 1;

    for (size_t j = 0; j < connections_list[i].size(); j++) {
      system.column[system.row_ptr[i] + j + 2] = connections_list[i][j] * 2;
      system.column[system.row_ptr[i] + j + 3] = connections_list[i][j] * 2 + 1;
    }
    // second row
    system.column[system.row_ptr[i + 1]] = i * 2 + 1;
    for (size_t j = 0; j < connections_list[i].size(); j++) {
      system.column[system.row_ptr[i + 1] + j + 1] = connections_list[i][j] * 2;
      system.column[system.row_ptr[i + 1] + j + 2] =
          connections_list[i][j] * 2 + 1;
    }
  }
};

/* This function parses the buffer and inits the domains
 * ======
 * IN : buffer -> the buffer that contains the data to be parsed
 *    : domain_list -> the vector that will hold the parsed domains (Output)
 *    : number_of_domains -> the number of domains that need to be parsed
 *    : elements_list -> list to the edges elements of the mesh, the domain
 * will pick the ones that belong to it
 */
void parse_domains(std::ifstream &stream, std::vector<fem::Domain> &domain_list,
                   int number_of_domains,
                   std::vector<fem::Element<fem::FEM_EDGE>> &elements_list) {
  for (int i = 0; i < number_of_domains; i++) {
    stream.ignore(29);
    domain_list.emplace_back(fem::Domain());
    char temp_trash_c_string[256]; // UTTER GARBAJE (256 chosen arbitrarily
                                   //  because I know it will work just fine
    stream.getline(temp_trash_c_string, 256);
    domain_list[i].name.append(temp_trash_c_string);
    stream.ignore(27);
    stream >> domain_list[i].number_of_elements;
    for (int j = 0; j < domain_list[i].number_of_elements; j++) {
      int element_number;
      stream >> element_number;
      domain_list[i].elements.push_back(&(elements_list[element_number]));
    }
  }
}

/* This function reads the in file and creates all the necessary structures
 * ======
 * IN : filename; the name of the file
 * OUT : fem::Problem; the filled Problem
 */

fem::Problem fem::read_mesh_file(const char *filename) {
  fem::Problem problem;
  std::ifstream mesh_file(filename); // open file stream

  mesh_file.ignore(16); // not very clean but intended way of discarding
  mesh_file >> problem.geometry.number_of_nodes; // parse number of nodes
  problem.geometry.nodes_list.resize(
      problem.geometry.number_of_nodes); // prealloc vector

  int lines_per_thread = problem.geometry.number_of_nodes / (MAX_THREADS);
  int remains = problem.geometry.number_of_nodes % (MAX_THREADS);
  std::array<char *, MAX_THREADS> buffer_list;
  std::array<std::thread, MAX_THREADS> thread_list;
  // first thread is special since it will pickup the remaining lines
  // init the buffer
  buffer_list[0] = (char *)malloc(sizeof(char) * NODES_LINE_LENGTH *
                                  (lines_per_thread + remains));
  // fill the first buffer
  mesh_file.read(buffer_list[0], sizeof(char) * NODES_LINE_LENGTH *
                                     (lines_per_thread + remains));
  // thread it
  thread_list[0] =
      std::thread(parse_nodes, buffer_list[0], (lines_per_thread + remains),
                  std::ref(problem.geometry.nodes_list));

  // all other threads
  for (int i = 1; i < MAX_THREADS; i++) {
    buffer_list[i] =
        (char *)malloc(sizeof(char) * NODES_LINE_LENGTH * lines_per_thread);
    mesh_file.read(buffer_list[i],
                   sizeof(char) * NODES_LINE_LENGTH * lines_per_thread);
    thread_list[i] = std::thread(parse_nodes, buffer_list[i], lines_per_thread,
                                 std::ref(problem.geometry.nodes_list));
  }

  // For the elements we need to make a list of "connections" for each node in
  // order to make the csr matrix. This array stores the nodes of greater ordrer
  // that the nodes are connected to
  std::vector<std::vector<int>> connections_list(
      problem.geometry.number_of_nodes);
  // mutexes to avoid changing the same vector at the same time
  std::vector<std::mutex> mutextes_list(problem.geometry.number_of_nodes);

  /*Now the edges*/
  int number_of_edges;
  mesh_file.ignore(17);
  mesh_file >> number_of_edges;
  problem.geometry.edge_mesh.elements_lists.resize(number_of_edges);

  lines_per_thread = number_of_edges / (MAX_THREADS);
  remains = number_of_edges % (MAX_THREADS);
  // first thread is special since it will pickup the remaining lines
  for (int i = 0; i < MAX_THREADS; i++) {
    thread_list[i].join();
  }
  // init the buffer
  buffer_list[0] =
      (char *)realloc(buffer_list[0], sizeof(char) * EDGES_LINE_LENGTH *
                                          (lines_per_thread + remains));
  // refill it
  mesh_file.read(buffer_list[0], sizeof(char) * EDGES_LINE_LENGTH *
                                     (lines_per_thread + remains));
  // thread it
  thread_list[0] = std::thread(
      parse_elements<fem::FEM_EDGE>, buffer_list[0],
      (lines_per_thread + remains),
      std::ref(problem.geometry.edge_mesh.elements_lists), problem.system,
      std::ref(connections_list), std::ref(mutextes_list));

  // all other threads
  for (int i = 1; i < MAX_THREADS; i++) {
    buffer_list[i] = (char *)realloc(
        buffer_list[i], sizeof(char) * EDGES_LINE_LENGTH * lines_per_thread);
    mesh_file.read(buffer_list[i],
                   sizeof(char) * EDGES_LINE_LENGTH * lines_per_thread);
    thread_list[i] = std::thread(
        parse_elements<fem::FEM_EDGE>, buffer_list[i], lines_per_thread,
        std::ref(problem.geometry.edge_mesh.elements_lists), problem.system,
        std::ref(connections_list), std::ref(mutextes_list));
  }
  /*Now the triangles*/
  int number_of_triangles;
  mesh_file.ignore(21);
  mesh_file >> number_of_triangles;
  mesh_file.ignore(2); // god knows why no space after triangle lines
  problem.geometry.full_mesh.elements_lists.resize(number_of_triangles);

  lines_per_thread = number_of_triangles / (MAX_THREADS);
  remains = number_of_triangles % (MAX_THREADS);
  // first thread is special since it will pickup the remaining lines
  for (int i = 0; i < MAX_THREADS; i++) {
    thread_list[i].join();
  }
  // init the buffer
  buffer_list[0] =
      (char *)realloc(buffer_list[0], sizeof(char) * TRIANGLES_LINE_LENGTH *
                                          (lines_per_thread + remains));
  // refill it
  mesh_file.read(buffer_list[0], sizeof(char) * TRIANGLES_LINE_LENGTH *
                                     (lines_per_thread + remains));
  // thread it
  thread_list[0] = std::thread(
      parse_elements<fem::FEM_TRIANGLE>, buffer_list[0],
      (lines_per_thread + remains),
      std::ref(problem.geometry.full_mesh.elements_lists), problem.system,
      std::ref(connections_list), std::ref(mutextes_list));

  // all other threads
  for (int i = 1; i < MAX_THREADS; i++) {
    buffer_list[i] =
        (char *)realloc(buffer_list[i], sizeof(char) * TRIANGLES_LINE_LENGTH *
                                            lines_per_thread);
    mesh_file.read(buffer_list[i],
                   sizeof(char) * TRIANGLES_LINE_LENGTH * lines_per_thread);
    thread_list[i] = std::thread(
        parse_elements<fem::FEM_TRIANGLE>, buffer_list[i], lines_per_thread,
        std::ref(problem.geometry.full_mesh.elements_lists), problem.system,
        std::ref(connections_list), std::ref(mutextes_list));
  }

  /*Now the domains*/
  mesh_file.ignore(18);
  mesh_file >> problem.geometry.number_of_domains;
  problem.geometry.domains_list.reserve(problem.geometry.number_of_domains);

  parse_domains(mesh_file, problem.geometry.domains_list,
                problem.geometry.number_of_domains,
                problem.geometry.edge_mesh.elements_lists);

  // Can join after bc vector does not  realloc
  for (int i = 0; i < MAX_THREADS; i++) {
    thread_list[i].join();
  }

  init_system(problem.system, connections_list,
              problem.geometry.number_of_nodes, number_of_edges,
              number_of_triangles);

  for (int i = 0; i < MAX_THREADS; i++) {
    free(buffer_list[i]);
  }
  return problem;
};
