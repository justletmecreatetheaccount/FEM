# Element finis : disque de frein en spirale sous pression

Ce projet s'intéresse aux déformations d'un disque de frein chaud sous pression et à l'arrêt après un premier échauffement.

Fais par le groupe 29.

## A propos du projet

Ce projet ce concentre sur un cas assez spécifique dans le fonctionnement d'un disque de frein, le cas suivant : 
- Un disque à l'arrêt
- Qui est a température élevée, on fais se base ici sur le fait que la température reste élévée beaucoup de temps apres utilisation
- Des déformations élastiques linéaires

Le but principal du projet est de trouver une solution qui puisse maximiser la largeur des canaux de refroidissement (trous) tout en gardant une résistance de la pièce assez élevée. Pour aller plus loin on pourrait également faire l'analyse pour maximiser le refroidissement par exemple mais cela sort du cadre de ce projet.
Plus concrêtement, le projet fournis un programme qui permet d'analyser les déformations pour une pièces dont les canaux de refroidissement sont déjà fournis (mais modifiable) sous une certaine conditions frontières (aussi fournie mais modifiable).

## Pour commencer

### Prérequis

- [`python`](https://www.python.org/downloads/)

- `matplotlib`
```sh
pip install matplotlib
```
- numpy
```sh
pip install numpy
```

- Cmake <= 3.13 
- X11 (Linux)
- tdm-GCC ou GCC
- GLFW et OpenGL

### Structure







.
├── benchmark.c
├── benchmark.h
├── build
│   ├── CMakeCache.txt
│   ├── CMakeFiles
│   │   ├── 3.28.3
│   │   │   ├── CMakeCCompiler.cmake
│   │   │   ├── CMakeCXXCompiler.cmake
│   │   │   ├── CMakeDetermineCompilerABI_C.bin
│   │   │   ├── CMakeDetermineCompilerABI_CXX.bin
│   │   │   ├── CMakeSystem.cmake
│   │   │   ├── CompilerIdC
│   │   │   │   ├── a.out
│   │   │   │   ├── CMakeCCompilerId.c
│   │   │   │   └── tmp
│   │   │   └── CompilerIdCXX
│   │   │       ├── a.out
│   │   │       ├── CMakeCXXCompilerId.cpp
│   │   │       └── tmp
│   │   ├── cmake.check_cache
│   │   ├── CMakeConfigureLog.yaml
│   │   ├── CMakeDirectoryInformation.cmake
│   │   ├── CMakeScratch
│   │   ├── Makefile2
│   │   ├── Makefile.cmake
│   │   ├── myFem.dir
│   │   │   ├── build.make
│   │   │   ├── cmake_clean.cmake
│   │   │   ├── compiler_depend.make
│   │   │   ├── compiler_depend.ts
│   │   │   ├── DependInfo.cmake
│   │   │   ├── depend.make
│   │   │   ├── flags.make
│   │   │   ├── link.txt
│   │   │   └── progress.make
│   │   ├── pkgRedirects
│   │   ├── progress.marks
│   │   └── TargetDirectories.txt
│   ├── cmake_install.cmake
│   ├── compile_commands.json
│   ├── Makefile
│   └── project
│       ├── CMakeFiles
│       │   ├── CMakeDirectoryInformation.cmake
│       │   ├── femsolver.dir
│       │   │   ├── build.make
│       │   │   ├── cmake_clean.cmake
│       │   │   ├── compiler_depend.make
│       │   │   ├── compiler_depend.ts
│       │   │   ├── DependInfo.cmake
│       │   │   ├── depend.make
│       │   │   ├── flags.make
│       │   │   ├── link.txt
│       │   │   └── progress.make
│       │   └── progress.marks
│       ├── cmake_install.cmake
│       └── Makefile
├── CMakeLists.txt
├── data
│   ├── fixMesh.py
│   ├── mesh.geo
│   ├── mesh.msh
│   ├── mesh_ref.txt
│   ├── mesh.txt
│   ├── UV_ref.txt
│   └── UV.txt
├── genmesh
│   ├── CMakeLists.txt
│   ├── data
│   │   └── mesh.txt
│   ├── glfw
│   │   ├── CMake
│   │   │   ├── amd64-mingw32msvc.cmake
│   │   │   ├── GenerateMappings.cmake
│   │   │   ├── i586-mingw32msvc.cmake
│   │   │   ├── i686-pc-mingw32.cmake
│   │   │   ├── i686-w64-mingw32.cmake
│   │   │   ├── MacOSXBundleInfo.plist.in
│   │   │   ├── modules
│   │   │   │   ├── FindMir.cmake
│   │   │   │   ├── FindOSMesa.cmake
│   │   │   │   ├── FindVulkan.cmake
│   │   │   │   ├── FindWaylandProtocols.cmake
│   │   │   │   └── FindXKBCommon.cmake
│   │   │   └── x86_64-w64-mingw32.cmake
│   │   ├── CMakeLists.txt
│   │   ├── include
│   │   │   └── GLFW
│   │   │       ├── glfw3.h
│   │   │       └── glfw3native.h
│   │   └── src
│   │       ├── CMakeLists.txt
│   │       ├── cocoa_init.m
│   │       ├── cocoa_joystick.h
│   │       ├── cocoa_joystick.m
│   │       ├── cocoa_monitor.m
│   │       ├── cocoa_platform.h
│   │       ├── cocoa_time.c
│   │       ├── cocoa_window.m
│   │       ├── context.c
│   │       ├── egl_context.c
│   │       ├── egl_context.h
│   │       ├── glfw3Config.cmake.in
│   │       ├── glfw3.pc.in
│   │       ├── glfw_config.h.in
│   │       ├── glx_context.c
│   │       ├── glx_context.h
│   │       ├── init.c
│   │       ├── input.c
│   │       ├── internal.h
│   │       ├── linux_joystick.c
│   │       ├── linux_joystick.h
│   │       ├── mappings.h
│   │       ├── mappings.h.in
│   │       ├── mir_init.c
│   │       ├── mir_monitor.c
│   │       ├── mir_platform.h
│   │       ├── mir_window.c
│   │       ├── monitor.c
│   │       ├── nsgl_context.h
│   │       ├── nsgl_context.m
│   │       ├── null_init.c
│   │       ├── null_joystick.c
│   │       ├── null_joystick.h
│   │       ├── null_monitor.c
│   │       ├── null_platform.h
│   │       ├── null_window.c
│   │       ├── osmesa_context.c
│   │       ├── osmesa_context.h
│   │       ├── posix_thread.c
│   │       ├── posix_thread.h
│   │       ├── posix_time.c
│   │       ├── posix_time.h
│   │       ├── vulkan.c
│   │       ├── wgl_context.c
│   │       ├── wgl_context.h
│   │       ├── win32_init.c
│   │       ├── win32_joystick.c
│   │       ├── win32_joystick.h
│   │       ├── win32_monitor.c
│   │       ├── win32_platform.h
│   │       ├── win32_thread.c
│   │       ├── win32_time.c
│   │       ├── win32_window.c
│   │       ├── window.c
│   │       ├── wl_init.c
│   │       ├── wl_monitor.c
│   │       ├── wl_platform.h
│   │       ├── wl_window.c
│   │       ├── x11_init.c
│   │       ├── x11_monitor.c
│   │       ├── x11_platform.h
│   │       ├── x11_window.c
│   │       ├── xkb_unicode.c
│   │       └── xkb_unicode.h
│   └── src
│       ├── fem.c
│       ├── fem.h
│       ├── glfem.c
│       ├── glfem.h
│       └── main.c
├── gmsh
│   └── gmsh-4.13.1-Linux64-sdk
│       ├── bin
│       │   ├── gmsh
│       │   └── onelab.py
│       ├── include
│       │   ├── gmshc.h
│       │   ├── gmsh.f90
│       │   ├── gmsh.h
│       │   └── gmsh.h_cwrap
│       ├── lib
│       │   ├── gmsh-4.13.1.dist-info
│       │   │   └── METADATA
│       │   ├── gmsh.jl
│       │   ├── gmsh.py
│       │   ├── libgmsh.so -> libgmsh.so.4.13
│       │   ├── libgmsh.so.4.13 -> libgmsh.so.4.13.1
│       │   └── libgmsh.so.4.13.1
│       ├── README.txt
│       └── share
│           ├── doc
│           │   └── gmsh
│           │       ├── CHANGELOG.txt
│           │       ├── CREDITS.txt
│           │       ├── examples
│           │       │   ├── api
│           │       │   │   ├── adapt_mesh.cpp
│           │       │   │   ├── adapt_mesh.py
│           │       │   │   ├── aneurysm_data.stl
│           │       │   │   ├── aneurysm.py
│           │       │   │   ├── as1-tu-203.stp
│           │       │   │   ├── bgmesh.pos
│           │       │   │   ├── boolean.cpp
│           │       │   │   ├── boolean.py
│           │       │   │   ├── bspline_bezier_patches.py
│           │       │   │   ├── bspline_bezier_trimmed.py
│           │       │   │   ├── bspline_filling.py
│           │       │   │   ├── circle_arc.py
│           │       │   │   ├── closest_point.py
│           │       │   │   ├── CMakeLists.txt
│           │       │   │   ├── copy_mesh.py
│           │       │   │   ├── crack3d.py
│           │       │   │   ├── crack.py
│           │       │   │   ├── custom_gui.cpp
│           │       │   │   ├── custom_gui.py
│           │       │   │   ├── cylinderFFD.py
│           │       │   │   ├── discrete.cpp
│           │       │   │   ├── discrete.jl
│           │       │   │   ├── discrete.py
│           │       │   │   ├── edges.cpp
│           │       │   │   ├── explore.cpp
│           │       │   │   ├── explore.jl
│           │       │   │   ├── explore.py
│           │       │   │   ├── extend_field.py
│           │       │   │   ├── faces.cpp
│           │       │   │   ├── flatten2.py
│           │       │   │   ├── flatten.py
│           │       │   │   ├── fragment_surfaces.cpp
│           │       │   │   ├── get_data_perf.cpp
│           │       │   │   ├── get_data_perf.py
│           │       │   │   ├── glue_and_remesh_stl.py
│           │       │   │   ├── gui.cpp
│           │       │   │   ├── gui.jl
│           │       │   │   ├── gui.py
│           │       │   │   ├── heal.py
│           │       │   │   ├── hex.py
│           │       │   │   ├── hybrid_order.py
│           │       │   │   ├── import_perf.c
│           │       │   │   ├── import_perf.cpp
│           │       │   │   ├── import_perf.jl
│           │       │   │   ├── import_perf.py
│           │       │   │   ├── mesh_from_discrete_curve.py
│           │       │   │   ├── mesh_quality.py
│           │       │   │   ├── mirror_mesh.py
│           │       │   │   ├── msh_attributes.py
│           │       │   │   ├── multi_process.py
│           │       │   │   ├── multi_thread.py
│           │       │   │   ├── naca_boundary_layer_2d.py
│           │       │   │   ├── naca_boundary_layer_3d.py
│           │       │   │   ├── neighbors.py
│           │       │   │   ├── normals.py
│           │       │   │   ├── object.stl
│           │       │   │   ├── ocean.py
│           │       │   │   ├── onelab_run_auto.c
│           │       │   │   ├── onelab_run_auto.cpp
│           │       │   │   ├── onelab_run_auto.py
│           │       │   │   ├── onelab_run.py
│           │       │   │   ├── onelab_test.jl
│           │       │   │   ├── onelab_test.py
│           │       │   │   ├── open.cpp
│           │       │   │   ├── open.py
│           │       │   │   ├── opt.py
│           │       │   │   ├── partition.cpp
│           │       │   │   ├── partition.py
│           │       │   │   ├── periodic.py
│           │       │   │   ├── pipe.py
│           │       │   │   ├── plugin.cpp
│           │       │   │   ├── plugin.py
│           │       │   │   ├── poisson.py
│           │       │   │   ├── prepro.py
│           │       │   │   ├── prim_axis.py
│           │       │   │   ├── raw_tetrahedralization.py
│           │       │   │   ├── raw_triangulation.py
│           │       │   │   ├── README.txt
│           │       │   │   ├── relocate_nodes.py
│           │       │   │   ├── remesh_stl.py
│           │       │   │   ├── remove_elements.py
│           │       │   │   ├── renumbering.py
│           │       │   │   ├── reparamOnFace.py
│           │       │   │   ├── select_elements.py
│           │       │   │   ├── simple.c
│           │       │   │   ├── simple.cpp
│           │       │   │   ├── simple.py
│           │       │   │   ├── spherical_surf.jl
│           │       │   │   ├── spherical_surf.py
│           │       │   │   ├── spline.cpp
│           │       │   │   ├── spline.py
│           │       │   │   ├── split_window.py
│           │       │   │   ├── square.cpp
│           │       │   │   ├── square.geo
│           │       │   │   ├── step_assembly.py
│           │       │   │   ├── step_boundary_colors.py
│           │       │   │   ├── step_boundary_colors.stp
│           │       │   │   ├── step_header_data.py
│           │       │   │   ├── step_header_data.stp
│           │       │   │   ├── stl_to_brep.py
│           │       │   │   ├── stl_to_mesh.py
│           │       │   │   ├── surface1.stl
│           │       │   │   ├── surface2.stl
│           │       │   │   ├── surface_filling.py
│           │       │   │   ├── terrain_bspline.jl
│           │       │   │   ├── terrain_bspline.py
│           │       │   │   ├── terrain.py
│           │       │   │   ├── terrain_stl_data.stl
│           │       │   │   ├── terrain_stl.py
│           │       │   │   ├── test.c
│           │       │   │   ├── test.py
│           │       │   │   ├── trimmed.py
│           │       │   │   ├── tube_boundary_layer.py
│           │       │   │   ├── view_adaptive_to_mesh.py
│           │       │   │   ├── view_combine.py
│           │       │   │   ├── view.cpp
│           │       │   │   ├── view_element_size.py
│           │       │   │   ├── viewlist.cpp
│           │       │   │   ├── viewlist.py
│           │       │   │   ├── view.py
│           │       │   │   ├── view_renumbering.py
│           │       │   │   ├── volume.py
│           │       │   │   └── x3d_export.py
│           │       │   ├── boolean
│           │       │   │   ├── as1-tu-203.stp
│           │       │   │   ├── baffles.geo
│           │       │   │   ├── boolean.geo
│           │       │   │   ├── chamfer.geo
│           │       │   │   ├── coherence.geo
│           │       │   │   ├── component8.step
│           │       │   │   ├── compsolid2.geo
│           │       │   │   ├── compsolid.geo
│           │       │   │   ├── extend_field.geo
│           │       │   │   ├── extrude2.geo
│           │       │   │   ├── extrude.geo
│           │       │   │   ├── fillet2.geo
│           │       │   │   ├── fillet3.geo
│           │       │   │   ├── fillet4.geo
│           │       │   │   ├── fillet_chamfer.geo
│           │       │   │   ├── fillet.geo
│           │       │   │   ├── fleur.geo
│           │       │   │   ├── fragment_numbering.geo
│           │       │   │   ├── hybrid_occ_builtin.geo
│           │       │   │   ├── hyperboloid.geo
│           │       │   │   ├── import2.geo
│           │       │   │   ├── import.geo
│           │       │   │   ├── intersect_line_volume.geo
│           │       │   │   ├── mesh_size_per_volume.geo
│           │       │   │   ├── neuron.geo
│           │       │   │   ├── number_of_tets.geo
│           │       │   │   ├── periodic_embedded.geo
│           │       │   │   ├── periodic.geo
│           │       │   │   ├── pipe.geo
│           │       │   │   ├── primitives.geo
│           │       │   │   ├── revolve2.geo
│           │       │   │   ├── revolve.geo
│           │       │   │   ├── shell_sewing.geo
│           │       │   │   ├── simple2.geo
│           │       │   │   ├── simple3.geo
│           │       │   │   ├── simple4.geo
│           │       │   │   ├── simple5.geo
│           │       │   │   ├── simple6.geo
│           │       │   │   ├── simple7.geo
│           │       │   │   ├── simple.geo
│           │       │   │   ├── slicer.geo
│           │       │   │   ├── slicer_surfaces.geo
│           │       │   │   ├── spherical_surf.geo
│           │       │   │   ├── spline.geo
│           │       │   │   ├── step_assembly.geo
│           │       │   │   ├── surface_filling.geo
│           │       │   │   ├── thicksolid.geo
│           │       │   │   ├── thrusections.geo
│           │       │   │   ├── transfinite.geo
│           │       │   │   ├── transform.geo
│           │       │   │   └── twist.geo
│           │       │   ├── post_processing
│           │       │   │   ├── anim.script
│           │       │   │   ├── compute_area_volume.geo
│           │       │   │   ├── encode.script
│           │       │   │   ├── isosurf.script
│           │       │   │   ├── lowmem-anim.geo
│           │       │   │   ├── multislice.script
│           │       │   │   ├── plot2d.geo
│           │       │   │   ├── primitives.pos
│           │       │   │   ├── right_scale_centered.geo
│           │       │   │   ├── rotate.script
│           │       │   │   ├── title.script
│           │       │   │   └── view_groups.geo
│           │       │   ├── simple_geo
│           │       │   │   ├── antenna.geo
│           │       │   │   ├── antenna.i1
│           │       │   │   ├── cone.geo
│           │       │   │   ├── cube.geo
│           │       │   │   ├── filter.geo
│           │       │   │   ├── hex.geo
│           │       │   │   ├── homology.geo
│           │       │   │   ├── indheat.geo
│           │       │   │   ├── machine.geo
│           │       │   │   ├── machine.i1
│           │       │   │   ├── machine.i2
│           │       │   │   ├── piece-extr.geo
│           │       │   │   ├── piece-extr-rec.geo
│           │       │   │   ├── piece.geo
│           │       │   │   ├── pripyrtet.geo
│           │       │   │   ├── sphere-discrete.geo
│           │       │   │   ├── sphere.geo
│           │       │   │   ├── sphere-surf.stl
│           │       │   │   ├── splines.geo
│           │       │   │   ├── square_regular.geo
│           │       │   │   ├── tower.geo
│           │       │   │   ├── tower.i1
│           │       │   │   ├── tower.i2
│           │       │   │   ├── tower.i3
│           │       │   │   ├── tower.i4
│           │       │   │   ├── tower.i5
│           │       │   │   └── transfinite.geo
│           │       │   └── struct
│           │       │       ├── Exists_GetForced.geo
│           │       │       └── struct.geo
│           │       ├── images
│           │       │   ├── gui_add_rectangle.png
│           │       │   ├── gui_mesh.png
│           │       │   ├── gui.png
│           │       │   ├── t10.png
│           │       │   ├── t11.png
│           │       │   ├── t12.png
│           │       │   ├── t13.png
│           │       │   ├── t14.png
│           │       │   ├── t15.png
│           │       │   ├── t16.png
│           │       │   ├── t17.png
│           │       │   ├── t18.png
│           │       │   ├── t19.png
│           │       │   ├── t1.png
│           │       │   ├── t20.png
│           │       │   ├── t21.png
│           │       │   ├── t2.png
│           │       │   ├── t3.png
│           │       │   ├── t4.png
│           │       │   ├── t5.png
│           │       │   ├── t6.png
│           │       │   ├── t7.png
│           │       │   ├── t8.png
│           │       │   ├── t9.png
│           │       │   ├── x1.png
│           │       │   ├── x2.png
│           │       │   ├── x3.png
│           │       │   ├── x4.png
│           │       │   ├── x5.png
│           │       │   └── x7.png
│           │       ├── LICENSE.txt
│           │       ├── README.txt
│           │       └── tutorials
│           │           ├── c
│           │           │   ├── README.txt
│           │           │   ├── t16.c
│           │           │   ├── t1.c
│           │           │   ├── t2.c
│           │           │   └── t6.c
│           │           ├── c++
│           │           │   ├── README.txt
│           │           │   ├── t10.cpp
│           │           │   ├── t11.cpp
│           │           │   ├── t12.cpp
│           │           │   ├── t13.cpp
│           │           │   ├── t14.cpp
│           │           │   ├── t15.cpp
│           │           │   ├── t16.cpp
│           │           │   ├── t17.cpp
│           │           │   ├── t18.cpp
│           │           │   ├── t19.cpp
│           │           │   ├── t1.cpp
│           │           │   ├── t20.cpp
│           │           │   ├── t21.cpp
│           │           │   ├── t2.cpp
│           │           │   ├── t3.cpp
│           │           │   ├── t4.cpp
│           │           │   ├── t5.cpp
│           │           │   ├── t6.cpp
│           │           │   ├── t7.cpp
│           │           │   ├── t8.cpp
│           │           │   ├── t9.cpp
│           │           │   ├── x1.cpp
│           │           │   ├── x2.cpp
│           │           │   ├── x3.cpp
│           │           │   ├── x4.cpp
│           │           │   ├── x5.cpp
│           │           │   ├── x6.cpp
│           │           │   └── x7.cpp
│           │           ├── fortran
│           │           │   ├── README.txt
│           │           │   ├── t10.f90
│           │           │   ├── t11.f90
│           │           │   ├── t12.f90
│           │           │   ├── t13.f90
│           │           │   ├── t14.f90
│           │           │   ├── t15.f90
│           │           │   ├── t16.f90
│           │           │   ├── t17.f90
│           │           │   ├── t18.f90
│           │           │   ├── t19.f90
│           │           │   ├── t1.f90
│           │           │   ├── t20.f90
│           │           │   ├── t21.f90
│           │           │   ├── t2.f90
│           │           │   ├── t3.f90
│           │           │   ├── t4.f90
│           │           │   ├── t5.f90
│           │           │   ├── t6.f90
│           │           │   ├── t7.f90
│           │           │   ├── t8.f90
│           │           │   └── t9.f90
│           │           ├── julia
│           │           │   ├── README.txt
│           │           │   ├── t10.jl
│           │           │   ├── t11.jl
│           │           │   ├── t12.jl
│           │           │   ├── t13.jl
│           │           │   ├── t14.jl
│           │           │   ├── t15.jl
│           │           │   ├── t16.jl
│           │           │   ├── t17.jl
│           │           │   ├── t18.jl
│           │           │   ├── t19.jl
│           │           │   ├── t1.jl
│           │           │   ├── t20.jl
│           │           │   ├── t21.jl
│           │           │   ├── t2.jl
│           │           │   ├── t3.jl
│           │           │   ├── t4.jl
│           │           │   ├── t5.jl
│           │           │   ├── t6.jl
│           │           │   ├── t7.jl
│           │           │   ├── t8.jl
│           │           │   ├── t9.jl
│           │           │   ├── x1.jl
│           │           │   ├── x2.jl
│           │           │   ├── x3.jl
│           │           │   ├── x4.jl
│           │           │   ├── x5.jl
│           │           │   ├── x6.jl
│           │           │   └── x7.jl
│           │           ├── python
│           │           │   ├── README.txt
│           │           │   ├── t10.py
│           │           │   ├── t11.py
│           │           │   ├── t12.py
│           │           │   ├── t13.py
│           │           │   ├── t14.py
│           │           │   ├── t15.py
│           │           │   ├── t16.py
│           │           │   ├── t17.py
│           │           │   ├── t18.py
│           │           │   ├── t19.py
│           │           │   ├── t1.py
│           │           │   ├── t20.py
│           │           │   ├── t21.py
│           │           │   ├── t2.py
│           │           │   ├── t3.py
│           │           │   ├── t4.py
│           │           │   ├── t5.py
│           │           │   ├── t6.py
│           │           │   ├── t7.py
│           │           │   ├── t8.py
│           │           │   ├── t9.py
│           │           │   ├── x1.py
│           │           │   ├── x2.py
│           │           │   ├── x3.py
│           │           │   ├── x4.py
│           │           │   ├── x5.py
│           │           │   ├── x6.py
│           │           │   └── x7.py
│           │           ├── README.txt
│           │           ├── t10.geo
│           │           ├── t11.geo
│           │           ├── t12.geo
│           │           ├── t13_data.stl
│           │           ├── t13.geo
│           │           ├── t14.geo
│           │           ├── t15.geo
│           │           ├── t16.geo
│           │           ├── t17_bgmesh.pos
│           │           ├── t17.geo
│           │           ├── t18.geo
│           │           ├── t19.geo
│           │           ├── t1.geo
│           │           ├── t20_data.step
│           │           ├── t20.geo
│           │           ├── t21.geo
│           │           ├── t2.geo
│           │           ├── t3.geo
│           │           ├── t4.geo
│           │           ├── t4_image.png
│           │           ├── t5.geo
│           │           ├── t6.geo
│           │           ├── t7_bgmesh.pos
│           │           ├── t7.geo
│           │           ├── t8.geo
│           │           ├── t9.geo
│           │           ├── view1.pos
│           │           ├── view2.pos
│           │           ├── view3.pos
│           │           ├── view4.pos
│           │           └── view5.msh
│           └── man
│               └── man1
│                   └── gmsh.1
├── plot.py
├── project
│   ├── CMakeLists.txt
│   ├── example.c
│   ├── fem.c
│   ├── fem.h
│   └── homework.c
├── README.md
└── validate.py