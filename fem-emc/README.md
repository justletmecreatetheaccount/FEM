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

Ci dessous la structure du projet pour un programme qui à déjà tourner et avec donc les fichiers de build et mesh
```
.
├── main.c                      <- the main program
├── main.h
├── build                       <- build folder (to create)
├── CMakeLists.txt              <- main cmake configuration
├── data                        <- meshes folder
│   ├── fixMesh.py              <- python program to fix meshes from gmsh
│   ├── mesh_ref.txt            <- a reference mesh for another geometry
│   ├── mesh.txt                <- the mesh fixed by fixMesh.py
│   ├── UV_ref.txt              <- the solution to the reference mesh
│   └── UV.txt                  <- the solution to mesh.txt problem
├── genmesh                     <- folder containing the mesh generator
│   ├── build                   <- build folder (to create)  
│   ├── CMakeLists.txt
│   ├── data                    
│   │   └── mesh.txt            <- the generated mesh
│   ├── glfw                    <- openGL interface
│   └── src                     
│       ├── fem.c               <- utils and mesh generator
│       ├── fem.h               
│       ├── glfem.c             <- utils for window creation
│       ├── glfem.h
│       └── main.c              <- core program to generate mesh
├── gmsh                        <- the gmsh folder used to generate mesh
├── plot.py                     <- python program to plot results
├── project                     <- folder containting the solver
│   ├── CMakeLists.txt
│   ├── solve.c                 <- the solver implementation
│   ├── fem.c                   <- utils
│   ├── fem.h
│   └── assemble.c              <- assemble the problem and matrices
├── README.md                   <- you are here
└── validate.py                 <- python program to compare meshes

```

### Initialisation

> ⚠️ **Attention**  

Pour initialiser le projet, il faut créer les dossiers de build dans `./build/` et `./genmesh/build/`.

Si vous êtes sur Linux, vous pouvez lancer `./_init.sh` (qui est l'équivalent) après avoir fait 
`chmod +x _init.sh` (vous pouvez le faire pour tout les fichiers .sh)

## Lancer le programme

### Générer un maillage

Pour générer un maillage il faut d'abords naviguer dans `./genmesh/build/` ensuite executer
```sh
cmake ..
```
ce qui va générer les fichiers de build, ensuite on compile le programme de génération avec
```sh
make
```
et enfin, on peut faire tourner l'exécutable avec 
```
./genmesh
```
L'exécutable peut également prendre deux paramètres en entrée (`quad`/`tri` et `plot`) qui nous disent si l'exécutable dois utiliser des quads ou non, le programme choisis les quads par défaut, exemple : `./genmesh tri` pour des triangles. Par défaut les quads sont choisis.

Le deuxième paramètre `plot` est pour visualiser ou non le champ de la taille des éléments par exemple : `./genmesh plot` pour ouvrir une fenetre avec des quads

Par défaut pas de plot est fait et les quads sont choisis, si `quad` et `tri` sont passés en argument, le programme prend le dernier choisis, `./genmesh quad tri` ferra des triangles

Si vous êtes sur Linux, vous pouvez lancer `./_mesh.sh` pour générer le maillage avec un plot du champ de taille des éléments ainsi que visualiser le maillage, c'est équivalent à utiliser les commandes ci-dessus


### Résoudre le problème


Les paramètres du problème sont définis dans `main.h` par
```
#define _DISPLACEMENT (2.0 * 10e-7)
#define _E (170.0 * 10e9)
#define _NU (0.29)
#define _RHO (7.5 * 10e3)
#define _G (9.81)
```
Ou `_DISPLACEMENT` est le deplacement par compression des bords du disque et les autres constantes associées sont les paramètres physiques trouver sur https://www.makeitfrom.com/material-properties/SAE-ASTM-Grade-G3500-F10007-Grey-Cast-Iron et https://en.wikipedia.org/wiki/Young%27s_modulus#Temperature_dependence

On peut donc facilement les changer pour adapter l'analyse a un autre matériaux par exemple un acier froid comparé a un acier chaud.

Pour résoudre le problème, et afficher les déformations subies, il faut naviguer dans `./build/`
puis exécuter
```
cmake ..
```
suivis par 
```
make
```
et finalement
```
./solve
```
le binaire solve peut prendre comme argument le chemin du fichier d'entrée et de sortie (dans cet ordre)

Si vous êtes sur Linux, vous pouvez lancer `./_run.sh` c'est équivalent à utiliser les commandes ci-dessus, cela va donc resoudre le probleme definis dans `solve.c` et faire une animation avec les déformations

Concernant la fenêtre ouverte, elle permet aussi de visualiser d'autre chose en appuyant sur la touche correspondante

| Touche appuyée     | Figure affichée                          |
|----------------------|--------------------------------------|
| `D`          | Permet d'iterer sur les domaines des frontières en appuyant sur `N`         |
| `V`  | Vue par défaut du maillage générer      |
| `X`           | Vue des forces sur les domaines frontières en X    |
| `Y`   | Vue des forces sur les domaines frontières en Y       |

En executant `./solve` un benchmark apparait dans les logs, c'est le temps utiliser pour resoudre le système
