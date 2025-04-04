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

- Cmake
- X11 (Linux)
- tdm-GCC ou GCC
- GLFW et OpenGL
- GMSH

### GMSH
Pour utiliser correctement gmsh, il faut télécharger la [version](https://gmsh.info/) correspondant à votre OS et la mettre dans `genmesh/gmsh` pour avoir 
`genmesh/gmsh/gmsh-4.13.1-.../...`

En cas de problème de linkage avec gmsh, veuillez vous référrez aux CMakeLists.txt ou le chemin utilisé est specifier. C'est primordial de bien installer gmsh, sans ça le programme de génération de maillage ne pourra pas fonctionner


### Structure

Ci dessous la structure du projet pour un programme qui à déjà tourner et avec donc les fichiers de build et mesh
```
.
├── build                       <- build folder (to create)
├── _clean.sh
├── CMakeLists.txt
├── data                        <- meshes folder
│   ├── mesh_ref.txt
│   ├── mesh.txt                <- the mesh (generated)
│   └── UV_ref.txt
├── genmesh                     <- folder containing the mesh generator
│   ├── build                   <- build folder (to create)  
│   ├── CMakeLists.txt
│   ├── glfw                    <- openGL interface for showing the mesh
│   ├── gmsh                    <- the gmsh folder used to generate mesh
│   └── src
│       ├── fem.c               <- utils and mesh generator
│       ├── fem.h
│       ├── glfem.c             <- utils for window creation (opengl)
│       ├── glfem.h
│       └── main.c              <- core program to generate mesh
├── _init.sh
├── main.c                      <- the main program
├── main.h
├── _mesh.sh
├── project                     <- folder containting the solver
│   ├── assemble.c              <- assemble the problem and matrices
│   ├── CMakeLists.txt
│   ├── fem.c                   <- utils and core of the solver
│   ├── fem.h
│   └── solve.c                 <- the problem definition and solver
├── README.md                   <- you are here
├── _run.sh
└── utils
    ├── fixMesh.py              <- python program to fix meshes from gmsh
    ├── plot.py                 <- python program to plot results
    └── validate.py             <- python program to compare solutions
```

### Initialisation

> ⚠️ **Attention**  

Pour initialiser le projet, il faut créer les dossiers de build dans `build/` et `genmesh/build/`.

Si vous êtes sur Linux, vous pouvez lancer `./_init.sh` (qui est l'équivalent) après avoir fait 
`chmod +x _init.sh` (vous pouvez le faire pour tout les fichiers .sh)

## Lancer le programme

### Générer un maillage

Pour générer un maillage il faut d'abords naviguer dans `genmesh/build/` ensuite executer
```sh
cmake ..
```
ce qui va générer les fichiers de build, ensuite on compile le programme de génération avec
```sh
make
```
et enfin, on peut faire tourner l'exécutable avec 
```sh
./genmesh
```
L'exécutable peut également prendre deux paramètres en entrée (`quad`/`tri` et `plot`) qui nous disent si l'exécutable dois utiliser des quads ou non, le programme choisis les quads par défaut, exemple : `./genmesh tri` pour des triangles. Par défaut les quads sont choisis.

Le deuxième paramètre `plot` est pour visualiser ou non le champ de la taille des éléments par exemple : `./genmesh plot` pour ouvrir une fenetre avec des quads

Par défaut pas de plot est fait et les quads sont choisis, si `quad` et `tri` sont passés en argument, le programme prend le dernier choisis, `./genmesh quad tri` ferra des triangles

Si vous êtes sur Linux, vous pouvez lancer `./_mesh.sh` pour générer le maillage avec un plot du champ de taille des éléments ainsi que visualiser le maillage, c'est équivalent à utiliser les commandes ci-dessus


### Résoudre le problème


Les paramètres du problème sont définis dans `main.h` par
```
#define _FACTOR (1e5)
#define _DISPLACEMENT (2.0 * 1e-6)
#define _E (170.0 * 1e9)
#define _NU (0.29)
#define _RHO (7.5 * 1e3)
#define _G (9.81)
```
Ou `_DISPLACEMENT` est le deplacement par compression des bords du disque et les autres constantes associées sont les paramètres physiques trouver sur https://www.makeitfrom.com/material-properties/SAE-ASTM-Grade-G3500-F10007-Grey-Cast-Iron et https://en.wikipedia.org/wiki/Young%27s_modulus#Temperature_dependence 

`_FACTOR` est le facteur de féformation utilisé pour amplifier les résultats dans les plots, si celui-ci est changé, faites attention a également le changer dans `utils/plot.py` pour que le plot python corresponde a la figure openGL


On peut donc facilement les changer pour adapter l'analyse a un autre matériaux par exemple un acier froid comparé a un acier chaud.

Pour résoudre le problème, et afficher les déformations subies, il faut naviguer dans `build/`
puis exécuter
```sh
cmake ..
```
suivis par 
```sh
make
```
et finalement
```sh
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

## Code

Cette section est pour détaillé les parties de code plus intéréssantes à voir

### Génération de Mesh

Pour générer le mesh, le code utiliser est inspiré du devoir 2 et 6, avec cependant une geométrie différente. La géométrie est définie dans `genmesh/src/fem.c` par la fonction `void geoMeshGenerate()`, avec la fonction  `double geoSize(double x, double y)` associée au problème

### Solveur

Pour résoudre le problème il y a deux options, soit le stockage de la stiffness matrix en matrice creuse, le code provient en grosse partie du devoir 6 et du template fournis pour le concours avec quelques ajustements.
Soit le stockage en compressed sparse row, les fonctions adoubées d'un csr sont fonctionnelement les exactes répliques de leurs contreparties "normales" si ce n'est qu'elles peuvent intéragir avec des matrices CSR.
Vous remarquerez aussi l'ajout d'un solver par Méthode du gradient conjugué avec ou sans préconditionnement.
Les modifications se trouvent dans fem.c.
Les conditions frontières sont définies dans `project/solve.c`
tandis que le solveur est défini dans `project/assemble.c` et est appelé par `double *femElasticitySolve(femProblem *theProblem)`

### Plots

Pour ce projet, plusieurs fenêtres sont créées, le code de ces dernières provient des devoirs et a été ajusté pour le probleme et les prérequis


## Notes générales

- Si en changeant le mesh vous tomber sur une erreur de pivot (ce n'est pas toujours le cas), c'est qu'il faut réparer le mesh en exécutant
```sh
python3 utils/fixMesh.py
```
- Il existe aussi un programme python pour visualiser les resultats avec son propre facteur de déformations, il suffit de run
```sh
python3 utils/plot.py
```
- Pour éviter le code inutile, les fonctions non-utilisées des fichiers fem.c et fem.h ont été retirées

- Tout le code a été testé sur Linux mais pas MacOs ou Windows si des erreurs de linkage arrive, il va peut être falloir trifouiller les CMakeLists, même si ils ont été faits pour ne pas que cela arrive
