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







