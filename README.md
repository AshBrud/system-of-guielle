# Système de Guielle ⚙️

> **Modélisation et simulation numérique d'un système mécanique élastique à rotation continue**
> *Modeling and numerical simulation of a continuous-rotation elastic mechanical system*

**Auteur / Author:** Herolde GUIELLE
**Date:** Mars 2025
**Langue / Language:** Français / English (bilingual)

---

## Table des matières / Table of Contents

- [Résumé / Abstract](#résumé--abstract)
- [Description du système / System Description](#description-du-système--system-description)
- [Modélisation mathématique / Mathematical Modeling](#modélisation-mathématique--mathematical-modeling)
- [Triple condition de mouvement / Triple Motion Condition](#triple-condition-de-mouvement--triple-motion-condition)
- [Simulation numérique / Numerical Simulation](#simulation-numérique--numerical-simulation)
- [Résultats / Results](#résultats--results)
- [Installation et utilisation / Installation & Usage](#installation-et-utilisation--installation--usage)
- [Structure du dépôt / Repository Structure](#structure-du-dépôt--repository-structure)
- [Nomenclature](#nomenclature)
- [Licence / License](#licence--license)

---

## Résumé / Abstract

**FR —** Le Système de Guielle est un mécanisme mécanique élastique constitué d'une roue de rayon $r$ et de masse $m$, entraînée par les forces de rappel de $n$ ressorts de raideur $k$, tous ancrés en un point unique $P$ de la roue, leurs extrémités fixes $A_i$ étant attachées à un support immobile extérieur. L'objectif est de démontrer que ce système peut générer un couple mécanique continu à partir de l'énergie élastique des ressorts, en vérifiant une triple condition sur le moment résultant adimensionnel $g(t)$.

**EN —** The Guielle System is an elastic mechanical device consisting of a wheel of radius $r$ and mass $m$, driven by the restoring forces of $n$ springs of stiffness $k$, all attached at a single point $P$ on the wheel, with their fixed ends $A_i$ anchored to an immobile external support. The goal is to demonstrate that this system can generate continuous mechanical torque from the elastic energy stored in the springs, by verifying a triple condition on the dimensionless resultant moment $g(t)$.

---

## Description du système / System Description

```
       Aᵢ (fixe / fixed)
        ●
         \   ← ressort i / spring i (raideur k / stiffness k)
          \
    O ────P  ← point d'accroche unique / single attachment point
    │
    └── roue de rayon r, masse m / wheel radius r, mass m
```

- $n$ ressorts identiques, extrémités fixes $A_i$ réparties angulairement
- Tous les ressorts sont attachés au **même point $P$** de la roue
- Les $A_i$ sont **strictement fixes** dans le référentiel du laboratoire
- La roue pivote autour de son centre $O$

---

## Modélisation mathématique / Mathematical Modeling

### 1. Géométrie du système / System Geometry

Les $n$ points fixes $A_i$ sont positionnés à distance $r \cdot \varepsilon$ du centre, avec l'angle fixe $\beta_i$ dans le référentiel du laboratoire :

$$
\begin{aligned}
\beta_i &= i \cdot \frac{2\pi}{n}, \quad i = 1, 2, \ldots, n
    && \text{(angle fixe dans le labo / fixed lab angle)} \\
\varepsilon &= 1 + \frac{l_0}{r}
    && \text{(rapport géométrique / geometric ratio)} \\
A_i &= O + r \cdot \varepsilon \cdot \begin{pmatrix} \cos \beta_i \\ \sin \beta_i \end{pmatrix}
    && \text{(position fixe de } A_i\text{)} \\
\tau_i &= \beta_i + \theta
    && \text{(angle caractéristique / characteristic angle)}
\end{aligned}
$$

Le point $P$, mobile sur la roue, est repéré par l'angle de rotation $\theta(t)$ :

$$
P(t) = O + r \begin{pmatrix} \cos\theta \\ \sin\theta \end{pmatrix}
$$

---

### 2. Élongation des ressorts / Spring Extension

Par le théorème de Pythagore appliqué au triangle $A_i O P$ :

$$
\left( x_i + l_0 \right)^2 = r^2 + (r + l_0)^2 + 2r(r + l_0)\cos\tau_i
$$

En posant $\varepsilon = 1 + l_0/r$, la forme adimensionnelle devient :

$$
\left( \frac{x_i + l_0}{r} \right)^2 = 1 + \varepsilon^2 + 2\varepsilon \cos\tau_i
$$

D'où l'**élongation adimensionnelle** :

$$
\boxed{
\frac{x_i}{r} = \sqrt{1 + \varepsilon^2 + 2\varepsilon\cos\tau_i} \;-\; \frac{l_0}{r}
}
$$

> Le ressort $i$ contribue **uniquement si** $x_i > 0$ (extension active).  
> Spring $i$ contributes **only if** $x_i > 0$ (under active extension).

---

### 3. Direction de la force de rappel / Spring Force Direction

L'angle $\alpha_i$ découle de la relation géométrique $HP = (x_i + l_0)\sin\alpha_i = r\sin\tau_i$ :

$$
\boxed{
\alpha_i = \arcsin\!\left( \frac{\sin\tau_i}{\sqrt{1 + \varepsilon^2 + 2\varepsilon\cos\tau_i}} \right)
}
$$

L'angle orienté entre $\hat{e}_\rho$ et $\hat{e}_{x_i}$ vaut alors :

$$
\left(\widehat{\hat{e}_\rho,\, \hat{e}_{x_i}}\right) = -\theta + \beta_i + \alpha_i
$$

---

### 4. Force de rappel et moment résultant / Restoring Force and Resultant Moment

La force de rappel du ressort $i$ est :

$$
\vec{F}_{r_i} = -k \, x_i \, \hat{e}_{x_i}
$$

Le moment de la force résultante par rapport à $O$ est :

$$
\mathcal{M}_{\vec{F}_r/O} = \vec{OP} \times \vec{F}_r
= -r\,\hat{e}_\rho \times \sum_{i=1}^{n} k\, x_i\, \hat{e}_{x_i}
$$

En développant le produit vectoriel $\hat{e}_\rho \times \hat{e}_{x_i} = \sin\!\left(\widehat{\hat{e}_\rho, \hat{e}_{x_i}}\right)\hat{e}_z$ :

$$
\boxed{
\mathcal{M}_{\vec{F}_r/O}(\theta)
= -r^2 k \sum_{i=1}^{n} \left[ \frac{x_i}{r} \cdot \sin\!\left(\theta - \beta_i - \alpha_i\right) \right]
}
$$

La **fonction adimensionnelle centrale** de l'étude est :

$$
\boxed{
g(t) = -\frac{\mathcal{M}}{r^2 k}
= \sum_{i=1}^{n} \left[ \frac{x_i}{r} \cdot \sin\!\left(\theta - \beta_i - \alpha_i\right) \right]
}
$$

---

### 5. Hypothèse de mouvement uniforme / Uniform Motion Assumption

Dans le document original, l'étude analytique de $g(t)$ est conduite en posant $\theta = \omega t$ :

$$
\begin{aligned}
g(t) &= \sum_{i=1}^{n} \left[ \frac{x_i}{r} \cdot \sin\!\left(\omega t - \mu_i\right) \right] \\
\frac{x_i}{r} &= \sqrt{1 + \varepsilon^2 + 2\varepsilon\cos(\omega t + \beta_i)} - \frac{l_0}{r} \\
\alpha_i &= \arcsin\!\left( \frac{\sin(\omega t + \beta_i)}{\sqrt{1 + \varepsilon^2 + 2\varepsilon\cos(\omega t + \beta_i)}} \right)
\end{aligned}
$$

avec $\mu_i = \beta_i + \alpha_i$.

> La simulation RK4 lève cette hypothèse en résolvant $\theta(t)$ dynamiquement.

---

### 6. Équation du mouvement / Equation of Motion

Le moment d'inertie est **paramétré librement** par le coefficient $\alpha$ :

$$
I = \alpha \cdot m \cdot r^2, \qquad \alpha \in [0.1,\; 1]
$$

Le **théorème du moment dynamique** donne l'équation différentielle du second ordre :

$$
\boxed{
I\,\ddot{\theta} = \mathcal{M}_{\vec{F}_r/O}(\theta) \;-\; b\,\dot{\theta} \;-\; C_s\,\mathrm{sgn}(\dot{\theta})
}
$$

où $b$ est le coefficient de frottement visqueux et $C_s$ le couple résistant externe.  
Dans le cas idéal du document : $b = 0$ et $C_s = 0$.

---

### 7. Bilan énergétique / Energy Balance

$$
\begin{aligned}
E_c &= \frac{1}{2} I \omega^2 = \frac{1}{2}\alpha m r^2 \omega^2
    && \text{(énergie cinétique / kinetic energy)} \
E_e &= \sum_{i=1}^{n} \frac{1}{2} k\, x_i^2 \quad (x_i > 0)
    && \text{(énergie élastique / elastic energy)} \\
P_f &= b\,\omega^2
    && \text{(puissance frottement / friction power)} \\
P_s &= C_s\,|\omega|
    && \text{(puissance résistante / load power)} \\
\dot{E}_c &= \mathcal{M}\cdot\omega - P_f - P_s
    && \text{(bilan / energy balance)}
\end{aligned}
$$

---

## Triple condition de mouvement / Triple Motion Condition

La roue pivote **indéfiniment** si et seulement si $g(t)$ vérifie **simultanément** pour tout $t \geq 0$ :

$$
\begin{cases}
\textbf{(C1)} \quad |g(t)| > 0
    & \text{moment moteur non nul en tout instant} \\
\textbf{(C2)} \quad g(t) \approx \text{constante}
    & \text{mouvement tend vers un régime uniforme} \\
\textbf{(C3)} \quad \mathrm{sgn}(g(t)) = \text{invariable}
    & \text{la roue ne peut pas s'inverser}
\end{cases}
$$

Formellement :

$$
\begin{cases}
\left| \mathcal{M}_{\vec{F}_r/O}(t) \right| > 0 \\
\mathcal{M}_{\vec{F}_r/O}(t + \Delta t) - \mathcal{M}_{\vec{F}_r/O}(t) \approx 0
\end{cases}
\quad \forall\, t \in \mathbb{R}^+
$$

### Effets des paramètres / Parameter Effects

| Paramètre | Effet sur la triple condition |
|---|---|
| $n \uparrow$ | $g(t)$ plus stable, oscillations réduites, (C2) mieux vérifiée |
| $k \uparrow$ | Moments plus grands, (C1) plus robuste |
| $l_0 > r$ | Assure plusieurs ressorts en extension simultanée |
| $\alpha \downarrow$ | Meilleur dépassement des positions d'équilibre |
| $b > 0$ | Décroissance de $\omega$ si $\mathcal{M} \ll b\,\omega$ |

---

## Simulation numérique / Numerical Simulation

### Méthode RK4 / RK4 Method

L'équation du second ordre est réduite à un système du premier ordre avec l'état $\mathbf{s} = (\theta,\, \omega)$ :

$$
\frac{d\mathbf{s}}{dt} = f(\mathbf{s})
= \begin{pmatrix} \omega \\ 
\dfrac{\mathcal{M}(\theta) - b\,\omega - C_s\,\mathrm{sgn}(\omega)}{I} \end{pmatrix}
$$

Le schéma **Runge-Kutta d'ordre 4** à chaque pas $\Delta t$ :

$$
\begin{aligned}
k_1 &= f(\mathbf{s}_n) \\
k_2 &= f\!\left(\mathbf{s}_n + \tfrac{\Delta t}{2}\,k_1\right) \\
k_3 &= f\!\left(\mathbf{s}_n + \tfrac{\Delta t}{2}\,k_2\right) \\
k_4 &= f\!\left(\mathbf{s}_n + \Delta t\, k_3\right)
\end{aligned}
$$

$$
\boxed{
\mathbf{s}_{n+1} = \mathbf{s}_n + \frac{\Delta t}{6}\left(k_1 + 2k_2 + 2k_3 + k_4\right)
}
$$

avec $\Delta t = 0.012\,\text{s}$.

### Configurations testées / Tested Configurations

| Config | $n$ | $r$ (m) | $k$ (N/m) | $l_0$ (m) | $m$ (kg) | $\alpha$ | Triple condition |
|---|---|---|---|---|---|---|---|
| C1 | 2 | 0.30 | 30 | 0.40 | 5 | 0.50 | (C1)✓ (C2)~ (C3)✓ |
| C2 | 4 | 0.30 | 30 | 0.40 | 5 | 0.50 | (C1)✓ (C2)✓ (C3)✓ |
| C3 | 6 | 0.30 | 50 | 0.50 | 8 | 0.40 | (C1)✓ (C2)✓ (C3)✓ |
| C4 | 8 | 0.40 | 80 | 0.60 | 10 | 0.50 | (C1)✓ (C2)✓ (C3)✓ |

> **Observation principale :** $n \geq 4$ avec $l_0 > r$ est la configuration minimale garantissant la triple condition simultanée.

---

## Résultats / Results

- Moment moteur non nul et de signe constant confirmé pour $n \geq 3$
- Constance de $g(t)$ croissante avec $n$
- La configuration (tous les ressorts en $P$) crée une **asymétrie dynamique favorable**
- $n \geq 4$ avec $l_0 > r$ est la configuration minimale recommandée

---

## Installation et utilisation / Installation & Usage

### Prérequis / Requirements

```bash
pip install streamlit numpy scipy matplotlib plotly pandas
```

### Lancer le simulateur / Run the simulator

```bash
streamlit run guielle_simulator.py
```

L'application s'ouvre dans le navigateur à `http://localhost:8501`.

### Fonctionnalités / Features

- Sliders interactifs pour tous les paramètres : $n$, $r$, $k$, $l_0$, $m$, $\alpha$, $b$, $C_s$, $\omega_0$
- Intégration RK4 avec pas $\Delta t$ configurable
- Graphiques : $\omega(t)$, $g(t)$, $E_c(t)$, $E_e(t)$, portrait de phase $\omega(\theta)$
- Triple condition (C1)(C2)(C3) évaluée quantitativement en temps réel
- Analyse comparative $n = 2$ à $8$ en un clic
- Export des données au format CSV

---

## Structure du dépôt / Repository Structure

```
systeme-de-guielle/
├── README.md                          # Ce fichier / This file
├── LICENSE                            # Licence CC BY 4.0
└── guielle_simulator.py               # Simulateur Streamlit principal
```

---

## Nomenclature

| Symbole | Définition |
|---|---|
| $n$ | Nombre de ressorts / Number of springs |
| $r$ | Rayon de la roue (m) |
| $m$ | Masse de la roue (kg) |
| $k$ | Constante de raideur (N/m) |
| $l_0$ | Longueur naturelle des ressorts (m) |
| $\alpha$ | Coefficient d'inertie libre |
| $I = \alpha m r^2$ | Moment d'inertie (kg·m²) |
| $\theta$ | Angle de rotation (rad) |
| $\omega = \dot{\theta}$ | Vitesse angulaire (rad/s) |
| $\beta_i$ | Angle fixe du point $A_i$ (rad) |
| $\tau_i = \beta_i + \theta$ | Angle caractéristique |
| $\varepsilon = 1 + l_0/r$ | Rapport géométrique |
| $x_i$ | Élongation du ressort $i$ (m) |
| $g(t)$ | Moment adimensionnel |
| $b$ | Coefficient de frottement visqueux |
| $C_s$ | Couple résistant externe (N·m) |

---

## Licence / License

Ce projet est distribué sous licence **Creative Commons Attribution 4.0 International (CC BY 4.0)**.

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

**© 2025 Herolde GUIELLE** — Le concept mécanique du Système de Guielle, sa modélisation et ses formules constituent une œuvre intellectuelle originale. Toute publication académique ou commerciale utilisant ce travail doit citer l'auteur.

> *Any academic or commercial publication using this work must cite the author: Herolde GUIELLE (2025).*

---

*Système de Guielle · Rapport de Recherche · H. GUIELLE · 2025*
