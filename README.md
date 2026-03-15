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
- [Licence / License](#licence--license)

---

## Résumé / Abstract

**FR —** Le Système de Guielle est un mécanisme mécanique élastique constitué d'une roue de rayon *r* et de masse *m*, entraînée par les forces de rappel de *n* ressorts de raideur *k*, tous ancrés en un point unique *P* de la roue, leurs extrémités fixes *Aᵢ* étant attachées à un support immobile extérieur. L'objectif est de démontrer que ce système peut générer un couple mécanique continu à partir de l'énergie élastique des ressorts, en vérifiant une triple condition sur le moment résultant adimensionnel *g(t)*.

**EN —** The Guielle System is an elastic mechanical device consisting of a wheel of radius *r* and mass *m*, driven by the restoring forces of *n* springs of stiffness *k*, all attached at a single point *P* on the wheel, with their fixed ends *Aᵢ* anchored to an immobile external support. The goal is to demonstrate that this system can generate continuous mechanical torque from the elastic energy stored in the springs, by verifying a triple condition on the dimensionless resultant moment *g(t)*.

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

- **n** ressorts identiques, extrémités fixes *Aᵢ* réparties angulairement : `βᵢ = i · 2π/n`
- Tous les ressorts sont attachés au **même point P** de la roue
- Les *Aᵢ* sont **strictement fixes** dans le référentiel du laboratoire
- La roue pivote autour de son centre *O*

---

## Modélisation mathématique / Mathematical Modeling

### Géométrie / Geometry

Les points fixes *Aᵢ* sont positionnés à distance `r·ε` du centre :

```
βᵢ = i · (2π/n)        (angle fixe dans le labo / fixed lab angle)
ε  = 1 + l₀/r          (rapport géométrique / geometric ratio)
Aᵢ = O + r·ε·(cos βᵢ,  sin βᵢ)
τᵢ = βᵢ + θ             (angle caractéristique / characteristic angle)
```

### Élongation des ressorts / Spring Extension

Par le théorème de Pythagore appliqué au triangle *AᵢOP* :

```
( (xᵢ + l₀)/r )² = 1 + ε² + 2ε·cos τᵢ

xᵢ/r = √(1 + ε² + 2ε·cos τᵢ)  −  l₀/r
```

> Le ressort contribue **uniquement si** `xᵢ > 0` (en extension).
> The spring contributes **only if** `xᵢ > 0` (under extension).

### Direction de la force de rappel / Spring Force Direction

```
αᵢ = arcsin( sin τᵢ / √(1 + ε² + 2ε·cos τᵢ) )

(ê_ρ, ê_{xᵢ}) = −θ + βᵢ + αᵢ
```

### Moment résultant / Resultant Moment

```
M(θ) = −r²·k · Σᵢ [ (xᵢ/r) · sin(θ − βᵢ − αᵢ) ]

g(t) = −M / (r²·k) = Σᵢ [ (xᵢ/r) · sin(θ − βᵢ − αᵢ) ]
```

### Équation du mouvement / Equation of Motion

Le moment d'inertie est paramétré librement (`α` dépend de la géométrie réelle de la roue) :

```
I = α · m · r²,    α ∈ [0.1, 1]
```

Équation différentielle du second ordre (théorème du moment dynamique) :

```
┌─────────────────────────────────────────────────────────┐
│  I · θ̈ = M(θ) − b · θ̇ − Cs · sgn(θ̇)                  │
└─────────────────────────────────────────────────────────┘
```

où `b` est le frottement visqueux et `Cs` le couple résistant externe.

### Bilan énergétique / Energy Balance

| Grandeur | Expression |
|---|---|
| Énergie cinétique / Kinetic energy | `Ec = ½ · I · ω²` |
| Énergie élastique / Elastic energy | `Ee = Σᵢ ½·k·xᵢ²  (xᵢ > 0)` |
| Puissance frottement / Friction power | `Pf = b · ω²` |
| Puissance résistante / Load power | `Ps = Cs · |ω|` |
| Bilan / Balance | `Ėc = M·ω − Pf − Ps` |

---

## Triple condition de mouvement / Triple Motion Condition

La roue pivote **indéfiniment** si et seulement si la fonction *g(t)* vérifie **simultanément** :

| # | Condition | Signification |
|---|---|---|
| (C1) | `\|g(t)\| > 0` | Moment moteur non nul en tout instant |
| (C2) | `g(t) ≈ constante` | Mouvement tend vers un régime uniforme |
| (C3) | `sgn(g(t))` invariable | La roue ne peut pas s'inverser |

> Dans le document original, l'étude analytique suppose `θ = ωt` (mouvement uniforme imposé).
> La simulation RK4 lève cette hypothèse en résolvant le régime transitoire complet.

---

## Simulation numérique / Numerical Simulation

### Méthode / Method

L'équation du second ordre est réduite à un système du premier ordre `s = (θ, ω)` :

```
ds/dt = f(s) = ( ω,  [M(θ) − b·ω − Cs·sgn(ω)] / I )
```

Intégration par **Runge-Kutta ordre 4 (RK4)**, pas `Δt = 0.012 s` :

```
s_{n+1} = sₙ + (Δt/6)·(k₁ + 2k₂ + 2k₃ + k₄)
```

### Configurations testées / Tested Configurations

| Config | n | r (m) | k (N/m) | l₀ (m) | m (kg) | α | Triple condition |
|---|---|---|---|---|---|---|---|
| C1 | 2 | 0.30 | 30 | 0.40 | 5 | 0.50 | (C1)✓ (C2)~ (C3)✓ |
| C2 | 4 | 0.30 | 30 | 0.40 | 5 | 0.50 | (C1)✓ (C2)✓ (C3)✓ |
| C3 | 6 | 0.30 | 50 | 0.50 | 8 | 0.40 | (C1)✓ (C2)✓ (C3)✓ |
| C4 | 8 | 0.40 | 80 | 0.60 | 10 | 0.50 | (C1)✓ (C2)✓ (C3)✓ |

> **Observation principale :** `n ≥ 4` avec `l₀ > r` est la configuration minimale garantissant la triple condition simultanée.

---

## Résultats / Results

- La simulation confirme qu'un moment moteur non nul et de signe constant est produit pour `n ≥ 3`
- La constance de `g(t)` s'améliore avec `n`
- La configuration géométrique particulière (tous les ressorts en *P*) crée une **asymétrie dynamique favorable**
- `n ≥ 4` avec `l₀ > r` est la configuration minimale recommandée

### Effets des paramètres / Parameter Effects

| Paramètre | Effet sur la triple condition |
|---|---|
| `n` ↑ | `g(t)` plus stable, oscillations réduites, (C2) mieux vérifiée |
| `k` ↑ | Moments plus grands, (C1) plus robuste |
| `l₀ > r` | Assure plusieurs ressorts en extension simultanée |
| `α` ↓ | Meilleur dépassement des positions d'équilibre |
| `b > 0` | Décroissance progressive de `ω` si `M ≪ b·ω` |

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

### Fonctionnalités du simulateur / Simulator Features

- Sliders interactifs pour tous les paramètres : `n`, `r`, `k`, `l₀`, `m`, `α`, `b`, `Cs`, `ω₀`
- Intégration RK4 avec pas configurable
- Graphiques en temps réel : `ω(t)`, `g(t)`, `Ec(t)`, `Ee(t)`, portrait de phase `ω(θ)`
- Triple condition évaluée quantitativement
- Analyse comparative `n = 2` à `8` en un clic
- Export des données au format CSV

---

## Structure du dépôt / Repository Structure

```
systeme-de-guielle/
├── README.md                          # Ce fichier / This file
├── LICENSE                            # Licence CC BY 4.0
├── guielle_simulator.py               # Simulateur Streamlit principal
├── docs/
│   └── Systeme_Guielle_Rapport.docx  # Rapport de recherche bilingue
└── notebooks/
    └── guielle_analysis.ipynb         # (optionnel) Jupyter Notebook d'analyse
```

---

## Nomenclature

| Symbole | Définition |
|---|---|
| `n` | Nombre de ressorts / Number of springs |
| `r` | Rayon de la roue (m) / Wheel radius (m) |
| `m` | Masse de la roue (kg) / Wheel mass (kg) |
| `k` | Constante de raideur (N/m) / Spring stiffness (N/m) |
| `l₀` | Longueur naturelle des ressorts (m) / Natural spring length (m) |
| `α` | Coefficient d'inertie (libre) / Free inertia coefficient |
| `I` | Moment d'inertie = α·m·r² (kg·m²) |
| `θ` | Angle de rotation (rad) |
| `ω` | Vitesse angulaire (rad/s) |
| `βᵢ` | Angle fixe du point Aᵢ (rad) |
| `τᵢ` | Angle caractéristique τᵢ = βᵢ + θ |
| `ε` | Rapport géométrique ε = 1 + l₀/r |
| `xᵢ` | Élongation du ressort i (m) |
| `g(t)` | Moment adimensionnel (sans unité) |
| `b` | Coefficient de frottement visqueux |
| `Cs` | Couple résistant externe (N·m) |

---

## Licence / License

Ce projet est distribué sous licence **Creative Commons Attribution 4.0 International (CC BY 4.0)**.

This project is distributed under the **Creative Commons Attribution 4.0 International (CC BY 4.0)** license.

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

Vous êtes libre de / You are free to:
- **Partager** — copier et redistribuer le matériel / **Share** — copy and redistribute
- **Adapter** — remixer, transformer et construire à partir du matériel / **Adapt** — remix, transform, and build upon

**Sous les conditions suivantes / Under the following terms:**
> Vous devez créditer l'œuvre, intégrer un lien vers la licence et indiquer si des modifications ont été effectuées.
> You must give appropriate credit, provide a link to the license, and indicate if changes were made.

**© 2025 Herolde GUIELLE — Tous droits réservés sur le concept original / All rights reserved on the original concept.**

---

*Système de Guielle · Rapport de Recherche · H. GUIELLE · 2025*
