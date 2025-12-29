# STATS 700: Phylodynamic Inference

> **Course**: STATS 700 - Topics in Applied Statistics: Phylodynamic Inference
> **Semester**: Fall 2025
> **Instructors**: Aaron A. King & Edward L. Ionides
> **Institution**: University of Michigan

This repository contains comprehensive reading notes, summaries, and Quarto documents for the STATS 700 course on phylodynamic inference. The course explores the mathematical and statistical foundations linking genealogies, epidemic dynamics, and evolutionary processes.

Official course materials: [ionides/700f25](https://github.com/ionides/700f25)

---

## Table of Contents

- [Weekly Content Overview](#weekly-content-overview)
- [Project Structure](#project-structure)
- [How to View the Notes](#how-to-view-the-notes)
- [Course Topics](#course-topics)
- [Final Project](#final-project)
- [License](#license)

---

## Weekly Content Overview

### Week 2: The Coalescent (Kingman 1982)
**Topic**: Foundations of coalescent theory
**Key Concepts**:
- Continuous-time Markov chain (CTMC) on partitions
- Theorem 1: Factorization into pure death process and jump chain
- Waiting times and merge probabilities
- Independence of time and structure in coalescent processes

**Files**: [`weeks/week-02/`](weeks/week-02/)

---

### Week 3: Phylodynamics of Infectious Disease Epidemics (Volz et al. 2009)
**Topic**: Linking epidemic dynamics to genealogies
**Key Concepts**:
- Forward-time SIR dynamics and backward-time coalescent
- Ancestor function $A(t;T)$ and coalescent rate derivation
- Cluster size distribution (CSD)
- Likelihood from coalescent times
- Effect of sampling fraction on tree shapes

**Files**: [`weeks/week-03/`](weeks/week-03/)

---

### Week 4: Sampling-through-time in Birth-Death Trees (Stadler 2010)
**Topic**: Birth-death models with heterochronous sampling
**Key Concepts**:
- Birth-death (BD) trees vs. coalescent: pros and cons
- Pruning and sampled (reconstructed) trees
- Master equations for $p_0(t)$ and $p_1(t)$
- Theorem 3.5: Exact likelihood formula
- Conditioning on $n$ samples and the $1/(n\lambda)$ weight

**Files**: [`weeks/week-04/`](weeks/week-04/)

---

### Week 5: Markov Genealogy Processes (King, Lin & Ionides 2022)
**Topic**: Structured genealogy processes and exact likelihoods
**Key Concepts**:
- From population dynamics to genealogies
- Inventory, node sequences, and colored-balls representation
- Pruning and visible genealogy $V_t$
- Theorem 1: Conditioning on history factorization
- Theorem 2: DMZ-style unnormalized nonlinear filter
- Specializations to linear birth-death-sampling, SIR, and SIRS

**Files**: [`weeks/week-05/`](weeks/week-05/)

---

### Week 6: Viral Phylodynamics Overview
**Topic**: Applied phylodynamics in virology
**Key Concepts**:
- Tree-to-process heuristics (population size change, structure, selection)
- Dating origins and estimating $R_0$
- Phylogeography and migration inference
- Molecular clocks and Bayesian phylogenetics
- Case studies: Influenza A/H3N2, HIV

**Files**: [`weeks/week-06/`](weeks/week-06/)

---

### Week 7: Complex Population Dynamics Under Neutrality (Volz 2012)
**Topic**: General coalescent framework for complex dynamics
**Key Concepts**:
- Coalescent rate: $\lambda_2(t) = 2f(t)/Y(t)^2$
- Faster-than-exponential (FTE) and slower-than-exponential (STE) growth
- Skyline bias under non-constant birth rates
- Structured populations with time-varying birth/migration matrices
- Invisible transmission states
- Exact genealogy likelihood and simulation algorithms

**Files**: [`weeks/week-07/`](weeks/week-07/)

---

### Week 8: Multitype Population Trajectories (Vaughan & Stadler 2025)
**Topic**: Bayesian inference of type-specific population sizes
**Key Concepts**:
- Linear multitype birth-death models
- Stochastic mapping of ancestral types
- Particle filtering for population trajectories
- BDMM-Prime implementation in BEAST 2
- MERS-CoV case study: camel-human spillovers
- Inferring events in the full population vs. sampled lineages

**Files**: [`weeks/week-08/`](weeks/week-08/)

---

### Week 9: Exact Phylodynamic Likelihood (King, Lin & Ionides 2025)
**Topic**: Structured Markov genealogy processes - complete framework
**Key Concepts**:
- Discretely structured, time-inhomogeneous Markov jump processes
- Exact likelihoods for pruned and obscured genealogies
- Forward filter equations and local multipliers $\varphi_u$
- Binomial ratios, production $r_u$, and saturation $s$
- Event-driven SMC (particle filter) implementation
- Specialization to Kingman coalescent and Stadler's birth-death model

**Files**: [`weeks/week-09/`](weeks/week-09/)

---

## Project Structure

```
stats-700-note/
├── README.md                    # This file
├── weeks/                       # Weekly course materials
│   ├── week-02/                # Kingman 1982
│   │   ├── class 2.qmd
│   │   ├── class 2.html
│   │   └── class 2_files/
│   ├── week-03/                # Volz et al. 2009
│   ├── week-04/                # Stadler 2010
│   ├── week-05/                # King, Lin & Ionides 2022
│   ├── week-06/                # Viral Phylodynamics
│   ├── week-07/                # Volz 2012
│   ├── week-08/                # Vaughan & Stadler 2025
│   └── week-09/                # King, Lin & Ionides 2025
└── project/                     # Final project materials
    ├── 700_proj_final.pdf
    ├── STATS700_Project.pdf
    └── code.R
```

---

## How to View the Notes

### Option 1: View HTML Files Directly

Each week contains pre-rendered HTML files that can be opened directly in your browser:

```bash
# Example: Open Week 3 notes
open weeks/week-03/class\ 3.html
```

### Option 2: Build from Quarto Source

If you want to modify or rebuild the notes:

1. Install [Quarto](https://quarto.org/)
2. Navigate to a week's directory:
   ```bash
   cd weeks/week-03
   ```
3. Render the document:
   ```bash
   # To HTML
   quarto render "class 3.qmd" --to html

   # To PDF (requires LaTeX)
   quarto render "class 3.qmd" --to pdf
   ```

---

## Course Topics

The course covers the following major themes:

1. **Coalescent Theory**: Mathematical foundations (Kingman 1982)
2. **Phylodynamics**: Linking genealogies to epidemic/evolutionary dynamics
3. **Birth-Death Models**: Alternative framework with sampling-through-time
4. **Structured Populations**: Migration, demes, and multitype processes
5. **Computational Methods**: Particle filtering, MCMC, stochastic mapping
6. **Applications**: Infectious diseases (HIV, MERS-CoV, Influenza), conservation biology

### Key Papers Covered

- Kingman (1982) - The Coalescent
- Volz et al. (2009) - Phylodynamics of infectious disease epidemics
- Stadler (2010) - Sampling-through-time in birth-death trees
- Volz (2012) - Complex population dynamics and the coalescent
- King, Lin & Ionides (2022/2025) - Markov genealogy processes
- Vaughan & Stadler (2025) - Bayesian phylodynamic inference of multitype population trajectories

---

## Final Project

The `project/` directory contains materials from the course final project, including:

- **700_proj_final.pdf**: Final project report
- **STATS700_Project.pdf**: Project presentation/additional materials
- **code.R**: Analysis code and simulations

---

## License

These notes are for personal study and academic discussion. Please consult the original papers and the [official course repository](https://github.com/ionides/700f25) for authoritative materials.

All course materials are subject to the University of Michigan's academic policies. Paper citations and references belong to their respective authors and publishers.

---

## Acknowledgments

- **Instructors**: [Aaron A. King](https://kinglab.eeb.lsa.umich.edu/) and [Edward L. Ionides](https://ionides.github.io/)
- **Course Repository**: [ionides/700f25](https://github.com/ionides/700f25)
- **Institution**: Department of Statistics, University of Michigan

---

*Last updated: December 2025*
