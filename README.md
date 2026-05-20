# A Unified Pore-Scale Multiphysics Model for the Integrated Soot Transport-Deposition-Oxidation in Catalytic Diesel Particulate Filters

[![arXiv](https://img.shields.io/badge/arXiv-2512.22230-b31b1b.svg)](https://arxiv.org/abs/2512.22230)
[![DOI](https://img.shields.io/badge/DOI-10.1063%2F5.0321009-red)](https://doi.org/10.1063/5.0321009)


This user-defined function (UDF) code implements a unified pore-scale multiphysics model based on the Eulerian-Lagrangian framework to comprehensively resolve the transport, deposition, and oxidation of soot within a catalytic diesel particulate filter (CDPF).  The soot oxidation model is derived based on the fundamental chemical kinetics and thermochemistry, enabling comprehensive coverage of all CDPF operating conditions. The particle-wall interaction is modeled based on the elastic deformation and surface adhesion theories. More details on our model and implementation can be found in [doi:10.1063/5.0321009](https://doi.org/10.1063/5.0321009).

- Authors: Yujing Zhang, Yunhua Zhang, Liang Fang, Diming Lou, Piqiang Tan, and Zhiyuan Hu

<img src=".\data\solver_strategy.png" width=70%>

## File structure

- **data**

  - `data/cdpf_geo`: This folder contains the CDPF porous structure used in this study. Details on its generation process can be referred to Sec. II.A in our publication.

  - `data/cdpf_sim`: This folder contains all post-processing data related to the simulation of a CDPF porous media (See Sec. IV in our paper).
  - `data/validation_and_verification`: This directory provides raw data of our benchmark simulations. (See Sec. III & IV.A in our paper).

- **src**

  - `src/support`: contains all `lib` files needed.
  - `src/CMakeList.txt`: top-level `CMake` file that controls other `CMake` files in `src/support`.
  - `src/chemicalReactions_38cb0c25.c`: UDF codes for the soot oxidation model.
  - `src/impactionModel_bdc776e3.c`: UDF codes for the particle deposit-rebound model.

---

## Requirements

This code is developed based on the `ANSYS Fluent 2024R2` (OS: Windows). It also requires a compiler to compile these UDF codes, which I recommend `Visual Studio 14.0`. For those who want to make further development, I recommend to use `CLion`. More details can be found in [![github](https://img.shields.io/badge/github-repo-blue?logo=github)](https://github.com/Izumiko/udf?tab=readme-ov-file).

---

## References

[1] Zhang, Y., Zhang, Y., Fang, L., Lou, D., Tan, P., & Hu, Z. (2026). A unified pore-scale multiphysics model for the integrated soot transport-deposition-oxidation in catalytic diesel particulate filters. Physics of Fluids, 38(3).

Feel free to copy and modify our code (`PR` is needed). If you find this useful, please cite our paper and consider give a little star :star: to this project. For any questions or comments, please contact [zhangyujing@tongji.edu.cn] or [zhangyunhua@tongji.edu.cn]. 





