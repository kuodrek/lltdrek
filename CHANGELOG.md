# LLTDREK CHANGELOG

All notable changes to this project will be documented in this file.
- Following [text](https://keepachangelog.com/en/1.0.0/)

## [Unreleased]

## 1.0.0

### TODO

- Transform post_processing methods into class methods
- Create SimulationResults object
- Create AirfoilDatabase object
  - Maybe use a dataframe library
- Create input data interface
  - Add airfoil data by yaml files
- move setup_airfoil_data method from Wing to WingPool
- Wing.generate_mesh() should automatically be called upon creation


## 0.2.1

### Added

- Add build tests for multiple python versions
    - Some users were reporting problems using python with these versions. Might be a dependency problem.
- Get dependencies from requirements.txt file in pyproject.toml

### Changed

- Also update requirements
