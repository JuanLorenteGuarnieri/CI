# Lab 4 Instructions

To execute this lab, follow these steps:

1. **Prepare the Data**:  
  Add all the datasets provided for the lab to the folder `/data`.

2. **Configuration**:  
  The MATLAB scripts use the `config.m` script. You can define the parameters here that all the scripts will use.

3. **Main Script**:  
  - We provide `lab4.m`, which performs the back-projection reconstruction.
  - Use `show_volume.m` to load the volume previously saved (with `lab4.m`) and visualize it. This allows you to see the different results without waiting for the main script to reconstruct volumes just to display a result.

4. **Results**:  
  - A `results.csv` file is provided, which contains all the execution times for the different parameters.
  - All the volume results are saved in the `/results` folder. These can be displayed using `show_volume.m` by modifying the parameters in `config.m`.

---

**Authors**:  
- Juan Lorente Guarnieri (816020)  
- Hugo Mateo Trejo (816678)