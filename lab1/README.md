**Execution**:

  You can execute the program with the image you want by changing the variable `imgName`. Assuming you have all the images in `./images_tiff/`, you can also change this by modifying the `imgPath` variable. The final result will be saved in `./images_output/` with the formats `.jpg` (with a lot of compression) and `.png` (with better quality).

  In the code, we automatically use the manual white balancing method. To change the method, you can modify the method on line 149.

  Also, you have to consider that we have used the Parallel Computing Toolbox for doing the calculations on the GPU, but the code will use the CPU if you don't have this Toolbox installed.

**Authors**:  

    Juan Lorente Guarnieri 816020

    Hugo Mateo Trejo 816678