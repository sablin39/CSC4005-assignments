# Description

In the project, you are required to calculate the Mandelbrot Set parallelly.

The user will provide the width `Y_RESN` and height `Y_RESN` of the resulting graph as well

as the max iteration value `max_iteration`.

The iteration of each point $z$ in complex plane is

$$
z_{n+1}=z_{n}^2+c

$$

Here $c$ is provided by  This iteration will end when $|z| \ge 2$ or $n+1=max\_iteration$, and the "color" of point $(z.real,z.imag)$ in the resulting graph is $n/max\_iteration$.

# Requirements

You need to implement two versions of the tasks, which are MPI version and a Pthread version. And hand in the codes for these two versions in two seperatecode files.

- In your submit code, it should display an image with size of 800 × 800.
- Include the results in your report. You need to specify the command line about how to compile and run your program.
- You need to compare the performance of different implementation andconfigurations in your report.
  - a. The number of processes or threads used in the program
  - b. MPI vs Sequential vs Pthread
  - c. Size of the output images (three different sizes ranging from small, medium to large)
- You need to include two figures describing the structure of your MPI program and Pthread program.
